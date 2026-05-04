from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
from shapely import affinity
from shapely.geometry import MultiPolygon, Polygon, mapping, shape
from shapely.ops import unary_union

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from histoseg.threed import (  # noqa: E402
    CodaImageRegistrationConfig,
    estimate_coda_image_registration,
    hard_align_geojson,
    hash_alignment_manifest,
)
from histoseg.threed.cell_cloud import canonical_json  # noqa: E402

try:
    from shapely import contains_xy as _contains_xy
except Exception:  # pragma: no cover - only used on older Shapely versions.
    _contains_xy = None


@dataclass(frozen=True)
class SimilaritySeed:
    origin_x: float
    origin_y: float
    rotation_degrees: float
    translate_x: float
    translate_y: float
    scale: float = 1.0


@dataclass(frozen=True)
class StressMetric:
    trial: int
    chaos_rotation_degrees: float
    chaos_translate_x: float
    chaos_translate_y: float
    raw_iou: float
    contour_seed_iou: float
    coda_seed_iou: float
    iou_gain_over_contour_seed: float
    coda_seed_success: bool
    contour_seed_success: bool
    coda_rotation_degrees: float
    expected_inverse_rotation_degrees: float
    angle_error_mod180_degrees: float
    angle_error_full_circle_degrees: float
    phase_shift_y: float
    phase_shift_x: float
    transform_hash_changed_on_angle_perturbation: bool


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Stress-test the CODA-inspired Radon/phase-correlation seed against "
            "large synthetic contour rotations and translations."
        )
    )
    parser.add_argument("--out-dir", default="output/coda_backend_stress_test")
    parser.add_argument("--trials", type=int, default=3)
    parser.add_argument("--random-seed", type=int, default=7)
    parser.add_argument("--rotation-degrees", type=float, default=None)
    parser.add_argument("--rotation-min", type=float, default=30.0)
    parser.add_argument("--rotation-max", type=float, default=180.0)
    parser.add_argument("--translation-fraction", type=float, default=0.35)
    parser.add_argument("--raster-size", type=int, default=512)
    parser.add_argument("--angle-step", type=float, default=1.0)
    parser.add_argument("--phase-upsample-factor", type=int, default=1)
    parser.add_argument("--seed-success-iou", type=float, default=0.5)
    parser.add_argument(
        "--baseline-maxiter",
        type=int,
        default=0,
        help="Use 0 for the naive contour PCA/centroid seed; increase to test full hard alignment.",
    )
    args = parser.parse_args()

    if args.trials < 1:
        raise ValueError("--trials must be at least 1.")
    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    fixed_geom = make_synthetic_tissue_geometry()
    fixed_geojson = out_dir / "fixed_ground_truth.geojson"
    write_geojson(fixed_geojson, fixed_geom)
    ground_truth_iou = iou(fixed_geom, fixed_geom)

    rng = np.random.default_rng(int(args.random_seed))
    bounds = fixed_geom.bounds
    span = max(bounds[2] - bounds[0], bounds[3] - bounds[1])
    metrics: list[StressMetric] = []
    for trial in range(1, int(args.trials) + 1):
        rotation = (
            float(args.rotation_degrees)
            if args.rotation_degrees is not None
            else float(rng.uniform(float(args.rotation_min), float(args.rotation_max)))
        )
        dx = float(rng.choice([-1.0, 1.0]) * rng.uniform(0.3, 1.0) * args.translation_fraction * span)
        dy = float(rng.choice([-1.0, 1.0]) * rng.uniform(0.3, 1.0) * args.translation_fraction * span)
        moving_geom = affinity.rotate(fixed_geom, rotation, origin="center")
        moving_geom = affinity.translate(moving_geom, xoff=dx, yoff=dy)
        moving_geojson = out_dir / f"trial_{trial:02d}_moving_chaos.geojson"
        write_geojson(moving_geojson, moving_geom)

        contour_iou = run_contour_seed_baseline(
            fixed_geojson=fixed_geojson,
            moving_geojson=moving_geojson,
            out_dir=out_dir / f"trial_{trial:02d}_contour_seed",
            maxiter=int(args.baseline_maxiter),
        )
        coda_seed, coda_payload = estimate_coda_seed(
            fixed_geom,
            moving_geom,
            raster_size=int(args.raster_size),
            angle_step=float(args.angle_step),
            phase_upsample_factor=int(args.phase_upsample_factor),
        )
        coda_aligned = apply_similarity(moving_geom, coda_seed)
        coda_iou = iou(fixed_geom, coda_aligned)
        write_geojson(out_dir / f"trial_{trial:02d}_coda_seed_aligned.geojson", coda_aligned)

        expected_inverse = _normalize_degrees_180(-rotation)
        metric = StressMetric(
            trial=trial,
            chaos_rotation_degrees=rotation,
            chaos_translate_x=dx,
            chaos_translate_y=dy,
            raw_iou=iou(fixed_geom, moving_geom),
            contour_seed_iou=contour_iou,
            coda_seed_iou=coda_iou,
            iou_gain_over_contour_seed=coda_iou - contour_iou,
            coda_seed_success=coda_iou >= float(args.seed_success_iou),
            contour_seed_success=contour_iou >= float(args.seed_success_iou),
            coda_rotation_degrees=coda_seed.rotation_degrees,
            expected_inverse_rotation_degrees=expected_inverse,
            angle_error_mod180_degrees=_angle_error_mod180(
                coda_seed.rotation_degrees,
                expected_inverse,
            ),
            angle_error_full_circle_degrees=_angle_error_full_circle(
                coda_seed.rotation_degrees,
                -rotation,
            ),
            phase_shift_y=float(coda_payload["phase_shift_y"]),
            phase_shift_x=float(coda_payload["phase_shift_x"]),
            transform_hash_changed_on_angle_perturbation=hash_changes_on_angle_perturbation(
                coda_seed,
                epsilon=0.0001,
            ),
        )
        metrics.append(metric)

    write_metrics(out_dir, metrics)
    summary = {
        "ground_truth_iou": ground_truth_iou,
        "trials": len(metrics),
        "seed_success_iou": float(args.seed_success_iou),
        "contour_seed_success_rate": float(np.mean([m.contour_seed_success for m in metrics])),
        "coda_seed_success_rate": float(np.mean([m.coda_seed_success for m in metrics])),
        "mean_iou_gain_over_contour_seed": float(
            np.mean([m.iou_gain_over_contour_seed for m in metrics])
        ),
        "mean_angle_error_mod180_degrees": float(
            np.mean([m.angle_error_mod180_degrees for m in metrics])
        ),
        "mean_angle_error_full_circle_degrees": float(
            np.mean([m.angle_error_full_circle_degrees for m in metrics])
        ),
        "hash_invalidation_all_passed": all(
            m.transform_hash_changed_on_angle_perturbation for m in metrics
        ),
        "radon_periodicity_note": (
            "Radon tissue-mask signatures are 180-degree periodic. Use the full-circle "
            "angle error and IoU together to catch half-turn ambiguity in extreme rotations."
        ),
    }
    (out_dir / "stress_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"Wrote metrics to {out_dir / 'stress_metrics.csv'}")


def make_synthetic_tissue_geometry() -> Polygon | MultiPolygon:
    outer = [
        (-70, -45),
        (-25, -70),
        (40, -62),
        (85, -30),
        (70, 20),
        (35, 58),
        (-18, 72),
        (-65, 42),
        (-92, -8),
    ]
    hole = [(-18, -16), (18, -12), (25, 18), (-12, 26), (-30, 6)]
    body = Polygon(outer, holes=[hole])
    lobe = affinity.rotate(Polygon([(35, 50), (95, 62), (105, 96), (45, 105)]), 12, origin=(45, 50))
    tail = Polygon([(-95, -8), (-135, -25), (-126, 18), (-88, 26)])
    notch = Polygon([(56, -46), (98, -50), (86, -18), (62, -24)])
    return unary_union([body, lobe, tail]).difference(notch).buffer(0)


def write_geojson(path: Path, geom: Any) -> None:
    payload = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "properties": {"structure": "Synthetic tissue"},
                "geometry": mapping(geom),
            }
        ],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")


def run_contour_seed_baseline(
    *,
    fixed_geojson: Path,
    moving_geojson: Path,
    out_dir: Path,
    maxiter: int,
) -> float:
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = hard_align_geojson(
        fixed_geojson=fixed_geojson,
        moving_geojson=moving_geojson,
        output_geojson=out_dir / "moving_contour_seed_aligned.geojson",
        summary_json=out_dir / "hard_similarity_alignment.json",
        maxiter=maxiter,
        multistart=False,
        overwrite=True,
    )
    return float(summary["union_iou_after_hard"])


def estimate_coda_seed(
    fixed_geom: Any,
    moving_geom: Any,
    *,
    raster_size: int,
    angle_step: float,
    phase_upsample_factor: int,
) -> tuple[SimilaritySeed, dict[str, float]]:
    fixed_mask, moving_mask, metadata = rasterize_pair(fixed_geom, moving_geom, raster_size=raster_size)
    result = estimate_coda_image_registration(
        fixed_mask,
        moving_mask,
        CodaImageRegistrationConfig(
            angle_step=angle_step,
            phase_upsample_factor=phase_upsample_factor,
        ),
    )
    bounds = metadata["square_bounds"]
    native_units_per_pixel = float(metadata["native_units_per_pixel"])
    seed = SimilaritySeed(
        origin_x=(bounds[0] + bounds[2]) / 2.0,
        origin_y=(bounds[1] + bounds[3]) / 2.0,
        rotation_degrees=float(result.rotation.rotation_degrees),
        translate_x=float(result.translation.shift_x) * native_units_per_pixel,
        translate_y=-float(result.translation.shift_y) * native_units_per_pixel,
    )
    payload = {
        "radon_rotation_degrees": float(result.rotation.rotation_degrees),
        "phase_shift_y": float(result.translation.shift_y),
        "phase_shift_x": float(result.translation.shift_x),
    }
    return seed, payload


def rasterize_pair(
    fixed_geom: Any,
    moving_geom: Any,
    *,
    raster_size: int,
    padding_fraction: float = 0.08,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    bounds = square_bounds(unary_union([fixed_geom, moving_geom]), padding_fraction=padding_fraction)
    x_min, y_min, x_max, y_max = bounds
    x_values = np.linspace(x_min, x_max, int(raster_size), dtype=float)
    y_values = np.linspace(y_max, y_min, int(raster_size), dtype=float)
    xx, yy = np.meshgrid(x_values, y_values)
    fixed_mask = geometry_contains_grid(fixed_geom, xx, yy)
    moving_mask = geometry_contains_grid(moving_geom, xx, yy)
    return fixed_mask, moving_mask, {
        "square_bounds": bounds,
        "native_units_per_pixel": (x_max - x_min) / max(int(raster_size) - 1, 1),
    }


def square_bounds(geom: Any, *, padding_fraction: float) -> tuple[float, float, float, float]:
    x_min, y_min, x_max, y_max = map(float, geom.bounds)
    side = max(x_max - x_min, y_max - y_min, 1.0)
    side *= 1.0 + 2.0 * float(padding_fraction)
    center_x = (x_min + x_max) / 2.0
    center_y = (y_min + y_max) / 2.0
    half = side / 2.0
    return (center_x - half, center_y - half, center_x + half, center_y + half)


def geometry_contains_grid(geom: Any, xx: np.ndarray, yy: np.ndarray) -> np.ndarray:
    if _contains_xy is not None:
        return np.asarray(_contains_xy(geom, xx, yy), dtype=bool)
    flat = [geom.contains(shape({"type": "Point", "coordinates": [x, y]})) for x, y in zip(xx.ravel(), yy.ravel())]
    return np.asarray(flat, dtype=bool).reshape(xx.shape)


def apply_similarity(geom: Any, seed: SimilaritySeed) -> Any:
    origin = (seed.origin_x, seed.origin_y)
    transformed = affinity.rotate(geom, seed.rotation_degrees, origin=origin)
    transformed = affinity.scale(transformed, xfact=seed.scale, yfact=seed.scale, origin=origin)
    return affinity.translate(transformed, xoff=seed.translate_x, yoff=seed.translate_y)


def iou(a: Any, b: Any) -> float:
    intersection = a.intersection(b).area
    union = a.union(b).area
    return 0.0 if union <= 0 else float(intersection / union)


def write_metrics(out_dir: Path, metrics: list[StressMetric]) -> None:
    path = out_dir / "stress_metrics.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(metrics[0]).keys()))
        writer.writeheader()
        for metric in metrics:
            writer.writerow(asdict(metric))


def hash_changes_on_angle_perturbation(seed: SimilaritySeed, *, epsilon: float) -> bool:
    manifest = {
        "schema_version": 2,
        "coordinate_contract": "histoseg_3d_cell_cloud",
        "slices": [{"hard_registration": {"registration_backend": "coda-image"}, "hard_transform": asdict(seed)}],
    }
    perturbed = json.loads(canonical_json(manifest))
    perturbed["slices"][0]["hard_transform"]["rotation_degrees"] += float(epsilon)
    return hash_alignment_manifest(manifest) != hash_alignment_manifest(perturbed)


def _normalize_degrees_180(angle: float) -> float:
    normalized = math.fmod(float(angle), 180.0)
    if normalized < 0:
        normalized += 180.0
    if normalized >= 90.0:
        normalized -= 180.0
    return normalized


def _angle_error_mod180(estimated: float, expected: float) -> float:
    return abs(_normalize_degrees_180(float(estimated) - float(expected)))


def _angle_error_full_circle(estimated: float, expected: float) -> float:
    delta = math.fmod(float(estimated) - float(expected), 360.0)
    if delta < -180.0:
        delta += 360.0
    if delta > 180.0:
        delta -= 360.0
    return abs(delta)


if __name__ == "__main__":
    main()
