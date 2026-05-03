"""Synthetic validation for HistoSeg SDF metric contracts.

The script is a consumer of the public HistoSeg API. It builds small synthetic
3D binary masks, evaluates ``compute_hotspot_sdf_metrics`` under controlled
anisotropic spacing, and writes CSV/PNG artifacts for manuscript Figure 2.
It does not modify or reimplement package internals.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from histoseg.threed import compute_hotspot_sdf_metrics


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT_DIR = REPO_ROOT / "reproducibility" / "validation" / "results"
DEFAULT_Z_FACTORS = (1.0, 2.0, 5.0, 10.0)


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    out_dir = Path(args.out_dir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)

    z_factors = tuple(float(value) for value in args.z_factors)
    truth_table = build_truth_table()
    empty_contract = build_empty_mask_contract()
    sweep = build_anisotropy_sweep(z_factors)

    truth_csv = out_dir / "sdf_truth_table.csv"
    empty_csv = out_dir / "empty_mask_contract.csv"
    sweep_csv = out_dir / "anisotropy_sensitivity.csv"
    fraction_csv = out_dir / "anisotropy_sensitivity_fraction_inside.csv"
    signed_csv = out_dir / "anisotropy_sensitivity_signed_distance.csv"
    png = out_dir / "anisotropy_sensitivity.png"
    manifest_json = out_dir / "validation_manifest.json"

    truth_table.to_csv(truth_csv, index=False)
    empty_contract.to_csv(empty_csv, index=False)
    sweep.to_csv(sweep_csv, index=False)
    sweep[
        [
            "z_spacing_factor",
            "z_spacing_um",
            "fraction_inside_structure",
            "fraction_inside_delta",
            "fraction_inside_percent_delta",
        ]
    ].to_csv(fraction_csv, index=False)
    sweep[
        [
            "z_spacing_factor",
            "z_spacing_um",
            "median_signed_distance_um",
            "median_signed_delta_um",
            "p05_signed_distance_um",
            "p05_signed_delta_um",
            "p95_signed_distance_um",
            "p95_signed_delta_um",
            "signed_summary_mae_um",
        ]
    ].to_csv(signed_csv, index=False)
    plot_anisotropy_sweep(sweep, png)

    manifest = {
        "schema_version": 1,
        "validation": "sdf_anisotropy_sweep",
        "description": (
            "Synthetic validation of compute_hotspot_sdf_metrics under "
            "anisotropic z sampling."
        ),
        "z_factors": list(z_factors),
        "outputs": {
            "sdf_truth_table_csv": str(truth_csv),
            "empty_mask_contract_csv": str(empty_csv),
            "anisotropy_sensitivity_csv": str(sweep_csv),
            "anisotropy_sensitivity_fraction_inside_csv": str(fraction_csv),
            "anisotropy_sensitivity_signed_distance_csv": str(signed_csv),
            "anisotropy_sensitivity_png": str(png),
        },
        "summary": summarize_sweep(sweep),
    }
    manifest_json.write_text(json.dumps(manifest, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, ensure_ascii=True))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run synthetic HistoSeg SDF anisotropy validation."
    )
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--z-factors", nargs="*", type=float, default=list(DEFAULT_Z_FACTORS))
    return parser


def build_truth_table() -> pd.DataFrame:
    structure_mask = np.zeros((5, 5, 5), dtype=bool)
    hotspot_mask = np.zeros_like(structure_mask)
    structure_mask[2, 2, 2] = True
    hotspot_mask[2, 2, 2] = True
    hotspot_mask[2, 2, 3] = True
    hotspot_mask[3, 2, 2] = True

    metrics = compute_hotspot_sdf_metrics(
        hotspot_mask,
        structure_mask,
        spacing_zyx_um=(5.0, 1.0, 1.0),
    )
    expected = {
        "n_hotspot_voxels": 3,
        "min_signed_distance_um": -1.0,
        "median_signed_distance_um": 1.0,
        "max_signed_distance_um": 5.0,
        "mean_unsigned_distance_um": 2.0,
        "fraction_inside_structure": 1.0 / 3.0,
        "fraction_touching_or_inside_structure": 1.0 / 3.0,
    }
    rows: list[dict[str, Any]] = []
    for key, expected_value in expected.items():
        observed = float(metrics[key])
        rows.append(
            {
                "case": "5x5x5_anisotropic_truth_table",
                "metric": key,
                "observed": observed,
                "expected": float(expected_value),
                "absolute_error": abs(observed - float(expected_value)),
                "passed": bool(math.isclose(observed, float(expected_value), rel_tol=1e-9, abs_tol=1e-9)),
            }
        )
    return pd.DataFrame(rows)


def build_empty_mask_contract() -> pd.DataFrame:
    structure_mask = np.zeros((5, 5, 5), dtype=bool)
    hotspot_mask = np.zeros_like(structure_mask)
    structure_mask[2, 2, 2] = True

    empty_hotspot = compute_hotspot_sdf_metrics(
        hotspot_mask,
        structure_mask,
        spacing_zyx_um=(5.0, 1.0, 1.0),
    )

    empty_structure_mask = np.zeros((5, 5, 5), dtype=bool)
    one_voxel_hotspot = np.zeros_like(empty_structure_mask)
    one_voxel_hotspot[2, 2, 2] = True
    empty_structure = compute_hotspot_sdf_metrics(
        one_voxel_hotspot,
        empty_structure_mask,
        spacing_zyx_um=(5.0, 1.0, 1.0),
    )

    rows = []
    for case, metrics in (
        ("empty_hotspot", empty_hotspot),
        ("empty_structure", empty_structure),
    ):
        distance_keys = [
            "min_signed_distance_um",
            "median_signed_distance_um",
            "mean_signed_distance_um",
            "max_signed_distance_um",
            "mean_unsigned_distance_um",
        ]
        rows.append(
            {
                "case": case,
                "n_hotspot_voxels": int(metrics["n_hotspot_voxels"]),
                "fraction_inside_structure": float(metrics["fraction_inside_structure"]),
                "fraction_touching_or_inside_structure": float(
                    metrics["fraction_touching_or_inside_structure"]
                ),
                "distance_summaries_are_nan": all(
                    math.isnan(float(metrics[key])) for key in distance_keys
                ),
                "contract_passed": bool(
                    float(metrics["fraction_inside_structure"]) == 0.0
                    and float(metrics["fraction_touching_or_inside_structure"]) == 0.0
                    and all(math.isnan(float(metrics[key])) for key in distance_keys)
                ),
            }
        )
    return pd.DataFrame(rows)


def build_anisotropy_sweep(z_factors: Sequence[float]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    baseline: dict[str, float] | None = None
    for factor in z_factors:
        structure_mask, hotspot_mask = synthetic_ellipsoid_scene(z_spacing_um=float(factor))
        metrics = compute_hotspot_sdf_metrics(
            hotspot_mask,
            structure_mask,
            spacing_zyx_um=(float(factor), 1.0, 1.0),
        )
        row = {
            "z_spacing_factor": float(factor),
            "z_spacing_um": float(factor),
            "grid_z_slices": int(structure_mask.shape[0]),
            "structure_voxels": int(structure_mask.sum()),
            "hotspot_voxels": int(hotspot_mask.sum()),
            "fraction_inside_structure": float(metrics["fraction_inside_structure"]),
            "median_signed_distance_um": float(metrics["median_signed_distance_um"]),
            "p05_signed_distance_um": float(metrics["p05_signed_distance_um"]),
            "p95_signed_distance_um": float(metrics["p95_signed_distance_um"]),
            "mean_signed_distance_um": float(metrics["mean_signed_distance_um"]),
            "mean_unsigned_distance_um": float(metrics["mean_unsigned_distance_um"]),
        }
        if baseline is None:
            baseline = dict(row)
        row["fraction_inside_delta"] = row["fraction_inside_structure"] - baseline["fraction_inside_structure"]
        row["fraction_inside_percent_delta"] = 100.0 * row["fraction_inside_delta"]
        row["median_signed_delta_um"] = row["median_signed_distance_um"] - baseline["median_signed_distance_um"]
        row["p05_signed_delta_um"] = row["p05_signed_distance_um"] - baseline["p05_signed_distance_um"]
        row["p95_signed_delta_um"] = row["p95_signed_distance_um"] - baseline["p95_signed_distance_um"]
        row["signed_summary_mae_um"] = float(
            np.mean(
                np.abs(
                    [
                        row["median_signed_delta_um"],
                        row["p05_signed_delta_um"],
                        row["p95_signed_delta_um"],
                    ]
                )
            )
        )
        rows.append(row)
    return pd.DataFrame(rows)


def synthetic_ellipsoid_scene(z_spacing_um: float) -> tuple[np.ndarray, np.ndarray]:
    z = np.arange(-30.0, 30.0 + 1e-9, float(z_spacing_um))
    y = np.arange(-42.0, 42.0 + 1e-9, 1.0)
    x = np.arange(-46.0, 46.0 + 1e-9, 1.0)
    z_grid, y_grid, x_grid = np.meshgrid(z, y, x, indexing="ij")

    tilted_x = x_grid - 0.22 * z_grid
    structure_mask = (
        (tilted_x / 31.0) ** 2
        + (y_grid / 25.0) ** 2
        + (z_grid / 19.0) ** 2
        <= 1.0
    )
    hotspot_mask = (
        ((tilted_x - 11.0) / 23.0) ** 2
        + ((y_grid + 3.0) / 16.0) ** 2
        + ((z_grid - 1.0) / 12.0) ** 2
        <= 1.0
    )
    return structure_mask, hotspot_mask


def plot_anisotropy_sweep(sweep: pd.DataFrame, out_png: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), dpi=180)
    x = sweep["z_spacing_factor"].to_numpy(dtype=float)

    axes[0].plot(x, sweep["fraction_inside_structure"], marker="o", color="#2563eb")
    axes[0].set_xscale("log")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([f"{value:g}x" for value in x])
    axes[0].set_xlabel("Z spacing factor")
    axes[0].set_ylabel("Fraction inside structure")
    axes[0].set_title("Mask-overlap stability")
    axes[0].grid(True, alpha=0.25)

    axes[1].plot(x, sweep["median_signed_distance_um"], marker="o", label="median", color="#0f766e")
    axes[1].plot(x, sweep["p05_signed_distance_um"], marker="^", label="5th percentile", color="#7c3aed")
    axes[1].plot(x, sweep["p95_signed_distance_um"], marker="s", label="95th percentile", color="#ea580c")
    axes[1].set_xscale("log")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([f"{value:g}x" for value in x])
    axes[1].set_xlabel("Z spacing factor")
    axes[1].set_ylabel("Signed distance (um)")
    axes[1].set_title("SDF summary drift")
    axes[1].grid(True, alpha=0.25)
    axes[1].legend(frameon=False, fontsize=8)

    fig.suptitle("Synthetic HistoSeg SDF anisotropy sweep", fontsize=13)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)


def summarize_sweep(sweep: pd.DataFrame) -> dict[str, float]:
    return {
        "min_fraction_inside_structure": float(sweep["fraction_inside_structure"].min()),
        "max_fraction_inside_structure": float(sweep["fraction_inside_structure"].max()),
        "max_abs_fraction_inside_delta": float(np.abs(sweep["fraction_inside_delta"]).max()),
        "max_abs_fraction_inside_percent_delta": float(
            np.abs(sweep["fraction_inside_percent_delta"]).max()
        ),
        "max_abs_median_signed_delta_um": float(np.abs(sweep["median_signed_delta_um"]).max()),
        "max_signed_summary_mae_um": float(sweep["signed_summary_mae_um"].max()),
    }


if __name__ == "__main__":
    raise SystemExit(main())
