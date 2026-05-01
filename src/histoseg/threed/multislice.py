from __future__ import annotations

import copy
import json
import math
import re
import shutil
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence, Union

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter, label as nd_label
from scipy.optimize import minimize
from shapely import affinity
from shapely.geometry import GeometryCollection, MultiPolygon, Polygon, mapping, shape
from shapely.ops import unary_union
from skimage import measure

from histoseg.contour import (
    MultiStructureContourConfig,
    MultiStructureSpec,
    run_multi_structure_contours,
)

from .soft_alignment import (
    ThreeDContourReconstructionConfig,
    _feature_group,
    _iou,
    _read_geojson,
    _records_from_geojson,
    run_3d_contour_reconstruction,
)

try:  # Shapely 2.x fast vectorized containment.
    from shapely import contains_xy as _contains_xy
except Exception:  # pragma: no cover - only used on older Shapely versions.
    _contains_xy = None

PathLike = Union[str, Path]


@dataclass(frozen=True)
class ThreeDStackReconstructionConfig:
    """Configuration for same-sample multi-slice Xenium contour reconstruction.

    This workflow reads Xenium output folders through GitHub pyXenium's
    ``XeniumSlide`` IO path, builds 2D structure contours for each slice, aligns
    slices in a fixed order, and exports a 3D contour stack plus surface meshes.
    """

    xenium_root: PathLike
    out_dir: PathLike = "outputs/3d_stack_reconstruction"
    segmentation_strategy: PathLike | None = None
    structures: Sequence[MultiStructureSpec | Mapping[str, Any]] | None = None
    slice_dirs: Sequence[PathLike] | None = None
    sample_glob: str = "*"
    reference_slice_index: int = 0
    z_spacing_um: float = 5.0
    xenium_pixel_size_um: float = 0.2125
    group_property: str = "structure"
    cluster_column_name: str = "cluster"
    clusters_relpath: str = "analysis/clustering/gene_expression_graphclust/clusters.csv"
    merged_h5ad: PathLike | None = None
    merged_cluster_column: str | None = None
    merged_sample_column: str = "sample_id"
    merged_barcode_column: str = "barcode"
    bins_x: int = 900
    bins_y: int = 700
    gaussian_sigma: float = 2.25
    density_scale_quantile: float = 0.98
    support_quantile: float = 0.18
    tissue_quantile: float = 0.06
    min_dominance: float = 0.34
    min_cells: int = 500
    min_component_pixels: int = 180
    save_slice_preview_png: bool = False
    hard_alignment_maxiter: int = 80
    run_soft_alignment: bool = True
    sampling_distance_um: float = 50.0
    max_landmark_distance_um: float = 180.0
    landmarks_per_structure: int | None = 260
    diagnostic_structure_landmarks: int | None = 620
    rbf_neighbors: int | None = 96
    rbf_smoothing: float = 1e-4
    diagnostic_structure: str | None = "Structure 5"
    save_alignment_preview_png: bool = True
    point_sample_distance_um: float = 80.0
    voxel_size_um: float = 80.0
    mesh_method: str = "marching_cubes"
    mesh_smoothing_sigma_um: float = 40.0
    mesh_level: float = 0.5
    mesh_export_formats: Sequence[str] = ("ply", "obj")
    mesh_cleanup: bool = True
    min_mesh_component_volume_um3: float | None = None
    mesh_max_faces_for_html: int = 25000
    overwrite: bool = False
    dpi: int = 180


@dataclass
class ThreeDStackReconstructionResult:
    out_dir: Path
    slice_manifest_csv: Path
    pairwise_metrics_csv: Path
    aligned_manifest_csv: Path
    contour_points_csv: Path
    summary_json: Path
    visualization_html: Path
    mesh_dir: Path


@dataclass(frozen=True)
class _SliceInput:
    order: int
    sample_id: str
    sample_dir: Path
    xenium_dir: Path


@dataclass(frozen=True)
class _SimilarityTransform:
    origin_x: float
    origin_y: float
    rotation_degrees: float
    scale: float
    translate_x: float
    translate_y: float


def run_3d_stack_reconstruction(
    cfg: ThreeDStackReconstructionConfig,
) -> ThreeDStackReconstructionResult:
    """Run 32-slice style Xenium contour alignment and 3D reconstruction."""

    _validate_stack_config(cfg)
    out_dir = Path(cfg.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    structures = list(cfg.structures) if cfg.structures is not None else _read_strategy_specs(cfg)
    slices = discover_xenium_slices(
        cfg.xenium_root,
        slice_dirs=cfg.slice_dirs,
        sample_glob=cfg.sample_glob,
    )
    if not slices:
        raise ValueError(f"No Xenium output folders were found under {cfg.xenium_root!s}.")
    if cfg.reference_slice_index != 0:
        raise ValueError("Only reference_slice_index=0 is supported in v1 to keep order fixed.")

    slice_manifest_path = out_dir / "xenium_slice_manifest.csv"
    pd.DataFrame(
        [
            {
                "order": item.order,
                "sample_id": item.sample_id,
                "sample_dir": str(item.sample_dir),
                "xenium_dir": str(item.xenium_dir),
            }
            for item in slices
        ]
    ).to_csv(slice_manifest_path, index=False)
    merged_cluster_tables = _load_merged_cluster_tables(cfg, slices)

    contour_dir = out_dir / "slice_contours"
    aligned_dir = out_dir / "aligned_contours"
    alignment_dir = out_dir / "pairwise_alignments"
    mesh_dir = out_dir / "meshes"
    for path in (contour_dir, aligned_dir, alignment_dir, mesh_dir):
        path.mkdir(parents=True, exist_ok=True)

    aligned_rows: list[dict[str, Any]] = []
    pairwise_rows: list[dict[str, Any]] = []
    aligned_paths: list[Path] = []

    for slice_input in slices:
        raw_geojson = _build_slice_contours(
            slice_input,
            structures,
            contour_dir,
            cfg,
            merged_clusters=merged_cluster_tables.get(slice_input.sample_id),
        )
        aligned_path = aligned_dir / f"{slice_input.order:03d}_{slice_input.sample_id}_aligned.geojson"

        if slice_input.order == 1:
            if cfg.overwrite or not aligned_path.exists():
                shutil.copy2(raw_geojson, aligned_path)
            aligned_paths.append(aligned_path)
            aligned_rows.append(
                {
                    "order": slice_input.order,
                    "sample_id": slice_input.sample_id,
                    "z_um": 0.0,
                    "raw_geojson": str(raw_geojson),
                    "aligned_geojson": str(aligned_path),
                    "alignment_role": "reference",
                }
            )
            continue

        fixed_path = aligned_paths[-1]
        pair_dir = alignment_dir / f"{slice_input.order - 1:03d}_to_{slice_input.order:03d}_{slice_input.sample_id}"
        pair_dir.mkdir(parents=True, exist_ok=True)
        hard_path = pair_dir / "moving_hard_aligned.geojson"
        hard_summary_path = pair_dir / "hard_similarity_alignment.json"

        hard_summary = hard_align_geojson(
            fixed_geojson=fixed_path,
            moving_geojson=raw_geojson,
            output_geojson=hard_path,
            summary_json=hard_summary_path,
            group_property=cfg.group_property,
            maxiter=cfg.hard_alignment_maxiter,
            overwrite=cfg.overwrite,
        )

        if cfg.run_soft_alignment:
            soft_result = run_3d_contour_reconstruction(
                ThreeDContourReconstructionConfig(
                    fixed_geojson=fixed_path,
                    moving_hard_aligned_geojson=hard_path,
                    out_dir=pair_dir / "soft_tps",
                    group_property=cfg.group_property,
                    sampling_distance_um=cfg.sampling_distance_um,
                    max_landmark_distance_um=cfg.max_landmark_distance_um,
                    landmarks_per_structure=cfg.landmarks_per_structure,
                    diagnostic_structure_landmarks=cfg.diagnostic_structure_landmarks,
                    rbf_neighbors=cfg.rbf_neighbors,
                    rbf_smoothing=cfg.rbf_smoothing,
                    diagnostic_structure=cfg.diagnostic_structure,
                    dpi=cfg.dpi,
                    save_preview_png=cfg.save_alignment_preview_png,
                )
            )
            soft_summary = json.loads(soft_result.summary_json.read_text(encoding="utf-8"))
            soft_improved = (
                soft_summary["qc"]["union_iou_soft_after"]
                >= soft_summary["qc"]["union_iou_hard_before_soft"]
            )
            if cfg.overwrite or not aligned_path.exists():
                shutil.copy2(
                    soft_result.soft_aligned_geojson if soft_improved else hard_path,
                    aligned_path,
                )
            pairwise_rows.append(
                _pairwise_row(
                    slice_input,
                    hard_summary,
                    soft_summary,
                    soft_result,
                    soft_accepted=soft_improved,
                )
            )
        else:
            if cfg.overwrite or not aligned_path.exists():
                shutil.copy2(hard_path, aligned_path)
            pairwise_rows.append(_pairwise_row(slice_input, hard_summary, None, None))

        aligned_paths.append(aligned_path)
        aligned_rows.append(
            {
                "order": slice_input.order,
                "sample_id": slice_input.sample_id,
                "z_um": (slice_input.order - 1) * cfg.z_spacing_um,
                "raw_geojson": str(raw_geojson),
                "aligned_geojson": str(aligned_path),
                "alignment_role": "moving",
            }
        )

    aligned_manifest_path = out_dir / "aligned_slice_manifest.csv"
    pd.DataFrame(aligned_rows).to_csv(aligned_manifest_path, index=False)

    pairwise_metrics_path = out_dir / "pairwise_alignment_metrics.csv"
    pd.DataFrame(pairwise_rows).to_csv(pairwise_metrics_path, index=False)

    contour_points_path = out_dir / "aligned_contour_3d_points.csv"
    contour_points = write_3d_contour_points(
        aligned_rows,
        contour_points_path,
        group_property=cfg.group_property,
        point_sample_distance_um=cfg.point_sample_distance_um,
        xenium_pixel_size_um=cfg.xenium_pixel_size_um,
    )

    mesh_payloads = reconstruct_3d_contour_meshes(
        aligned_rows,
        mesh_dir,
        group_property=cfg.group_property,
        voxel_size_um=cfg.voxel_size_um,
        z_spacing_um=cfg.z_spacing_um,
        xenium_pixel_size_um=cfg.xenium_pixel_size_um,
        mesh_method=cfg.mesh_method,
        mesh_smoothing_sigma_um=cfg.mesh_smoothing_sigma_um,
        mesh_level=cfg.mesh_level,
        mesh_export_formats=cfg.mesh_export_formats,
        mesh_cleanup=cfg.mesh_cleanup,
        min_mesh_component_volume_um3=cfg.min_mesh_component_volume_um3,
        max_faces_for_html=cfg.mesh_max_faces_for_html,
    )

    visualization_path = out_dir / "histoseg_3d_contour_stack.html"
    write_3d_visualization_html(
        contour_points,
        mesh_payloads,
        visualization_path,
        title="HistoSeg 3D contour reconstruction",
    )

    summary_path = out_dir / "3d_stack_reconstruction_summary.json"
    summary = {
        "config": _jsonable_config(cfg),
        "slice_count": len(slices),
        "structure_count": len(structures),
        "outputs": {
            "slice_manifest_csv": str(slice_manifest_path),
            "aligned_manifest_csv": str(aligned_manifest_path),
            "pairwise_metrics_csv": str(pairwise_metrics_path),
            "contour_points_csv": str(contour_points_path),
            "visualization_html": str(visualization_path),
            "mesh_dir": str(mesh_dir),
            "mesh_qc_summary_json": str(mesh_dir / "mesh_qc_summary.json"),
        },
    }
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    return ThreeDStackReconstructionResult(
        out_dir=out_dir,
        slice_manifest_csv=slice_manifest_path,
        pairwise_metrics_csv=pairwise_metrics_path,
        aligned_manifest_csv=aligned_manifest_path,
        contour_points_csv=contour_points_path,
        summary_json=summary_path,
        visualization_html=visualization_path,
        mesh_dir=mesh_dir,
    )


def discover_xenium_slices(
    xenium_root: PathLike,
    *,
    slice_dirs: Sequence[PathLike] | None = None,
    sample_glob: str = "*",
) -> list[_SliceInput]:
    """Discover Xenium outs directories and return them in numeric slice order."""

    root = Path(xenium_root).expanduser().resolve()
    if slice_dirs is not None:
        candidates = [Path(path).expanduser().resolve() for path in slice_dirs]
    elif _looks_like_xenium_dir(root):
        candidates = [root]
    else:
        candidates = [
            child
            for child in root.glob(sample_glob)
            if child.is_dir() and _find_xenium_output_dir(child) is not None
        ]

    slices: list[_SliceInput] = []
    for candidate in candidates:
        xenium_dir = _find_xenium_output_dir(candidate)
        if xenium_dir is None:
            continue
        slices.append(
            _SliceInput(
                order=0,
                sample_id=candidate.name,
                sample_dir=candidate,
                xenium_dir=xenium_dir,
            )
        )

    ordered = sorted(slices, key=lambda item: _sample_sort_key(item.sample_id))
    return [
        _SliceInput(
            order=index,
            sample_id=item.sample_id,
            sample_dir=item.sample_dir,
            xenium_dir=item.xenium_dir,
        )
        for index, item in enumerate(ordered, start=1)
    ]


def hard_align_geojson(
    *,
    fixed_geojson: PathLike,
    moving_geojson: PathLike,
    output_geojson: PathLike,
    summary_json: PathLike | None = None,
    group_property: str = "structure",
    maxiter: int = 80,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Hard-align one contour GeoJSON to another with a similarity transform."""

    output_path = Path(output_geojson)
    summary_path = Path(summary_json) if summary_json is not None else None
    if output_path.exists() and not overwrite and summary_path is not None and summary_path.exists():
        return json.loads(summary_path.read_text(encoding="utf-8"))

    fixed_payload = _read_geojson(Path(fixed_geojson))
    moving_payload = _read_geojson(Path(moving_geojson))
    fixed_records = _records_from_geojson(fixed_payload, group_property, "fixed")
    moving_records = _records_from_geojson(moving_payload, group_property, "moving")
    fixed_union = unary_union([record.geometry for record in fixed_records])
    moving_union = unary_union([record.geometry for record in moving_records])

    transform, optimization = _estimate_similarity_transform(
        fixed_union,
        moving_union,
        maxiter=maxiter,
    )
    aligned_union = _apply_similarity_to_geometry(moving_union, transform)
    before_iou = _iou(fixed_union, moving_union)
    after_iou = _iou(fixed_union, aligned_union)
    accepted = after_iou >= before_iou
    if not accepted:
        transform = _SimilarityTransform(
            origin_x=float(moving_union.centroid.x),
            origin_y=float(moving_union.centroid.y),
            rotation_degrees=0.0,
            scale=1.0,
            translate_x=0.0,
            translate_y=0.0,
        )
        optimization = {
            **optimization,
            "accepted": False,
            "accepted_reason": "identity_kept_because_similarity_alignment_reduced_iou",
        }
        aligned_union = moving_union
        after_iou = before_iou
    else:
        optimization = {**optimization, "accepted": True}

    aligned_payload = _apply_similarity_to_geojson(moving_payload, transform)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(aligned_payload, ensure_ascii=False), encoding="utf-8")

    per_structure = _per_structure_iou(
        fixed_records,
        _records_from_geojson(aligned_payload, group_property, "aligned"),
    )
    summary = {
        "fixed_geojson": str(Path(fixed_geojson)),
        "moving_geojson": str(Path(moving_geojson)),
        "output_geojson": str(output_path),
        "transform": asdict(transform),
        "optimization": optimization,
        "union_iou_before_hard": before_iou,
        "union_iou_after_hard": after_iou,
        "hard_alignment_accepted": accepted,
        "per_structure_iou_after_hard": per_structure,
    }
    if summary_path is not None:
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text(
            json.dumps(summary, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
    return summary


def write_3d_contour_points(
    aligned_rows: Sequence[Mapping[str, Any]],
    output_csv: PathLike,
    *,
    group_property: str,
    point_sample_distance_um: float,
    xenium_pixel_size_um: float,
) -> pd.DataFrame:
    """Sample aligned slice boundaries into a 3D contour point table."""

    rows: list[dict[str, Any]] = []
    native_distance = max(point_sample_distance_um / xenium_pixel_size_um, 1e-6)
    for row in aligned_rows:
        payload = _read_geojson(Path(row["aligned_geojson"]))
        z_um = float(row["z_um"])
        for feature_index, feature in enumerate(payload["features"]):
            group = _feature_group(feature.get("properties") or {}, group_property)
            geom = _valid_polygonal_geometry(shape(feature["geometry"]))
            if geom is None:
                continue
            for poly_index, polygon in enumerate(_iter_polygons(geom)):
                boundary = polygon.exterior
                if boundary.length <= 0:
                    continue
                distances = np.arange(0.0, boundary.length, native_distance)
                if len(distances) == 0:
                    distances = np.array([0.0])
                polyline_id = (
                    f"{int(row['order']):03d}_{feature_index:05d}_{poly_index:03d}"
                )
                for point_index, distance in enumerate(distances):
                    point = boundary.interpolate(float(distance))
                    rows.append(
                        {
                            "slice_order": int(row["order"]),
                            "sample_id": row["sample_id"],
                            "z_um": z_um,
                            "structure": str(group),
                            "polyline_id": polyline_id,
                            "point_index": point_index,
                            "x_um": float(point.x) * xenium_pixel_size_um,
                            "y_um": float(point.y) * xenium_pixel_size_um,
                        }
                    )

    df = pd.DataFrame(
        rows,
        columns=[
            "slice_order",
            "sample_id",
            "z_um",
            "structure",
            "polyline_id",
            "point_index",
            "x_um",
            "y_um",
        ],
    )
    output_path = Path(output_csv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    return df


def reconstruct_3d_contour_meshes(
    aligned_rows: Sequence[Mapping[str, Any]],
    mesh_dir: PathLike,
    *,
    group_property: str,
    voxel_size_um: float,
    z_spacing_um: float,
    xenium_pixel_size_um: float,
    mesh_method: str = "marching_cubes",
    mesh_smoothing_sigma_um: float = 40.0,
    mesh_level: float = 0.5,
    mesh_export_formats: Sequence[str] = ("ply", "obj"),
    mesh_cleanup: bool = True,
    min_mesh_component_volume_um3: float | None = None,
    max_faces_for_html: int = 25000,
) -> list[dict[str, Any]]:
    """Voxelize aligned structure masks and export continuous surface meshes."""

    if mesh_method != "marching_cubes":
        raise ValueError("mesh_method currently supports only 'marching_cubes'.")
    export_formats = _normalize_mesh_export_formats(mesh_export_formats)
    trimesh = _import_trimesh()
    mesh_path = Path(mesh_dir)
    mesh_path.mkdir(parents=True, exist_ok=True)
    slice_records = _load_aligned_slice_geometries(
        aligned_rows,
        group_property=group_property,
        xenium_pixel_size_um=xenium_pixel_size_um,
    )
    structures = sorted({record["structure"] for record in slice_records})
    if not structures:
        return []

    bounds = _combined_scaled_bounds(slice_records)
    x_min, y_min, x_max, y_max = bounds
    x_values = np.arange(x_min, x_max + voxel_size_um, voxel_size_um, dtype=float)
    y_values = np.arange(y_min, y_max + voxel_size_um, voxel_size_um, dtype=float)
    xx, yy = np.meshgrid(x_values, y_values)
    volume_metadata = {
        "x_min_um": float(x_min),
        "y_min_um": float(y_min),
        "x_max_um": float(x_max),
        "y_max_um": float(y_max),
        "voxel_size_um": float(voxel_size_um),
        "z_spacing_um": float(z_spacing_um),
        "mesh_method": mesh_method,
        "mesh_smoothing_sigma_um": float(mesh_smoothing_sigma_um),
        "mesh_level": float(mesh_level),
        "mesh_cleanup": bool(mesh_cleanup),
        "min_mesh_component_volume_um3": min_mesh_component_volume_um3,
    }

    mesh_payloads: list[dict[str, Any]] = []
    for structure in structures:
        volume = np.zeros((len(aligned_rows), len(y_values), len(x_values)), dtype=bool)
        for record in slice_records:
            if record["structure"] != structure:
                continue
            slice_index = int(record["slice_order"]) - 1
            mask = _geometry_contains_grid(record["geometry_um"], xx, yy)
            volume[slice_index] |= mask
        if not volume.any():
            continue

        component_stats_before = _volume_component_stats(volume, voxel_size_um, z_spacing_um)
        filtered_volume = _filter_volume_components(
            volume,
            min_mesh_component_volume_um3=min_mesh_component_volume_um3,
            voxel_size_um=voxel_size_um,
            z_spacing_um=z_spacing_um,
        )
        if not filtered_volume.any():
            continue

        component_stats_after = _volume_component_stats(
            filtered_volume,
            voxel_size_um,
            z_spacing_um,
        )
        padded = np.pad(filtered_volume.astype(np.float32), ((1, 1), (1, 1), (1, 1)))
        field = _smooth_volume_field(
            padded,
            mesh_smoothing_sigma_um=mesh_smoothing_sigma_um,
            voxel_size_um=voxel_size_um,
            z_spacing_um=z_spacing_um,
        )
        level = float(mesh_level)
        smoothing_applied = bool(mesh_smoothing_sigma_um > 0)
        if field.max() <= level and smoothing_applied:
            field = padded
            level = 0.5
            smoothing_applied = False

        try:
            vertices, faces, _, _ = measure.marching_cubes(
                field,
                level=level,
                spacing=(z_spacing_um, voxel_size_um, voxel_size_um),
            )
        except ValueError:
            continue

        xyz = np.column_stack(
            [
                x_min + (vertices[:, 2] - voxel_size_um),
                y_min + (vertices[:, 1] - voxel_size_um),
                vertices[:, 0] - z_spacing_um,
            ]
        )

        mesh = trimesh.Trimesh(vertices=xyz, faces=faces, process=False)
        if mesh_cleanup:
            mesh = _cleanup_trimesh_mesh(mesh)

        base_path = mesh_path / _safe_name(structure)
        export_paths = _export_mesh(mesh, base_path, export_formats)
        component_count = _mesh_component_count(mesh)
        payload = {
            "structure": structure,
            "ply": str(export_paths.get("ply")) if "ply" in export_paths else None,
            "obj": str(export_paths.get("obj")) if "obj" in export_paths else None,
            "vertex_count": int(len(mesh.vertices)),
            "face_count": int(len(mesh.faces)),
            "surface_area_um2": float(mesh.area),
            "volume_um3": float(abs(mesh.volume)) if len(mesh.faces) else 0.0,
            "is_watertight": bool(mesh.is_watertight),
            "euler_number": int(mesh.euler_number),
            "component_count": int(component_count),
            "volume_component_count_before_filter": int(
                component_stats_before["component_count"]
            ),
            "volume_component_count_after_filter": int(component_stats_after["component_count"]),
            "volume_voxel_count_after_filter": int(component_stats_after["voxel_count"]),
            "smoothing_applied": smoothing_applied,
            "mesh_level_used": float(level),
        }
        if len(mesh.faces) <= max_faces_for_html:
            payload.update(
                {
                    "x": np.asarray(mesh.vertices)[:, 0].round(3).tolist(),
                    "y": np.asarray(mesh.vertices)[:, 1].round(3).tolist(),
                    "z": np.asarray(mesh.vertices)[:, 2].round(3).tolist(),
                    "i": np.asarray(mesh.faces)[:, 0].astype(int).tolist(),
                    "j": np.asarray(mesh.faces)[:, 1].astype(int).tolist(),
                    "k": np.asarray(mesh.faces)[:, 2].astype(int).tolist(),
                }
            )
        mesh_payloads.append(payload)

    manifest_columns = [
        "structure",
        "ply",
        "obj",
        "vertex_count",
        "face_count",
        "surface_area_um2",
        "volume_um3",
        "is_watertight",
        "euler_number",
        "component_count",
        "volume_component_count_before_filter",
        "volume_component_count_after_filter",
        "volume_voxel_count_after_filter",
        "smoothing_applied",
        "mesh_level_used",
    ]
    pd.DataFrame(
        [{column: item.get(column) for column in manifest_columns} for item in mesh_payloads],
        columns=manifest_columns,
    ).to_csv(mesh_path / "mesh_manifest.csv", index=False)
    _write_mesh_qc_summary(mesh_path / "mesh_qc_summary.json", mesh_payloads, volume_metadata)
    return mesh_payloads


def write_3d_visualization_html(
    contour_points: pd.DataFrame,
    mesh_payloads: Sequence[Mapping[str, Any]],
    output_html: PathLike,
    *,
    title: str,
) -> Path:
    """Write a self-contained Plotly HTML view for the 3D contour stack."""

    traces: list[dict[str, Any]] = []
    palette = [
        "#3ddbc7",
        "#5aa9ff",
        "#ffb156",
        "#b891ff",
        "#ff7188",
        "#8be563",
        "#ffd15c",
    ]
    has_mesh_surface = any({"x", "y", "z", "i", "j", "k"}.issubset(mesh.keys()) for mesh in mesh_payloads)
    for color_index, mesh in enumerate(mesh_payloads):
        if not {"x", "y", "z", "i", "j", "k"}.issubset(mesh.keys()):
            continue
        traces.append(
            {
                "type": "mesh3d",
                "name": f"{mesh['structure']} mesh",
                "x": mesh["x"],
                "y": mesh["y"],
                "z": mesh["z"],
                "i": mesh["i"],
                "j": mesh["j"],
                "k": mesh["k"],
                "color": palette[color_index % len(palette)],
                "opacity": 0.36,
                "flatshading": True,
            }
        )
    for color_index, (structure, group) in enumerate(contour_points.groupby("structure")):
        x_values: list[float | None] = []
        y_values: list[float | None] = []
        z_values: list[float | None] = []
        for _, polyline in group.sort_values(
            ["slice_order", "polyline_id", "point_index"]
        ).groupby("polyline_id", sort=False):
            x_values.extend(polyline["x_um"].round(3).tolist())
            y_values.extend(polyline["y_um"].round(3).tolist())
            z_values.extend(polyline["z_um"].round(3).tolist())
            x_values.append(None)
            y_values.append(None)
            z_values.append(None)
        trace = {
            "type": "scatter3d",
            "mode": "lines",
            "name": f"{structure} contours",
            "x": x_values,
            "y": y_values,
            "z": z_values,
            "line": {"width": 3, "color": palette[color_index % len(palette)]},
            "opacity": 0.9,
        }
        if has_mesh_surface:
            trace["visible"] = "legendonly"
        traces.append(trace)

    html = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8" />
  <title>{title}</title>
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <style>
    body {{ margin: 0; background: #07111d; color: #d8e2ef; font-family: Arial, sans-serif; }}
    #plot {{ width: 100vw; height: 100vh; }}
  </style>
</head>
<body>
  <div id="plot"></div>
  <script>
    const data = {json.dumps(traces)};
    const layout = {{
      title: {{text: {json.dumps(title)}, font: {{color: "#d8e2ef", size: 22}}}},
      paper_bgcolor: "#07111d",
      plot_bgcolor: "#07111d",
      scene: {{
        xaxis: {{title: "X (um)", color: "#b9c7d6", gridcolor: "#27364a"}},
        yaxis: {{title: "Y (um)", color: "#b9c7d6", gridcolor: "#27364a"}},
        zaxis: {{title: "Z (um)", color: "#b9c7d6", gridcolor: "#27364a"}},
        aspectmode: "data"
      }},
      legend: {{font: {{color: "#d8e2ef"}}}},
      margin: {{l: 0, r: 0, t: 50, b: 0}}
    }};
    Plotly.newPlot("plot", data, layout, {{responsive: true}});
  </script>
</body>
</html>
"""
    output_path = Path(output_html)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html, encoding="utf-8")
    return output_path


def _validate_stack_config(cfg: ThreeDStackReconstructionConfig) -> None:
    if cfg.z_spacing_um <= 0:
        raise ValueError("z_spacing_um must be greater than 0.")
    if cfg.xenium_pixel_size_um <= 0:
        raise ValueError("xenium_pixel_size_um must be greater than 0.")
    if cfg.voxel_size_um <= 0:
        raise ValueError("voxel_size_um must be greater than 0.")
    if cfg.point_sample_distance_um <= 0:
        raise ValueError("point_sample_distance_um must be greater than 0.")
    if cfg.mesh_method != "marching_cubes":
        raise ValueError("mesh_method currently supports only 'marching_cubes'.")
    if cfg.mesh_smoothing_sigma_um < 0:
        raise ValueError("mesh_smoothing_sigma_um must be non-negative.")
    if cfg.mesh_level <= 0:
        raise ValueError("mesh_level must be greater than 0.")
    if (
        cfg.min_mesh_component_volume_um3 is not None
        and cfg.min_mesh_component_volume_um3 < 0
    ):
        raise ValueError("min_mesh_component_volume_um3 must be non-negative.")
    _normalize_mesh_export_formats(cfg.mesh_export_formats)


def _read_strategy_specs(cfg: ThreeDStackReconstructionConfig) -> list[MultiStructureSpec]:
    strategy_path = Path(cfg.segmentation_strategy) if cfg.segmentation_strategy else None
    if strategy_path is None:
        candidate = Path(cfg.xenium_root) / "contour for alignment" / "segmentationstrategy.txt"
        if candidate.exists():
            strategy_path = candidate
    if strategy_path is None or not strategy_path.exists():
        raise ValueError(
            "Provide structures=... or --segmentation-strategy. "
            "No contour for alignment/segmentationstrategy.txt was found."
        )
    specs: list[MultiStructureSpec] = []
    for index, line in enumerate(strategy_path.read_text(encoding="utf-8").splitlines(), start=1):
        text = line.strip()
        if not text:
            continue
        clusters = [value.strip() for value in text.split(",") if value.strip()]
        specs.append(
            MultiStructureSpec(
                structure_name=f"Structure {len(specs) + 1}",
                structure_id=len(specs) + 1,
                cluster_ids=clusters,
            )
        )
    if not specs:
        raise ValueError(f"No structures were parsed from {strategy_path}.")
    return specs


def _build_slice_contours(
    slice_input: _SliceInput,
    structures: Sequence[MultiStructureSpec | Mapping[str, Any]],
    contour_root: Path,
    cfg: ThreeDStackReconstructionConfig,
    *,
    merged_clusters: pd.DataFrame | None = None,
) -> Path:
    slice_out = contour_root / f"{slice_input.order:03d}_{slice_input.sample_id}"
    contour_geojson = slice_out / "xenium_explorer_annotations.geojson"
    if contour_geojson.exists() and not cfg.overwrite:
        return contour_geojson

    slice_out.mkdir(parents=True, exist_ok=True)
    slide = _read_xenium_slide(slice_input.xenium_dir, cfg)
    inputs = _write_slide_contour_inputs(
        slide,
        slice_out,
        cfg,
        merged_clusters=merged_clusters,
    )
    result = run_multi_structure_contours(
        MultiStructureContourConfig(
            clusters_csv=inputs["clusters_csv"],
            cells_parquet=inputs["cells_parquet"],
            out_dir=slice_out,
            structures=structures,
            bins_x=cfg.bins_x,
            bins_y=cfg.bins_y,
            gaussian_sigma=cfg.gaussian_sigma,
            density_scale_quantile=cfg.density_scale_quantile,
            support_quantile=cfg.support_quantile,
            tissue_quantile=cfg.tissue_quantile,
            min_dominance=cfg.min_dominance,
            min_cells=cfg.min_cells,
            min_component_pixels=cfg.min_component_pixels,
            xenium_pixel_size_um=cfg.xenium_pixel_size_um,
            save_preview_png=cfg.save_slice_preview_png,
        )
    )
    return result.geojson


def _read_xenium_slide(xenium_dir: Path, cfg: ThreeDStackReconstructionConfig) -> Any:
    read_xenium = _import_pyxenium_read_xenium()
    try:
        return read_xenium(
            str(xenium_dir),
            as_="slide",
            prefer="auto",
            include_transcripts=False,
            stream_transcripts=True,
            include_boundaries=False,
            include_images=False,
            clusters_relpath=cfg.clusters_relpath,
            cluster_column_name=cfg.cluster_column_name,
        )
    except TypeError:
        return read_xenium(str(xenium_dir), as_="slide")


def _import_pyxenium_read_xenium():
    try:
        from pyXenium.io import read_xenium

        return read_xenium
    except Exception as first_error:
        try:
            from pyxenium.io import read_xenium

            return read_xenium
        except Exception as second_error:
            raise ImportError(
                "HistoSeg 3D stack reconstruction requires GitHub pyXenium. "
                "Install it with: pip install "
                "'pyXenium @ git+https://github.com/hutaobo/pyXenium.git'"
            ) from second_error or first_error


def _write_slide_contour_inputs(
    slide: Any,
    out_dir: Path,
    cfg: ThreeDStackReconstructionConfig,
    *,
    merged_clusters: pd.DataFrame | None = None,
) -> dict[str, Path]:
    table = getattr(slide, "table", None)
    if table is None or not hasattr(table, "obs"):
        raise ValueError("pyXenium did not return a XeniumSlide with a table.obs annotation table.")
    obs = table.obs.copy()
    obs["Barcode"] = obs.index.astype(str)

    x_values, y_values = _extract_slide_xy(obs, table)

    cells = pd.DataFrame(
        {
            "Barcode": obs["Barcode"].astype(str).to_numpy(),
            "cell_id": obs["Barcode"].astype(str).to_numpy(),
            "x_centroid": np.asarray(x_values, dtype=float),
            "y_centroid": np.asarray(y_values, dtype=float),
        }
    )
    if merged_clusters is not None:
        clusters = merged_clusters.loc[:, ["Barcode", "Cluster"]].copy()
    else:
        cluster_col = _find_cluster_column(obs, cfg.cluster_column_name)
        if cluster_col is None:
            raise ValueError(
                f"XeniumSlide.table.obs does not contain cluster column "
                f"{cfg.cluster_column_name!r}, and no merged_h5ad clusters were provided."
            )
        clusters = pd.DataFrame(
            {
                "Barcode": obs["Barcode"].astype(str).to_numpy(),
                "Cluster": obs[cluster_col].astype(str).to_numpy(),
            }
        )
    cell_keep = (
        np.isfinite(cells["x_centroid"].to_numpy(dtype=float))
        & np.isfinite(cells["y_centroid"].to_numpy(dtype=float))
    )
    cells = cells.loc[cell_keep].reset_index(drop=True)
    clusters = clusters.loc[clusters["Cluster"].astype(str).str.strip() != ""].reset_index(
        drop=True
    )

    cells_path = out_dir / "pyxenium_cells_for_contours.parquet"
    clusters_path = out_dir / "pyxenium_clusters_for_contours.csv"
    cells.to_parquet(cells_path, index=False)
    clusters.to_csv(clusters_path, index=False)
    return {"cells_parquet": cells_path, "clusters_csv": clusters_path}


def _load_merged_cluster_tables(
    cfg: ThreeDStackReconstructionConfig,
    slices: Sequence[_SliceInput],
) -> dict[str, pd.DataFrame]:
    h5ad_path = _find_merged_h5ad(cfg)
    if h5ad_path is None:
        return {}
    try:
        import scanpy as sc
    except Exception as exc:  # pragma: no cover - scanpy is a package dependency.
        raise ImportError("Reading merged_h5ad cluster labels requires scanpy.") from exc

    adata = sc.read_h5ad(h5ad_path, backed="r")
    try:
        obs = adata.obs
        sample_col = cfg.merged_sample_column
        barcode_col = cfg.merged_barcode_column
        if sample_col not in obs.columns:
            raise ValueError(f"merged_h5ad obs is missing sample column {sample_col!r}.")
        if barcode_col not in obs.columns:
            raise ValueError(f"merged_h5ad obs is missing barcode column {barcode_col!r}.")
        cluster_col = _choose_merged_cluster_column(obs, cfg.merged_cluster_column)
        sample_ids = {item.sample_id for item in slices}
        selected = obs.loc[
            obs[sample_col].astype(str).isin(sample_ids),
            [sample_col, barcode_col, cluster_col],
        ].copy()
    finally:
        if getattr(adata, "isbacked", False):
            adata.file.close()

    selected.columns = ["sample_id", "Barcode", "Cluster"]
    selected["sample_id"] = selected["sample_id"].astype(str)
    selected["Barcode"] = selected["Barcode"].astype(str)
    selected["Cluster"] = selected["Cluster"].astype(str)
    tables: dict[str, pd.DataFrame] = {}
    for sample_id, group in selected.groupby("sample_id", sort=False):
        tables[str(sample_id)] = group.loc[:, ["Barcode", "Cluster"]].reset_index(drop=True)
    return tables


def _find_merged_h5ad(cfg: ThreeDStackReconstructionConfig) -> Path | None:
    if cfg.merged_h5ad is not None:
        path = Path(cfg.merged_h5ad).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(path)
        return path
    root = Path(cfg.xenium_root).expanduser().resolve()
    merge_dir = root / "pdc_merge_leiden"
    if not merge_dir.exists():
        return None
    candidates = sorted(merge_dir.glob("*processed_leiden*.h5ad"))
    if candidates:
        return candidates[0]
    candidates = sorted(merge_dir.glob("*.h5ad"))
    return candidates[0] if candidates else None


def _choose_merged_cluster_column(obs: pd.DataFrame, requested: str | None) -> str:
    if requested is not None:
        if requested not in obs.columns:
            raise ValueError(f"merged_h5ad obs is missing cluster column {requested!r}.")
        return requested
    if "leiden_1_0" in obs.columns:
        return "leiden_1_0"
    candidates = [column for column in obs.columns if str(column).startswith("leiden")]
    if candidates:
        return sorted(candidates)[0]
    raise ValueError("Could not infer a merged_h5ad Leiden cluster column.")


def _extract_slide_xy(obs: pd.DataFrame, table: Any) -> tuple[np.ndarray, np.ndarray]:
    x_candidates = [
        "x_centroid",
        "x_centroid_um",
        "cell_x_centroid",
        "x",
        "x_location",
        "x_coord",
    ]
    y_candidates = [
        "y_centroid",
        "y_centroid_um",
        "cell_y_centroid",
        "y",
        "y_location",
        "y_coord",
    ]
    for x_col in x_candidates:
        for y_col in y_candidates:
            if x_col in obs.columns and y_col in obs.columns:
                return obs[x_col].to_numpy(dtype=float), obs[y_col].to_numpy(dtype=float)
    if hasattr(table, "obsm") and "spatial" in table.obsm:
        spatial = np.asarray(table.obsm["spatial"], dtype=float)
        if spatial.ndim == 2 and spatial.shape[1] >= 2:
            return spatial[:, 0], spatial[:, 1]
    raise ValueError("XeniumSlide.table does not expose x/y centroid coordinates.")


def _find_cluster_column(obs: pd.DataFrame, preferred: str) -> str | None:
    if preferred in obs.columns:
        return preferred
    for candidate in ("Cluster", "cluster", "graphclust", "leiden", "leiden_cluster"):
        if candidate in obs.columns:
            return candidate
    return None


def _find_xenium_output_dir(path: Path) -> Path | None:
    if _looks_like_xenium_dir(path):
        return path
    for child in sorted(path.iterdir()) if path.exists() else []:
        if child.is_dir() and _looks_like_xenium_dir(child):
            return child
    return None


def _looks_like_xenium_dir(path: Path) -> bool:
    return (path / "cells.parquet").exists() and (
        (path / "cell_feature_matrix.h5").exists()
        or (path / "cell_feature_matrix").exists()
    )


def _sample_sort_key(sample_id: str) -> tuple[int, str]:
    matches = re.findall(r"\d+", sample_id)
    if matches:
        return int(matches[-1]), sample_id
    return 10**9, sample_id


def _estimate_similarity_transform(
    fixed_union: Any,
    moving_union: Any,
    *,
    maxiter: int,
) -> tuple[_SimilarityTransform, dict[str, Any]]:
    origin = moving_union.centroid
    scale0 = math.sqrt(max(fixed_union.area, 1e-9) / max(moving_union.area, 1e-9))
    rotation0 = _principal_axis_angle(fixed_union) - _principal_axis_angle(moving_union)
    initial = _SimilarityTransform(
        origin_x=float(origin.x),
        origin_y=float(origin.y),
        rotation_degrees=float(rotation0),
        scale=float(scale0),
        translate_x=0.0,
        translate_y=0.0,
    )
    initial_aligned = _apply_similarity_to_geometry(moving_union, initial)
    centroid_dx = float(fixed_union.centroid.x - initial_aligned.centroid.x)
    centroid_dy = float(fixed_union.centroid.y - initial_aligned.centroid.y)
    initial = _SimilarityTransform(
        origin_x=initial.origin_x,
        origin_y=initial.origin_y,
        rotation_degrees=initial.rotation_degrees,
        scale=initial.scale,
        translate_x=centroid_dx,
        translate_y=centroid_dy,
    )
    before = _iou(fixed_union, moving_union)
    initial_iou = _iou(fixed_union, _apply_similarity_to_geometry(moving_union, initial))

    if maxiter <= 0:
        return initial, {
            "method": "pca_centroid",
            "success": True,
            "union_iou_before": before,
            "union_iou_initial": initial_iou,
            "union_iou_final": initial_iou,
        }

    def objective(params: np.ndarray) -> float:
        transform = _SimilarityTransform(
            origin_x=initial.origin_x,
            origin_y=initial.origin_y,
            rotation_degrees=float(params[0]),
            scale=float(math.exp(params[1])),
            translate_x=float(params[2]),
            translate_y=float(params[3]),
        )
        return -_iou(fixed_union, _apply_similarity_to_geometry(moving_union, transform))

    params0 = np.array(
        [
            initial.rotation_degrees,
            math.log(max(initial.scale, 1e-9)),
            initial.translate_x,
            initial.translate_y,
        ],
        dtype=float,
    )
    result = minimize(
        objective,
        params0,
        method="Nelder-Mead",
        options={"maxiter": int(maxiter), "xatol": 1e-3, "fatol": 1e-5},
    )
    params = result.x if result.fun <= objective(params0) else params0
    transform = _SimilarityTransform(
        origin_x=initial.origin_x,
        origin_y=initial.origin_y,
        rotation_degrees=float(params[0]),
        scale=float(math.exp(params[1])),
        translate_x=float(params[2]),
        translate_y=float(params[3]),
    )
    final_iou = _iou(fixed_union, _apply_similarity_to_geometry(moving_union, transform))
    return transform, {
        "method": "pca_centroid_nelder_mead",
        "success": bool(result.success) or final_iou >= initial_iou,
        "message": str(result.message),
        "iterations": int(getattr(result, "nit", 0)),
        "union_iou_before": before,
        "union_iou_initial": initial_iou,
        "union_iou_final": final_iou,
    }


def _principal_axis_angle(geom: Any) -> float:
    hull = geom.convex_hull
    if hull.is_empty or not hasattr(hull, "exterior"):
        return 0.0
    coords = np.asarray(hull.exterior.coords, dtype=float)
    if len(coords) < 3:
        return 0.0
    centered = coords - coords.mean(axis=0)
    cov = np.cov(centered.T)
    eigvals, eigvecs = np.linalg.eigh(cov)
    axis = eigvecs[:, int(np.argmax(eigvals))]
    return float(math.degrees(math.atan2(axis[1], axis[0])))


def _apply_similarity_to_geojson(
    geojson: dict[str, Any],
    transform: _SimilarityTransform,
) -> dict[str, Any]:
    payload = copy.deepcopy(geojson)
    for feature in payload["features"]:
        geom = _apply_similarity_to_geometry(shape(feature["geometry"]), transform)
        feature["geometry"] = mapping(geom)
    return payload


def _apply_similarity_to_geometry(geom: Any, transform: _SimilarityTransform) -> Any:
    origin = (transform.origin_x, transform.origin_y)
    transformed = affinity.rotate(geom, transform.rotation_degrees, origin=origin)
    transformed = affinity.scale(
        transformed,
        xfact=transform.scale,
        yfact=transform.scale,
        origin=origin,
    )
    transformed = affinity.translate(
        transformed,
        xoff=transform.translate_x,
        yoff=transform.translate_y,
    )
    return transformed


def _per_structure_iou(fixed_records: Sequence[Any], moving_records: Sequence[Any]) -> dict[str, float]:
    fixed_by_group = _group_union(fixed_records)
    moving_by_group = _group_union(moving_records)
    return {
        group: _iou(fixed_by_group[group], moving_by_group[group])
        for group in sorted(set(fixed_by_group) & set(moving_by_group))
    }


def _group_union(records: Sequence[Any]) -> dict[str, Any]:
    grouped: dict[str, list[Any]] = {}
    for record in records:
        grouped.setdefault(record.group, []).append(record.geometry)
    return {group: unary_union(geoms) for group, geoms in grouped.items()}


def _pairwise_row(
    slice_input: _SliceInput,
    hard_summary: Mapping[str, Any],
    soft_summary: Mapping[str, Any] | None,
    soft_result: Any,
    soft_accepted: bool | None = None,
) -> dict[str, Any]:
    row = {
        "moving_order": slice_input.order,
        "moving_sample_id": slice_input.sample_id,
        "hard_union_iou_before": hard_summary.get("union_iou_before_hard"),
        "hard_union_iou_after": hard_summary.get("union_iou_after_hard"),
        "hard_transform_rotation_degrees": hard_summary["transform"]["rotation_degrees"],
        "hard_transform_scale": hard_summary["transform"]["scale"],
        "hard_transform_translate_x": hard_summary["transform"]["translate_x"],
        "hard_transform_translate_y": hard_summary["transform"]["translate_y"],
        "hard_accepted": hard_summary.get("hard_alignment_accepted"),
        "soft_union_iou_before": None,
        "soft_union_iou_after": None,
        "soft_accepted": soft_accepted,
        "soft_boundary_landmarks": None,
        "soft_summary_json": None,
    }
    if soft_summary is not None and soft_result is not None:
        row.update(
            {
                "soft_union_iou_before": soft_summary["qc"]["union_iou_hard_before_soft"],
                "soft_union_iou_after": soft_summary["qc"]["union_iou_soft_after"],
                "soft_boundary_landmarks": soft_summary["landmarks"]["boundary_landmark_count"],
                "soft_summary_json": str(soft_result.summary_json),
            }
        )
    return row


def _import_trimesh():
    try:
        import trimesh

        return trimesh
    except Exception as exc:
        raise ImportError(
            "3D surface mesh export requires trimesh. Install with "
            "`pip install -e '.[threed]'` or `pip install trimesh`."
        ) from exc


def _normalize_mesh_export_formats(mesh_export_formats: Sequence[str] | str) -> tuple[str, ...]:
    if isinstance(mesh_export_formats, str):
        raw_formats = [
            value.strip().lower()
            for value in mesh_export_formats.split(",")
            if value.strip()
        ]
    else:
        raw_formats = [str(value).strip().lower() for value in mesh_export_formats]
    formats = tuple(dict.fromkeys(raw_formats))
    if not formats:
        raise ValueError("mesh_export_formats must include at least one format.")
    unsupported = sorted(set(formats) - {"ply", "obj"})
    if unsupported:
        raise ValueError(f"Unsupported mesh export format(s): {', '.join(unsupported)}.")
    return formats


def _smooth_volume_field(
    volume: np.ndarray,
    *,
    mesh_smoothing_sigma_um: float,
    voxel_size_um: float,
    z_spacing_um: float,
) -> np.ndarray:
    field = np.asarray(volume, dtype=np.float32)
    if mesh_smoothing_sigma_um <= 0:
        return field
    sigma = (
        float(mesh_smoothing_sigma_um) / float(z_spacing_um),
        float(mesh_smoothing_sigma_um) / float(voxel_size_um),
        float(mesh_smoothing_sigma_um) / float(voxel_size_um),
    )
    return gaussian_filter(field, sigma=sigma, mode="constant", cval=0.0)


def _volume_component_stats(
    volume: np.ndarray,
    voxel_size_um: float,
    z_spacing_um: float,
) -> dict[str, float | int]:
    labels, component_count = nd_label(volume, structure=np.ones((3, 3, 3), dtype=bool))
    counts = np.bincount(labels.ravel())
    voxel_count = int(volume.sum())
    voxel_volume = float(voxel_size_um) * float(voxel_size_um) * float(z_spacing_um)
    if component_count == 0:
        largest_component_voxels = 0
    else:
        largest_component_voxels = int(counts[1:].max())
    return {
        "component_count": int(component_count),
        "voxel_count": voxel_count,
        "volume_um3": float(voxel_count * voxel_volume),
        "largest_component_voxels": largest_component_voxels,
    }


def _filter_volume_components(
    volume: np.ndarray,
    *,
    min_mesh_component_volume_um3: float | None,
    voxel_size_um: float,
    z_spacing_um: float,
) -> np.ndarray:
    if min_mesh_component_volume_um3 is None or min_mesh_component_volume_um3 <= 0:
        return np.asarray(volume, dtype=bool)

    labels, component_count = nd_label(volume, structure=np.ones((3, 3, 3), dtype=bool))
    if component_count == 0:
        return np.asarray(volume, dtype=bool)
    voxel_volume = float(voxel_size_um) * float(voxel_size_um) * float(z_spacing_um)
    min_voxels = max(int(math.ceil(float(min_mesh_component_volume_um3) / voxel_volume)), 1)
    counts = np.bincount(labels.ravel())
    keep_labels = np.flatnonzero(counts >= min_voxels)
    keep_labels = keep_labels[keep_labels != 0]
    if len(keep_labels) == 0:
        return np.zeros_like(volume, dtype=bool)
    return np.isin(labels, keep_labels)


def _cleanup_trimesh_mesh(mesh):
    mesh = mesh.copy()
    for method_name in (
        "remove_infinite_values",
        "remove_degenerate_faces",
        "remove_duplicate_faces",
    ):
        method = getattr(mesh, method_name, None)
        if callable(method):
            try:
                method()
            except Exception:
                pass
    try:
        mesh.remove_unreferenced_vertices()
    except Exception:
        pass
    try:
        mesh.merge_vertices()
    except Exception:
        pass
    try:
        mesh.fill_holes()
    except Exception:
        pass
    try:
        processed = mesh.process(validate=True)
        if processed is not None:
            mesh = processed
    except Exception:
        pass
    return mesh


def _mesh_component_count(mesh) -> int:
    try:
        return len(mesh.split(only_watertight=False))
    except Exception:
        return 1 if len(mesh.faces) else 0


def _export_mesh(mesh, base_path: Path, export_formats: Sequence[str]) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    for export_format in export_formats:
        path = base_path.with_suffix(f".{export_format}")
        mesh.export(path)
        paths[export_format] = path
    return paths


def _write_mesh_qc_summary(
    path: Path,
    mesh_payloads: Sequence[Mapping[str, Any]],
    volume_metadata: Mapping[str, Any],
) -> None:
    manifest = [
        {
            "structure": item.get("structure"),
            "vertex_count": item.get("vertex_count"),
            "face_count": item.get("face_count"),
            "surface_area_um2": item.get("surface_area_um2"),
            "volume_um3": item.get("volume_um3"),
            "is_watertight": item.get("is_watertight"),
            "euler_number": item.get("euler_number"),
            "component_count": item.get("component_count"),
            "ply": item.get("ply"),
            "obj": item.get("obj"),
        }
        for item in mesh_payloads
    ]
    summary = {
        "mesh_count": len(mesh_payloads),
        "watertight_count": int(sum(1 for item in mesh_payloads if item.get("is_watertight"))),
        "total_vertex_count": int(sum(int(item.get("vertex_count", 0)) for item in mesh_payloads)),
        "total_face_count": int(sum(int(item.get("face_count", 0)) for item in mesh_payloads)),
        "total_surface_area_um2": float(
            sum(float(item.get("surface_area_um2", 0.0)) for item in mesh_payloads)
        ),
        "total_volume_um3": float(
            sum(float(item.get("volume_um3", 0.0)) for item in mesh_payloads)
        ),
        "volume_metadata": dict(volume_metadata),
        "meshes": manifest,
    }
    path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")


def _valid_polygonal_geometry(geom: Any) -> Any | None:
    if geom.is_empty:
        return None
    if not geom.is_valid:
        geom = geom.buffer(0)
    if geom.is_empty:
        return None
    if isinstance(geom, (Polygon, MultiPolygon)):
        return geom
    if isinstance(geom, GeometryCollection):
        polygons = [part for part in geom.geoms if isinstance(part, (Polygon, MultiPolygon))]
        if not polygons:
            return None
        return unary_union(polygons)
    return None


def _iter_polygons(geom: Any) -> Iterable[Polygon]:
    if isinstance(geom, Polygon):
        yield geom
    elif isinstance(geom, MultiPolygon):
        yield from geom.geoms


def _load_aligned_slice_geometries(
    aligned_rows: Sequence[Mapping[str, Any]],
    *,
    group_property: str,
    xenium_pixel_size_um: float,
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for row in aligned_rows:
        payload = _read_geojson(Path(row["aligned_geojson"]))
        for feature in payload["features"]:
            group = _feature_group(feature.get("properties") or {}, group_property)
            geom = _valid_polygonal_geometry(shape(feature["geometry"]))
            if geom is None:
                continue
            scaled = affinity.scale(
                geom,
                xfact=xenium_pixel_size_um,
                yfact=xenium_pixel_size_um,
                origin=(0.0, 0.0),
            )
            records.append(
                {
                    "slice_order": int(row["order"]),
                    "structure": str(group),
                    "geometry_um": scaled,
                }
            )
    return records


def _combined_scaled_bounds(records: Sequence[Mapping[str, Any]]) -> tuple[float, float, float, float]:
    union = unary_union([record["geometry_um"] for record in records])
    return tuple(map(float, union.bounds))


def _geometry_contains_grid(geom: Any, xx: np.ndarray, yy: np.ndarray) -> np.ndarray:
    if _contains_xy is not None:
        return np.asarray(_contains_xy(geom, xx, yy), dtype=bool)
    flat = [geom.contains(shape({"type": "Point", "coordinates": [x, y]})) for x, y in zip(xx.ravel(), yy.ravel())]
    return np.asarray(flat, dtype=bool).reshape(xx.shape)


def _write_ascii_ply(path: Path, vertices: np.ndarray, faces: np.ndarray) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("ply\n")
        handle.write("format ascii 1.0\n")
        handle.write(f"element vertex {len(vertices)}\n")
        handle.write("property float x\n")
        handle.write("property float y\n")
        handle.write("property float z\n")
        handle.write(f"element face {len(faces)}\n")
        handle.write("property list uchar int vertex_indices\n")
        handle.write("end_header\n")
        for vertex in vertices:
            handle.write(f"{vertex[0]:.6f} {vertex[1]:.6f} {vertex[2]:.6f}\n")
        for face in faces:
            handle.write(f"3 {int(face[0])} {int(face[1])} {int(face[2])}\n")


def _safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_") or "structure"


def _jsonable_config(cfg: ThreeDStackReconstructionConfig) -> dict[str, Any]:
    payload = asdict(cfg)
    for key in ("xenium_root", "out_dir", "segmentation_strategy", "merged_h5ad"):
        if payload.get(key) is not None:
            payload[key] = str(payload[key])
    if payload.get("slice_dirs") is not None:
        payload["slice_dirs"] = [str(path) for path in payload["slice_dirs"]]
    return payload
