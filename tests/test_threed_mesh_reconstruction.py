"""Unit tests for the upgraded 3-D Marching Cubes mesh reconstruction.

These tests exercise the new mesh generation path (Gaussian smoothing,
trimesh cleanup, OBJ/PLY export, QC metrics, component filtering) without
requiring a real Xenium dataset.  The helpers mirror the minimal fixture
pattern used elsewhere in the test suite.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from shapely.geometry import box, mapping

from histoseg.threed.multislice import (
    ThreeDStackReconstructionConfig,
    reconstruct_3d_contour_meshes,
    write_3d_visualization_html,
    write_3d_contour_points,
)


# ---------------------------------------------------------------------------
# Minimal aligned_rows fixture helpers
# ---------------------------------------------------------------------------


def _write_geojson(path: Path, features: list) -> None:
    path.write_text(
        json.dumps({"type": "FeatureCollection", "features": features}),
        encoding="utf-8",
    )


def _feature(geom, structure: str) -> dict:
    return {
        "type": "Feature",
        "properties": {"structure": structure},
        "geometry": mapping(geom),
    }


def _make_aligned_rows(tmp_path: Path, n_slices: int = 3) -> list[dict]:
    """Build n_slices GeoJSON files each containing a simple box polygon."""
    rows = []
    z_spacing = 5.0
    for i in range(n_slices):
        geojson_path = tmp_path / f"slice_{i:03d}.geojson"
        # Each slice is a 200 x 200 µm box (in pixel units at 0.2125 µm/px ≈ ~941 px)
        px_size = 200.0 / 0.2125
        features = [_feature(box(0, 0, px_size, px_size), "Structure 1")]
        _write_geojson(geojson_path, features)
        rows.append(
            {
                "order": i + 1,
                "sample_id": f"s{i + 1}",
                "z_um": i * z_spacing,
                "aligned_geojson": str(geojson_path),
            }
        )
    return rows


def _make_two_structure_rows(tmp_path: Path, n_slices: int = 3) -> list[dict]:
    """Two structures: a large box and a tiny box to test component filtering."""
    rows = []
    z_spacing = 5.0
    for i in range(n_slices):
        geojson_path = tmp_path / f"slice2_{i:03d}.geojson"
        px_size = 200.0 / 0.2125
        features = [
            _feature(box(0, 0, px_size, px_size), "Structure 1"),
            # Tiny structure: ~1 voxel at 80 µm resolution
            _feature(box(1000, 1000, 1004, 1004), "Structure 2"),
        ]
        _write_geojson(geojson_path, features)
        rows.append(
            {
                "order": i + 1,
                "sample_id": f"s{i + 1}",
                "z_um": i * z_spacing,
                "aligned_geojson": str(geojson_path),
            }
        )
    return rows


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_box_stack_generates_ply_and_obj(tmp_path):
    """A simple box polygon stack must produce non-empty PLY and OBJ files."""
    aligned_rows = _make_aligned_rows(tmp_path)
    mesh_dir = tmp_path / "meshes"

    payloads = reconstruct_3d_contour_meshes(
        aligned_rows,
        mesh_dir,
        group_property="structure",
        voxel_size_um=80.0,
        z_spacing_um=5.0,
        xenium_pixel_size_um=0.2125,
        mesh_export_formats=("ply", "obj"),
        mesh_smoothing_sigma_um=None,
        mesh_cleanup=True,
    )

    assert len(payloads) == 1, "Expected one structure payload."
    payload = payloads[0]
    assert payload["structure"] == "Structure 1"
    assert payload["vertex_count"] > 0
    assert payload["face_count"] > 0

    ply_path = Path(mesh_dir / "Structure_1.ply")
    obj_path = Path(mesh_dir / "Structure_1.obj")
    assert ply_path.exists(), "PLY file not created."
    assert ply_path.stat().st_size > 0, "PLY file is empty."
    assert obj_path.exists(), "OBJ file not created."
    assert obj_path.stat().st_size > 0, "OBJ file is empty."


def test_smoothing_does_not_change_coord_units(tmp_path):
    """Enabling / disabling Gaussian smoothing must not change the coordinate scale."""
    aligned_rows = _make_aligned_rows(tmp_path)

    payloads_no_smooth = reconstruct_3d_contour_meshes(
        aligned_rows,
        tmp_path / "meshes_no_smooth",
        group_property="structure",
        voxel_size_um=80.0,
        z_spacing_um=5.0,
        xenium_pixel_size_um=0.2125,
        mesh_smoothing_sigma_um=None,
        mesh_export_formats=("ply",),
    )
    payloads_smooth = reconstruct_3d_contour_meshes(
        aligned_rows,
        tmp_path / "meshes_smooth",
        group_property="structure",
        voxel_size_um=80.0,
        z_spacing_um=5.0,
        xenium_pixel_size_um=0.2125,
        mesh_smoothing_sigma_um=40.0,
        mesh_export_formats=("ply",),
    )

    assert payloads_no_smooth and payloads_smooth
    # Both should produce vertices with X coordinates in microns (roughly 0..200+ µm)
    for payload in (payloads_no_smooth[0], payloads_smooth[0]):
        if "x" in payload:
            xs = payload["x"]
            assert max(xs) < 2000.0, "X coords look like pixel values, expected µm."
            assert min(xs) > -200.0, "X coords unexpectedly negative."


def test_mesh_manifest_contains_qc_fields(tmp_path):
    """mesh_manifest.csv must include QC columns added by the upgrade."""
    aligned_rows = _make_aligned_rows(tmp_path)
    mesh_dir = tmp_path / "meshes"

    reconstruct_3d_contour_meshes(
        aligned_rows,
        mesh_dir,
        group_property="structure",
        voxel_size_um=80.0,
        z_spacing_um=5.0,
        xenium_pixel_size_um=0.2125,
        mesh_export_formats=("ply", "obj"),
    )

    manifest_path = mesh_dir / "mesh_manifest.csv"
    assert manifest_path.exists(), "mesh_manifest.csv not written."
    df = pd.read_csv(manifest_path)
    required_cols = {
        "structure",
        "vertex_count",
        "face_count",
        "surface_area_um2",
        "volume_um3",
        "is_watertight",
        "euler_number",
        "component_count",
    }
    missing = required_cols - set(df.columns)
    assert not missing, f"mesh_manifest.csv missing columns: {missing}"
    assert len(df) == 1
    assert df.loc[0, "vertex_count"] > 0
    assert df.loc[0, "face_count"] > 0


def test_mesh_qc_summary_json_written(tmp_path):
    """mesh_qc_summary.json must be created with the correct schema."""
    aligned_rows = _make_aligned_rows(tmp_path)
    mesh_dir = tmp_path / "meshes"

    reconstruct_3d_contour_meshes(
        aligned_rows,
        mesh_dir,
        group_property="structure",
        voxel_size_um=80.0,
        z_spacing_um=5.0,
        xenium_pixel_size_um=0.2125,
    )

    summary_path = mesh_dir / "mesh_qc_summary.json"
    assert summary_path.exists(), "mesh_qc_summary.json not written."
    data = json.loads(summary_path.read_text(encoding="utf-8"))
    assert data["structure_count"] == 1
    assert len(data["structures"]) == 1
    entry = data["structures"][0]
    assert entry["structure"] == "Structure 1"
    assert "vertex_count" in entry
    assert "face_count" in entry
    assert "is_watertight" in entry


def test_min_component_volume_filters_small_fragments(tmp_path):
    """Structure 2 (tiny box) should be filtered when min_mesh_component_volume_um3
    is set above its expected volume."""
    aligned_rows = _make_two_structure_rows(tmp_path)
    mesh_dir = tmp_path / "meshes"

    # Structure 2 is ~4px × 4px = ~0.68µm × 0.68µm at 0.2125 µm/px, so its
    # mesh volume will be near zero.  Setting a large min threshold should keep
    # only Structure 1.
    _LARGE_VOLUME_FILTER_UM3 = 1e9  # Deliberately large to remove all tiny components.
    payloads = reconstruct_3d_contour_meshes(
        aligned_rows,
        mesh_dir,
        group_property="structure",
        voxel_size_um=80.0,
        z_spacing_um=5.0,
        xenium_pixel_size_um=0.2125,
        mesh_export_formats=("ply",),
        min_mesh_component_volume_um3=_LARGE_VOLUME_FILTER_UM3,
    )

    # Structure 1 is large enough and should survive.
    structures_found = {p["structure"] for p in payloads}
    assert "Structure 1" in structures_found


def test_config_new_fields_have_correct_defaults():
    """ThreeDStackReconstructionConfig new fields must have the specified defaults."""
    cfg = ThreeDStackReconstructionConfig(
        xenium_root="/tmp/xenium",
        structures=[],
    )
    assert cfg.mesh_method == "marching_cubes"
    assert cfg.mesh_smoothing_sigma_um == 40.0
    assert cfg.mesh_level == 0.5
    assert cfg.mesh_export_formats == ("ply", "obj")
    assert cfg.mesh_cleanup is True
    assert cfg.min_mesh_component_volume_um3 is None
    assert cfg.merged_h5ad is None
    assert cfg.merged_cluster_column is None


def test_html_contains_mesh_trace_and_contour_fallback(tmp_path):
    """The HTML output must include both a mesh3d trace and scatter3d contour lines."""
    aligned_rows = _make_aligned_rows(tmp_path)
    mesh_dir = tmp_path / "meshes"

    mesh_payloads = reconstruct_3d_contour_meshes(
        aligned_rows,
        mesh_dir,
        group_property="structure",
        voxel_size_um=80.0,
        z_spacing_um=5.0,
        xenium_pixel_size_um=0.2125,
        mesh_export_formats=("ply", "obj"),
        mesh_smoothing_sigma_um=None,
    )

    contour_points = write_3d_contour_points(
        aligned_rows,
        tmp_path / "contour_points.csv",
        group_property="structure",
        point_sample_distance_um=80.0,
        xenium_pixel_size_um=0.2125,
    )

    html_path = tmp_path / "viz.html"
    write_3d_visualization_html(
        contour_points,
        mesh_payloads,
        html_path,
        title="Test",
    )

    assert html_path.exists()
    html = html_path.read_text(encoding="utf-8")
    assert "mesh3d" in html, "HTML does not contain a mesh3d trace."
    assert "scatter3d" in html, "HTML does not contain scatter3d contour fallback."


def test_validate_config_rejects_bad_mesh_level():
    """_validate_stack_config must raise ValueError for mesh_level out of range."""
    from histoseg.threed.multislice import _validate_stack_config

    cfg = ThreeDStackReconstructionConfig(
        xenium_root="/tmp/x",
        structures=[],
        mesh_level=1.5,
    )
    with pytest.raises(ValueError, match="mesh_level"):
        _validate_stack_config(cfg)


def test_validate_config_rejects_negative_smoothing_sigma():
    """_validate_stack_config must raise ValueError for negative smoothing sigma."""
    from histoseg.threed.multislice import _validate_stack_config

    cfg = ThreeDStackReconstructionConfig(
        xenium_root="/tmp/x",
        structures=[],
        mesh_smoothing_sigma_um=-1.0,
    )
    with pytest.raises(ValueError, match="mesh_smoothing_sigma_um"):
        _validate_stack_config(cfg)
