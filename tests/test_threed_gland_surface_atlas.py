from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from shapely.geometry import Polygon, mapping

from histoseg.threed import GlandSurfaceAtlasConfig, render_gland_surface_atlas
from histoseg.threed.gland_surface_atlas import (
    _build_gland_volume,
    _laplacian_smooth,
    _mesh_from_volume,
    _simplify_mesh,
)


def test_stacked_gland_polygons_produce_shell_and_lumen_mesh(tmp_path: Path):
    gland_dir = _write_surface_fixture(tmp_path, gland_ids=("gland_0001",))
    result = render_gland_surface_atlas(
        GlandSurfaceAtlasConfig(
            gland_instance_dir=gland_dir,
            out_dir=tmp_path / "surface",
            gland_ids=("gland_0001",),
            voxel_size_xy_um=8.0,
            z_interpolation_factor=2,
            max_faces_per_mesh=4000,
        )
    )

    manifest = pd.read_csv(result.gland_surface_mesh_manifest_csv)
    assert result.rendered_gland_count == 1
    assert int(manifest.loc[0, "shell_vertex_count"]) > 0
    assert int(manifest.loc[0, "shell_face_count"]) > 0
    assert bool(manifest.loc[0, "lumen_available"])
    html = (tmp_path / "surface" / "glands" / "gland_0001_surface.html").read_text(encoding="utf-8")
    assert "mesh3d" in html
    assert "epithelial shell" in html


def test_z_interpolation_increases_volume_depth(tmp_path: Path):
    gland_dir = _write_surface_fixture(tmp_path, gland_ids=("gland_0001",))
    tracks = pd.read_csv(gland_dir / "gland_instance_tracks.csv")
    lookup = _feature_lookup(gland_dir)
    cfg_no_interp = GlandSurfaceAtlasConfig(
        gland_instance_dir=gland_dir,
        out_dir=tmp_path / "unused1",
        gland_ids=("gland_0001",),
        z_interpolation_factor=1,
    )
    cfg_interp = GlandSurfaceAtlasConfig(
        gland_instance_dir=gland_dir,
        out_dir=tmp_path / "unused2",
        gland_ids=("gland_0001",),
        z_interpolation_factor=3,
    )

    volume_a = _build_gland_volume(tracks, lookup, cfg_no_interp)
    volume_b = _build_gland_volume(tracks, lookup, cfg_interp)

    assert volume_b.shell_volume.shape[0] > volume_a.shell_volume.shape[0]
    assert volume_b.spacing_zyx_um[0] < volume_a.spacing_zyx_um[0]


def test_smoothing_preserves_vertex_count_and_finite_coordinates():
    vertices = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    )
    faces = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int32)

    smoothed = _laplacian_smooth(vertices, faces, iterations=3, lam=0.25)

    assert smoothed.shape == vertices.shape
    assert np.isfinite(smoothed).all()


def test_mesh_simplification_respects_face_limit():
    vertices = np.random.default_rng(0).normal(size=(120, 3))
    faces = np.array([[i, i + 1, i + 2] for i in range(0, 100)], dtype=np.int32)

    simplified_vertices, simplified_faces = _simplify_mesh(vertices, faces, max_faces=25, max_vertices=80)

    assert len(simplified_faces) <= 25
    assert len(simplified_vertices) <= 80


def test_explicit_gland_ids_limit_output_and_ply_export(tmp_path: Path):
    gland_dir = _write_surface_fixture(tmp_path, gland_ids=("gland_0001", "gland_0002"))
    result = render_gland_surface_atlas(
        GlandSurfaceAtlasConfig(
            gland_instance_dir=gland_dir,
            out_dir=tmp_path / "surface",
            gland_ids=("gland_0002",),
            voxel_size_xy_um=8.0,
            max_faces_per_mesh=5000,
            export_meshes=True,
        )
    )

    manifest = pd.read_csv(result.gland_surface_mesh_manifest_csv)
    assert result.rendered_gland_count == 1
    assert manifest["gland_id"].tolist() == ["gland_0002"]
    assert (tmp_path / "surface" / "meshes" / "gland_0002_shell.ply").exists()
    assert (tmp_path / "surface" / "gland_surface_atlas.html").exists()


def test_surface_atlas_cli_smoke(tmp_path: Path):
    from histoseg.threed.cli import main as threed_cli_main

    gland_dir = _write_surface_fixture(tmp_path, gland_ids=("gland_0001",))
    out_dir = tmp_path / "surface_cli"

    threed_cli_main(
        [
            "render-gland-surface-atlas",
            "--gland-instance-dir",
            str(gland_dir),
            "--out-dir",
            str(out_dir),
            "--gland-ids",
            "gland_0001",
            "--voxel-size-xy-um",
            "8",
            "--max-faces-per-mesh",
            "5000",
        ]
    )

    assert (out_dir / "gland_surface_atlas.html").exists()
    assert (out_dir / "glands" / "gland_0001_surface.html").exists()


def _write_surface_fixture(tmp_path: Path, *, gland_ids: tuple[str, ...]) -> Path:
    gland_dir = tmp_path / "glands"
    gland_dir.mkdir()
    features = []
    rows = []
    index_rows = []
    for gland_number, gland_id in enumerate(gland_ids, start=1):
        offset = gland_number * 240.0
        gland_rows = []
        for slice_order, radius in [(1, 46.0), (2, 52.0), (3, 48.0)]:
            outer = _circle(offset, 100.0, radius)
            lumen = _circle(offset, 100.0, radius * 0.35)
            instance_id = f"{slice_order:03d}_Structure_3_{gland_number:04d}"
            props = {
                "slice_instance_id": instance_id,
                "slice_order": slice_order,
                "sample_id": f"slice_{slice_order}",
                "z_um": float((slice_order - 1) * 5.0),
                "semantic_structure": "Structure 3",
                "area_um2": float(outer.area),
                "centroid_x_um": float(outer.centroid.x),
                "centroid_y_um": float(outer.centroid.y),
                "lumen_area_um2": float(lumen.area),
                "lumen_centroid_x_um": float(lumen.centroid.x),
                "lumen_centroid_y_um": float(lumen.centroid.y),
                "ring_support_score": 1.0,
                "epithelial_marker_score": 2.0,
                "stromal_immune_contamination_score": 0.0,
                "cell_count": 20,
                "qc_flags": "",
                "flag_no_lumen_seed": False,
                "flag_weak_ring": False,
                "flag_merged_candidate": False,
                "flag_small_fragment": False,
                "marker_profile_json": "{}",
                "lumen_polygon_coordinates": json.dumps([[float(x), float(y)] for x, y in lumen.exterior.coords]),
                "gland_id": gland_id,
                "prev_slice_instance_id": "",
                "next_slice_instance_id": "",
                "prev_link_score": "",
                "next_link_score": "",
                "link_score": "",
                "branch_merge_candidates": "",
            }
            rows.append(props)
            gland_rows.append(props)
            features.append(
                {
                    "type": "Feature",
                    "geometry": mapping(outer),
                    "properties": props,
                }
            )
        index_rows.append(
            {
                "gland_id": gland_id,
                "semantic_structure": "Structure 3",
                "slice_count": len(gland_rows),
                "component_count": len(gland_rows),
                "first_slice_order": 1,
                "last_slice_order": 3,
                "z_min_um": 0.0,
                "z_max_um": 10.0,
                "z_span_um": 10.0,
                "area_median_um2": float(np.median([row["area_um2"] for row in gland_rows])),
                "area_cv": 0.0,
                "median_ring_support_score": 1.0,
                "median_epithelial_marker_score": 2.0,
                "max_stromal_immune_contamination_score": 0.0,
                "branch_merge_candidate_count": 0,
                "qc_flags": "",
                "page": f"glands/{gland_id}.html",
                "page_rendered": True,
            }
        )
    (gland_dir / "slice_gland_instances.geojson").write_text(
        json.dumps({"type": "FeatureCollection", "features": features}),
        encoding="utf-8",
    )
    pd.DataFrame(rows).to_csv(gland_dir / "gland_instance_tracks.csv", index=False)
    pd.DataFrame(index_rows).to_csv(gland_dir / "gland_instance_qc_index.csv", index=False)
    return gland_dir


def _circle(cx: float, cy: float, radius: float, *, n: int = 48) -> Polygon:
    angles = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    coords = [(cx + radius * np.cos(a), cy + radius * np.sin(a)) for a in angles]
    coords.append(coords[0])
    return Polygon(coords)


def _feature_lookup(gland_dir: Path) -> dict[str, dict]:
    payload = json.loads((gland_dir / "slice_gland_instances.geojson").read_text(encoding="utf-8"))
    return {
        feature["properties"]["slice_instance_id"]: feature
        for feature in payload["features"]
    }
