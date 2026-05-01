from __future__ import annotations

import json

import pandas as pd
from shapely.geometry import box, mapping, shape

from histoseg.threed import (
    discover_xenium_slices,
    hard_align_geojson,
    reconstruct_3d_contour_meshes,
    write_3d_contour_points,
    write_3d_visualization_html,
)
from histoseg.threed.multislice import _read_strategy_specs
from histoseg.threed import ThreeDStackReconstructionConfig


def test_discover_xenium_slices_uses_numeric_order(tmp_path):
    for name in ["A079-C-008_10", "A079-C-008_1", "A079-C-008_2"]:
        xenium = tmp_path / name / "output-XETG000"
        xenium.mkdir(parents=True)
        (xenium / "cells.parquet").write_text("", encoding="utf-8")
        (xenium / "cell_feature_matrix.h5").write_text("", encoding="utf-8")

    slices = discover_xenium_slices(tmp_path, sample_glob="A079-C-008_*")

    assert [item.sample_id for item in slices] == [
        "A079-C-008_1",
        "A079-C-008_2",
        "A079-C-008_10",
    ]
    assert [item.order for item in slices] == [1, 2, 3]


def test_segmentation_strategy_parses_one_structure_per_line(tmp_path):
    strategy = tmp_path / "segmentationstrategy.txt"
    strategy.write_text("18\n31,3,30\n", encoding="utf-8")

    specs = _read_strategy_specs(
        ThreeDStackReconstructionConfig(
            xenium_root=tmp_path,
            segmentation_strategy=strategy,
        )
    )

    assert [spec.structure_name for spec in specs] == ["Structure 1", "Structure 2"]
    assert list(specs[1].cluster_ids) == ["31", "3", "30"]


def test_hard_alignment_similarity_improves_union_iou(tmp_path):
    fixed = tmp_path / "fixed.geojson"
    moving = tmp_path / "moving.geojson"
    out = tmp_path / "hard.geojson"
    _write_geojson(
        fixed,
        [
            _feature(box(0, 0, 20, 10), "Structure 1"),
            _feature(box(30, 0, 40, 12), "Structure 2"),
        ],
    )
    _write_geojson(
        moving,
        [
            _feature(box(5, 2, 25, 12), "Structure 1"),
            _feature(box(35, 2, 45, 14), "Structure 2"),
        ],
    )

    summary = hard_align_geojson(
        fixed_geojson=fixed,
        moving_geojson=moving,
        output_geojson=out,
        group_property="structure",
        maxiter=20,
        overwrite=True,
    )

    assert out.exists()
    assert summary["union_iou_after_hard"] > summary["union_iou_before_hard"]
    aligned = json.loads(out.read_text(encoding="utf-8"))
    assert aligned["features"][0]["properties"]["structure"] == "Structure 1"


def test_3d_contour_points_mesh_and_html_outputs(tmp_path):
    slice1 = tmp_path / "slice1.geojson"
    slice2 = tmp_path / "slice2.geojson"
    _write_geojson(slice1, [_feature(box(0, 0, 20, 20), "Structure 1")])
    _write_geojson(slice2, [_feature(box(2, 0, 22, 20), "Structure 1")])
    aligned_rows = [
        {
            "order": 1,
            "sample_id": "s1",
            "z_um": 0.0,
            "aligned_geojson": str(slice1),
        },
        {
            "order": 2,
            "sample_id": "s2",
            "z_um": 5.0,
            "aligned_geojson": str(slice2),
        },
    ]

    points_csv = tmp_path / "points.csv"
    points = write_3d_contour_points(
        aligned_rows,
        points_csv,
        group_property="structure",
        point_sample_distance_um=5.0,
        xenium_pixel_size_um=1.0,
    )
    meshes = reconstruct_3d_contour_meshes(
        aligned_rows,
        tmp_path / "meshes",
        group_property="structure",
        voxel_size_um=5.0,
        z_spacing_um=5.0,
        xenium_pixel_size_um=1.0,
        mesh_smoothing_sigma_um=None,
    )
    html = write_3d_visualization_html(points, meshes, tmp_path / "view.html", title="test")

    assert points_csv.exists()
    assert not pd.read_csv(points_csv).empty
    assert (tmp_path / "meshes" / "Structure_1.ply").exists()
    assert html.exists()
    assert "Plotly.newPlot" in html.read_text(encoding="utf-8")


def _feature(geom, structure: str):
    return {
        "type": "Feature",
        "properties": {"structure": structure},
        "geometry": mapping(geom),
    }


def _write_geojson(path, features):
    path.write_text(
        json.dumps({"type": "FeatureCollection", "features": features}),
        encoding="utf-8",
    )


def _union_iou(features_a, features_b):
    union_a = shape(features_a[0]["geometry"]).union(shape(features_a[1]["geometry"]))
    union_b = shape(features_b[0]["geometry"]).union(shape(features_b[1]["geometry"]))
    return union_a.intersection(union_b).area / union_a.union(union_b).area
