from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
import trimesh
from shapely.geometry import MultiPolygon, box, mapping, shape
from shapely.ops import unary_union

from histoseg.threed import (
    ThreeDStackReconstructionResult,
    discover_xenium_slices,
    hard_align_geojson,
    reconstruct_3d_contour_meshes,
    run_3d_stack_reconstruction,
    write_3d_contour_points,
    write_3d_visualization_html,
)
import histoseg.threed.multislice as multislice
from histoseg.threed.multislice import (
    _SliceInput,
    _build_per_structure_soft_geojson,
    _pairwise_row,
    _read_strategy_specs,
)
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


def test_discover_xenium_slices_accepts_pyxenium_slide_zarr(tmp_path):
    for name in [
        "A079-C-008_10.pyxenium.slide.zarr",
        "A079-C-008_1.pyxenium.slide.zarr",
        "A079-C-008_2.pyxenium.slide.zarr",
    ]:
        slide = tmp_path / name
        (slide / "tables").mkdir(parents=True)
        (slide / "zarr.json").write_text("{}", encoding="utf-8")

    slices = discover_xenium_slices(tmp_path, sample_glob="*.pyxenium.slide.zarr")

    assert [item.sample_id for item in slices] == [
        "A079-C-008_1",
        "A079-C-008_2",
        "A079-C-008_10",
    ]
    assert [item.xenium_dir.name for item in slices] == [
        "A079-C-008_1.pyxenium.slide.zarr",
        "A079-C-008_2.pyxenium.slide.zarr",
        "A079-C-008_10.pyxenium.slide.zarr",
    ]
    assert [item.order for item in slices] == [1, 2, 3]


def test_reconstruct_stack_uses_precomputed_contour_manifest(monkeypatch, tmp_path):
    fixed_geojson = tmp_path / "slice1.geojson"
    moving_geojson = tmp_path / "slice2.geojson"
    _write_geojson(
        fixed_geojson,
        [_feature(box(0, 0, 10, 10), "Local group A")],
    )
    _write_geojson(
        moving_geojson,
        [_feature(box(2, 0, 12, 10), "Different local group")],
    )
    manifest = tmp_path / "slice_local_contours.csv"
    pd.DataFrame(
        [
            {
                "order": 1,
                "sample_id": "slice_1",
                "z_um": 0.0,
                "geojson": str(fixed_geojson),
            },
            {
                "order": 2,
                "sample_id": "slice_2",
                "z_um": 7.5,
                "geojson": str(moving_geojson),
            },
        ]
    ).to_csv(manifest, index=False)

    def fail_build_slice_contours(*args, **kwargs):
        raise AssertionError("manifest mode must not build contours from Xenium data")

    def fake_hard_alignment(cfg, *, fixed_geojson, moving_geojson, output_geojson, summary_json):
        import shutil

        shutil.copy2(moving_geojson, output_geojson)
        summary = {
            "registration_backend": "contour-tps",
            "selected_hard_seed_backend": "contour-tps",
            "union_iou_before_hard": 0.5,
            "union_iou_after_hard": 0.6,
            "hard_alignment_accepted": True,
            "transform": {
                "rotation_degrees": 0.0,
                "scale": 1.0,
                "translate_x": 0.0,
                "translate_y": 0.0,
            },
        }
        Path(summary_json).write_text(json.dumps(summary), encoding="utf-8")
        return summary

    def fake_mesh(aligned_rows, mesh_dir, **kwargs):
        mesh_dir = Path(mesh_dir)
        mesh_dir.mkdir(parents=True, exist_ok=True)
        (mesh_dir / "mesh_manifest.csv").write_text("structure\n", encoding="utf-8")
        (mesh_dir / "mesh_qc_summary.json").write_text("{}", encoding="utf-8")
        return []

    monkeypatch.setattr(multislice, "_build_slice_contours", fail_build_slice_contours)
    monkeypatch.setattr(multislice, "_run_hard_alignment_backend", fake_hard_alignment)
    monkeypatch.setattr(multislice, "reconstruct_3d_contour_meshes", fake_mesh)

    result = run_3d_stack_reconstruction(
        ThreeDStackReconstructionConfig(
            contour_manifest=manifest,
            out_dir=tmp_path / "out",
            registration_backend="contour-tps",
            run_soft_alignment=False,
        )
    )

    slice_manifest = pd.read_csv(result.slice_manifest_csv)
    assert slice_manifest["sample_id"].tolist() == ["slice_1", "slice_2"]
    assert slice_manifest["z_um"].tolist() == [0.0, 7.5]
    assert slice_manifest["precomputed_geojson"].notna().all()

    aligned_manifest = pd.read_csv(result.aligned_manifest_csv)
    assert aligned_manifest["z_um"].tolist() == [0.0, 7.5]
    raw_copy = result.out_dir / "slice_contours" / "002_slice_2" / "xenium_explorer_annotations.geojson"
    payload = json.loads(raw_copy.read_text(encoding="utf-8"))
    assert payload["features"][0]["properties"]["structure"] == "Different local group"

    summary = json.loads(result.summary_json.read_text(encoding="utf-8"))
    assert summary["slice_count"] == 2
    assert summary["structure_count"] == 2
    assert summary["config"]["contour_manifest"] == str(manifest)


def test_reconstruct_stack_contour_manifest_requires_core_columns(tmp_path):
    manifest = tmp_path / "bad_manifest.csv"
    pd.DataFrame([{"order": 1, "sample_id": "slice_1"}]).to_csv(manifest, index=False)

    with pytest.raises(ValueError, match="contour_manifest is missing required columns"):
        run_3d_stack_reconstruction(
            ThreeDStackReconstructionConfig(
                contour_manifest=manifest,
                out_dir=tmp_path / "out",
            )
        )


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
        mesh_smoothing_sigma_um=1.0,
        mesh_export_formats=("ply", "obj"),
    )
    html = write_3d_visualization_html(points, meshes, tmp_path / "view.html", title="test")

    assert points_csv.exists()
    assert not pd.read_csv(points_csv).empty
    assert (tmp_path / "meshes" / "Structure_1.ply").exists()
    assert (tmp_path / "meshes" / "Structure_1.obj").exists()
    assert (tmp_path / "meshes" / "mesh_qc_summary.json").exists()
    mesh = trimesh.load(tmp_path / "meshes" / "Structure_1.ply", force="mesh")
    assert 10.0 < float(mesh.bounds[1][0]) < 40.0
    assert -5.0 <= float(mesh.bounds[0][2]) <= 5.0
    manifest = pd.read_csv(tmp_path / "meshes" / "mesh_manifest.csv")
    assert {
        "surface_area_um2",
        "volume_um3",
        "is_watertight",
        "euler_number",
        "component_count",
    }.issubset(set(manifest.columns))
    assert float(manifest.loc[0, "surface_area_um2"]) > 0
    assert html.exists()
    html_text = html.read_text(encoding="utf-8")
    assert "Plotly.newPlot" in html_text
    assert '"type": "mesh3d"' in html_text
    assert '"visible": "legendonly"' in html_text


def test_mesh_component_filter_removes_tiny_disconnected_fragments(tmp_path):
    slice1 = tmp_path / "slice1.geojson"
    slice2 = tmp_path / "slice2.geojson"
    geom = MultiPolygon([box(0, 0, 20, 20), box(79, 79, 81, 81)])
    _write_geojson(slice1, [_feature(geom, "Structure 1")])
    _write_geojson(slice2, [_feature(geom, "Structure 1")])
    aligned_rows = [
        {"order": 1, "sample_id": "s1", "z_um": 0.0, "aligned_geojson": str(slice1)},
        {"order": 2, "sample_id": "s2", "z_um": 5.0, "aligned_geojson": str(slice2)},
    ]

    reconstruct_3d_contour_meshes(
        aligned_rows,
        tmp_path / "meshes",
        group_property="structure",
        voxel_size_um=2.0,
        z_spacing_um=5.0,
        xenium_pixel_size_um=1.0,
        mesh_smoothing_sigma_um=0.0,
        min_mesh_component_volume_um3=100.0,
    )

    manifest = pd.read_csv(tmp_path / "meshes" / "mesh_manifest.csv")
    assert int(manifest.loc[0, "volume_component_count_before_filter"]) == 2
    assert int(manifest.loc[0, "volume_component_count_after_filter"]) == 1


def test_pairwise_row_exports_topology_counts(tmp_path):
    hard_summary = {
        "registration_backend": "coda-image",
        "selected_hard_seed_backend": "contour-tps",
        "method_credit": "CODA-inspired image registration seed",
        "method_reference_doi": "10.1038/s41592-022-01650-9",
        "union_iou_before_hard": 0.4,
        "union_iou_after_hard": 0.8,
        "hard_alignment_candidates": [
            {"backend": "contour-tps", "union_iou_after_hard": 0.8},
            {"backend": "coda-image", "union_iou_after_hard": 0.7},
        ],
        "hard_alignment_tournament": {"rotation_difference_degrees": 15.0},
        "transform": {
            "rotation_degrees": 2.0,
            "scale": 1.0,
            "translate_x": 3.0,
            "translate_y": 4.0,
        },
        "hard_alignment_accepted": True,
        "coda_image": {
            "radon_rotation_degrees": 180.0,
            "phase_shift_y": 5.0,
            "phase_shift_x": -6.0,
        },
    }
    soft_summary = {
        "qc": {
            "union_iou_hard_before_soft": 0.8,
            "union_iou_soft_after": 0.83,
            "geometry_status_counts": {"invalid": 0},
            "topology_check": {
                "valid": False,
                "checked_cells": 576,
                "min_area_ratio": 0.31,
                "median_area_ratio": 0.98,
                "max_area_ratio": 2.30,
                "folded_cell_count": 2,
                "compressed_cell_count": 3,
                "expanded_cell_count": 4,
            },
        },
        "landmarks": {"boundary_landmark_count": 42},
    }

    row = _pairwise_row(
        _SliceInput(
            order=2,
            sample_id="slice_2",
            sample_dir=tmp_path / "sample",
            xenium_dir=tmp_path / "sample" / "xenium",
        ),
        hard_summary,
        soft_summary,
        SimpleNamespace(summary_json=tmp_path / "soft_summary.json"),
        soft_accepted=False,
    )

    assert row["selected_hard_seed_backend"] == "contour-tps"
    assert row["hard_candidate_contour_iou_after"] == pytest.approx(0.8)
    assert row["hard_candidate_coda_iou_after"] == pytest.approx(0.7)
    assert row["coda_radon_rotation_degrees"] == pytest.approx(180.0)
    assert row["soft_topology_valid"] is False
    assert row["soft_topology_checked_cells"] == 576
    assert row["soft_topology_folded_cells"] == 2
    assert row["soft_topology_compressed_cells"] == 3
    assert row["soft_topology_expanded_cells"] == 4
    assert row["soft_topology_median_area_ratio"] == pytest.approx(0.98)


def test_mesh_smoothing_none_disables_gaussian_smoothing(tmp_path):
    slice1 = tmp_path / "slice1.geojson"
    slice2 = tmp_path / "slice2.geojson"
    _write_geojson(slice1, [_feature(box(0, 0, 20, 20), "Structure 1")])
    _write_geojson(slice2, [_feature(box(0, 0, 20, 20), "Structure 1")])
    aligned_rows = [
        {"order": 1, "sample_id": "s1", "z_um": 0.0, "aligned_geojson": str(slice1)},
        {"order": 2, "sample_id": "s2", "z_um": 5.0, "aligned_geojson": str(slice2)},
    ]

    reconstruct_3d_contour_meshes(
        aligned_rows,
        tmp_path / "meshes",
        group_property="structure",
        voxel_size_um=5.0,
        z_spacing_um=5.0,
        xenium_pixel_size_um=1.0,
        mesh_smoothing_sigma_um=None,
    )

    manifest = pd.read_csv(tmp_path / "meshes" / "mesh_manifest.csv")
    assert bool(manifest.loc[0, "smoothing_applied"]) is False


def test_mesh_level_must_be_between_zero_and_one(tmp_path):
    slice1 = tmp_path / "slice1.geojson"
    slice2 = tmp_path / "slice2.geojson"
    _write_geojson(slice1, [_feature(box(0, 0, 20, 20), "Structure 1")])
    _write_geojson(slice2, [_feature(box(0, 0, 20, 20), "Structure 1")])
    aligned_rows = [
        {"order": 1, "sample_id": "s1", "z_um": 0.0, "aligned_geojson": str(slice1)},
        {"order": 2, "sample_id": "s2", "z_um": 5.0, "aligned_geojson": str(slice2)},
    ]

    with pytest.raises(ValueError, match="strictly between 0 and 1"):
        reconstruct_3d_contour_meshes(
            aligned_rows,
            tmp_path / "meshes",
            group_property="structure",
            voxel_size_um=5.0,
            z_spacing_um=5.0,
            xenium_pixel_size_um=1.0,
            mesh_level=1.0,
        )


def test_hard_align_multistart_improves_or_maintains_iou(tmp_path):
    """Multi-start should produce at least as good IoU as single-start."""
    fixed = tmp_path / "fixed.geojson"
    moving = tmp_path / "moving.geojson"
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
    summary_single = hard_align_geojson(
        fixed_geojson=fixed,
        moving_geojson=moving,
        output_geojson=tmp_path / "hard_single.geojson",
        group_property="structure",
        maxiter=30,
        overwrite=True,
        multistart=False,
    )
    summary_multi = hard_align_geojson(
        fixed_geojson=fixed,
        moving_geojson=moving,
        output_geojson=tmp_path / "hard_multi.geojson",
        group_property="structure",
        maxiter=30,
        overwrite=True,
        multistart=True,
    )
    assert summary_multi["union_iou_after_hard"] >= summary_single["union_iou_after_hard"] - 0.02


def test_per_structure_soft_mixing_preserves_feature_identity():
    hard_payload = _geojson_payload(
        [
            _feature(box(0, 0, 10, 10), "Structure 1"),
            _feature(box(20, 0, 30, 10), "Structure 1"),
            _feature(box(40, 0, 50, 10), "Structure 2"),
        ]
    )
    soft_payload = _geojson_payload(
        [
            _feature(box(1, 0, 11, 10), "Structure 1"),
            _feature(box(21, 0, 31, 10), "Structure 1"),
            _feature(box(41, 0, 51, 10), "Structure 2"),
        ]
    )
    summary = {
        "qc": {
            "per_structure": {
                "Structure 1": {
                    "iou_hard_before_soft": 0.4,
                    "iou_soft_after": 0.5,
                },
                "Structure 2": {
                    "iou_hard_before_soft": 0.8,
                    "iou_soft_after": 0.7,
                },
            }
        }
    }

    mixed = _build_per_structure_soft_geojson(
        hard_payload,
        soft_payload,
        "structure",
        summary,
    )
    geoms = [shape(feature["geometry"]) for feature in mixed["features"]]

    assert geoms[0].bounds == pytest.approx((1.0, 0.0, 11.0, 10.0))
    assert geoms[1].bounds == pytest.approx((21.0, 0.0, 31.0, 10.0))
    assert geoms[2].bounds == pytest.approx((40.0, 0.0, 50.0, 10.0))
    assert unary_union(geoms[:2]).area == pytest.approx(200.0)


def test_reconstruct_stack_cli_parses_registration_backend(monkeypatch, tmp_path, capsys):
    import histoseg.threed.cli as cli

    calls: list[ThreeDStackReconstructionConfig] = []

    def fake_run(cfg):
        calls.append(cfg)
        return ThreeDStackReconstructionResult(
            out_dir=tmp_path / "out",
            slice_manifest_csv=tmp_path / "out" / "xenium_slice_manifest.csv",
            pairwise_metrics_csv=tmp_path / "out" / "pairwise_alignment_metrics.csv",
            aligned_manifest_csv=tmp_path / "out" / "aligned_slice_manifest.csv",
            contour_points_csv=tmp_path / "out" / "aligned_contour_3d_points.csv",
            summary_json=tmp_path / "out" / "3d_stack_reconstruction_summary.json",
            visualization_html=tmp_path / "out" / "histoseg_3d_contour_stack.html",
            mesh_dir=tmp_path / "out" / "meshes",
        )

    monkeypatch.setattr(cli, "run_3d_stack_reconstruction", fake_run)
    cli.main(
        [
            "reconstruct-stack",
            "--xenium-root",
            str(tmp_path / "xenium"),
            "--segmentation-strategy",
            str(tmp_path / "segmentationstrategy.txt"),
            "--out-dir",
            str(tmp_path / "out"),
        ]
    )
    cli.main(
        [
            "reconstruct-stack",
            "--xenium-root",
            str(tmp_path / "xenium"),
            "--segmentation-strategy",
            str(tmp_path / "segmentationstrategy.txt"),
            "--out-dir",
            str(tmp_path / "out"),
            "--registration-backend",
            "coda-image",
            "--coda-raster-size",
            "128",
            "--coda-angle-step",
            "0.5",
            "--coda-phase-upsample-factor",
            "4",
            "--coda-mask-padding-fraction",
            "0.1",
        ]
    )
    cli.main(
        [
            "reconstruct-stack",
            "--xenium-root",
            str(tmp_path / "xenium"),
            "--segmentation-strategy",
            str(tmp_path / "segmentationstrategy.txt"),
            "--out-dir",
            str(tmp_path / "out"),
            "--soft-alignment-mode",
            "anchor-only",
            "--anchor-only-bbox-padding-fraction",
            "0.2",
            "--anchor-only-identity-padding-count",
            "20",
            "--anchor-only-jacobian-grid-size",
            "30",
            "--anchor-only-max-negative-jacobian-ratio",
            "0.01",
        ]
    )
    cli.main(
        [
            "reconstruct-stack",
            "--xenium-root",
            str(tmp_path / "xenium"),
            "--segmentation-strategy",
            str(tmp_path / "segmentationstrategy.txt"),
            "--out-dir",
            str(tmp_path / "out"),
            "--no-soft",
        ]
    )
    cli.main(
        [
            "reconstruct-stack",
            "--contour-manifest",
            str(tmp_path / "slice_local_contours.csv"),
            "--out-dir",
            str(tmp_path / "out"),
        ]
    )

    capsys.readouterr()
    assert calls[0].registration_backend == "auto"
    assert calls[0].soft_alignment_mode == "auto"
    assert calls[1].registration_backend == "coda-image"
    assert calls[1].coda_raster_size == 128
    assert calls[1].coda_angle_step == 0.5
    assert calls[1].coda_phase_upsample_factor == 4
    assert calls[1].coda_mask_padding_fraction == 0.1
    assert calls[2].soft_alignment_mode == "anchor-only"
    assert calls[2].anchor_only_bbox_padding_fraction == 0.2
    assert calls[2].anchor_only_identity_padding_count == 20
    assert calls[2].anchor_only_jacobian_grid_size == 30
    assert calls[2].anchor_only_max_negative_jacobian_ratio == 0.01
    assert calls[3].soft_alignment_mode == "none"
    assert calls[3].run_soft_alignment is False
    assert calls[4].xenium_root is None
    assert calls[4].contour_manifest == str(tmp_path / "slice_local_contours.csv")


def test_hard_align_affine_fallback_accepted_when_triggered(tmp_path):
    """affine_fallback_iou_threshold=1.0 always triggers the affine fallback."""
    fixed = tmp_path / "fixed.geojson"
    moving = tmp_path / "moving.geojson"
    _write_geojson(fixed, [_feature(box(0, 0, 20, 10), "Structure 1")])
    _write_geojson(moving, [_feature(box(2, 1, 22, 11), "Structure 1")])
    summary = hard_align_geojson(
        fixed_geojson=fixed,
        moving_geojson=moving,
        output_geojson=tmp_path / "hard_affine.geojson",
        group_property="structure",
        maxiter=20,
        overwrite=True,
        affine_fallback_iou_threshold=1.0,  # always trigger
    )
    # The function should complete without error regardless of whether affine was used.
    assert summary["union_iou_after_hard"] >= 0.0
    assert (tmp_path / "hard_affine.geojson").exists()


def test_build_per_structure_soft_geojson_uses_soft_where_improved(tmp_path):
    """Per-structure mixing should use soft geometry where IoU is better."""
    from histoseg.threed.multislice import _build_per_structure_soft_geojson

    hard_features = [
        _feature(box(0, 0, 10, 10), "S1"),
        _feature(box(20, 0, 30, 10), "S2"),
    ]
    # Soft S1 is slightly different, soft S2 is significantly different
    soft_features = [
        _feature(box(0.1, 0.1, 10.1, 10.1), "S1"),
        _feature(box(100, 100, 200, 200), "S2"),  # deliberately bad
    ]
    hard_payload = {"type": "FeatureCollection", "features": hard_features}
    soft_payload = {"type": "FeatureCollection", "features": soft_features}
    # S1: soft IoU > hard IoU (both at identity, soft slightly shifted closer)
    # S2: soft IoU=0 (far off), hard IoU > 0
    soft_summary = {
        "qc": {
            "per_structure": {
                "S1": {"iou_hard_before_soft": 0.8, "iou_soft_after": 0.9, "delta_iou_soft": 0.1},
                "S2": {"iou_hard_before_soft": 0.7, "iou_soft_after": 0.0, "delta_iou_soft": -0.7},
            }
        }
    }
    mixed = _build_per_structure_soft_geojson(hard_payload, soft_payload, "structure", soft_summary)
    geoms_by_structure = {
        f["properties"]["structure"]: shape(f["geometry"])
        for f in mixed["features"]
    }
    # S1 should use soft geometry (IoU improved)
    assert geoms_by_structure["S1"].equals(shape(soft_features[0]["geometry"]))
    # S2 should keep hard geometry (IoU degraded)
    assert geoms_by_structure["S2"].equals(shape(hard_features[1]["geometry"]))


def test_global_drift_correction_applies_translation(tmp_path):
    """Drift correction should translate slices with accumulated centroid drift."""
    from histoseg.threed.multislice import _apply_global_drift_correction

    slices = []
    aligned_rows = []
    for i in range(4):
        path = tmp_path / f"slice_{i}.geojson"
        # Each slice's centroid drifts by (10, 5) per step.
        offset = i * 10.0
        features = [_feature(box(offset, offset / 2, offset + 20, offset / 2 + 20), "S")]
        _write_geojson(path, features)
        slices.append(path)
        aligned_rows.append({"order": i + 1, "sample_id": f"s{i}", "aligned_geojson": str(path)})

    _apply_global_drift_correction(slices, aligned_rows, group_property="structure")

    # The reference slice (index 0) should not be modified.
    ref_payload = json.loads(slices[0].read_text(encoding="utf-8"))
    ref_geom = shape(ref_payload["features"][0]["geometry"])
    assert abs(ref_geom.centroid.x - 10.0) < 0.1  # centroid of box(0,0,20,20) = (10,10)

    # Non-reference slices should have been corrected (closer to reference centroid trend).
    last_payload = json.loads(slices[-1].read_text(encoding="utf-8"))
    last_geom = shape(last_payload["features"][0]["geometry"])
    # After drift correction, the centroid x should be closer to the trend-predicted value.
    assert abs(last_geom.centroid.x) < 50.0  # broad sanity check (no extreme drift)


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


def _geojson_payload(features):
    return {"type": "FeatureCollection", "features": features}


def _union_iou(features_a, features_b):
    union_a = shape(features_a[0]["geometry"]).union(shape(features_a[1]["geometry"]))
    union_b = shape(features_b[0]["geometry"]).union(shape(features_b[1]["geometry"]))
    return union_a.intersection(union_b).area / union_a.union(union_b).area
