from __future__ import annotations

import json

import numpy as np
import pandas as pd
from shapely.geometry import box, mapping, shape

from histoseg.threed import (
    ThreeDContourReconstructionConfig,
    ThreeDContourReconstructionResult,
    run_3d_contour_reconstruction,
)
from histoseg.threed.cli import main
from histoseg.threed.soft_alignment import (
    _FeatureRecord,
    _check_tps_topology,
    _compute_boundary_curvature,
    _filter_landmark_outliers,
    _filter_mutual_nn_landmarks,
    _generate_curvature_landmarks,
    _select_rbf_smoothing_cv,
    _zero_anchor_landmarks,
)


def test_3d_soft_alignment_preserves_geojson_properties_and_improves_iou(tmp_path):
    fixed_path, moving_path = _write_pair(tmp_path)

    result = run_3d_contour_reconstruction(
        ThreeDContourReconstructionConfig(
            fixed_geojson=fixed_path,
            moving_hard_aligned_geojson=moving_path,
            out_dir=tmp_path / "out",
            sampling_distance_um=5.0,
            max_landmark_distance_um=3.0,
            rbf_neighbors=None,
            rbf_smoothing=1e-4,
            diagnostic_structure="Structure 2",
            dpi=80,
        )
    )

    assert result.soft_aligned_geojson.exists()
    assert result.metrics_csv.exists()
    assert result.landmarks_csv.exists()
    assert result.summary_json.exists()
    assert result.diagnostic_report_png and result.diagnostic_report_png.exists()
    assert result.residuals_csv and result.residuals_csv.exists()

    soft = json.loads(result.soft_aligned_geojson.read_text(encoding="utf-8"))
    assert len(soft["features"]) == 2
    assert soft["features"][0]["properties"]["structure"] == "Structure 1"
    assert soft["features"][0]["properties"]["note"] == "preserve-me"
    assert soft["features"][0]["geometry"]["type"] == "Polygon"

    summary = json.loads(result.summary_json.read_text(encoding="utf-8"))
    assert summary["qc"]["geometry_status_counts"] == {"valid": 2}
    assert (
        summary["qc"]["union_iou_soft_after"]
        >= summary["qc"]["union_iou_hard_before_soft"]
    )
    assert summary["landmarks"]["boundary_landmark_count"] > 0
    assert summary["landmarks"]["zero_anchor_count"] == 16
    assert summary["qc"]["topology_check"]["valid"] is True
    assert summary["qc"]["topology_check"]["checked_cells"] > 0


def test_3d_soft_alignment_can_use_normal_aware_landmark_matching(tmp_path):
    fixed_path, moving_path = _write_pair(tmp_path)

    result = run_3d_contour_reconstruction(
        ThreeDContourReconstructionConfig(
            fixed_geojson=fixed_path,
            moving_hard_aligned_geojson=moving_path,
            out_dir=tmp_path / "normal_aware",
            sampling_distance_um=5.0,
            max_landmark_distance_um=3.0,
            landmark_normal_weight_um=10.0,
            landmark_candidate_count=4,
            rbf_neighbors=None,
            save_preview_png=False,
        )
    )

    landmarks = result.landmarks_csv.read_text(encoding="utf-8")
    summary = json.loads(result.summary_json.read_text(encoding="utf-8"))
    assert "normal_aware_kdtree" in landmarks
    assert summary["method"]["landmark_matching"] == "normal_aware_kdtree"
    assert summary["qc"]["topology_check"]["valid"] is True


def test_tps_topology_check_rejects_grid_fold():
    geom = box(0, 0, 10, 10)
    record = _FeatureRecord(
        feature={"type": "Feature", "properties": {"structure": "Structure 1"}},
        group="Structure 1",
        geometry=geom,
    )

    class FlipModel:
        def warp(self, xy):
            arr = xy.copy()
            arr[:, 0] *= -1.0
            return arr

    topology = _check_tps_topology(
        [record],
        FlipModel(),
        ThreeDContourReconstructionConfig(
            fixed_geojson="fixed.geojson",
            moving_hard_aligned_geojson="moving.geojson",
            topology_grid_size=6,
        ),
    )

    assert topology["valid"] is False
    assert topology["folded_cell_count"] > 0


def test_3d_cli_help_exits_successfully(capsys):
    try:
        main(["--help"])
    except SystemExit as exc:
        assert exc.code == 0

    assert "Run HistoSeg 3D Analysis workflows" in capsys.readouterr().out


def test_3d_cli_reconstruct_parses_landmark_diagnostic_parameters(monkeypatch, tmp_path, capsys):
    import histoseg.threed.cli as cli

    calls: list[ThreeDContourReconstructionConfig] = []

    def fake_run(cfg):
        calls.append(cfg)
        return ThreeDContourReconstructionResult(
            out_dir=tmp_path / "out",
            soft_aligned_geojson=tmp_path / "out" / "soft_aligned_contours.geojson",
            metrics_csv=tmp_path / "out" / "soft_tps_alignment_metrics.csv",
            landmarks_csv=tmp_path / "out" / "soft_tps_landmarks.csv",
            summary_json=tmp_path / "out" / "soft_tps_alignment_summary.json",
        )

    monkeypatch.setattr(cli, "run_3d_contour_reconstruction", fake_run)
    base_args = [
        "reconstruct",
        "--fixed-geojson",
        str(tmp_path / "fixed.geojson"),
        "--moving-hard-aligned-geojson",
        str(tmp_path / "moving.geojson"),
        "--out-dir",
        str(tmp_path / "out"),
    ]

    cli.main(base_args)
    cli.main(
        base_args
        + [
            "--landmarks-per-structure",
            "420",
            "--diagnostic-structure-landmarks",
            "900",
            "--diagnostic-structure",
            "Structure 3",
        ]
    )
    cli.main(
        base_args
        + [
            "--landmarks-per-structure",
            "0",
            "--diagnostic-structure-landmarks",
            "0",
            "--diagnostic-structure",
            "",
        ]
    )

    capsys.readouterr()
    assert calls[0].landmarks_per_structure == 260
    assert calls[0].diagnostic_structure_landmarks == 620
    assert calls[0].diagnostic_structure == "Structure 5"
    assert calls[1].landmarks_per_structure == 420
    assert calls[1].diagnostic_structure_landmarks == 900
    assert calls[1].diagnostic_structure == "Structure 3"
    assert calls[2].landmarks_per_structure is None
    assert calls[2].diagnostic_structure_landmarks is None
    assert calls[2].diagnostic_structure is None


def test_3d_cli_reconstruct_writes_expected_outputs(tmp_path):
    fixed_path, moving_path = _write_pair(tmp_path)
    out_dir = tmp_path / "cli_out"

    main(
        [
            "reconstruct",
            "--fixed-geojson",
            str(fixed_path),
            "--moving-hard-aligned-geojson",
            str(moving_path),
            "--out-dir",
            str(out_dir),
            "--sampling-distance-um",
            "5",
            "--max-landmark-distance-um",
            "3",
            "--rbf-neighbors",
            "0",
            "--diagnostic-structure",
            "Structure 2",
            "--dpi",
            "80",
        ]
    )

    assert (out_dir / "soft_aligned_contours.geojson").exists()
    assert (out_dir / "soft_tps_alignment_metrics.csv").exists()
    assert (out_dir / "soft_tps_landmarks.csv").exists()
    assert (out_dir / "soft_tps_alignment_summary.json").exists()
    assert (out_dir / "soft_tps_overlay_hard_before.png").exists()
    assert (out_dir / "soft_tps_overlay_soft_after.png").exists()
    assert (out_dir / "soft_tps_diagnostic_report.png").exists()


def test_zero_anchor_landmarks_uses_configured_count():
    records = [
        _FeatureRecord(
            feature={"type": "Feature", "properties": {}},
            group="S",
            geometry=box(0, 0, 10, 10),
        )
    ]
    for count in (4, 8, 16, 32):
        anchors = _zero_anchor_landmarks(records, zero_anchor_count=count)
        assert len(anchors) == count
        assert (anchors["kind"] == "anchor").all()
        assert (anchors["dx"] == 0.0).all()
        assert (anchors["dy"] == 0.0).all()


def test_compute_boundary_curvature_returns_high_values_at_corners():
    """Box corners should have significantly higher curvature than flat edges."""
    boundary = box(0, 0, 10, 10).boundary
    pts, curvatures = _compute_boundary_curvature(boundary, spacing=0.5)
    assert len(pts) > 0
    assert len(curvatures) == len(pts)
    # Corners should stand out; not all curvatures should be near zero.
    assert curvatures.max() > 0.1


def test_curvature_landmarks_do_not_flood_for_box_geometry():
    """curvature_landmark_weight=0.5 should generate only corner-region extras for a box."""
    fixed_records = [
        _FeatureRecord(
            feature={"type": "Feature", "properties": {}},
            group="S",
            geometry=box(0, 0, 100, 100),
        )
    ]
    moving_records = [
        _FeatureRecord(
            feature={"type": "Feature", "properties": {}},
            group="S",
            geometry=box(1, 0.5, 101, 100.5),
        )
    ]
    cfg = ThreeDContourReconstructionConfig(
        fixed_geojson="f.geojson",
        moving_hard_aligned_geojson="m.geojson",
        curvature_landmark_weight=0.5,
        max_landmark_distance_um=5.0,
        sampling_distance_um=10.0,
    )
    extras = _generate_curvature_landmarks(fixed_records, moving_records, cfg)
    # For a box, only corner-region points should pass the strict-positive curvature filter.
    # We should get some extras (near corners) but not dozens of flat-edge points.
    assert len(extras) < 20 if not extras.empty else True


def test_filter_mutual_nn_landmarks_removes_non_mutual():
    """Create a landmark whose reverse projection is far from the source; it should be removed."""
    moving_boundary = box(0, 0, 10, 10)
    moving_records = [
        _FeatureRecord(
            feature={"type": "Feature", "properties": {}},
            group="S",
            geometry=moving_boundary,
        )
    ]
    cfg = ThreeDContourReconstructionConfig(
        fixed_geojson="f.geojson",
        moving_hard_aligned_geojson="m.geojson",
        sampling_distance_um=5.0,
        mutual_nn_check=True,
    )
    # Valid landmark: src on left edge, dst very close → reverse maps back near src.
    # abs_tolerance = 2 * 5 = 10.0
    # Invalid landmark: src on left edge, dst far away (50,50).
    #   Reverse projects (50,50) onto box boundary → corner (10,10).
    #   rev_dist = ||(0,5)-(10,10)|| ≈ 11.18 > 10 → removed.
    landmarks = pd.DataFrame(
        [
            {
                "kind": "boundary",
                "structure": "S",
                "src_x": 0.0,  # on left edge of moving boundary
                "src_y": 5.0,
                "dst_x": 0.2,  # very close – reverse projects back near src
                "dst_y": 5.0,
                "dx": 0.2,
                "dy": 0.0,
                "source_distance_um": 0.2,
            },
            {
                "kind": "boundary",
                "structure": "S",
                "src_x": 0.0,  # on left edge of moving boundary
                "src_y": 5.0,
                "dst_x": 50.0,  # very far – reverse projects to (10,10), rev_dist ≈ 11.18
                "dst_y": 50.0,
                "dx": 50.0,
                "dy": 50.0,
                "source_distance_um": 70.7,
            },
        ]
    )
    filtered = _filter_mutual_nn_landmarks(landmarks, moving_records, cfg)
    # The far-target landmark should be dropped (rev_dist ≈ 11.18 > 2*5=10).
    assert len(filtered) < len(landmarks)
    assert float(filtered.iloc[0]["src_y"]) == 5.0
    assert float(filtered.iloc[0]["dst_x"]) < 5.0


def test_filter_landmark_outliers_removes_anomalous_displacements():
    """An outlier with wildly different displacement should be removed by MAD filter."""
    normal_rows = [
        {"kind": "boundary", "structure": "S", "src_x": float(i), "src_y": 0.0,
         "dst_x": float(i) + 1.0, "dst_y": 0.1, "dx": 1.0, "dy": 0.1,
         "source_distance_um": 1.0}
        for i in range(20)
    ]
    # Add one wild outlier
    normal_rows.append(
        {"kind": "boundary", "structure": "S", "src_x": 100.0, "src_y": 0.0,
         "dst_x": 200.0, "dst_y": 100.0, "dx": 100.0, "dy": 100.0,
         "source_distance_um": 141.0}
    )
    lm = pd.DataFrame(normal_rows)
    filtered = _filter_landmark_outliers(lm, mad_threshold=3.5)
    assert len(filtered) < len(lm)
    # The outlier should be gone
    assert all(filtered["dx"].abs() < 10.0)


def test_select_rbf_smoothing_cv_returns_valid_candidate():
    """Auto smoothing should return one of the candidates."""
    rng = np.random.default_rng(0)
    n = 50
    src = rng.random((n, 2))
    dst = rng.random((n, 2)) * 0.01  # small displacement
    cfg = ThreeDContourReconstructionConfig(
        fixed_geojson="f.geojson",
        moving_hard_aligned_geojson="m.geojson",
        rbf_smoothing_candidates=(1e-5, 1e-4, 1e-3, 1e-2),
    )
    chosen = _select_rbf_smoothing_cv(src, dst, cfg, neighbors=None)
    assert chosen in (1e-5, 1e-4, 1e-3, 1e-2)


def test_icp_iterations_run_in_soft_alignment(tmp_path):
    """icp_iterations=2 should still produce valid output and improve or maintain IoU."""
    fixed_path, moving_path = _write_pair(tmp_path)
    result = run_3d_contour_reconstruction(
        ThreeDContourReconstructionConfig(
            fixed_geojson=fixed_path,
            moving_hard_aligned_geojson=moving_path,
            out_dir=tmp_path / "icp_out",
            sampling_distance_um=5.0,
            max_landmark_distance_um=3.0,
            rbf_neighbors=None,
            rbf_smoothing=1e-4,
            icp_iterations=2,
            save_preview_png=False,
        )
    )
    summary = json.loads(result.summary_json.read_text(encoding="utf-8"))
    assert summary["method"]["icp_iterations_configured"] == 2
    assert summary["method"]["icp_iterations_executed"] >= 1
    assert result.soft_aligned_geojson.exists()


def test_rbf_smoothing_auto_runs_without_error(tmp_path):
    """rbf_smoothing='auto' should select smoothing via CV and produce output."""
    fixed_path, moving_path = _write_pair(tmp_path)
    result = run_3d_contour_reconstruction(
        ThreeDContourReconstructionConfig(
            fixed_geojson=fixed_path,
            moving_hard_aligned_geojson=moving_path,
            out_dir=tmp_path / "auto_smooth",
            sampling_distance_um=5.0,
            max_landmark_distance_um=3.0,
            rbf_neighbors=None,
            rbf_smoothing="auto",
            rbf_smoothing_candidates=(1e-4, 1e-3),
            save_preview_png=False,
        )
    )
    assert result.soft_aligned_geojson.exists()
    summary = json.loads(result.summary_json.read_text(encoding="utf-8"))
    assert summary["method"]["rbf_smoothing"] == "auto"


def _write_pair(tmp_path):
    fixed_features = [
        _feature(box(0, 0, 10, 10), "Structure 1", note="preserve-me"),
        _feature(box(20, 0, 30, 10), "Structure 2"),
    ]
    moving_features = [
        _feature(box(0.6, 0.2, 10.6, 10.2), "Structure 1", note="preserve-me"),
        _feature(box(20.5, -0.1, 30.5, 9.9), "Structure 2"),
    ]
    fixed_path = tmp_path / "fixed.geojson"
    moving_path = tmp_path / "moving_hard.geojson"
    fixed_path.write_text(
        json.dumps({"type": "FeatureCollection", "features": fixed_features}),
        encoding="utf-8",
    )
    moving_path.write_text(
        json.dumps({"type": "FeatureCollection", "features": moving_features}),
        encoding="utf-8",
    )
    assert _union_iou(fixed_features, moving_features) < 1.0
    return fixed_path, moving_path


def _feature(geom, structure: str, **properties):
    return {
        "type": "Feature",
        "properties": {"structure": structure, **properties},
        "geometry": mapping(geom),
    }


def _union_iou(features_a, features_b):
    union_a = shape(features_a[0]["geometry"]).union(shape(features_a[1]["geometry"]))
    union_b = shape(features_b[0]["geometry"]).union(shape(features_b[1]["geometry"]))
    return union_a.intersection(union_b).area / union_a.union(union_b).area
