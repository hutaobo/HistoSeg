from __future__ import annotations

import json

from shapely.geometry import box, mapping, shape

from histoseg.threed import ThreeDContourReconstructionConfig, run_3d_contour_reconstruction
from histoseg.threed.cli import main


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
    assert summary["landmarks"]["zero_anchor_count"] == 8


def test_3d_cli_help_exits_successfully(capsys):
    try:
        main(["--help"])
    except SystemExit as exc:
        assert exc.code == 0

    assert "Run HistoSeg 3D Analysis workflows" in capsys.readouterr().out


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
