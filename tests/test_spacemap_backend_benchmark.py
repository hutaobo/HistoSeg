from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from scripts.compare_registration_backends import (
    BackendComparisonConfig,
    compare_registration_backends,
)
from scripts.import_spacemap_aligned_coordinates import (
    SpaceMapImportConfig,
    import_spacemap_aligned_coordinates,
)
from scripts.harmonize_spacemap_frame import (
    SpaceMapHarmonizationConfig,
    harmonize_spacemap_frame,
)
from scripts.qc_spacemap_coordinate_frame import (
    SpaceMapCoordinateQcConfig,
    run_spacemap_coordinate_frame_qc,
)


def test_import_spacemap_aligned_coordinates_matches_by_barcode_and_converts_pixels(tmp_path):
    manifest = _write_manifest(tmp_path)
    reference = tmp_path / "histoseg_cells.parquet"
    pd.DataFrame(
        {
            "sample_id": ["s2", "s1", "s1", "s2"],
            "barcode": ["b4", "b1", "b2", "b3"],
            "leiden_1_0": ["1", "0", "0", "1"],
            "slice_order": [2, 1, 1, 2],
            "x_3d_um": [40.0, 10.0, 20.0, 30.0],
            "y_3d_um": [4.0, 1.0, 2.0, 3.0],
            "z_um": [5.0, 0.0, 0.0, 5.0],
        }
    ).to_parquet(reference, index=False)
    spacemap_csv = tmp_path / "aligned.csv.gz"
    pd.DataFrame(
        {
            "layer": [0, 0, 1, 1],
            "barcode": ["b1", "b2", "b3", "b4"],
            "x": [100.0, 200.0, 300.0, 400.0],
            "y": [10.0, 20.0, 30.0, 40.0],
        }
    ).to_csv(spacemap_csv, index=False, compression="gzip")

    out_parquet = tmp_path / "space_map_aligned_cells.parquet"
    result = import_spacemap_aligned_coordinates(
        SpaceMapImportConfig(
            spacemap_aligned_csv=spacemap_csv,
            slice_manifest=manifest,
            histoseg_aligned_cells_parquet=reference,
            out_parquet=out_parquet,
            coordinate_unit="pixel",
            pixel_size_um=0.5,
        )
    )

    table = pd.read_parquet(out_parquet)
    assert result.cell_count == 4
    assert result.match_strategy == "sample_barcode"
    assert result.layer_mapping_mode == "zero_based_order"
    assert table["barcode"].tolist() == ["b4", "b1", "b2", "b3"]
    assert table["slice_order"].tolist() == [2, 1, 1, 2]
    np.testing.assert_allclose(table["x_3d_um"], [200.0, 50.0, 100.0, 150.0])
    np.testing.assert_allclose(table["y_3d_um"], [20.0, 5.0, 10.0, 15.0])
    assert table["coordinate_backend"].unique().tolist() == ["spacemap"]
    assert {"sample_id", "barcode", "slice_order", "x_3d_um", "y_3d_um", "z_um"}.issubset(table.columns)
    assert result.summary_json.exists()


def test_import_spacemap_aligned_coordinates_requires_explicit_row_order_matching(tmp_path):
    manifest = _write_manifest(tmp_path)
    reference = tmp_path / "histoseg_cells.parquet"
    pd.DataFrame(
        {
            "sample_id": ["s1", "s1", "s2", "s2"],
            "barcode": ["b1", "b2", "b3", "b4"],
        }
    ).to_parquet(reference, index=False)
    spacemap_csv = tmp_path / "aligned.csv.gz"
    pd.DataFrame(
        {
            "layer": [1, 1, 2, 2],
            "x": [1.0, 2.0, 3.0, 4.0],
            "y": [5.0, 6.0, 7.0, 8.0],
        }
    ).to_csv(spacemap_csv, index=False, compression="gzip")

    with pytest.raises(ValueError, match="barcode"):
        import_spacemap_aligned_coordinates(
            SpaceMapImportConfig(
                spacemap_aligned_csv=spacemap_csv,
                slice_manifest=manifest,
                histoseg_aligned_cells_parquet=reference,
                out_parquet=tmp_path / "fail.parquet",
            )
        )

    result = import_spacemap_aligned_coordinates(
        SpaceMapImportConfig(
            spacemap_aligned_csv=spacemap_csv,
            slice_manifest=manifest,
            histoseg_aligned_cells_parquet=reference,
            out_parquet=tmp_path / "space_map_row_order.parquet",
            allow_row_order_match=True,
        )
    )
    assert result.match_strategy == "within_slice_row_order"
    table = pd.read_parquet(result.out_parquet)
    assert table["x_3d_um"].tolist() == [1.0, 2.0, 3.0, 4.0]


def test_import_spacemap_aligned_coordinates_fails_on_row_order_count_mismatch(tmp_path):
    manifest = _write_manifest(tmp_path)
    reference = tmp_path / "histoseg_cells.parquet"
    pd.DataFrame({"sample_id": ["s1", "s1", "s2"], "barcode": ["b1", "b2", "b3"]}).to_parquet(
        reference,
        index=False,
    )
    spacemap_csv = tmp_path / "aligned.csv.gz"
    pd.DataFrame({"layer": [1, 1, 2, 2], "x": [1, 2, 3, 4], "y": [1, 2, 3, 4]}).to_csv(
        spacemap_csv,
        index=False,
        compression="gzip",
    )

    with pytest.raises(ValueError, match="identical per-slice counts"):
        import_spacemap_aligned_coordinates(
            SpaceMapImportConfig(
                spacemap_aligned_csv=spacemap_csv,
                slice_manifest=manifest,
                histoseg_aligned_cells_parquet=reference,
                out_parquet=tmp_path / "fail.parquet",
                allow_row_order_match=True,
            )
        )


def test_compare_registration_backends_writes_panel_d_outputs(tmp_path):
    histoseg_dir = tmp_path / "histoseg_batch"
    spacemap_dir = tmp_path / "spacemap_batch"
    _write_backend_matrices(histoseg_dir, offset=0.0)
    _write_backend_matrices(spacemap_dir, offset=0.05)

    result = compare_registration_backends(
        BackendComparisonConfig(
            histoseg_batch_dir=histoseg_dir,
            spacemap_batch_dir=spacemap_dir,
            out_dir=tmp_path / "panel_d",
            hotspot="top05",
            focus_structure="Structure 5",
        )
    )

    for path in (
        result.metric_correlations_csv,
        result.sdf_distribution_metrics_csv,
        result.rank_stability_csv,
        result.delta_heatmaps_csv,
        result.delta_heatmap_png,
        result.sdf_ridge_plot_png,
        result.summary_json,
    ):
        assert path.exists()
        assert path.stat().st_size > 0
    delta = pd.read_csv(result.delta_heatmaps_csv)
    assert {"metric", "gene", "structure", "delta_spacemap_minus_histoseg"}.issubset(delta.columns)
    summary = json.loads(result.summary_json.read_text(encoding="utf-8"))
    assert "signed_distance" in summary["metrics_compared"]


def test_qc_spacemap_coordinate_frame_writes_audit_outputs(tmp_path):
    stack_root = tmp_path / "stack"
    stack_root.mkdir()
    histoseg = tmp_path / "histoseg.parquet"
    spacemap = tmp_path / "spacemap.parquet"
    pd.DataFrame(
        {
            "sample_id": ["s1", "s1", "s2", "s2"],
            "slice_order": [1, 1, 2, 2],
            "x_3d_um": [0.0, 10.0, 0.0, 10.0],
            "y_3d_um": [0.0, 10.0, 0.0, 10.0],
            "z_um": [0.0, 0.0, 5.0, 5.0],
        }
    ).to_parquet(histoseg, index=False)
    pd.DataFrame(
        {
            "sample_id": ["s1", "s1", "s2", "s2"],
            "slice_order": [1, 1, 2, 2],
            "x_3d_um": [0.0, 10.0, 0.0, 10.0],
            "y_3d_um": [0.0, 10.0, 0.0, 10.0],
            "z_um": [0.0, 0.0, 5.0, 5.0],
        }
    ).to_parquet(spacemap, index=False)
    pd.DataFrame(
        {
            "slice_order": [1, 1, 1, 2, 2, 2],
            "structure": ["Structure 1"] * 6,
            "polyline_id": ["a", "a", "a", "b", "b", "b"],
            "point_index": [0, 1, 2, 0, 1, 2],
            "x_um": [0.0, 10.0, 10.0, 0.0, 10.0, 10.0],
            "y_um": [0.0, 0.0, 10.0, 0.0, 0.0, 10.0],
        }
    ).to_csv(stack_root / "aligned_contour_3d_points.csv", index=False)

    result = run_spacemap_coordinate_frame_qc(
        SpaceMapCoordinateQcConfig(
            spacemap_aligned_cells_parquet=spacemap,
            histoseg_aligned_cells_parquet=histoseg,
            stack_root=stack_root,
            out_dir=tmp_path / "qc",
            overlay_slices=(1, 2),
        )
    )

    for path in (
        result.summary_json,
        result.by_slice_csv,
        result.bbox_overlay_png,
        result.centroid_delta_png,
        *result.overlay_pngs,
    ):
        assert path.exists()
        assert path.stat().st_size > 0
    summary = json.loads(result.summary_json.read_text(encoding="utf-8"))
    assert summary["cell_counts_match_by_slice"] is True
    assert summary["global_scale_ratio"]["width"] == pytest.approx(1.0)
    assert "coordinate_frame_appears_compatible" in summary["diagnosis"]
    by_slice = pd.read_csv(result.by_slice_csv)
    assert by_slice["centroid_delta_x_um"].tolist() == [0.0, 0.0]


def test_harmonize_spacemap_frame_recovers_global_similarity(tmp_path):
    histoseg = tmp_path / "histoseg.parquet"
    spacemap = tmp_path / "spacemap.parquet"
    dst = np.array(
        [
            [0.0, 0.0],
            [10.0, 0.0],
            [0.0, 10.0],
            [10.0, 10.0],
            [5.0, 20.0],
        ]
    )
    theta = np.deg2rad(20.0)
    rotation = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    scale = 1.5
    translation = np.array([100.0, -30.0])
    src = ((rotation.T @ ((dst - translation).T / scale))).T
    pd.DataFrame(
        {
            "sample_id": ["s1"] * len(dst),
            "barcode": [f"c{i}" for i in range(len(dst))],
            "x_3d_um": dst[:, 0],
            "y_3d_um": dst[:, 1],
            "z_um": 0.0,
        }
    ).to_parquet(histoseg, index=False)
    pd.DataFrame(
        {
            "sample_id": ["s1"] * len(src),
            "barcode": [f"c{i}" for i in range(len(src))],
            "x_3d_um": src[:, 0],
            "y_3d_um": src[:, 1],
            "z_um": 0.0,
        }
    ).to_parquet(spacemap, index=False)

    result = harmonize_spacemap_frame(
        SpaceMapHarmonizationConfig(
            spacemap_aligned_cells_parquet=spacemap,
            histoseg_aligned_cells_parquet=histoseg,
            out_parquet=tmp_path / "spacemap_harmonized.parquet",
            max_points=5,
        )
    )

    assert result.scale == pytest.approx(scale)
    assert result.rotation_degrees == pytest.approx(20.0)
    assert result.rmse_used_after_um < 1e-10
    table = pd.read_parquet(result.out_parquet)
    np.testing.assert_allclose(table[["x_3d_um", "y_3d_um"]].to_numpy(), dst, atol=1e-10)
    assert table["coordinate_backend"].unique().tolist() == ["spacemap_global_similarity_harmonized"]
    summary = json.loads(result.summary_json.read_text(encoding="utf-8"))
    assert summary["local_nonrigid_warp_applied"] is False


def _write_manifest(tmp_path):
    manifest = tmp_path / "aligned_slice_manifest.csv"
    pd.DataFrame(
        [
            {"order": 1, "sample_id": "s1", "z_um": 0.0},
            {"order": 2, "sample_id": "s2", "z_um": 5.0},
        ]
    ).to_csv(manifest, index=False)
    return manifest


def _write_backend_matrices(path, offset: float) -> None:
    path.mkdir(parents=True)
    columns = ["top05|Structure 3", "top05|Structure 5"]
    pd.DataFrame(
        [[0.2 + offset, 0.8 + offset], [0.6 + offset, 0.1 + offset], [0.3 + offset, 0.4 + offset]],
        index=["GREM1", "EPCAM", "FAP"],
        columns=columns,
    ).to_csv(path / "gene_structure_fraction_inside_matrix.csv")
    pd.DataFrame(
        [[0.1 + offset, 0.7 + offset], [0.5 + offset, 0.1 + offset], [0.3 + offset, 0.5 + offset]],
        index=["GREM1", "EPCAM", "FAP"],
        columns=columns,
    ).to_csv(path / "gene_structure_overlap_fraction_matrix.csv")
    pd.DataFrame(
        [[50.0 + offset, -35.0 + offset], [-20.0 + offset, 120.0 + offset], [10.0 + offset, -12.0 + offset]],
        index=["GREM1", "EPCAM", "FAP"],
        columns=columns,
    ).to_csv(path / "gene_structure_signed_distance_matrix.csv")
