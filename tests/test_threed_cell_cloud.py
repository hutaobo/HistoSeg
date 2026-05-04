from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from histoseg.threed import (
    CODA_METHOD_CREDIT,
    CODA_METHOD_REFERENCE_DOI,
    CELL_CLOUD_ALIGNED_XY_OBSM_KEY,
    CELL_CLOUD_OBS_SLICE_KEY,
    CELL_CLOUD_OBSM_KEY,
    CELL_CLOUD_UNS_KEY,
    CellCloudProjectionConfig,
    CellCloudProjectionResult,
    CellCloudRenderConfig,
    CellCloudRenderResult,
    build_alignment_manifest,
    cell_cloud_cache_status,
    cell_cloud_dataframe_from_coordinates,
    hash_alignment_manifest,
    load_cell_alignment_transforms,
    project_cell_coordinates,
    render_cell_cloud_html,
    write_cell_cloud_cache,
)


def test_project_cell_coordinates_applies_hard_transform_in_physical_units(tmp_path):
    stack_root = _write_minimal_stack(tmp_path)
    obs = pd.DataFrame(
        {
            "sample_id": ["s1", "s2", "unknown"],
            "barcode": ["c1", "c2", "c3"],
            "x_centroid": [4.0, 4.0, 100.0],
            "y_centroid": [6.0, 6.0, 100.0],
            "leiden": ["0", "1", "2"],
        }
    )

    transforms = load_cell_alignment_transforms(
        stack_root,
        pixel_size_um=2.0,
    )
    coords, slice_order = project_cell_coordinates(
        obs,
        transforms,
        sample_column="sample_id",
        x_column="x_centroid",
        y_column="y_centroid",
        pixel_size_um=2.0,
    )

    np.testing.assert_allclose(coords[0], [4.0, 6.0, 0.0])
    np.testing.assert_allclose(coords[1], [24.0, 2.0, 5.0])
    assert slice_order.tolist() == [1, 2, -1]

    table = cell_cloud_dataframe_from_coordinates(
        obs,
        coords,
        slice_order,
        CellCloudProjectionConfig(
            h5ad="unused.h5ad",
            stack_root=stack_root,
            out_parquet=tmp_path / "cells.parquet",
            label_columns=("leiden",),
        ),
    )
    assert table["barcode"].tolist() == ["c1", "c2"]
    assert table["leiden"].tolist() == ["0", "1"]
    assert table.loc[1, "x_3d_um"] == pytest.approx(24.0)
    assert table.loc[1, "y_3d_um"] == pytest.approx(2.0)


def test_alignment_manifest_hash_changes_only_with_geometry_state(tmp_path):
    stack_root = _write_minimal_stack(tmp_path)

    manifest = build_alignment_manifest(stack_root, pixel_size_um=2.0)
    first_hash = hash_alignment_manifest(manifest)
    assert first_hash == hash_alignment_manifest(json.loads(json.dumps(manifest)))

    (stack_root / "preview_only.html").write_text("side effect", encoding="utf-8")
    assert hash_alignment_manifest(build_alignment_manifest(stack_root, pixel_size_um=2.0)) == first_hash

    hard_summary = stack_root / "pairwise_alignments" / "001_to_002_s2" / "hard_similarity_alignment.json"
    payload = json.loads(hard_summary.read_text(encoding="utf-8"))
    payload["union_iou_after_hard"] = 0.95
    hard_summary.write_text(json.dumps(payload), encoding="utf-8")
    assert hash_alignment_manifest(build_alignment_manifest(stack_root, pixel_size_um=2.0)) == first_hash

    payload["transform"]["translate_x"] = 11.0
    hard_summary.write_text(json.dumps(payload), encoding="utf-8")

    changed_hash = hash_alignment_manifest(build_alignment_manifest(stack_root, pixel_size_um=2.0))
    assert changed_hash != first_hash


def test_alignment_manifest_tracks_coda_backend_credit_and_parameters(tmp_path):
    stack_root = _write_minimal_stack(tmp_path)
    hard_summary = stack_root / "pairwise_alignments" / "001_to_002_s2" / "hard_similarity_alignment.json"
    payload = json.loads(hard_summary.read_text(encoding="utf-8"))
    payload.update(
        {
            "registration_backend": "coda-image",
            "method_credit": CODA_METHOD_CREDIT,
            "method_reference_doi": CODA_METHOD_REFERENCE_DOI,
            "coda_image": {
                "radon_rotation_degrees": -12.0,
                "radon_score": 0.9,
                "radon_angle_range": [0.0, 180.0],
                "radon_angle_step": 1.0,
                "phase_shift_y": -3.0,
                "phase_shift_x": 4.0,
                "phase_error": 0.1,
                "phase_difference": 0.0,
                "phase_upsample_factor": 1,
                "preprocessing": {
                    "input_source": "contour_union_raster_proxy",
                    "raster_size": 512,
                    "square_bounds": [0.0, 0.0, 100.0, 100.0],
                    "native_units_per_pixel": 0.2,
                    "fixed_positive_pixels": 10,
                    "moving_positive_pixels": 11,
                    "mask_padding_fraction": 0.05,
                },
            },
        }
    )
    hard_summary.write_text(json.dumps(payload), encoding="utf-8")

    manifest = build_alignment_manifest(stack_root, pixel_size_um=2.0)
    registration = manifest["slices"][1]["hard_registration"]
    first_hash = hash_alignment_manifest(manifest)

    assert registration["registration_backend"] == "coda-image"
    assert registration["method_credit"] == CODA_METHOD_CREDIT
    assert registration["method_reference_doi"] == CODA_METHOD_REFERENCE_DOI
    assert registration["coda_image"]["radon_angle_step"] == 1.0

    payload["coda_image"]["phase_error"] = 0.8
    hard_summary.write_text(json.dumps(payload), encoding="utf-8")
    assert hash_alignment_manifest(build_alignment_manifest(stack_root, pixel_size_um=2.0)) == first_hash

    payload["coda_image"]["radon_angle_step"] = 0.5
    hard_summary.write_text(json.dumps(payload), encoding="utf-8")
    assert hash_alignment_manifest(build_alignment_manifest(stack_root, pixel_size_um=2.0)) != first_hash


def test_write_cell_cloud_cache_uses_histoseg_keys_without_touching_spatial(tmp_path):
    adata = _MiniAnnData(
        pd.DataFrame(
            {"sample_id": ["s1", "s2"]},
            index=["c1", "c2"],
        )
    )
    coords = np.array([[1.0, 2.0, 0.0], [3.0, 4.0, 5.0]])
    slice_order = np.array([1, 2])
    provenance = {
        "alignment_hash": "abc123",
        "schema_version": 1,
    }

    assert cell_cloud_cache_status(adata, "abc123") == "missing"
    write_cell_cloud_cache(adata, coords, slice_order, provenance)

    assert cell_cloud_cache_status(adata, "abc123") == "valid"
    assert cell_cloud_cache_status(adata, "different") == "stale"
    np.testing.assert_allclose(adata.obsm[CELL_CLOUD_OBSM_KEY], coords)
    np.testing.assert_allclose(adata.obsm[CELL_CLOUD_ALIGNED_XY_OBSM_KEY], coords[:, :2])
    assert adata.obs[CELL_CLOUD_OBS_SLICE_KEY].tolist() == [1, 2]
    assert adata.uns[CELL_CLOUD_UNS_KEY]["alignment_hash"] == "abc123"
    assert "spatial" not in adata.obsm


def test_write_cell_cloud_cache_can_opt_into_scanpy_spatial(tmp_path):
    adata = _MiniAnnData(pd.DataFrame(index=["c1"]))
    coords = np.array([[1.0, 2.0, 0.0]])

    write_cell_cloud_cache(
        adata,
        coords,
        np.array([1]),
        {"alignment_hash": "abc123"},
        write_scanpy_spatial=True,
    )

    np.testing.assert_allclose(adata.obsm["spatial"], [[1.0, 2.0]])


def test_project_cells_cli_parses_cache_flags(monkeypatch, tmp_path, capsys):
    import histoseg.threed.cli as cli

    calls: list[CellCloudProjectionConfig] = []

    def fake_run(cfg):
        calls.append(cfg)
        return CellCloudProjectionResult(
            out_parquet=Path(cfg.out_parquet),
            cell_count=2,
            projected_cell_count=2,
            stack_root=Path(cfg.stack_root),
            alignment_hash="abc123",
            cache_status="ignored",
        )

    monkeypatch.setattr(cli, "run_cell_cloud_projection", fake_run)
    cli.main(
        [
            "project-cells",
            "--h5ad",
            "cells.h5ad",
            "--stack-root",
            str(tmp_path / "stack"),
            "--out-parquet",
            str(tmp_path / "cells.parquet"),
            "--label-columns",
            "leiden,cell_type",
            "sample_batch",
            "--pixel-size-um",
            "2.0",
            "--chunk-size",
            "42",
            "--ignore-cache",
            "--fail-on-stale-cache",
            "--write-cache",
            "--cache-h5ad",
            str(tmp_path / "cached.h5ad"),
            "--write-scanpy-spatial",
        ]
    )

    capsys.readouterr()
    assert len(calls) == 1
    cfg = calls[0]
    assert cfg.label_columns == ("leiden", "cell_type", "sample_batch")
    assert cfg.pixel_size_um == 2.0
    assert cfg.chunk_size == 42
    assert cfg.ignore_cache is True
    assert cfg.fail_on_stale_cache is True
    assert cfg.write_cache is True
    assert cfg.write_scanpy_spatial is True
    assert Path(cfg.cache_h5ad).name == "cached.h5ad"


def test_render_cell_cloud_html_writes_label_and_contour_traces(tmp_path, caplog):
    stack_root = tmp_path / "stack"
    stack_root.mkdir()
    pd.DataFrame(
        [
            {
                "structure": "Structure 1",
                "slice_order": 1,
                "polyline_id": 0,
                "point_index": 0,
                "x_um": 0.0,
                "y_um": 0.0,
                "z_um": 0.0,
            },
            {
                "structure": "Structure 1",
                "slice_order": 1,
                "polyline_id": 0,
                "point_index": 1,
                "x_um": 1.0,
                "y_um": 1.0,
                "z_um": 0.0,
            },
        ]
    ).to_csv(stack_root / "aligned_contour_3d_points.csv", index=False)

    aligned = tmp_path / "cells.parquet"
    pd.DataFrame(
        {
            "sample_id": ["s1", "s1", "s2", "s2"],
            "barcode": ["c1", "c2", "c3", "c4"],
            "leiden_1_0": ["0", "0", "1", "1"],
            "slice_order": [1, 1, 2, 2],
            "x_3d_um": [0.0, 1.0, 2.0, 3.0],
            "y_3d_um": [0.0, 1.0, 2.0, 3.0],
            "z_um": [0.0, 0.0, 5.0, 5.0],
        }
    ).to_parquet(aligned, index=False)

    caplog.set_level(logging.WARNING)
    result = render_cell_cloud_html(
        CellCloudRenderConfig(
            aligned_cells_parquet=aligned,
            stack_root=stack_root,
            out_html=tmp_path / "cells.html",
            max_points=3,
            performance_warning_threshold=3,
        )
    )

    assert result.total_cells == 4
    assert result.rendered_cells == 3
    assert result.label_count == 2
    assert result.contour_trace_count == 1
    html = result.out_html.read_text(encoding="utf-8")
    assert "Leiden 0" in html
    assert "Leiden 1" in html
    assert "Structure 1 contour" in html
    assert "Use --max-points to downsample" in caplog.text


def test_render_cell_cloud_cli_parses_existing_parquet(monkeypatch, tmp_path, capsys):
    import histoseg.threed.cli as cli

    calls: list[CellCloudRenderConfig] = []

    def fake_render(cfg):
        calls.append(cfg)
        return CellCloudRenderResult(
            out_html=Path(cfg.out_html),
            aligned_cells_parquet=Path(cfg.aligned_cells_parquet),
            stack_root=Path(cfg.stack_root),
            total_cells=4,
            rendered_cells=2,
            label_column=cfg.label_column,
            label_count=2,
            contour_trace_count=0,
            z_visual_scale=cfg.z_visual_scale,
        )

    monkeypatch.setattr(cli, "render_cell_cloud_html", fake_render)
    cli.main(
        [
            "render-cell-cloud",
            "--stack-root",
            str(tmp_path / "stack"),
            "--aligned-cells-parquet",
            str(tmp_path / "cells.parquet"),
            "--out-html",
            str(tmp_path / "cells.html"),
            "--label-column",
            "cell_type",
            "--max-points",
            "42",
            "--random-state",
            "7",
            "--z-visual-scale",
            "3.5",
            "--marker-size",
            "2.0",
            "--opacity",
            "0.4",
            "--no-contours",
            "--title",
            "Custom cell cloud",
            "--performance-warning-threshold",
            "123",
        ]
    )

    capsys.readouterr()
    assert len(calls) == 1
    cfg = calls[0]
    assert cfg.label_column == "cell_type"
    assert cfg.max_points == 42
    assert cfg.random_state == 7
    assert cfg.z_visual_scale == 3.5
    assert cfg.marker_size == 2.0
    assert cfg.opacity == 0.4
    assert cfg.include_contours is False
    assert cfg.title == "Custom cell cloud"
    assert cfg.performance_warning_threshold == 123


def test_render_cell_cloud_cli_projects_h5ad_before_rendering(monkeypatch, tmp_path, capsys):
    import histoseg.threed.cli as cli

    projection_calls: list[CellCloudProjectionConfig] = []
    render_calls: list[CellCloudRenderConfig] = []

    def fake_project(cfg):
        projection_calls.append(cfg)
        return CellCloudProjectionResult(
            out_parquet=Path(cfg.out_parquet),
            cell_count=4,
            projected_cell_count=4,
            stack_root=Path(cfg.stack_root),
            alignment_hash="abc123",
            cache_status="missing",
        )

    def fake_render(cfg):
        render_calls.append(cfg)
        return CellCloudRenderResult(
            out_html=Path(cfg.out_html),
            aligned_cells_parquet=Path(cfg.aligned_cells_parquet),
            stack_root=Path(cfg.stack_root),
            total_cells=4,
            rendered_cells=4,
            label_column=cfg.label_column,
            label_count=2,
            contour_trace_count=1,
            z_visual_scale=cfg.z_visual_scale,
        )

    monkeypatch.setattr(cli, "run_cell_cloud_projection", fake_project)
    monkeypatch.setattr(cli, "render_cell_cloud_html", fake_render)
    cli.main(
        [
            "render-cell-cloud",
            "--stack-root",
            str(tmp_path / "stack"),
            "--h5ad",
            str(tmp_path / "cells.h5ad"),
            "--out-parquet",
            str(tmp_path / "cells.parquet"),
            "--out-html",
            str(tmp_path / "cells.html"),
            "--label-column",
            "leiden_1_0",
            "--label-columns",
            "cell_type,batch",
            "--ignore-cache",
            "--write-cache",
        ]
    )

    capsys.readouterr()
    assert len(projection_calls) == 1
    assert projection_calls[0].label_columns == ("leiden_1_0", "cell_type", "batch")
    assert projection_calls[0].ignore_cache is True
    assert projection_calls[0].write_cache is True
    assert len(render_calls) == 1
    assert Path(render_calls[0].aligned_cells_parquet).name == "cells.parquet"
    assert render_calls[0].label_column == "leiden_1_0"


class _MiniAnnData:
    def __init__(self, obs: pd.DataFrame):
        self.obs = obs.copy()
        self.obsm = {}
        self.uns = {}

    @property
    def n_obs(self) -> int:
        return len(self.obs)


def _write_minimal_stack(root: Path) -> Path:
    stack_root = root / "stack"
    stack_root.mkdir()
    pd.DataFrame(
        [
            {"order": 1, "sample_id": "s1", "z_um": 0.0, "alignment_role": "reference"},
            {"order": 2, "sample_id": "s2", "z_um": 5.0, "alignment_role": "moving"},
        ]
    ).to_csv(stack_root / "aligned_slice_manifest.csv", index=False)
    pd.DataFrame(
        [
            {
                "moving_order": 2,
                "moving_sample_id": "s2",
                "soft_accepted": False,
                "soft_summary_json": "",
            }
        ]
    ).to_csv(stack_root / "pairwise_alignment_metrics.csv", index=False)

    pair_dir = stack_root / "pairwise_alignments" / "001_to_002_s2"
    pair_dir.mkdir(parents=True)
    hard_summary = {
        "transform": {
            "origin_x": 0.0,
            "origin_y": 0.0,
            "rotation_degrees": 0.0,
            "scale": 1.0,
            "translate_x": 10.0,
            "translate_y": -2.0,
        },
        "union_iou_before_hard": 0.1,
        "union_iou_after_hard": 0.9,
        "hard_alignment_accepted": True,
    }
    (pair_dir / "hard_similarity_alignment.json").write_text(
        json.dumps(hard_summary),
        encoding="utf-8",
    )
    return stack_root
