from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from histoseg.threed import (
    CELL_CLOUD_ALIGNED_XY_OBSM_KEY,
    CELL_CLOUD_OBS_SLICE_KEY,
    CELL_CLOUD_OBSM_KEY,
    CELL_CLOUD_UNS_KEY,
    CellCloudProjectionConfig,
    CellCloudProjectionResult,
    build_alignment_manifest,
    cell_cloud_cache_status,
    cell_cloud_dataframe_from_coordinates,
    hash_alignment_manifest,
    load_cell_alignment_transforms,
    project_cell_coordinates,
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
