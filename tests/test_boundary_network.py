from __future__ import annotations

import json

import pandas as pd
import pytest

from histoseg.contour import (
    BoundaryNetworkConfig,
    build_group_boundary_graph,
    normalize_group_boundary_overlap,
    run_group_boundary_network,
)
from histoseg.contour.cli import main


def _legacy_edges() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "group_a": ["Structure 6", "Structure 2", "Structure 5"],
            "group_b": ["Structure 7", "Structure 6", "Structure 7"],
            "n_boundary_pairs": [112, 38, 10],
            "total_shared_boundary": [51297.10488215108, 21518.52879491277, 693.2867248572746],
            "mean_shared_boundary": [458.0098650192061, 566.277073550336, 69.32867248572747],
        }
    )


def test_normalize_group_boundary_overlap_accepts_legacy_columns():
    normalized = normalize_group_boundary_overlap(_legacy_edges())

    assert list(normalized.columns) == [
        "group_a",
        "group_b",
        "n_boundary_pairs",
        "shared_boundary_length_um",
        "mean_shared_boundary_length_um",
    ]
    assert normalized.iloc[0]["group_a"] == "Structure 6"
    assert normalized.iloc[0]["group_b"] == "Structure 7"
    assert normalized.iloc[0]["shared_boundary_length_um"] == pytest.approx(
        51297.10488215108
    )


def test_normalize_group_boundary_overlap_accepts_topology_columns_and_derives_mean():
    normalized = normalize_group_boundary_overlap(
        pd.DataFrame(
            {
                "group_a": ["B"],
                "group_b": ["A"],
                "n_boundary_pairs": [5],
                "shared_boundary_length_um": [20.0],
            }
        )
    )

    row = normalized.iloc[0]
    assert row["group_a"] == "A"
    assert row["group_b"] == "B"
    assert row["mean_shared_boundary_length_um"] == pytest.approx(4.0)


def test_build_group_boundary_graph_derives_node_attributes():
    graph = build_group_boundary_graph(_legacy_edges())

    assert graph.number_of_nodes() == 4
    assert graph.number_of_edges() == 3
    assert graph.nodes["Structure 6"]["degree"] == 2
    assert graph.nodes["Structure 6"]["total_boundary_pairs"] == 150
    assert graph.nodes["Structure 6"]["weighted_degree_shared_boundary_um"] == pytest.approx(
        72815.63367706384
    )


def test_run_group_boundary_network_filters_drop_structures_and_writes_outputs(tmp_path):
    boundary_csv = tmp_path / "group_boundary_overlap_filtered.csv"
    _legacy_edges().to_csv(boundary_csv, index=False)

    result = run_group_boundary_network(
        BoundaryNetworkConfig(
            boundary_csv=boundary_csv,
            out_dir=tmp_path / "network",
            drop_structures=("6",),
        )
    )

    filtered = pd.read_csv(result.filtered_edges_csv)
    assert set(filtered["group_a"]).isdisjoint({"Structure 6"})
    assert set(filtered["group_b"]).isdisjoint({"Structure 6"})
    assert result.preview_png is not None
    assert result.preview_png.exists()
    assert result.nodes_csv.exists()
    metrics = json.loads(result.metrics_json.read_text(encoding="utf-8"))
    assert metrics["n_edges"] == 1
    assert metrics["n_nodes"] == 2


def test_run_group_boundary_network_writes_empty_graph_outputs(tmp_path):
    boundary_csv = tmp_path / "group_boundary_overlap_filtered.csv"
    _legacy_edges().to_csv(boundary_csv, index=False)

    result = run_group_boundary_network(
        BoundaryNetworkConfig(
            boundary_csv=boundary_csv,
            out_dir=tmp_path / "empty_network",
            drop_structures=("2", "5", "6", "7"),
        )
    )

    assert result.preview_png is not None
    assert result.preview_png.exists()
    assert pd.read_csv(result.filtered_edges_csv).empty
    assert pd.read_csv(result.nodes_csv).empty
    metrics = json.loads(result.metrics_json.read_text(encoding="utf-8"))
    assert metrics["n_edges"] == 0
    assert metrics["n_nodes"] == 0


def test_boundary_network_cli_writes_expected_outputs(tmp_path):
    boundary_csv = tmp_path / "group_boundary_overlap_filtered.csv"
    _legacy_edges().to_csv(boundary_csv, index=False)
    out_dir = tmp_path / "cli_network"

    main(
        [
            "boundary-network",
            "--boundary-csv",
            str(boundary_csv),
            "--out-dir",
            str(out_dir),
        ]
    )

    assert (out_dir / "group_boundary_network.png").exists()
    assert (out_dir / "group_boundary_overlap_filtered.csv").exists()
    assert (out_dir / "group_boundary_network_nodes.csv").exists()
    assert (out_dir / "group_boundary_network_metrics.json").exists()
