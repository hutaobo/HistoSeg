from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest
from shapely.geometry import box

from histoseg.contour import (
    ContourAdjacencyConfig,
    build_contour_adjacency_edges,
    run_contour_adjacency,
    summarize_contour_topology,
)
from histoseg.contour.cli import main


def _adjacency_contours() -> pd.DataFrame:
    return pd.DataFrame.from_records(
        [
            {"contour_id": "a1", "structure": "A", "geometry": box(0, 0, 10, 10)},
            {"contour_id": "b1", "structure": "B", "geometry": box(10, 0, 20, 10)},
            {"contour_id": "a2", "structure": "A", "geometry": box(0, 10, 10, 20)},
            {"contour_id": "outer", "structure": "Outer", "geometry": box(30, 0, 50, 20)},
            {"contour_id": "inner", "structure": "Inner", "geometry": box(35, 5, 45, 15)},
        ]
    )


def test_contour_adjacency_edges_include_boundary_and_enclosure_but_not_self_pairs():
    contours = _adjacency_contours()
    topology = summarize_contour_topology(
        contours,
        groupby="structure",
        boundary_tolerance=0.0,
    )

    edges = build_contour_adjacency_edges(
        topology,
        groups=["A", "B", "Inner", "Outer"],
        groupby="structure",
    )

    assert len(edges) == 6
    assert not (edges["group_a"] == edges["group_b"]).any()

    ab = edges.loc[(edges["group_a"] == "A") & (edges["group_b"] == "B")].iloc[0]
    assert ab["boundary_pair_count"] == 1
    assert ab["boundary_shared_length_um"] == pytest.approx(10.0)
    assert ab["enclosure_pair_count"] == 0
    assert ab["total_adjacency_length_um"] == pytest.approx(10.0)

    nested = edges.loc[
        (edges["group_a"] == "Inner") & (edges["group_b"] == "Outer")
    ].iloc[0]
    assert nested["boundary_pair_count"] == 0
    assert nested["enclosure_pair_count"] == 1
    assert nested["enclosure_inner_boundary_length_um"] == pytest.approx(40.0)
    assert nested["total_adjacency_length_um"] == pytest.approx(40.0)


def test_run_contour_adjacency_writes_matrix_and_plots(tmp_path):
    result = run_contour_adjacency(
        ContourAdjacencyConfig(
            contours=_adjacency_contours(),
            groupby="structure",
            out_dir=tmp_path / "adjacency",
            boundary_tolerance=0.0,
        )
    )

    edges = pd.read_csv(result.edges_csv)
    matrix = pd.read_csv(result.matrix_csv, index_col=0)
    metrics = json.loads(result.metrics_json.read_text(encoding="utf-8"))

    assert len(edges) == 6
    assert matrix.loc["A", "B"] == pytest.approx(10.0)
    assert matrix.loc["B", "A"] == pytest.approx(10.0)
    assert matrix.loc["Inner", "Outer"] == pytest.approx(40.0)
    assert np.isnan(matrix.loc["A", "A"])
    assert result.network_png is not None and result.network_png.exists()
    assert result.heatmap_png is not None and result.heatmap_png.exists()
    assert result.combined_png is not None and result.combined_png.exists()
    assert metrics["boundary_pair_count"] == 1
    assert metrics["enclosure_pair_count"] == 1


def test_contour_adjacency_cli_accepts_wkt_contours_csv(tmp_path):
    contours = _adjacency_contours().copy()
    contours["geometry"] = contours["geometry"].map(lambda geometry: geometry.wkt)
    contours_csv = tmp_path / "contours.csv"
    contours.to_csv(contours_csv, index=False)
    out_dir = tmp_path / "cli_adjacency"

    main(
        [
            "adjacency",
            "--contours-csv",
            str(contours_csv),
            "--groupby",
            "structure",
            "--out-dir",
            str(out_dir),
            "--boundary-tolerance",
            "0",
        ]
    )

    assert (out_dir / "contour_adjacency_edges.csv").exists()
    assert (out_dir / "contour_adjacency_matrix.csv").exists()
    assert (out_dir / "contour_adjacency_network.png").exists()
    assert (out_dir / "contour_adjacency_heatmap.png").exists()
    assert (out_dir / "contour_adjacency_overview.png").exists()
