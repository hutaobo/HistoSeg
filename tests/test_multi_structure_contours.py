from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from histoseg import (
    MultiStructureContourConfig,
    MultiStructureSpec,
    run_multi_structure_contours,
)


def _write_synthetic_inputs(tmp_path: Path) -> tuple[Path, Path]:
    cells = pd.DataFrame(
        {
            "cell_id": [f"cell_{idx}" for idx in range(10)],
            "x_centroid": [5.0, 8.0, 11.0, 14.0, 17.0, 65.0, 68.0, 71.0, 74.0, 77.0],
            "y_centroid": [10.0, 14.0, 18.0, 22.0, 26.0, 10.0, 14.0, 18.0, 22.0, 26.0],
        }
    )
    clusters = pd.DataFrame(
        {
            "Barcode": cells["cell_id"],
            "Cluster": [1, 1, 1, 1, 1, 2, 2, 2, 2, 2],
        }
    )

    cells_path = tmp_path / "cells.parquet"
    clusters_path = tmp_path / "clusters.csv"
    cells.to_parquet(cells_path, index=False)
    clusters.to_csv(clusters_path, index=False)
    return clusters_path, cells_path


def test_multi_structure_contours_create_non_overlapping_exports(tmp_path):
    clusters_path, cells_path = _write_synthetic_inputs(tmp_path)
    out_dir = tmp_path / "outputs"

    result = run_multi_structure_contours(
        MultiStructureContourConfig(
            clusters_csv=clusters_path,
            cells_parquet=cells_path,
            out_dir=out_dir,
            structures=[
                MultiStructureSpec(
                    structure_name="Left niche",
                    cluster_ids=[1],
                    structure_color="#6EF0D4",
                    structure_id=1,
                ),
                MultiStructureSpec(
                    structure_name="Right niche",
                    cluster_ids=[2],
                    structure_color="#78B9FF",
                    structure_id=2,
                ),
            ],
            bins_x=60,
            bins_y=40,
            gaussian_sigma=1.2,
            support_quantile=0.1,
            tissue_quantile=0.01,
            min_dominance=0.0,
            min_cells=2,
            min_component_pixels=1,
        )
    )

    assert result.geojson.exists()
    assert result.csv.exists()
    assert result.summary.exists()
    assert result.partition_table.exists()
    assert result.structure_count_csv.exists()
    assert result.metrics_json.exists()
    assert result.preview_png is not None and result.preview_png.exists()

    if result.partition_table.suffix == ".parquet":
        partition_df = pd.read_parquet(result.partition_table)
    else:
        partition_df = pd.read_csv(result.partition_table)
    assert set(partition_df["isoline_structure_id"].astype(int)) == {1, 2}
    assert partition_df["isoline_structure_name"].isin({"Left niche", "Right niche"}).all()

    geojson_payload = json.loads(result.geojson.read_text(encoding="utf-8"))
    assert geojson_payload["type"] == "FeatureCollection"
    assert len(geojson_payload["features"]) >= 2

    summary_df = pd.read_csv(result.summary)
    assert len(summary_df) == len(geojson_payload["features"])
    assert set(summary_df["AssignedStructure"]) == {"Left niche", "Right niche"}

    structure_count_df = pd.read_csv(result.structure_count_csv)
    assert set(structure_count_df["structure_name"]) == {"Left niche", "Right niche"}
    assert np.all(structure_count_df["assigned_cell_count"].to_numpy(dtype=int) > 0)
    assert np.all(structure_count_df["contour_count"].to_numpy(dtype=int) > 0)


def test_xenium_explorer_feature_schema_matches_expected_fields(tmp_path):
    clusters_path, cells_path = _write_synthetic_inputs(tmp_path)
    result = run_multi_structure_contours(
        MultiStructureContourConfig(
            clusters_csv=clusters_path,
            cells_parquet=cells_path,
            out_dir=tmp_path / "schema_outputs",
            structures=[
                {"structure_name": "Left niche", "cluster_ids": [1], "structure_color": "#6EF0D4"},
                {"structure_name": "Right niche", "cluster_ids": [2], "structure_color": "#78B9FF"},
            ],
            bins_x=60,
            bins_y=40,
            gaussian_sigma=1.2,
            support_quantile=0.1,
            tissue_quantile=0.01,
            min_dominance=0.0,
            min_cells=2,
            min_component_pixels=1,
        )
    )

    payload = json.loads(result.geojson.read_text(encoding="utf-8"))
    feature = payload["features"][0]
    props = feature["properties"]

    assert feature["type"] == "Feature"
    assert feature["geometry"]["type"] == "Polygon"
    assert props["objectType"] == "annotation"
    assert isinstance(props["name"], str) and props["name"]
    assert props["classification"]["name"] == props["name"]
    assert len(props["classification"]["color"]) == 3
    assert "assigned_structure" in props
    assert "component_index" in props
    assert "polygon_index" in props

    summary_df = pd.read_csv(result.summary)
    assert summary_df["Selection"].is_unique
    assert set(summary_df.columns) == {
        "Selection",
        "StructureID",
        "AssignedStructure",
        "ComponentIndex",
        "PolygonIndex",
        "VertexCount",
    }
