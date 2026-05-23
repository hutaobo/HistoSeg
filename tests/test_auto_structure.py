from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from shapely.geometry import Polygon, mapping

from histoseg.contour import AutoStructureDiscoveryConfig, discover_auto_structures, resolve_xenium_output_folder
from histoseg.threed import LabelFreeBeforeAfterFigureConfig, render_label_free_before_after_panel


def _write_graphclust_case(root: Path, *, cluster_count: int = 4, cells_per_cluster: int = 12) -> tuple[Path, Path]:
    root.mkdir(parents=True, exist_ok=True)
    rows = []
    clusters = []
    centers = [(120 * (idx % 6), 120 * (idx // 6)) for idx in range(max(cluster_count, 6))]
    cell_id = 1
    for cluster_id in range(1, cluster_count + 1):
        cx, cy = centers[cluster_id - 1]
        for offset in range(cells_per_cluster):
            rows.append(
                {
                    "cell_id": str(cell_id),
                    "x_centroid": cx + (offset % 4) * 2.0,
                    "y_centroid": cy + (offset // 4) * 2.0,
                }
            )
            clusters.append({"Barcode": str(cell_id), "Cluster": cluster_id})
            cell_id += 1
    cells_path = root / "cells.parquet"
    clusters_path = root / "analysis" / "clustering" / "gene_expression_graphclust" / "clusters.csv"
    clusters_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_parquet(cells_path, index=False)
    pd.DataFrame(clusters).to_csv(clusters_path, index=False)
    return clusters_path, cells_path


def test_auto_structure_discovers_valid_structures_json(tmp_path):
    clusters_path, cells_path = _write_graphclust_case(tmp_path / "outs")

    result = discover_auto_structures(
        AutoStructureDiscoveryConfig(
            clusters_csv=clusters_path,
            cells_parquet=cells_path,
            out_dir=tmp_path / "auto",
            cluster_count="auto",
            min_structure_count=3,
            max_structure_count=4,
            min_structure_cell_fraction=0.05,
            use_cophenetic=False,
        )
    )

    assert result.structures_json.exists()
    structures = json.loads(result.structures_json.read_text(encoding="utf-8"))
    assert 3 <= len(structures) <= 4
    assert all(item["structure_name"].startswith("Structure ") for item in structures)
    assert all(item["cluster_ids"] for item in structures)
    assert result.cluster_structure_csv.exists()
    table = pd.read_csv(result.cluster_structure_csv)
    assert set(table["cluster"].astype(str)) == {"1", "2", "3", "4"}


def test_auto_structure_leaf_balanced_prevents_giant_root_groups(tmp_path):
    clusters_path, cells_path = _write_graphclust_case(
        tmp_path / "outs",
        cluster_count=12,
        cells_per_cluster=8,
    )

    result = discover_auto_structures(
        AutoStructureDiscoveryConfig(
            clusters_csv=clusters_path,
            cells_parquet=cells_path,
            out_dir=tmp_path / "leaf_balanced",
            cluster_count="leaf-balanced",
            max_leaf_clusters_per_structure=5,
            min_leaf_clusters_per_structure=2,
            min_structure_cell_fraction=0.0,
            use_cophenetic=False,
        )
    )

    table = pd.read_csv(result.cluster_structure_csv)
    leaves_per_structure = table.groupby("structure_name")["cluster"].nunique()
    assert leaves_per_structure.max() <= 5
    assert leaves_per_structure.min() >= 2
    assert len(leaves_per_structure) >= 3
    summary = json.loads(result.summary_json.read_text(encoding="utf-8"))
    assert summary["cluster_count_mode"] == "leaf_balanced"


def test_auto_structure_ignores_unassigned_cluster_labels(tmp_path):
    clusters_path, cells_path = _write_graphclust_case(
        tmp_path / "outs",
        cluster_count=6,
        cells_per_cluster=8,
    )
    cells = pd.read_parquet(cells_path)
    clusters = pd.read_csv(clusters_path)
    extra_cells = pd.DataFrame(
        [
            {
                "cell_id": f"u{idx}",
                "x_centroid": 1000.0 + idx,
                "y_centroid": 1000.0 + idx,
            }
            for idx in range(6)
        ]
    )
    extra_clusters = pd.DataFrame(
        {"Barcode": extra_cells["cell_id"], "Cluster": ["Unassigned", "nan", "None", "null", "", np.nan]}
    )
    pd.concat([cells, extra_cells], ignore_index=True).to_parquet(cells_path, index=False)
    pd.concat([clusters, extra_clusters], ignore_index=True).to_csv(clusters_path, index=False)

    result = discover_auto_structures(
        AutoStructureDiscoveryConfig(
            clusters_csv=clusters_path,
            cells_parquet=cells_path,
            out_dir=tmp_path / "auto_no_unassigned",
            cluster_count="leaf-balanced",
            max_leaf_clusters_per_structure=3,
            min_leaf_clusters_per_structure=1,
            min_structure_cell_fraction=0.0,
            use_cophenetic=False,
        )
    )

    table = pd.read_csv(result.cluster_structure_csv)
    assert set(table["cluster"].astype(str)) == {"1", "2", "3", "4", "5", "6"}
    structures = json.loads(result.structures_json.read_text(encoding="utf-8"))
    assert all("Unassigned" not in item["cluster_ids"] for item in structures)


def test_resolve_xenium_output_folder_supports_parent_with_nested_outs(tmp_path):
    parent = tmp_path / "case"
    nested = parent / "case_outs"
    _write_graphclust_case(nested, cluster_count=2)

    assert resolve_xenium_output_folder(parent) == nested.resolve()
    assert resolve_xenium_output_folder(nested) == nested.resolve()


def test_auto_structure_cli_dry_run_does_not_write_outputs(tmp_path):
    outs = tmp_path / "case_outs"
    _write_graphclust_case(outs, cluster_count=3)
    out_dir = tmp_path / "dry"

    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "histoseg.contour.cli",
            "auto-structure",
            "--xenium-output",
            str(outs),
            "--out-dir",
            str(out_dir),
            "--dry-run",
            "--no-cophenetic",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(completed.stdout)
    assert payload["dry_run"] is True
    assert not out_dir.exists()
    assert payload["selected_structure_count"] >= 2


def test_label_free_before_after_renderer_writes_clean_panel(tmp_path):
    fixed = tmp_path / "fixed.geojson"
    moving = tmp_path / "moving.geojson"
    aligned = tmp_path / "aligned.geojson"
    anchors = tmp_path / "anchors.csv"
    _write_geojson(fixed, [Polygon([(0, 0), (20, 0), (20, 20), (0, 20)])])
    _write_geojson(moving, [Polygon([(40, 5), (60, 5), (60, 25), (40, 25)])])
    _write_geojson(aligned, [Polygon([(2, 1), (22, 1), (22, 21), (2, 21)])])
    pd.DataFrame(
        {
            "fixed_centroid_x": [10.0],
            "fixed_centroid_y": [10.0],
            "moving_centroid_x": [50.0],
            "moving_centroid_y": [15.0],
            "aligned_moving_centroid_x": [12.0],
            "aligned_moving_centroid_y": [11.0],
            "used_for_transform": [True],
        }
    ).to_csv(anchors, index=False)

    result = render_label_free_before_after_panel(
        LabelFreeBeforeAfterFigureConfig(
            fixed_geojson=fixed,
            moving_geojson=moving,
            aligned_geojson=aligned,
            anchors_csv=anchors,
            out_png=tmp_path / "figure.png",
            show_titles=False,
        )
    )

    assert result.out_png.exists()
    assert result.out_png.stat().st_size > 0


def _write_geojson(path: Path, geoms: list[Polygon]) -> None:
    payload = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "geometry": mapping(geom),
                "properties": {"assigned_structure": "Structure 1"},
            }
            for geom in geoms
        ],
    }
    path.write_text(json.dumps(payload), encoding="utf-8")
