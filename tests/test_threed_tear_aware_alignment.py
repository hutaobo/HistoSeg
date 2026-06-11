from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from shapely.geometry import Polygon, mapping, shape
from shapely.ops import unary_union

from histoseg.threed import ContourTearClosureConfig, run_contour_tear_closure
from histoseg.threed.cli import main


def test_contour_tear_closure_moves_shifted_block_without_tps(tmp_path):
    fixed = tmp_path / "fixed.geojson"
    moving = tmp_path / "moving.geojson"
    cells = tmp_path / "moving_cells.csv"

    fixed_geoms = [
        Polygon([(0, 0), (50, 0), (50, 50), (0, 50)]),
        Polygon([(90, 0), (140, 0), (140, 50), (90, 50)]),
    ]
    moving_geoms = [
        fixed_geoms[0],
        Polygon([(150, 0), (200, 0), (200, 50), (150, 50)]),
    ]
    _write_geojson(fixed, fixed_geoms, ["Epithelial", "Epithelial"])
    _write_geojson(moving, moving_geoms, ["Epithelial", "Epithelial"])
    pd.DataFrame(
        {
            "cell_id": ["left", "right"],
            "cluster": ["Epithelial", "Epithelial"],
            "x_aligned_um": [25.0, 175.0],
            "y_aligned_um": [25.0, 25.0],
        }
    ).to_csv(cells, index=False)

    result = run_contour_tear_closure(
        ContourTearClosureConfig(
            fixed_geojson=fixed,
            moving_aligned_geojson=moving,
            moving_cells_csv=cells,
            out_dir=tmp_path / "closed",
            centroid_search_radius_um=120.0,
            influence_radius_um=80.0,
            max_neighbors=1,
            max_displacement_um=80.0,
            save_preview_png=False,
            overwrite=True,
        )
    )

    closed_cells = pd.read_csv(result.closed_cells_csv)
    right = closed_cells.loc[closed_cells["cell_id"] == "right"].iloc[0]
    left = closed_cells.loc[closed_cells["cell_id"] == "left"].iloc[0]
    assert right["x_contour_tear_closure_um"] == 115.0
    assert right["dx_contour_tear_closure_um"] == -60.0
    assert left["x_contour_tear_closure_um"] == 25.0

    before_payload = json.loads(moving.read_text(encoding="utf-8"))
    after_payload = json.loads(result.closed_contours_geojson.read_text(encoding="utf-8"))
    fixed_union = unary_union(fixed_geoms)
    before_union = unary_union([shape(feature["geometry"]) for feature in before_payload["features"]])
    after_union = unary_union([shape(feature["geometry"]) for feature in after_payload["features"]])
    assert _iou(fixed_union, after_union) > _iou(fixed_union, before_union)

    summary = json.loads(result.summary_json.read_text(encoding="utf-8"))
    assert summary["method"] == "contour_tear_closure"
    assert summary["method_family"] == "contour_guided_local_translation_no_tps"
    assert summary["contours"]["reciprocal_correspondence_pairs"] == 2
    assert summary["landmarks"]["used_local_translation_anchor_count"] == 2


def test_contour_tear_closure_cli_smoke(tmp_path):
    fixed = tmp_path / "fixed.geojson"
    moving = tmp_path / "moving.geojson"
    cells = tmp_path / "moving_cells.csv"
    square = Polygon([(0, 0), (40, 0), (40, 40), (0, 40)])
    shifted = Polygon([(20, 0), (60, 0), (60, 40), (20, 40)])
    _write_geojson(fixed, [square], ["A"])
    _write_geojson(moving, [shifted], ["A"])
    pd.DataFrame(
        {
            "cell_id": ["cell-1"],
            "cluster": ["A"],
            "x_aligned_um": [40.0],
            "y_aligned_um": [20.0],
        }
    ).to_csv(cells, index=False)

    main(
        [
            "close-contour-tear",
            "--fixed-geojson",
            str(fixed),
            "--moving-aligned-geojson",
            str(moving),
            "--moving-cells-csv",
            str(cells),
            "--out-dir",
            str(tmp_path / "cli"),
            "--centroid-search-radius-um",
            "80",
            "--influence-radius-um",
            "80",
            "--max-neighbors",
            "1",
            "--max-displacement-um",
            "50",
            "--no-preview",
            "--overwrite",
        ]
    )

    assert (tmp_path / "cli" / "moving_contour_tear_closed_cells.csv").exists()
    assert (tmp_path / "cli" / "moving_contour_tear_closed_contours.geojson").exists()
    assert (tmp_path / "cli" / "contour_tear_closure_landmarks.csv").exists()
    assert (tmp_path / "cli" / "contour_tear_closure_summary.json").exists()


def _write_geojson(path: Path, geometries: list[Polygon], labels: list[str]) -> None:
    payload = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "properties": {
                    "assigned_structure": label,
                    "structure_id": index + 1,
                    "name": f"{label} #{index + 1}",
                },
                "geometry": mapping(geom),
            }
            for index, (geom, label) in enumerate(zip(geometries, labels))
        ],
    }
    path.write_text(json.dumps(payload), encoding="utf-8")


def _iou(a, b) -> float:
    return float(a.intersection(b).area / a.union(b).area)
