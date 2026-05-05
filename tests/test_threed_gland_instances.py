from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from shapely.geometry import Polygon, mapping

from histoseg.threed import (
    GlandInstanceDetectionResult,
    GlandInstanceSegmentationConfig,
    segment_gland_instances,
)
from histoseg.threed.cli import main as threed_cli_main


def test_two_lumens_split_one_semantic_component(tmp_path: Path):
    stack_root, cells_path = _write_segmentation_fixture(
        tmp_path,
        semantic_polygon=_box(0, 0, 220, 110),
        cells=_ring_cells([(60, 55), (160, 55)], radius=24, n=28),
    )

    result = segment_gland_instances(
        _seg_cfg(stack_root, cells_path, tmp_path / "out")
    )

    qc = pd.read_csv(result.slice_qc_csv)
    assert result.total_gland_instances == 2
    assert qc["slice_instance_id"].nunique() == 2
    assert not qc["flag_no_lumen_seed"].any()


def test_adjacent_glands_do_not_merge(tmp_path: Path):
    stack_root, cells_path = _write_segmentation_fixture(
        tmp_path,
        semantic_polygon=_box(0, 0, 170, 100),
        cells=_ring_cells([(62, 50), (108, 50)], radius=18, n=24),
    )

    result = segment_gland_instances(
        _seg_cfg(stack_root, cells_path, tmp_path / "out", lumen_min_area_um2=25.0)
    )

    qc = pd.read_csv(result.slice_qc_csv)
    assert len(qc) == 2
    assert qc["centroid_x_um"].max() - qc["centroid_x_um"].min() > 20


def test_no_lumen_epithelial_fragment(tmp_path: Path):
    filled_cells = pd.DataFrame(
        [
            {
                "sample_id": "slice_1",
                "barcode": f"cell_{i}_{j}",
                "slice_order": 1,
                "x_3d_um": float(20 + i * 8),
                "y_3d_um": float(20 + j * 8),
                "EPCAM": 2.0,
            }
            for i in range(8)
            for j in range(8)
        ]
    )
    stack_root, cells_path = _write_segmentation_fixture(
        tmp_path,
        semantic_polygon=_box(0, 0, 100, 100),
        cells=filled_cells,
    )

    result = segment_gland_instances(
        _seg_cfg(stack_root, cells_path, tmp_path / "out")
    )

    qc = pd.read_csv(result.slice_qc_csv)
    assert len(qc) == 1
    assert bool(qc.loc[0, "flag_no_lumen_seed"])
    assert "no_lumen_seed" in str(qc.loc[0, "qc_flags"])


def test_qc_flags_additive(tmp_path: Path):
    stack_root, cells_path = _write_segmentation_fixture(
        tmp_path,
        semantic_polygon=_box(0, 0, 35, 35),
        cells=pd.DataFrame(
            [
                {
                    "sample_id": "slice_1",
                    "barcode": "cell_1",
                    "slice_order": 1,
                    "x_3d_um": 15.0,
                    "y_3d_um": 15.0,
                    "EPCAM": 1.0,
                }
            ]
        ),
    )

    result = segment_gland_instances(
        _seg_cfg(
            stack_root,
            cells_path,
            tmp_path / "out",
            min_gland_area_um2=5000.0,
            min_ring_support_score=1.01,
        )
    )

    qc = pd.read_csv(result.slice_qc_csv)
    flags = set(str(qc.loc[0, "qc_flags"]).split(";"))
    assert {"no_lumen_seed", "weak_ring", "small_fragment"}.issubset(flags)
    assert bool(qc.loc[0, "flag_no_lumen_seed"])
    assert bool(qc.loc[0, "flag_weak_ring"])
    assert bool(qc.loc[0, "flag_small_fragment"])


def test_detect_gland_instances_cli_smoke(tmp_path: Path):
    stack_root, cells_path = _write_segmentation_fixture(
        tmp_path,
        semantic_polygon=_box(0, 0, 220, 110),
        cells=_ring_cells([(60, 55), (160, 55)], radius=24, n=28),
    )
    out_dir = tmp_path / "cli_out"

    threed_cli_main(
        [
            "detect-gland-instances",
            "--stack-root",
            str(stack_root),
            "--aligned-cells-parquet",
            str(cells_path),
            "--out-dir",
            str(out_dir),
            "--pixel-size-um",
            "1",
            "--raster-pixel-size-um",
            "2",
            "--lumen-min-area-um2",
            "50",
            "--cell-density-sigma-um",
            "5",
            "--support-closing-um",
            "8",
            "--max-gland-pages",
            "2",
        ]
    )

    assert (out_dir / "slice_gland_instances.geojson").exists()
    assert (out_dir / "slice_gland_instances_qc.csv").exists()
    assert (out_dir / "gland_instance_tracks.csv").exists()
    assert (out_dir / "gland_instance_atlas.html").exists()


def test_detect_gland_instances_cli_wires_many_to_many(monkeypatch, tmp_path: Path, capsys):
    import histoseg.threed.cli as cli

    calls = {}

    def fake_detection(*, segmentation_config, tracking_config, max_gland_pages, z_visual_scale):
        calls["segmentation"] = segmentation_config
        calls["tracking"] = tracking_config
        calls["max_gland_pages"] = max_gland_pages
        calls["z_visual_scale"] = z_visual_scale
        return GlandInstanceDetectionResult(
            out_dir=tmp_path / "out",
            slice_gland_instances_geojson=tmp_path / "out" / "slice_gland_instances.geojson",
            slice_qc_csv=tmp_path / "out" / "slice_gland_instances_qc.csv",
            gland_instance_tracks_csv=tmp_path / "out" / "gland_instance_tracks.csv",
            gland_instance_qc_index_csv=tmp_path / "out" / "gland_instance_qc_index.csv",
            gland_instance_atlas_html=tmp_path / "out" / "gland_instance_atlas.html",
            slice_count=0,
            total_gland_instances=0,
            gland_count=0,
        )

    monkeypatch.setattr(cli, "run_gland_instance_detection", fake_detection)

    cli.main(
        [
            "detect-gland-instances",
            "--stack-root",
            str(tmp_path / "stack"),
            "--aligned-cells-parquet",
            str(tmp_path / "aligned_cells.parquet"),
            "--out-dir",
            str(tmp_path / "out"),
            "--epithelial-structures",
            "Structure 3,Structure 4",
            "--markers",
            "EPCAM,MUC2",
            "--allow-many-to-many",
            "--max-gland-pages",
            "0",
            "--z-visual-scale",
            "6",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert payload["gland_count"] == 0
    assert calls["segmentation"].structures == ("Structure 3", "Structure 4")
    assert calls["segmentation"].epithelial_markers == ("EPCAM", "MUC2")
    assert calls["tracking"].use_one_to_one is False
    assert calls["max_gland_pages"] is None
    assert calls["z_visual_scale"] == 6.0


def _seg_cfg(
    stack_root: Path,
    cells_path: Path,
    out_dir: Path,
    **overrides,
) -> GlandInstanceSegmentationConfig:
    params = {
        "stack_root": stack_root,
        "aligned_cells_parquet": cells_path,
        "out_dir": out_dir,
        "structures": ("Structure 3",),
        "epithelial_markers": ("EPCAM",),
        "pixel_size_um": 1.0,
        "raster_pixel_size_um": 2.0,
        "lumen_min_area_um2": 50.0,
        "lumen_cell_density_threshold": 0.08,
        "cell_density_sigma_um": 5.0,
        "epithelial_density_threshold": 0.05,
        "support_closing_um": 8.0,
        "min_gland_area_um2": 100.0,
        "min_ring_support_score": 0.2,
    }
    params.update(overrides)
    return GlandInstanceSegmentationConfig(**params)


def _write_segmentation_fixture(
    tmp_path: Path,
    *,
    semantic_polygon: Polygon,
    cells: pd.DataFrame,
) -> tuple[Path, Path]:
    stack_root = tmp_path / "stack"
    stack_root.mkdir()
    xenium_dir = tmp_path / "xenium" / "slice_1"
    xenium_dir.mkdir(parents=True)
    geojson_path = stack_root / "slice_1_aligned.geojson"
    geojson_path.write_text(
        json.dumps(
            {
                "type": "FeatureCollection",
                "features": [
                    {
                        "type": "Feature",
                        "geometry": mapping(semantic_polygon),
                        "properties": {
                            "structure": "Structure 3",
                            "assigned_structure": "Structure 3",
                            "structure_id": 3,
                        },
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    pd.DataFrame(
        [
            {
                "order": 1,
                "sample_id": "slice_1",
                "z_um": 0.0,
                "aligned_geojson": str(geojson_path),
            }
        ]
    ).to_csv(stack_root / "aligned_slice_manifest.csv", index=False)
    pd.DataFrame(
        [
            {
                "order": 1,
                "sample_id": "slice_1",
                "sample_dir": str(xenium_dir.parent),
                "xenium_dir": str(xenium_dir),
            }
        ]
    ).to_csv(stack_root / "xenium_slice_manifest.csv", index=False)
    cells_path = tmp_path / "aligned_cells.parquet"
    cells.to_parquet(cells_path, index=False)
    return stack_root, cells_path


def _ring_cells(centers: list[tuple[float, float]], *, radius: float, n: int) -> pd.DataFrame:
    rows = []
    for center_index, (cx, cy) in enumerate(centers):
        for point_index, theta in enumerate(np.linspace(0, 2 * np.pi, n, endpoint=False)):
            rows.append(
                {
                    "sample_id": "slice_1",
                    "barcode": f"cell_{center_index}_{point_index}",
                    "slice_order": 1,
                    "x_3d_um": float(cx + radius * np.cos(theta)),
                    "y_3d_um": float(cy + radius * np.sin(theta)),
                    "EPCAM": 2.0,
                }
            )
    return pd.DataFrame(rows)


def _box(x0: float, y0: float, x1: float, y1: float) -> Polygon:
    return Polygon([(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)])
