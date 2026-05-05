from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from shapely.geometry import box, mapping

from histoseg.threed import (
    GlandQCAtlasConfig,
    GlandQCAtlasResult,
    render_gland_qc_atlas,
)


def test_gland_qc_tracks_clean_three_slice_gland(tmp_path):
    stack_root = _write_stack(
        tmp_path,
        [
            [_feature(box(0, 0, 20, 20), "Structure 3", 1)],
            [_feature(box(2, 0, 22, 20), "Structure 3", 1)],
            [_feature(box(4, 0, 24, 20), "Structure 3", 1)],
        ],
    )

    result = render_gland_qc_atlas(
        GlandQCAtlasConfig(
            stack_root=stack_root,
            out_dir=tmp_path / "qc",
            structures=("Structure 3",),
            pixel_size_um=1.0,
        )
    )

    assert isinstance(result, GlandQCAtlasResult)
    index = pd.read_csv(result.gland_qc_index_csv)
    tracks = pd.read_csv(result.gland_tracks_csv)
    assert len(index) == 1
    assert index.loc[0, "quality_flag"] == "good"
    assert int(index.loc[0, "slice_count"]) == 3
    assert tracks["gland_id"].nunique() == 1
    assert (tmp_path / "qc" / "glands" / "gland_0001.html").exists()
    assert "Plotly.newPlot" in (tmp_path / "qc" / "glands" / "gland_0001.html").read_text(encoding="utf-8")


def test_gland_qc_keeps_nearby_glands_separate(tmp_path):
    stack_root = _write_stack(
        tmp_path,
        [
            [
                _feature(box(0, 0, 20, 20), "Structure 3", 1),
                _feature(box(45, 0, 65, 20), "Structure 3", 2),
            ],
            [
                _feature(box(2, 0, 22, 20), "Structure 3", 1),
                _feature(box(47, 0, 67, 20), "Structure 3", 2),
            ],
            [
                _feature(box(4, 0, 24, 20), "Structure 3", 1),
                _feature(box(49, 0, 69, 20), "Structure 3", 2),
            ],
        ],
    )

    result = render_gland_qc_atlas(
        GlandQCAtlasConfig(
            stack_root=stack_root,
            out_dir=tmp_path / "qc",
            structures=("Structure 3",),
            pixel_size_um=1.0,
        )
    )

    index = pd.read_csv(result.gland_qc_index_csv)
    assert len(index) == 2
    assert sorted(index["slice_count"].astype(int).tolist()) == [3, 3]
    assert set(index["quality_flag"]) == {"good"}


def test_gland_qc_flags_branch_and_isolated_components(tmp_path):
    stack_root = _write_stack(
        tmp_path,
        [
            [
                _feature(box(0, 0, 30, 20), "Structure 3", 1),
                _feature(box(100, 0, 120, 20), "Structure 3", 9),
            ],
            [
                _feature(box(0, 0, 14, 20), "Structure 3", 1),
                _feature(box(16, 0, 30, 20), "Structure 3", 2),
            ],
            [
                _feature(box(1, 0, 15, 20), "Structure 3", 1),
                _feature(box(17, 0, 31, 20), "Structure 3", 2),
            ],
        ],
    )

    result = render_gland_qc_atlas(
        GlandQCAtlasConfig(
            stack_root=stack_root,
            out_dir=tmp_path / "qc",
            structures=("Structure 3",),
            pixel_size_um=1.0,
        )
    )

    index = pd.read_csv(result.gland_qc_index_csv)
    assert "isolated" in set(index["quality_flag"])
    branch_rows = index.loc[index["branch_merge_count"].astype(int) > 0]
    assert not branch_rows.empty
    assert set(branch_rows["quality_flag"]) == {"review"}


def test_gland_qc_cli_writes_atlas_outputs(tmp_path, capsys):
    import histoseg.threed.cli as cli

    stack_root = _write_stack(
        tmp_path,
        [
            [_feature(box(0, 0, 20, 20), "Structure 4", 1)],
            [_feature(box(2, 0, 22, 20), "Structure 4", 1)],
        ],
    )
    out_dir = tmp_path / "qc"

    cli.main(
        [
            "render-gland-qc-atlas",
            "--stack-root",
            str(stack_root),
            "--out-dir",
            str(out_dir),
            "--structures",
            "Structure 4",
            "--pixel-size-um",
            "1.0",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert payload["gland_count"] == 1
    assert (out_dir / "gland_tracks.csv").exists()
    assert (out_dir / "gland_qc_index.csv").exists()
    assert (out_dir / "gland_qc_atlas.html").exists()
    assert (out_dir / "glands" / "gland_0001.html").exists()


def _write_stack(tmp_path: Path, slice_features: list[list[dict]]) -> Path:
    stack_root = tmp_path / "stack"
    aligned_dir = stack_root / "aligned_contours"
    aligned_dir.mkdir(parents=True)
    rows = []
    for index, features in enumerate(slice_features, start=1):
        sample_id = f"s{index}"
        geojson = aligned_dir / f"{index:03d}_{sample_id}_aligned.geojson"
        geojson.write_text(
            json.dumps({"type": "FeatureCollection", "features": features}),
            encoding="utf-8",
        )
        rows.append(
            {
                "order": index,
                "sample_id": sample_id,
                "z_um": float(index - 1) * 5.0,
                "aligned_geojson": str(geojson),
            }
        )
    pd.DataFrame(rows).to_csv(stack_root / "aligned_slice_manifest.csv", index=False)
    return stack_root


def _feature(geom, structure: str, component_index: int) -> dict:
    return {
        "type": "Feature",
        "properties": {
            "assigned_structure": structure,
            "structure_id": int(structure.split()[-1]),
            "component_index": component_index,
            "polygon_index": 1,
            "name": f"{structure} #{component_index}.1",
        },
        "geometry": mapping(geom),
    }
