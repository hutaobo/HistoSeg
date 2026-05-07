from __future__ import annotations

import builtins
import json
from pathlib import Path
import sys
import types

import numpy as np
import pandas as pd
import pytest

from histoseg.threed import (
    LocalZOrientationConfig,
    LocalZOrientationResult,
    run_local_z_orientation_correction,
)


def test_local_z_orientation_corrects_reversed_transcript_slice(tmp_path):
    stack_root = _write_stack(tmp_path, ["s1", "s2", "s3"])
    xenium_root = tmp_path / "xenium"
    _write_transcripts(xenium_root / "s1", low_gene="A", high_gene="B")
    # s2 is experimentally flipped: its raw low band is the global high-side marker.
    _write_transcripts(xenium_root / "s2", low_gene="C", high_gene="B")
    _write_transcripts(xenium_root / "s3", low_gene="C", high_gene="D")

    result = run_local_z_orientation_correction(
        LocalZOrientationConfig(
            xenium_root=xenium_root,
            stack_root=stack_root,
            vertical_qc_backend="none",
            apply_local_z_flip=True,
            z_band_fraction=0.2,
        )
    )

    manifest = pd.read_csv(result.manifest_csv)
    assert manifest["selected_orientation"].tolist() == ["preserve", "reverse", "preserve"]
    assert manifest["scoring_backend"].tolist() == ["global", "global", "global"]
    assert bool(manifest["low_confidence"].iloc[0]) is False

    aligned = pd.read_parquet(result.aligned_transcripts_parquet)
    s2_c = aligned[(aligned["sample_id"] == "s2") & (aligned["gene"] == "C")]
    s2_b = aligned[(aligned["sample_id"] == "s2") & (aligned["gene"] == "B")]
    assert float(s2_c["z_raw_um"].iloc[0]) == pytest.approx(1.0)
    assert float(s2_c["local_z_corrected_um"].iloc[0]) == pytest.approx(9.0)
    assert float(s2_c["z_3d_um"].iloc[0]) == pytest.approx(24.0)
    assert float(s2_b["z_raw_um"].iloc[0]) == pytest.approx(9.0)
    assert float(s2_b["local_z_corrected_um"].iloc[0]) == pytest.approx(1.0)
    assert float(s2_b["z_3d_um"].iloc[0]) == pytest.approx(16.0)
    np.testing.assert_allclose(aligned["x_3d_um"], aligned["x_raw_um"])
    np.testing.assert_allclose(aligned["y_3d_um"], aligned["y_raw_um"])
    assert result.biological_report_html.exists()


def test_doublet_regions_are_masked_from_orientation_signature(monkeypatch, tmp_path):
    import histoseg.threed.local_z_orientation as lzo

    stack_root = _write_stack(tmp_path, ["s1", "s2", "s3"])
    xenium_root = tmp_path / "xenium"
    _write_transcripts(xenium_root / "s1", low_gene="A", high_gene="B")
    _write_transcripts(
        xenium_root / "s2",
        low_gene="C",
        high_gene="B",
        high_doublet_gene="A",
        high_doublet_count=60,
    )
    _write_transcripts(xenium_root / "s3", low_gene="C", high_gene="D")

    def fake_qc(transcripts, sample_qc_dir, cfg):
        doublet_rows = transcripts.loc[
            transcripts["x_raw_um"] > 90.0,
            ["x_raw_um", "y_raw_um"],
        ]
        doublets = doublet_rows.rename(columns={"x_raw_um": "x", "y_raw_um": "y"})
        path = sample_qc_dir / "doublets.csv"
        doublets.to_csv(path, index=False)
        return lzo._OvrlpyQC(status="ok", doublets_csv=path, doublet_count=len(doublets))

    monkeypatch.setattr(lzo, "_run_vertical_qc", fake_qc)

    result = run_local_z_orientation_correction(
        LocalZOrientationConfig(
            xenium_root=xenium_root,
            stack_root=stack_root,
            vertical_qc_backend="none",
            apply_local_z_flip=True,
            z_band_fraction=0.2,
            doublet_exclusion_radius_um=0.5,
        )
    )

    manifest = pd.read_csv(result.manifest_csv)
    assert manifest["selected_orientation"].tolist() == ["preserve", "reverse", "preserve"]
    s2 = manifest.loc[manifest["sample_id"] == "s2"].iloc[0]
    assert int(s2["excluded_doublet_transcripts"]) >= 60


def test_contour_aware_orientation_resolves_global_signature_tie(tmp_path):
    stack_root = _write_contour_stack(tmp_path, ["s1", "s2", "s3"])
    xenium_root = tmp_path / "xenium"
    _write_contour_transcripts(
        xenium_root / "s1",
        low_genes={"left": "L1", "right": "L2"},
        high_genes={"left": "X", "right": "Y"},
    )
    _write_contour_transcripts(
        xenium_root / "s2",
        low_genes={"left": "Y", "right": "X"},
        high_genes={"left": "X", "right": "Y"},
    )
    _write_contour_transcripts(
        xenium_root / "s3",
        low_genes={"left": "Y", "right": "X"},
        high_genes={"left": "H1", "right": "H2"},
    )

    result = run_local_z_orientation_correction(
        LocalZOrientationConfig(
            xenium_root=xenium_root,
            stack_root=stack_root,
            vertical_qc_backend="none",
            apply_local_z_flip=True,
            pixel_size_um=1.0,
            z_band_fraction=0.2,
            contour_min_transcripts=20,
            orientation_spatial_unit="contour",
            orientation_bootstrap_iterations=20,
        )
    )

    manifest = pd.read_csv(result.manifest_csv)
    assert manifest["selected_orientation"].tolist() == ["preserve", "reverse", "preserve"]
    assert manifest["scoring_backend"].tolist() == ["contour", "contour", "contour"]
    assert int(manifest.loc[manifest["sample_id"] == "s2", "contour_pairs_to_previous"].iloc[0]) == 2
    assert int(manifest.loc[manifest["sample_id"] == "s2", "contour_pairs_to_next"].iloc[0]) == 2

    aligned = pd.read_parquet(result.aligned_transcripts_parquet)
    assert {"contour_id", "contour_structure", "contour_feature_index"}.issubset(aligned.columns)
    assert aligned["contour_id"].notna().all()


def test_contour_matching_allows_split_merge_pairs(tmp_path):
    stack_root = _write_contour_stack(
        tmp_path,
        ["s1", "s2"],
        layouts={
            "s1": [
                ("left-a", "epithelium", (0.0, 0.0, 80.0, 80.0)),
                ("left-b", "epithelium", (100.0, 0.0, 180.0, 80.0)),
            ],
            "s2": [
                ("merged", "epithelium", (0.0, 0.0, 180.0, 80.0)),
            ],
        },
    )
    xenium_root = tmp_path / "xenium"
    _write_contour_transcripts(
        xenium_root / "s1",
        low_genes={"left-a": "A", "left-b": "B"},
        high_genes={"left-a": "C", "left-b": "D"},
        centers={"left-a": (20.0, 20.0), "left-b": (130.0, 20.0)},
    )
    _write_contour_transcripts(
        xenium_root / "s2",
        low_genes={"merged": "C"},
        high_genes={"merged": "E"},
        centers={"merged": (80.0, 20.0)},
        per_band=80,
    )

    result = run_local_z_orientation_correction(
        LocalZOrientationConfig(
            xenium_root=xenium_root,
            stack_root=stack_root,
            vertical_qc_backend="none",
            pixel_size_um=1.0,
            contour_min_transcripts=20,
            orientation_spatial_unit="contour",
            orientation_bootstrap_iterations=0,
        )
    )

    pairs = pd.read_csv(result.contour_pair_continuity_csv)
    assert len(pairs) >= 2
    assert pairs["next_contour_id"].nunique() == 1
    assert pairs["prev_contour_id"].nunique() == 2


def test_contour_doublet_burden_masks_affected_profile(monkeypatch, tmp_path):
    import histoseg.threed.local_z_orientation as lzo

    stack_root = _write_contour_stack(tmp_path, ["s1", "s2"])
    xenium_root = tmp_path / "xenium"
    _write_contour_transcripts(
        xenium_root / "s1",
        low_genes={"left": "A", "right": "B"},
        high_genes={"left": "C", "right": "D"},
    )
    _write_contour_transcripts(
        xenium_root / "s2",
        low_genes={"left": "C", "right": "D"},
        high_genes={"left": "E", "right": "F"},
    )

    def fake_qc(transcripts, sample_qc_dir, cfg):
        doublets = pd.DataFrame({"x": [10.0], "y": [10.0]})
        path = sample_qc_dir / "doublets.csv"
        doublets.to_csv(path, index=False)
        return lzo._OvrlpyQC(status="ok", doublets_csv=path, doublet_count=1)

    monkeypatch.setattr(lzo, "_run_vertical_qc", fake_qc)

    result = run_local_z_orientation_correction(
        LocalZOrientationConfig(
            xenium_root=xenium_root,
            stack_root=stack_root,
            vertical_qc_backend="none",
            pixel_size_um=1.0,
            contour_min_transcripts=20,
            doublet_exclusion_radius_um=5.0,
            orientation_spatial_unit="contour",
            orientation_bootstrap_iterations=0,
        )
    )

    profiles = pd.read_parquet(result.vertical_qc_dir / "s1" / "contour_layer_profiles.parquet")
    affected = profiles.sort_values("doublet_burden", ascending=False).iloc[0]
    assert float(affected["doublet_burden"]) == pytest.approx(1.0)
    assert int(affected["excluded_doublet_transcripts"]) >= 60
    assert bool(affected["usable"]) is False


def test_auto_falls_back_to_global_without_aligned_geojson(tmp_path):
    stack_root = _write_stack(tmp_path, ["s1", "s2"])
    xenium_root = tmp_path / "xenium"
    _write_transcripts(xenium_root / "s1", low_gene="A", high_gene="B")
    _write_transcripts(xenium_root / "s2", low_gene="B", high_gene="C")

    result = run_local_z_orientation_correction(
        LocalZOrientationConfig(
            xenium_root=xenium_root,
            stack_root=stack_root,
            vertical_qc_backend="none",
            orientation_spatial_unit="auto",
        )
    )

    manifest = pd.read_csv(result.manifest_csv)
    assert manifest["scoring_backend"].tolist() == ["global", "global"]
    assert manifest["contour_status"].tolist() == ["missing_aligned_geojson", "missing_aligned_geojson"]


def test_forced_contour_mode_requires_some_contour_data(tmp_path):
    stack_root = _write_stack(tmp_path, ["s1", "s2"])
    xenium_root = tmp_path / "xenium"
    _write_transcripts(xenium_root / "s1", low_gene="A", high_gene="B")
    _write_transcripts(xenium_root / "s2", low_gene="B", high_gene="C")

    with pytest.raises(ValueError, match="no aligned contour"):
        run_local_z_orientation_correction(
            LocalZOrientationConfig(
                xenium_root=xenium_root,
                stack_root=stack_root,
                vertical_qc_backend="none",
                orientation_spatial_unit="contour",
            )
        )


def test_local_z_orientation_requires_transcript_z_column(tmp_path):
    stack_root = _write_stack(tmp_path, ["s1"])
    sample_dir = tmp_path / "xenium" / "s1"
    sample_dir.mkdir(parents=True)
    pd.DataFrame({"gene": ["A"], "x": [0.0], "y": [0.0]}).to_parquet(
        sample_dir / "transcripts.parquet",
        index=False,
    )

    with pytest.raises(ValueError, match="missing a z column"):
        run_local_z_orientation_correction(
            LocalZOrientationConfig(
                xenium_root=tmp_path / "xenium",
                stack_root=stack_root,
                vertical_qc_backend="none",
            )
        )


def test_local_z_orientation_reads_pyxenium_slide_zarr(monkeypatch, tmp_path):
    stack_root = _write_stack(tmp_path, ["s1"])
    slide_dir = tmp_path / "slides" / "s1.pyxenium.slide.zarr"
    slide_dir.mkdir(parents=True)

    fake_module = types.ModuleType("pyXenium")

    class FakeSlide:
        points = {
            "transcripts": pd.DataFrame(
                {
                    "gene_name": ["A", "B", "drop"],
                    "x": [1.0, 2.0, 3.0],
                    "y": [4.0, 5.0, 6.0],
                    "z": [7.0, 8.0, 9.0],
                    "quality_score": [20.0, 21.0, 22.0],
                    "cell_id": ["c1", "UNASSIGNED", "c3"],
                    "valid": [True, True, False],
                }
            )
        }

    def fake_read_slide(path):
        assert Path(path) == slide_dir
        return FakeSlide()

    fake_module.read_slide = fake_read_slide
    monkeypatch.setitem(sys.modules, "pyXenium", fake_module)

    result = run_local_z_orientation_correction(
        LocalZOrientationConfig(
            xenium_root=tmp_path / "slides",
            stack_root=stack_root,
            sample_glob="*.pyxenium.slide.zarr",
            vertical_qc_backend="none",
        )
    )

    manifest = pd.read_csv(result.manifest_csv)
    assert manifest["sample_id"].tolist() == ["s1"]
    assert int(manifest["transcript_count"].iloc[0]) == 2

    aligned = pd.read_parquet(result.aligned_transcripts_parquet)
    assert aligned["gene"].tolist() == ["A", "B"]
    assert aligned["qv"].tolist() == [20.0, 21.0]
    assert aligned["cell_id"].tolist() == ["c1", "UNASSIGNED"]


def test_pyxenium_slide_zarr_can_use_sidecar_transcript_table(monkeypatch, tmp_path):
    stack_root = _write_stack(tmp_path, ["s1"])
    slides_root = tmp_path / "slides"
    slide_dir = slides_root / "s1.pyxenium.slide.zarr"
    slide_dir.mkdir(parents=True)
    sidecar_dir = slides_root / "s1"
    sidecar_dir.mkdir()
    pd.DataFrame(
        {
            "gene": ["A", "B"],
            "x": [1.0, 2.0],
            "y": [3.0, 4.0],
            "z": [5.0, 6.0],
        }
    ).to_parquet(sidecar_dir / "transcripts.parquet", index=False)

    fake_module = types.ModuleType("pyXenium")

    class FakeSlide:
        points = {}

    fake_module.read_slide = lambda path: FakeSlide()
    monkeypatch.setitem(sys.modules, "pyXenium", fake_module)

    result = run_local_z_orientation_correction(
        LocalZOrientationConfig(
            xenium_root=slides_root,
            stack_root=stack_root,
            sample_glob="*.pyxenium.slide.zarr",
            vertical_qc_backend="none",
        )
    )

    manifest = pd.read_csv(result.manifest_csv)
    assert manifest["sample_id"].tolist() == ["s1"]
    assert int(manifest["transcript_count"].iloc[0]) == 2
    aligned = pd.read_parquet(result.aligned_transcripts_parquet)
    assert aligned["gene"].tolist() == ["A", "B"]


def test_aligned_transcripts_keep_stable_optional_schema(tmp_path):
    stack_root = _write_stack(tmp_path, ["s1", "s2"])
    xenium_root = tmp_path / "xenium"
    _write_transcripts(xenium_root / "s1", low_gene="A", high_gene="B")

    s2 = xenium_root / "s2"
    s2.mkdir(parents=True)
    rows = []
    for index in range(30):
        rows.append(
            {
                "gene": "B",
                "x": float(index),
                "y": 0.0,
                "z": 1.0,
                "transcript_id": f"tx-low-{index}",
                "cell_id": "c1",
                "overlaps_nucleus": 0,
                "qv": 20.0,
                "structure": "low",
            }
        )
        rows.append(
            {
                "gene": "C",
                "x": float(index),
                "y": 10.0,
                "z": 9.0,
                "transcript_id": f"tx-high-{index}",
                "cell_id": "c2",
                "overlaps_nucleus": 1,
                "qv": 21.0,
                "structure": "high",
            }
        )
    pd.DataFrame(rows).to_parquet(s2 / "transcripts.parquet", index=False)

    result = run_local_z_orientation_correction(
        LocalZOrientationConfig(
            xenium_root=xenium_root,
            stack_root=stack_root,
            vertical_qc_backend="none",
            apply_local_z_flip=True,
            z_band_fraction=0.2,
        )
    )

    aligned = pd.read_parquet(result.aligned_transcripts_parquet)
    for column in ["transcript_id", "cell_id", "overlaps_nucleus", "qv", "structure"]:
        assert column in aligned.columns
    assert len(aligned) == 120
    assert aligned.loc[aligned["sample_id"] == "s1", "transcript_id"].isna().all()
    assert set(aligned.loc[aligned["sample_id"] == "s2", "overlaps_nucleus"]) == {0, 1}


def test_local_z_orientation_reports_missing_ovrlpy(monkeypatch, tmp_path):
    stack_root = _write_stack(tmp_path, ["s1"])
    xenium_root = tmp_path / "xenium"
    _write_transcripts(xenium_root / "s1", low_gene="A", high_gene="B")
    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "ovrlpy":
            raise ImportError("blocked for test")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError, match="ovrlpy is not installed"):
        run_local_z_orientation_correction(
            LocalZOrientationConfig(
                xenium_root=xenium_root,
                stack_root=stack_root,
                vertical_qc_backend="ovrlpy",
            )
        )


def test_infer_local_z_orientation_cli_parses_config(monkeypatch, tmp_path, capsys):
    import histoseg.threed.cli as cli

    calls: list[LocalZOrientationConfig] = []

    def fake_run(cfg):
        calls.append(cfg)
        return LocalZOrientationResult(
            out_dir=Path(cfg.out_dir or cfg.stack_root),
            manifest_csv=tmp_path / "local_z_orientation_manifest.csv",
            aligned_transcripts_parquet=tmp_path / "aligned_transcripts_3d.parquet",
            biological_report_html=tmp_path / "biological_z_report.html",
            marker_gradients_csv=tmp_path / "top_bottom_marker_gradients.csv",
            vertical_qc_dir=tmp_path / "vertical_qc",
            slice_count=2,
            transcript_count=10,
            best_score=1.0,
            second_score=0.5,
            confidence_margin=0.5,
            low_confidence=False,
        )

    monkeypatch.setattr(cli, "run_local_z_orientation_correction", fake_run)
    cli.main(
        [
            "infer-local-z-orientation",
            "--stack-root",
            str(tmp_path / "stack"),
            "--xenium-root",
            str(tmp_path / "xenium"),
            "--out-dir",
            str(tmp_path / "out"),
            "--vertical-qc-backend",
            "ovrlpy",
            "--apply-local-z-flip",
            "--z-band-fraction",
            "0.2",
            "--doublet-exclusion-radius-um",
            "12.5",
            "--no-ovrlpy-fit-umap",
            "--ovrlpy-min-transcripts",
            "3",
            "--orientation-spatial-unit",
            "contour",
            "--contour-group-property",
            "structure_label",
            "--contour-min-transcripts",
            "25",
            "--contour-match-min-iou",
            "0.05",
            "--contour-match-max-distance-um",
            "90",
            "--orientation-bootstrap-iterations",
            "7",
            "--orientation-bootstrap-seed",
            "11",
        ]
    )

    capsys.readouterr()
    assert len(calls) == 1
    cfg = calls[0]
    assert cfg.vertical_qc_backend == "ovrlpy"
    assert cfg.apply_local_z_flip is True
    assert cfg.z_band_fraction == pytest.approx(0.2)
    assert cfg.doublet_exclusion_radius_um == pytest.approx(12.5)
    assert cfg.ovrlpy_fit_umap is False
    assert cfg.ovrlpy_min_transcripts == pytest.approx(3.0)
    assert cfg.orientation_spatial_unit == "contour"
    assert cfg.contour_group_property == "structure_label"
    assert cfg.contour_min_transcripts == 25
    assert cfg.contour_match_min_iou == pytest.approx(0.05)
    assert cfg.contour_match_max_distance_um == pytest.approx(90.0)
    assert cfg.orientation_bootstrap_iterations == 7
    assert cfg.orientation_bootstrap_seed == 11


def _write_stack(tmp_path: Path, sample_ids: list[str]) -> Path:
    stack_root = tmp_path / "stack"
    stack_root.mkdir()
    pd.DataFrame(
        [
            {
                "order": index,
                "sample_id": sample_id,
                "z_um": float((index - 1) * 20.0),
                "alignment_role": "reference" if index == 1 else "moving",
            }
            for index, sample_id in enumerate(sample_ids, start=1)
        ]
    ).to_csv(stack_root / "aligned_slice_manifest.csv", index=False)
    pairwise_rows = [
        {
            "moving_order": index,
            "moving_sample_id": sample_id,
            "soft_accepted": False,
            "soft_summary_json": "",
        }
        for index, sample_id in enumerate(sample_ids, start=1)
        if index > 1
    ]
    pd.DataFrame(
        pairwise_rows,
        columns=["moving_order", "moving_sample_id", "soft_accepted", "soft_summary_json"],
    ).to_csv(stack_root / "pairwise_alignment_metrics.csv", index=False)
    for index, sample_id in enumerate(sample_ids, start=1):
        if index == 1:
            continue
        pair_dir = stack_root / "pairwise_alignments" / f"{index - 1:03d}_to_{index:03d}_{sample_id}"
        pair_dir.mkdir(parents=True)
        hard_summary = {
            "transform": {
                "origin_x": 0.0,
                "origin_y": 0.0,
                "rotation_degrees": 0.0,
                "scale": 1.0,
                "translate_x": 0.0,
                "translate_y": 0.0,
            },
            "union_iou_before_hard": 1.0,
            "union_iou_after_hard": 1.0,
            "hard_alignment_accepted": True,
        }
        (pair_dir / "hard_similarity_alignment.json").write_text(
            json.dumps(hard_summary),
            encoding="utf-8",
        )
    return stack_root


def _write_contour_stack(
    tmp_path: Path,
    sample_ids: list[str],
    *,
    layouts: dict[str, list[tuple[str, str, tuple[float, float, float, float]]]] | None = None,
) -> Path:
    stack_root = _write_stack(tmp_path, sample_ids)
    aligned_dir = stack_root / "aligned_contours"
    aligned_dir.mkdir()
    default_layout = [
        ("left", "epithelium", (0.0, 0.0, 80.0, 80.0)),
        ("right", "epithelium", (200.0, 0.0, 280.0, 80.0)),
    ]
    manifest = pd.read_csv(stack_root / "aligned_slice_manifest.csv")
    aligned_paths = []
    for row in manifest.itertuples(index=False):
        sample_id = str(row.sample_id)
        path = aligned_dir / f"{int(row.order):03d}_{sample_id}_aligned.geojson"
        features = []
        for name, structure, bounds in (layouts or {}).get(sample_id, default_layout):
            minx, miny, maxx, maxy = bounds
            features.append(
                {
                    "type": "Feature",
                    "properties": {"name": name, "structure": structure},
                    "geometry": {
                        "type": "Polygon",
                        "coordinates": [
                            [
                                [minx, miny],
                                [maxx, miny],
                                [maxx, maxy],
                                [minx, maxy],
                                [minx, miny],
                            ]
                        ],
                    },
                }
            )
        path.write_text(
            json.dumps({"type": "FeatureCollection", "features": features}),
            encoding="utf-8",
        )
        aligned_paths.append(str(path))
    manifest["aligned_geojson"] = aligned_paths
    manifest.to_csv(stack_root / "aligned_slice_manifest.csv", index=False)
    return stack_root


def _write_contour_transcripts(
    sample_dir: Path,
    *,
    low_genes: dict[str, str],
    high_genes: dict[str, str],
    centers: dict[str, tuple[float, float]] | None = None,
    per_band: int = 30,
) -> None:
    sample_dir.mkdir(parents=True, exist_ok=True)
    default_centers = {
        "left": (10.0, 10.0),
        "right": (210.0, 10.0),
    }
    centers = {**default_centers, **(centers or {})}
    rows = []
    for contour_name in sorted(set(low_genes) | set(high_genes)):
        cx, cy = centers[contour_name]
        for index in range(per_band):
            dx = float(index % 5) * 0.1
            dy = float(index // 5) * 0.1
            if contour_name in low_genes:
                rows.append(
                    {
                        "gene": low_genes[contour_name],
                        "x": cx + dx,
                        "y": cy + dy,
                        "z": 1.0,
                    }
                )
            if contour_name in high_genes:
                rows.append(
                    {
                        "gene": high_genes[contour_name],
                        "x": cx + dx,
                        "y": cy + dy,
                        "z": 9.0,
                    }
                )
    pd.DataFrame(rows).to_parquet(sample_dir / "transcripts.parquet", index=False)


def _write_transcripts(
    sample_dir: Path,
    *,
    low_gene: str,
    high_gene: str,
    high_doublet_gene: str | None = None,
    high_doublet_count: int = 0,
) -> None:
    sample_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for index in range(30):
        rows.append({"gene": low_gene, "x": float(index), "y": 0.0, "z": 1.0})
        rows.append({"gene": high_gene, "x": float(index), "y": 10.0, "z": 9.0})
    for index in range(high_doublet_count):
        rows.append(
            {
                "gene": high_doublet_gene,
                "x": 100.0 + float(index) * 0.01,
                "y": 100.0,
                "z": 9.0,
            }
        )
    pd.DataFrame(rows).to_parquet(sample_dir / "transcripts.parquet", index=False)
