from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from histoseg.contour import (
    GeneTranscriptIsolineConfig,
    GeneTranscriptIsolineResult,
    run_gene_transcript_isoline,
)
from histoseg.contour.gene_isoline import _discover_xenium_transcript_samples
from histoseg.contour.pattern1_isoline import Pattern1IsolineResult


def test_gene_isoline_sample_discovery_uses_numeric_order(tmp_path):
    for name in ["A079-C-008_10", "A079-C-008_1", "A079-C-008_2"]:
        xenium = tmp_path / name / "output-XETG000"
        xenium.mkdir(parents=True)
        (xenium / "cells.parquet").write_text("", encoding="utf-8")
        (xenium / "transcripts.parquet").write_text("", encoding="utf-8")

    samples = _discover_xenium_transcript_samples(tmp_path, sample_glob="A079-C-008_*")

    assert [item.sample_id for item in samples] == [
        "A079-C-008_1",
        "A079-C-008_2",
        "A079-C-008_10",
    ]
    assert [item.order for item in samples] == [1, 2, 3]


def test_gene_isoline_sample_discovery_accepts_single_sample_root(tmp_path):
    xenium = tmp_path / "A079-C-008_1" / "output-XETG000"
    xenium.mkdir(parents=True)
    (xenium / "cells.parquet").write_text("", encoding="utf-8")
    (xenium / "transcripts.parquet").write_text("", encoding="utf-8")

    samples = _discover_xenium_transcript_samples(tmp_path / "A079-C-008_1")

    assert len(samples) == 1
    assert samples[0].sample_id == "A079-C-008_1"
    assert samples[0].xenium_dir == xenium.resolve()


def test_gene_isoline_multigene_outputs_geojson_and_run_log(tmp_path):
    _write_xenium_sample(tmp_path, "A079-C-008_1")
    out_dir = tmp_path / "gene_outputs"

    result = run_gene_transcript_isoline(
        GeneTranscriptIsolineConfig(
            xenium_root=tmp_path,
            out_dir=out_dir,
            genes=["GREM1", "COL1A1"],
            sample_glob="A079-C-008_*",
            qv_min=20,
            min_transcripts=10,
            grid_n=70,
            knn_k=5,
            smooth_sigma=1.0,
            min_cells_inside=3,
            use_synth_bg=False,
            xenium_pixel_size_um=1.0,
        )
    )

    assert result.run_log_csv.exists()
    log_df = pd.read_csv(result.run_log_csv)
    assert set(log_df["gene"]) == {"GREM1", "COL1A1"}

    grem1_row = log_df.loc[log_df["gene"] == "GREM1"].iloc[0]
    assert grem1_row["status"] == "ok"
    assert int(grem1_row["n_transcripts"]) == 12
    assert int(grem1_row["n_contours"]) >= 1
    assert pd.isna(grem1_row["confidence"])

    col1a1_row = log_df.loc[log_df["gene"] == "COL1A1"].iloc[0]
    assert col1a1_row["status"] == "skip_too_few_transcripts"
    assert int(col1a1_row["n_transcripts"]) == 4

    geojson_path = out_dir / "A079-C-008_1" / "xenium_explorer_annotations.geojson"
    csv_path = out_dir / "A079-C-008_1" / "xenium_explorer_annotations.csv"
    summary_path = out_dir / "A079-C-008_1" / "xenium_explorer_annotations_summary.csv"
    assert geojson_path.exists()
    assert csv_path.exists()
    assert summary_path.exists()

    payload = json.loads(geojson_path.read_text(encoding="utf-8"))
    assert payload["type"] == "FeatureCollection"
    assert len(payload["features"]) >= 1
    feature = payload["features"][0]
    props = feature["properties"]
    assert feature["geometry"]["type"] == "Polygon"
    assert props["objectType"] == "annotation"
    assert props["gene"] == "GREM1"
    assert props["structure"] == "GREM1"
    assert props["sample_id"] == "A079-C-008_1"
    assert props["source"] == "gene_transcript_isoline"
    assert len(props["classification"]["color"]) == 3

    summary_df = pd.read_csv(summary_path)
    assert set(summary_df["Gene"]) == {"GREM1"}
    assert set(summary_df["Source"]) == {"gene_transcript_isoline"}
    assert not (out_dir / "A079-C-008_1" / "GREM1" / "prepared_inputs").exists()


def test_gene_isoline_confidence_default_and_flag_are_passed_through(tmp_path, monkeypatch):
    _write_xenium_sample(tmp_path, "A079-C-008_1")
    captured: list[bool] = []

    def fake_run_pattern1(cfg):
        captured.append(bool(cfg.compute_confidence_score))
        return Pattern1IsolineResult(
            out_dir=Path(cfg.out_dir),
            id_col_used="cell_id",
            x_col="x_centroid",
            y_col="y_centroid",
            n_target_cells=12,
            n_bg0_points=20,
            contours=[
                np.array(
                    [
                        [40.0, 40.0],
                        [60.0, 40.0],
                        [60.0, 60.0],
                        [40.0, 60.0],
                        [40.0, 40.0],
                    ]
                )
            ],
            label_scheme="p1_is_one",
            segmentation_confidence_score=0.42 if cfg.compute_confidence_score else None,
        )

    import histoseg.contour.gene_isoline as gene_isoline

    monkeypatch.setattr(gene_isoline, "run_pattern1_isoline", fake_run_pattern1)

    run_gene_transcript_isoline(
        GeneTranscriptIsolineConfig(
            xenium_root=tmp_path,
            out_dir=tmp_path / "default_confidence",
            genes=["GREM1"],
            sample_glob="A079-C-008_*",
            use_synth_bg=False,
        )
    )
    run_gene_transcript_isoline(
        GeneTranscriptIsolineConfig(
            xenium_root=tmp_path,
            out_dir=tmp_path / "enabled_confidence",
            genes=["GREM1"],
            sample_glob="A079-C-008_*",
            use_synth_bg=False,
            compute_confidence_score=True,
        )
    )

    assert captured == [False, True]


def test_gene_isoline_cli_parses_arguments(monkeypatch, tmp_path, capsys):
    import histoseg.contour.cli as cli

    calls: list[GeneTranscriptIsolineConfig] = []

    def fake_run(cfg):
        calls.append(cfg)
        return GeneTranscriptIsolineResult(
            out_dir=Path(cfg.out_dir),
            run_log_csv=Path(cfg.out_dir) / "gene_isoline_run_log.csv",
            records=[],
        )

    monkeypatch.setattr(cli, "run_gene_transcript_isoline", fake_run)
    cli.main(
        [
            "gene-isoline",
            "--xenium-root",
            str(tmp_path),
            "--out-dir",
            str(tmp_path / "out"),
            "--genes",
            "GREM1,COL1A1",
            "--sample-glob",
            "A079-C-008_*",
            "--qv-min",
            "25",
            "--min-transcripts",
            "12",
            "--grid-n",
            "90",
            "--knn-k",
            "7",
            "--smooth-sigma",
            "2.5",
            "--min-cells-inside",
            "4",
            "--alpha",
            "0.04",
            "--xenium-pixel-size-um",
            "1.0",
            "--compute-confidence-score",
            "--keep-prepared-inputs",
            "--no-synth-bg",
            "--fail-fast",
        ]
    )

    capsys.readouterr()
    assert len(calls) == 1
    cfg = calls[0]
    assert list(cfg.genes) == ["GREM1", "COL1A1"]
    assert cfg.sample_glob == "A079-C-008_*"
    assert cfg.qv_min == 25
    assert cfg.min_transcripts == 12
    assert cfg.grid_n == 90
    assert cfg.knn_k == 7
    assert cfg.smooth_sigma == 2.5
    assert cfg.min_cells_inside == 4
    assert cfg.alpha == 0.04
    assert cfg.xenium_pixel_size_um == 1.0
    assert cfg.compute_confidence_score is True
    assert cfg.keep_prepared_inputs is True
    assert cfg.use_synth_bg is False
    assert cfg.fail_fast is True


def _write_xenium_sample(root: Path, sample_name: str) -> Path:
    xenium = root / sample_name / "output-XETG000"
    xenium.mkdir(parents=True)

    angles = np.linspace(0.0, 2.0 * np.pi, 96, endpoint=False)
    cell_rows = []
    for index, angle in enumerate(angles):
        radius = 38.0 if index % 2 == 0 else 55.0
        cell_rows.append(
            {
                "cell_id": f"cell_{index}",
                "x_centroid": 50.0 + radius * np.cos(angle),
                "y_centroid": 50.0 + radius * np.sin(angle),
            }
        )
    pd.DataFrame(cell_rows).to_parquet(xenium / "cells.parquet", index=False)

    transcript_rows = []
    target_angles = np.linspace(0.0, 2.0 * np.pi, 15, endpoint=False)
    for index, angle in enumerate(target_angles):
        transcript_rows.append(
            {
                "transcript_id": f"grem1_{index}",
                "feature_name": "GREM1",
                "x_location": 50.0 + 5.0 * np.cos(angle),
                "y_location": 50.0 + 5.0 * np.sin(angle),
                "qv": 30 if index < 12 else 10,
            }
        )
    for index, angle in enumerate(target_angles[:4]):
        transcript_rows.append(
            {
                "transcript_id": f"col1a1_{index}",
                "feature_name": "COL1A1",
                "x_location": 52.0 + 4.0 * np.cos(angle),
                "y_location": 52.0 + 4.0 * np.sin(angle),
                "qv": 35,
            }
        )
    pd.DataFrame(transcript_rows).to_parquet(xenium / "transcripts.parquet", index=False)
    return xenium
