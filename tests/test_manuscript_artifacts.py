from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_render_manuscript_figures_writes_expected_outputs(tmp_path):
    out_root = tmp_path / "static"
    subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "scripts" / "render_manuscript_figures.py"),
            "--out-root",
            str(out_root),
        ],
        cwd=REPO_ROOT,
        check=True,
    )

    expected = [
        out_root / "manuscript" / "figure1_workflow_schematic.png",
        out_root / "manuscript" / "figure1_workflow_schematic.svg",
        out_root / "manuscript" / "figure1_histology_contour_overlay_template.png",
        out_root / "manuscript" / "figure1_histology_contour_overlay_template.svg",
        out_root / "threed" / "polyp" / "24gene" / "signed_distance_marker_group_violin.png",
        out_root / "threed" / "polyp" / "24gene" / "compartment_summary_table.csv",
        out_root / "threed" / "polyp" / "24gene" / "compartment_summary_table.png",
    ]
    for path in expected:
        assert path.exists()
        assert path.stat().st_size > 0


def test_benchmark_smoke_writes_full_method_matrix(tmp_path):
    out_dir = tmp_path / "benchmarks"
    subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "reproducibility" / "benchmarks" / "run_method_benchmark.py"),
            "--out-dir",
            str(out_dir),
            "--n-replicates",
            "5",
        ],
        cwd=REPO_ROOT,
        check=True,
    )

    metrics_path = out_dir / "method_benchmark_metrics.csv"
    summary_path = out_dir / "benchmark_summary.csv"
    manifest_path = out_dir / "benchmark_manifest.json"
    figure_path = out_dir / "figure4_benchmark_matrix.png"
    for path in (metrics_path, summary_path, manifest_path, figure_path):
        assert path.exists()
        assert path.stat().st_size > 0

    metrics = pd.read_csv(metrics_path)
    summary = pd.read_csv(summary_path)
    assert set(metrics["method"]) == {
        "naive_stack",
        "histoseg_hard_only",
        "histoseg_full",
        "paste",
        "paste2",
        "gpsa",
        "staligner",
        "spacel_scube",
    }
    mean_iou = summary.set_index("method")["union_iou"]
    assert mean_iou["histoseg_full"] > mean_iou["naive_stack"]
    assert "synthetic_proxy" in set(metrics["execution_status"])
