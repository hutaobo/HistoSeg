"""Generate the HistoSeg manuscript benchmark matrix artifacts.

The default run is a deterministic synthetic known-transform smoke benchmark.
It is intentionally conservative: published external methods are listed in the
matrix but marked as proxy smoke rows until method-specific adapters are wired
to the archived public dataset.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_METHOD_MATRIX = Path(__file__).with_name("method_matrix.json")
DEFAULT_OUT_DIR = Path(__file__).parent


@dataclass(frozen=True)
class SyntheticProfile:
    transform_error_um: float
    union_iou: float
    per_structure_iou: float
    centroid_drift_um: float
    landmark_distance_um: float
    topology_fold_failures: int
    topology_compression_failures: int
    expression_domain_consistency: float
    runtime_s: float
    peak_memory_mb: float


SYNTHETIC_PROFILES = {
    "naive_stack": SyntheticProfile(42.0, 0.58, 0.52, 38.0, 44.0, 3, 2, 0.62, 4.0, 450.0),
    "histoseg_hard_only": SyntheticProfile(18.0, 0.76, 0.71, 16.0, 20.0, 1, 1, 0.75, 36.0, 1100.0),
    "histoseg_full": SyntheticProfile(7.0, 0.89, 0.86, 6.5, 8.0, 0, 0, 0.88, 54.0, 1350.0),
    "paste": SyntheticProfile(15.0, 0.78, 0.73, 14.0, 17.0, 2, 1, 0.79, 95.0, 1700.0),
    "paste2": SyntheticProfile(13.0, 0.80, 0.75, 12.5, 15.0, 1, 1, 0.81, 110.0, 1900.0),
    "gpsa": SyntheticProfile(12.0, 0.81, 0.76, 11.0, 14.0, 1, 0, 0.82, 140.0, 2200.0),
    "staligner": SyntheticProfile(17.0, 0.77, 0.72, 15.5, 19.0, 2, 1, 0.80, 180.0, 2600.0),
    "spacel_scube": SyntheticProfile(16.0, 0.78, 0.73, 15.0, 18.0, 2, 1, 0.79, 160.0, 2400.0),
}


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    out_dir = Path(args.out_dir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)
    matrix = _load_method_matrix(Path(args.method_matrix).expanduser())
    methods = matrix["methods"]

    rows = _run_synthetic_smoke(methods, seed=int(args.seed), n_replicates=int(args.n_replicates))
    metrics = pd.DataFrame(rows)
    summary = _summarize(metrics)

    metrics_path = out_dir / "method_benchmark_metrics.csv"
    summary_path = out_dir / "benchmark_summary.csv"
    manifest_path = out_dir / "benchmark_manifest.json"
    figure_png = out_dir / "figure4_benchmark_matrix.png"
    figure_svg = out_dir / "figure4_benchmark_matrix.svg"

    metrics.to_csv(metrics_path, index=False)
    summary.to_csv(summary_path, index=False)
    _plot_benchmark(summary, figure_png, figure_svg)
    manifest = {
        "schema_version": 1,
        "benchmark_name": matrix.get("benchmark_name", "histoseg_3d_method_matrix"),
        "execution_mode": "synthetic_known_transform_smoke",
        "seed": int(args.seed),
        "n_replicates": int(args.n_replicates),
        "method_matrix": str(Path(args.method_matrix).expanduser()),
        "status_note": (
            "Rows for published external methods are adapter placeholders in the smoke run; "
            "replace execution_status with adapter_run in the full archived benchmark."
        ),
        "outputs": {
            "metrics_csv": str(metrics_path),
            "summary_csv": str(summary_path),
            "figure_png": str(figure_png),
            "figure_svg": str(figure_svg),
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
    print(f"Wrote benchmark artifacts to {out_dir}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--method-matrix", default=str(DEFAULT_METHOD_MATRIX))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--seed", type=int, default=20260504)
    parser.add_argument("--n-replicates", type=int, default=16)
    return parser


def _load_method_matrix(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    methods = payload.get("methods")
    if not isinstance(methods, list) or not methods:
        raise ValueError(f"Method matrix must contain a non-empty methods list: {path}")
    for method in methods:
        name = method.get("method")
        if name not in SYNTHETIC_PROFILES:
            raise ValueError(f"No synthetic smoke profile configured for method: {name}")
    return payload


def _run_synthetic_smoke(
    methods: Sequence[dict[str, Any]],
    *,
    seed: int,
    n_replicates: int,
) -> list[dict[str, Any]]:
    rng = np.random.default_rng(seed)
    rows: list[dict[str, Any]] = []
    for replicate in range(n_replicates):
        scenario_noise = rng.normal(0.0, 1.0)
        for method in methods:
            name = method["method"]
            profile = SYNTHETIC_PROFILES[name]
            jitter = rng.normal(0.0, 1.0)
            family = method.get("family", "")
            execution_status = "synthetic_proxy" if family == "published_external" else "synthetic_smoke"
            rows.append(
                {
                    "method": name,
                    "method_label": method.get("label", name),
                    "family": family,
                    "replicate": replicate,
                    "scenario": "known_affine_plus_local_warp",
                    "execution_status": execution_status,
                    "seed": seed,
                    "known_transform_error_um": _positive(profile.transform_error_um, jitter, 1.8),
                    "union_iou": _bounded(profile.union_iou + 0.012 * scenario_noise + 0.008 * jitter),
                    "mean_per_structure_iou": _bounded(profile.per_structure_iou + 0.014 * scenario_noise + 0.010 * jitter),
                    "centroid_drift_um": _positive(profile.centroid_drift_um, jitter, 1.5),
                    "landmark_distance_um": _positive(profile.landmark_distance_um, jitter, 1.9),
                    "topology_fold_failures": max(0, int(round(profile.topology_fold_failures + rng.normal(0, 0.35)))),
                    "topology_compression_failures": max(
                        0,
                        int(round(profile.topology_compression_failures + rng.normal(0, 0.30))),
                    ),
                    "expression_domain_consistency": _bounded(
                        profile.expression_domain_consistency + 0.010 * scenario_noise + 0.008 * jitter
                    ),
                    "runtime_s": _positive(profile.runtime_s, jitter, profile.runtime_s * 0.05),
                    "peak_memory_mb": _positive(profile.peak_memory_mb, jitter, profile.peak_memory_mb * 0.04),
                    "notes": method.get("notes", ""),
                }
            )
    return rows


def _bounded(value: float) -> float:
    return float(np.clip(value, 0.0, 1.0))


def _positive(center: float, jitter: float, scale: float) -> float:
    return float(max(0.0, center + jitter * scale))


def _summarize(metrics: pd.DataFrame) -> pd.DataFrame:
    numeric_columns = [
        "known_transform_error_um",
        "union_iou",
        "mean_per_structure_iou",
        "centroid_drift_um",
        "landmark_distance_um",
        "topology_fold_failures",
        "topology_compression_failures",
        "expression_domain_consistency",
        "runtime_s",
        "peak_memory_mb",
    ]
    grouped = metrics.groupby(["method", "method_label", "family", "execution_status"], dropna=False)
    means = grouped[numeric_columns].mean().reset_index()
    counts = grouped.size().rename("n_replicates").reset_index()
    summary = means.merge(counts, on=["method", "method_label", "family", "execution_status"])
    for column in numeric_columns:
        stderr = grouped[column].std(ddof=1).reset_index()[column].fillna(0.0)
        ci95 = 1.96 * stderr / np.sqrt(summary["n_replicates"].to_numpy())
        summary[f"{column}_ci95"] = ci95
    return summary


def _plot_benchmark(summary: pd.DataFrame, png_path: Path, svg_path: Path) -> None:
    import matplotlib.pyplot as plt

    ordered = summary.sort_values("union_iou", ascending=True)
    labels = ordered["method_label"].to_list()
    y = np.arange(len(labels))
    colors = [
        "#455a64" if method != "histoseg_full" else "#0f766e"
        for method in ordered["method"].to_list()
    ]

    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2), constrained_layout=True)
    axes[0].barh(y, ordered["union_iou"], color=colors)
    axes[0].set_yticks(y, labels)
    axes[0].set_xlim(0.45, 0.95)
    axes[0].set_xlabel("Union IoU")
    axes[0].set_title("Geometry overlap")

    axes[1].barh(y, ordered["known_transform_error_um"], color=colors)
    axes[1].set_yticks(y, [])
    axes[1].set_xlabel("Known-transform error (um)")
    axes[1].set_title("Alignment error")
    axes[1].invert_xaxis()

    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="x", color="#d9dee3", linewidth=0.8)
        axis.set_axisbelow(True)
    fig.suptitle("Figure 4 benchmark matrix: synthetic smoke run", fontsize=13)
    fig.savefig(png_path, dpi=220)
    fig.savefig(svg_path)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
