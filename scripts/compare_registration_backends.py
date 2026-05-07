from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial.distance import jensenshannon
from scipy.stats import kendalltau, pearsonr, spearmanr, wasserstein_distance


MATRIX_FILES = {
    "fraction_inside": "gene_structure_fraction_inside_matrix.csv",
    "overlap_fraction": "gene_structure_overlap_fraction_matrix.csv",
    "signed_distance": "gene_structure_signed_distance_matrix.csv",
}


@dataclass(frozen=True)
class BackendComparisonConfig:
    histoseg_batch_dir: str | Path
    spacemap_batch_dir: str | Path
    out_dir: str | Path
    hotspot: str = "top05"
    genes: Sequence[str] = ()
    structures: Sequence[str] = ()
    focus_structure: str = "Structure 5"


@dataclass(frozen=True)
class BackendComparisonResult:
    out_dir: Path
    summary_json: Path
    metric_correlations_csv: Path
    sdf_distribution_metrics_csv: Path
    rank_stability_csv: Path
    delta_heatmaps_csv: Path
    delta_heatmap_png: Path
    sdf_ridge_plot_png: Path


def compare_registration_backends(cfg: BackendComparisonConfig) -> BackendComparisonResult:
    """Compare HistoSeg and Space-map coordinate backends using spatial-module matrices."""

    histoseg_dir = Path(cfg.histoseg_batch_dir).expanduser()
    spacemap_dir = Path(cfg.spacemap_batch_dir).expanduser()
    out_dir = Path(cfg.out_dir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)

    metric_rows: list[dict] = []
    delta_rows: list[dict] = []
    matrices: dict[str, tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]] = {}

    for metric, filename in MATRIX_FILES.items():
        hs_path = histoseg_dir / filename
        sm_path = spacemap_dir / filename
        if not hs_path.exists() or not sm_path.exists():
            continue
        hs = _load_hotspot_matrix(hs_path, cfg.hotspot, cfg.genes, cfg.structures)
        sm = _load_hotspot_matrix(sm_path, cfg.hotspot, cfg.genes, cfg.structures)
        hs, sm = _align_matrices(hs, sm)
        delta = sm - hs
        matrices[metric] = (hs, sm, delta)
        delta.to_csv(out_dir / f"backend_delta_{metric}_matrix.csv")

        metric_rows.extend(_matrix_similarity_rows(metric, hs, sm))
        for gene, row in delta.iterrows():
            for structure, value in row.items():
                delta_rows.append(
                    {
                        "metric": metric,
                        "hotspot": cfg.hotspot,
                        "gene": gene,
                        "structure": structure,
                        "histoseg": _finite_or_nan(hs.loc[gene, structure]),
                        "spacemap": _finite_or_nan(sm.loc[gene, structure]),
                        "delta_spacemap_minus_histoseg": _finite_or_nan(value),
                    }
                )

    if not matrices:
        raise ValueError("No matching spatial-module matrices were found in the two backend directories.")

    metric_correlations_csv = out_dir / "backend_metric_correlations.csv"
    pd.DataFrame(metric_rows).to_csv(metric_correlations_csv, index=False)

    delta_heatmaps_csv = out_dir / "backend_delta_heatmaps.csv"
    pd.DataFrame(delta_rows).to_csv(delta_heatmaps_csv, index=False)

    signed_tuple = matrices.get("signed_distance")
    if signed_tuple is None:
        rank_rows: list[dict] = []
        sdf_rows: list[dict] = []
    else:
        hs_signed, sm_signed, delta_signed = signed_tuple
        rank_rows = _rank_stability_rows(hs_signed, sm_signed)
        sdf_rows = _sdf_profile_rows(cfg.hotspot, hs_signed, sm_signed, delta_signed)

    rank_stability_csv = out_dir / "backend_rank_stability.csv"
    pd.DataFrame(rank_rows).to_csv(rank_stability_csv, index=False)

    sdf_distribution_metrics_csv = out_dir / "backend_sdf_distribution_metrics.csv"
    pd.DataFrame(sdf_rows).to_csv(sdf_distribution_metrics_csv, index=False)

    delta_heatmap_png = out_dir / "backend_delta_heatmap.png"
    _plot_delta_heatmap(matrices, delta_heatmap_png)

    sdf_ridge_plot_png = out_dir / "backend_sdf_ridge_plot.png"
    _plot_sdf_backend_profile(matrices.get("signed_distance"), cfg.focus_structure, sdf_ridge_plot_png)

    summary_json = out_dir / "backend_comparison_summary.json"
    _write_json(
        summary_json,
        {
            "histoseg_batch_dir": str(histoseg_dir),
            "spacemap_batch_dir": str(spacemap_dir),
            "out_dir": str(out_dir),
            "hotspot": cfg.hotspot,
            "genes": list(cfg.genes),
            "structures": list(cfg.structures),
            "focus_structure": cfg.focus_structure,
            "metrics_compared": sorted(matrices),
            "outputs": {
                "metric_correlations_csv": str(metric_correlations_csv),
                "sdf_distribution_metrics_csv": str(sdf_distribution_metrics_csv),
                "rank_stability_csv": str(rank_stability_csv),
                "delta_heatmaps_csv": str(delta_heatmaps_csv),
                "delta_heatmap_png": str(delta_heatmap_png),
                "sdf_ridge_plot_png": str(sdf_ridge_plot_png),
            },
        },
    )

    return BackendComparisonResult(
        out_dir=out_dir,
        summary_json=summary_json,
        metric_correlations_csv=metric_correlations_csv,
        sdf_distribution_metrics_csv=sdf_distribution_metrics_csv,
        rank_stability_csv=rank_stability_csv,
        delta_heatmaps_csv=delta_heatmaps_csv,
        delta_heatmap_png=delta_heatmap_png,
        sdf_ridge_plot_png=sdf_ridge_plot_png,
    )


def _load_hotspot_matrix(
    path: Path,
    hotspot: str,
    genes: Sequence[str],
    structures: Sequence[str],
) -> pd.DataFrame:
    matrix = pd.read_csv(path, index_col=0)
    columns = [column for column in matrix.columns if str(column).startswith(f"{hotspot}|")]
    if not columns:
        raise ValueError(f"No {hotspot!r} columns found in {path}.")
    matrix = matrix.loc[:, columns].copy()
    matrix.columns = [str(column).split("|", 1)[1] for column in matrix.columns]
    if genes:
        selected = [gene for gene in genes if gene in matrix.index]
        matrix = matrix.loc[selected]
    if structures:
        selected_cols = [structure for structure in structures if structure in matrix.columns]
        matrix = matrix.loc[:, selected_cols]
    return matrix.apply(pd.to_numeric, errors="coerce")


def _align_matrices(left: pd.DataFrame, right: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    genes = left.index.intersection(right.index)
    structures = left.columns.intersection(right.columns)
    if genes.empty or structures.empty:
        raise ValueError("Backend matrices do not share any genes or structures.")
    return left.loc[genes, structures].sort_index().sort_index(axis=1), right.loc[genes, structures].sort_index().sort_index(axis=1)


def _matrix_similarity_rows(metric: str, histoseg: pd.DataFrame, spacemap: pd.DataFrame) -> list[dict]:
    rows = [_similarity_row(metric, "all_values", histoseg.to_numpy().ravel(), spacemap.to_numpy().ravel())]
    for structure in histoseg.columns:
        rows.append(_similarity_row(metric, f"structure:{structure}", histoseg[structure], spacemap[structure]))
    for gene in histoseg.index:
        rows.append(_similarity_row(metric, f"gene:{gene}", histoseg.loc[gene], spacemap.loc[gene]))
    return rows


def _similarity_row(metric: str, scope: str, histoseg_values, spacemap_values) -> dict:
    hs, sm = _finite_pair(histoseg_values, spacemap_values)
    return {
        "metric": metric,
        "scope": scope,
        "n_values": int(len(hs)),
        "pearson_r": _safe_corr(pearsonr, hs, sm),
        "spearman_r": _safe_corr(spearmanr, hs, sm),
        "wasserstein_distance": _safe_wasserstein(hs, sm),
        "jensen_shannon_distance": _safe_js(hs, sm),
        "mean_abs_delta": float(np.mean(np.abs(sm - hs))) if len(hs) else math.nan,
        "max_abs_delta": float(np.max(np.abs(sm - hs))) if len(hs) else math.nan,
    }


def _rank_stability_rows(histoseg: pd.DataFrame, spacemap: pd.DataFrame) -> list[dict]:
    rows: list[dict] = []
    for structure in histoseg.columns:
        hs, sm = _finite_pair(histoseg[structure], spacemap[structure])
        rows.append(
            {
                "axis": "genes_ranked_by_structure",
                "structure": structure,
                "gene": None,
                "n_values": int(len(hs)),
                "spearman_r": _safe_corr(spearmanr, hs, sm),
                "kendall_tau": _safe_corr(kendalltau, hs, sm),
            }
        )
    for gene in histoseg.index:
        hs, sm = _finite_pair(histoseg.loc[gene], spacemap.loc[gene])
        rows.append(
            {
                "axis": "structures_ranked_by_gene",
                "structure": None,
                "gene": gene,
                "n_values": int(len(hs)),
                "spearman_r": _safe_corr(spearmanr, hs, sm),
                "kendall_tau": _safe_corr(kendalltau, hs, sm),
            }
        )
    return rows


def _sdf_profile_rows(hotspot: str, histoseg: pd.DataFrame, spacemap: pd.DataFrame, delta: pd.DataFrame) -> list[dict]:
    rows: list[dict] = []
    for gene in histoseg.index:
        for structure in histoseg.columns:
            rows.append(
                {
                    "hotspot": hotspot,
                    "gene": gene,
                    "structure": structure,
                    "histoseg_median_signed_distance_um": _finite_or_nan(histoseg.loc[gene, structure]),
                    "spacemap_median_signed_distance_um": _finite_or_nan(spacemap.loc[gene, structure]),
                    "delta_spacemap_minus_histoseg_um": _finite_or_nan(delta.loc[gene, structure]),
                }
            )
    return rows


def _plot_delta_heatmap(
    matrices: dict[str, tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]],
    out_png: Path,
) -> None:
    metrics = [metric for metric in ("fraction_inside", "signed_distance") if metric in matrices]
    if not metrics:
        metrics = list(matrices)[:1]
    fig, axes = plt.subplots(1, len(metrics), figsize=(7 * len(metrics), 5), squeeze=False)
    for ax, metric in zip(axes[0], metrics):
        delta = matrices[metric][2]
        cmap = "vlag" if metric == "signed_distance" else "coolwarm"
        sns.heatmap(delta, ax=ax, cmap=cmap, center=0.0, cbar_kws={"label": "Space-map - HistoSeg"})
        ax.set_title(metric.replace("_", " "))
        ax.set_xlabel("Structure")
        ax.set_ylabel("Gene")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=220)
    plt.close(fig)


def _plot_sdf_backend_profile(
    signed_tuple: tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame] | None,
    focus_structure: str,
    out_png: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    if signed_tuple is None:
        ax.text(0.5, 0.5, "No signed-distance matrix found", ha="center", va="center")
        ax.axis("off")
    else:
        histoseg, spacemap, _ = signed_tuple
        if focus_structure not in histoseg.columns:
            focus_structure = str(histoseg.columns[0])
        frame = pd.DataFrame(
            {
                "gene": histoseg.index,
                "HistoSeg TPS": histoseg[focus_structure].to_numpy(dtype=float),
                "Space-map LDDMM": spacemap[focus_structure].to_numpy(dtype=float),
            }
        ).dropna()
        frame = frame.assign(sort_value=frame["HistoSeg TPS"]).sort_values("sort_value")
        y = np.arange(len(frame))
        ax.axvline(0.0, color="0.55", linestyle="--", linewidth=1)
        ax.hlines(y, frame["HistoSeg TPS"], frame["Space-map LDDMM"], color="0.75", linewidth=1)
        ax.scatter(frame["HistoSeg TPS"], y, label="HistoSeg TPS", s=42)
        ax.scatter(frame["Space-map LDDMM"], y, label="Space-map LDDMM", marker="^", s=42)
        ax.set_yticks(y)
        ax.set_yticklabels(frame["gene"])
        ax.set_xlabel(f"Median signed distance to {focus_structure} (µm)")
        ax.set_title("Backend SDF profile stability")
        ax.legend(frameon=False)
        ax.grid(axis="x", alpha=0.2)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=220)
    plt.close(fig)


def _finite_pair(left, right) -> tuple[np.ndarray, np.ndarray]:
    lhs = np.asarray(left, dtype=float)
    rhs = np.asarray(right, dtype=float)
    mask = np.isfinite(lhs) & np.isfinite(rhs)
    return lhs[mask], rhs[mask]


def _safe_corr(func, left: np.ndarray, right: np.ndarray) -> float:
    if len(left) < 2 or np.nanstd(left) == 0 or np.nanstd(right) == 0:
        return math.nan
    result = func(left, right)
    value = result.statistic if hasattr(result, "statistic") else result[0]
    return float(value)


def _safe_wasserstein(left: np.ndarray, right: np.ndarray) -> float:
    if len(left) == 0:
        return math.nan
    return float(wasserstein_distance(left, right))


def _safe_js(left: np.ndarray, right: np.ndarray, bins: int = 24) -> float:
    if len(left) < 2:
        return math.nan
    combined = np.concatenate([left, right])
    if np.nanmin(combined) == np.nanmax(combined):
        return 0.0
    edges = np.linspace(float(np.nanmin(combined)), float(np.nanmax(combined)), bins + 1)
    left_hist, _ = np.histogram(left, bins=edges)
    right_hist, _ = np.histogram(right, bins=edges)
    left_prob = left_hist.astype(float) + 1e-12
    right_prob = right_hist.astype(float) + 1e-12
    left_prob /= left_prob.sum()
    right_prob /= right_prob.sum()
    return float(jensenshannon(left_prob, right_prob))


def _finite_or_nan(value) -> float:
    number = float(value)
    return number if np.isfinite(number) else math.nan


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=str), encoding="utf-8")


def _parse_list(values: Sequence[str]) -> tuple[str, ...]:
    parsed: list[str] = []
    for value in values:
        parsed.extend(part.strip() for part in str(value).split(",") if part.strip())
    return tuple(parsed)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Compare HistoSeg and Space-map 3D spatial module outputs.")
    parser.add_argument("--histoseg-batch-dir", required=True)
    parser.add_argument("--spacemap-batch-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--hotspot", default="top05", choices=["top15", "top10", "top05"])
    parser.add_argument("--genes", nargs="*", default=())
    parser.add_argument("--structures", nargs="*", default=())
    parser.add_argument("--focus-structure", default="Structure 5")
    args = parser.parse_args(argv)
    result = compare_registration_backends(
        BackendComparisonConfig(
            histoseg_batch_dir=args.histoseg_batch_dir,
            spacemap_batch_dir=args.spacemap_batch_dir,
            out_dir=args.out_dir,
            hotspot=args.hotspot,
            genes=_parse_list(args.genes),
            structures=_parse_list(args.structures),
            focus_structure=args.focus_structure,
        )
    )
    print(json.dumps(asdict(result), indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
