"""Render manuscript-ready schematic and evidence-summary panels.

This script creates lightweight figure panels that do not require private raw
data. Panels derived from ``polyp_24_gene_evidence.tsv`` are suitable as figure
source-data views; morphology overlays are emitted as templates until the
public H&E/morphology images are added to the data bundle.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EVIDENCE = REPO_ROOT / "docs" / "manuscripts" / "polyp_24_gene_evidence.tsv"
DEFAULT_OUT_ROOT = REPO_ROOT / "docs" / "_static"


STRUCTURE_COLORS = {
    "Structure 1": "#159895",
    "Structure 2": "#3f7cac",
    "Structure 3": "#e0a21a",
    "Structure 4": "#8b6bb8",
    "Structure 5": "#c84c61",
}


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    evidence = _load_evidence(Path(args.evidence_tsv).expanduser())
    out_root = Path(args.out_root).expanduser()
    manuscript_dir = out_root / "manuscript"
    polyp_dir = out_root / "threed" / "polyp" / "24gene"
    manuscript_dir.mkdir(parents=True, exist_ok=True)
    polyp_dir.mkdir(parents=True, exist_ok=True)

    outputs = [
        *_render_workflow_schematic(manuscript_dir),
        *_render_histology_overlay_template(manuscript_dir),
        _render_signed_distance_violin(evidence, polyp_dir),
        *_render_compartment_summary(evidence, polyp_dir),
    ]
    for path in outputs:
        print(path)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence-tsv", default=str(DEFAULT_EVIDENCE))
    parser.add_argument("--out-root", default=str(DEFAULT_OUT_ROOT))
    return parser


def _load_evidence(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t")
    required = {
        "gene",
        "histoseg_dominant_structure",
        "top05_fraction_inside",
        "signed_distance_um",
        "nonzero_cell_count",
        "marker_class",
        "cell_tissue_interpretation",
    }
    missing = sorted(required.difference(table.columns))
    if missing:
        raise ValueError(f"Evidence table is missing required column(s): {', '.join(missing)}")
    table = table.copy()
    table["top05_fraction_inside"] = pd.to_numeric(table["top05_fraction_inside"], errors="coerce")
    table["signed_distance_um"] = pd.to_numeric(table["signed_distance_um"], errors="coerce")
    table["marker_group"] = table["marker_class"].map(_marker_group)
    return table


def _marker_group(marker_class: str) -> str:
    text = str(marker_class).lower()
    if any(token in text for token in ("immune", "t-cell", "b-cell", "macrophage", "myeloid", "leukocyte")):
        return "Immune"
    if any(token in text for token in ("epithelial", "goblet", "mucin", "stem-like", "adhesion")):
        return "Epithelial"
    if any(token in text for token in ("endothelial", "pericyte", "vascular")):
        return "Vascular"
    if any(token in text for token in ("stromal", "collagen", "fibroblast", "contractile", "matrix", "mesenchymal")):
        return "Stromal"
    return "Other"


def _render_workflow_schematic(out_dir: Path) -> list[Path]:
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    out_png = out_dir / "figure1_workflow_schematic.png"
    out_svg = out_dir / "figure1_workflow_schematic.svg"
    fig, ax = plt.subplots(figsize=(11.5, 4.6))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    steps = [
        ("A", "Serial Xenium\nslices", "GeoJSON contours\ncell coordinates"),
        ("B", "Topology-aware\nalignment", "hard contours\nsoft TPS/RBF correction"),
        ("C", "Voxelized 3D\nstructures", "SDF volumes\nmesh export"),
        ("D", "Gene-structure\nmetrics", "fraction inside\nsigned distance"),
    ]
    x_positions = [0.08, 0.33, 0.58, 0.83]
    box_w, box_h = 0.19, 0.42
    for index, ((panel, title, body), x) in enumerate(zip(steps, x_positions)):
        y = 0.48
        box = FancyBboxPatch(
            (x - box_w / 2, y - box_h / 2),
            box_w,
            box_h,
            boxstyle="round,pad=0.012,rounding_size=0.018",
            linewidth=1.1,
            edgecolor="#253238",
            facecolor="#f7f9fb",
        )
        ax.add_patch(box)
        ax.text(x - box_w / 2 + 0.015, y + box_h / 2 - 0.055, panel, weight="bold", fontsize=12)
        ax.text(x, y + 0.08, title, ha="center", va="center", fontsize=9.6, weight="bold", linespacing=1.15)
        ax.text(x, y - 0.08, body, ha="center", va="center", fontsize=9.0, linespacing=1.25)
        _draw_icon(ax, x, y - 0.185, index)
    for left, right in zip(x_positions[:-1], x_positions[1:]):
        ax.add_patch(
            FancyArrowPatch(
                (left + box_w / 2 + 0.01, 0.48),
                (right - box_w / 2 - 0.01, 0.48),
                arrowstyle="-|>",
                mutation_scale=14,
                linewidth=1.4,
                color="#455a64",
            )
        )
    ax.text(0.5, 0.92, "HistoSeg 2D-to-3D spatial morphology workflow", ha="center", fontsize=14, weight="bold")
    ax.text(
        0.5,
        0.08,
        "Benchmark and sensitivity panels report geometry, topology, expression consistency, runtime and uncertainty.",
        ha="center",
        fontsize=9.5,
        color="#455a64",
    )
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)
    return [out_png, out_svg]


def _draw_icon(ax, x: float, y: float, index: int) -> None:
    theta = np.linspace(0, 2 * np.pi, 80)
    if index == 0:
        for offset in (-0.03, 0.0, 0.03):
            ax.plot(x + 0.045 * np.cos(theta), y + offset + 0.018 * np.sin(theta), color="#3f7cac", lw=1)
    elif index == 1:
        ax.plot([x - 0.05, x, x + 0.05], [y - 0.025, y + 0.025, y - 0.015], color="#159895", lw=2)
        ax.scatter([x - 0.05, x, x + 0.05], [y - 0.025, y + 0.025, y - 0.015], s=16, color="#c84c61")
    elif index == 2:
        for radius, color in ((0.055, "#159895"), (0.037, "#e0a21a"), (0.020, "#c84c61")):
            ax.plot(x + radius * np.cos(theta), y + 0.65 * radius * np.sin(theta), color=color, lw=1.5)
    else:
        ax.bar(
            [x - 0.045, x - 0.015, x + 0.015, x + 0.045],
            [0.05, 0.025, 0.04, 0.065],
            bottom=y - 0.035,
            width=0.018,
            color="#8b6bb8",
        )
        ax.plot([x - 0.055, x + 0.055], [y - 0.035, y - 0.035], color="#455a64", lw=1)


def _render_histology_overlay_template(out_dir: Path) -> list[Path]:
    import matplotlib.pyplot as plt

    out_png = out_dir / "figure1_histology_contour_overlay_template.png"
    out_svg = out_dir / "figure1_histology_contour_overlay_template.svg"
    fig, ax = plt.subplots(figsize=(6.5, 5.2))
    ax.set_facecolor("#f4e7eb")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 8)
    ax.set_xticks([])
    ax.set_yticks([])
    rng = np.random.default_rng(20260504)
    for _ in range(90):
        cx, cy = rng.uniform(0.4, 9.6), rng.uniform(0.4, 7.6)
        radius = rng.uniform(0.035, 0.085)
        ax.add_patch(plt.Circle((cx, cy), radius, color="#b77aa4", alpha=0.22, lw=0))
    theta = np.linspace(0, 2 * np.pi, 240)
    centers = [(2.0, 5.9), (5.0, 5.6), (7.4, 5.5), (3.8, 2.7), (6.8, 2.4)]
    scales = [(1.6, 0.72), (1.9, 0.78), (1.35, 0.60), (1.8, 0.84), (2.3, 0.95)]
    for index, ((cx, cy), (sx, sy)) in enumerate(zip(centers, scales), start=1):
        wiggle = 1 + 0.08 * np.sin(3 * theta + index)
        x = cx + sx * wiggle * np.cos(theta)
        y = cy + sy * wiggle * np.sin(theta)
        structure = f"Structure {index}"
        ax.plot(x, y, color=STRUCTURE_COLORS[structure], lw=2.2, label=structure)
        ax.fill(x, y, color=STRUCTURE_COLORS[structure], alpha=0.10)
    ax.legend(loc="upper right", frameon=True, fontsize=8)
    ax.set_title("Morphology/contour overlay template", fontsize=12, weight="bold")
    ax.text(
        0.02,
        0.02,
        "Replace with public H&E or morphology image before submission.",
        transform=ax.transAxes,
        fontsize=8.5,
        color="#253238",
    )
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)
    return [out_png, out_svg]


def _render_signed_distance_violin(evidence: pd.DataFrame, out_dir: Path) -> Path:
    import matplotlib.pyplot as plt

    out_png = out_dir / "signed_distance_marker_group_violin.png"
    groups = ["Stromal", "Epithelial", "Immune", "Vascular", "Other"]
    data = [
        evidence.loc[evidence["marker_group"] == group, "signed_distance_um"].dropna().to_numpy()
        for group in groups
    ]
    keep = [index for index, values in enumerate(data) if len(values)]
    data = [data[index] for index in keep]
    labels = [groups[index] for index in keep]
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    parts = ax.violinplot(data, showmeans=False, showmedians=True, widths=0.72)
    for body in parts["bodies"]:
        body.set_facecolor("#6b8f71")
        body.set_edgecolor("#253238")
        body.set_alpha(0.40)
    parts["cmedians"].set_color("#253238")
    parts["cmedians"].set_linewidth(1.4)
    rng = np.random.default_rng(20260504)
    for position, values in enumerate(data, start=1):
        jitter = rng.uniform(-0.08, 0.08, size=len(values))
        ax.scatter(position + jitter, values, s=32, color="#253238", alpha=0.78, zorder=3)
    ax.axhline(0, color="#c84c61", lw=1.2, linestyle="--")
    ax.set_xticks(np.arange(1, len(labels) + 1), labels)
    ax.set_ylabel("Signed distance to dominant structure (um)")
    ax.set_title("Figure 3E marker-group signed-distance distributions", fontsize=12, weight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#d9dee3", linewidth=0.8)
    ax.set_axisbelow(True)
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return out_png


def _render_compartment_summary(evidence: pd.DataFrame, out_dir: Path) -> list[Path]:
    import matplotlib.pyplot as plt

    summary = (
        evidence.groupby(["marker_group", "histoseg_dominant_structure"], dropna=False)
        .agg(
            genes=("gene", lambda values: ", ".join(values)),
            n_genes=("gene", "size"),
            median_fraction_inside=("top05_fraction_inside", "median"),
            median_signed_distance_um=("signed_distance_um", "median"),
            total_nonzero_cells=("nonzero_cell_count", "sum"),
        )
        .reset_index()
        .sort_values(["marker_group", "histoseg_dominant_structure"])
    )
    csv_path = out_dir / "compartment_summary_table.csv"
    png_path = out_dir / "compartment_summary_table.png"
    summary.to_csv(csv_path, index=False)

    display = summary.copy()
    display["median_fraction_inside"] = display["median_fraction_inside"].map(lambda value: f"{value:.2f}")
    display["median_signed_distance_um"] = display["median_signed_distance_um"].map(lambda value: f"{value:.1f}")
    display["genes"] = display["genes"].map(_shorten_gene_list)
    display = display[
        [
            "marker_group",
            "histoseg_dominant_structure",
            "n_genes",
            "median_fraction_inside",
            "median_signed_distance_um",
            "genes",
        ]
    ]
    fig, ax = plt.subplots(figsize=(10.8, 0.62 * len(display) + 1.4))
    ax.axis("off")
    table = ax.table(
        cellText=display.values,
        colLabels=["Group", "Dominant structure", "n", "Median inside", "Median SDF", "Genes"],
        loc="center",
        cellLoc="left",
        colLoc="left",
        colWidths=[0.13, 0.17, 0.06, 0.12, 0.12, 0.40],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1.0, 1.35)
    for (row, _col), cell in table.get_celld().items():
        cell.set_edgecolor("#d9dee3")
        if row == 0:
            cell.set_facecolor("#253238")
            cell.get_text().set_color("white")
            cell.get_text().set_weight("bold")
        else:
            cell.set_facecolor("#ffffff" if row % 2 else "#f7f9fb")
    ax.set_title("Compartment-level marker summary", fontsize=12, weight="bold", pad=10)
    fig.savefig(png_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return [csv_path, png_path]


def _shorten_gene_list(value: str, max_chars: int = 58) -> str:
    if len(value) <= max_chars:
        return value
    return value[: max_chars - 3].rstrip(", ") + "..."


if __name__ == "__main__":
    raise SystemExit(main())
