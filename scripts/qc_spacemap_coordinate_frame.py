from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


@dataclass(frozen=True)
class SpaceMapCoordinateQcConfig:
    spacemap_aligned_cells_parquet: str | Path
    histoseg_aligned_cells_parquet: str | Path
    stack_root: str | Path
    out_dir: str | Path
    sample_column: str = "sample_id"
    slice_order_column: str = "slice_order"
    x_column: str = "x_3d_um"
    y_column: str = "y_3d_um"
    z_column: str = "z_um"
    overlay_slices: Sequence[int] = ()
    max_points_per_overlay: int = 80000
    random_state: int = 0
    contour_csv: str | Path | None = None
    scale_ratio_ok_min: float = 0.95
    scale_ratio_ok_max: float = 1.05


@dataclass(frozen=True)
class SpaceMapCoordinateQcResult:
    out_dir: Path
    summary_json: Path
    by_slice_csv: Path
    bbox_overlay_png: Path
    centroid_delta_png: Path
    overlay_pngs: tuple[Path, ...]


def run_spacemap_coordinate_frame_qc(cfg: SpaceMapCoordinateQcConfig) -> SpaceMapCoordinateQcResult:
    """Audit Space-map aligned cells against the HistoSeg coordinate frame without modifying coordinates."""

    out_dir = Path(cfg.out_dir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)
    hs = pd.read_parquet(Path(cfg.histoseg_aligned_cells_parquet).expanduser())
    sm = pd.read_parquet(Path(cfg.spacemap_aligned_cells_parquet).expanduser())
    _require_columns(hs, [cfg.sample_column, cfg.slice_order_column, cfg.x_column, cfg.y_column, cfg.z_column])
    _require_columns(sm, [cfg.sample_column, cfg.slice_order_column, cfg.x_column, cfg.y_column, cfg.z_column])

    by_slice = _coordinate_qc_by_slice(hs, sm, cfg)
    by_slice_csv = out_dir / "coordinate_frame_qc_by_slice.csv"
    by_slice.to_csv(by_slice_csv, index=False)

    global_summary = _global_coordinate_summary(hs, sm, by_slice, cfg)
    summary_json = out_dir / "coordinate_frame_qc_summary.json"
    _write_json(
        summary_json,
        {
            "histoseg_aligned_cells_parquet": str(Path(cfg.histoseg_aligned_cells_parquet).expanduser()),
            "spacemap_aligned_cells_parquet": str(Path(cfg.spacemap_aligned_cells_parquet).expanduser()),
            "stack_root": str(Path(cfg.stack_root).expanduser()),
            "out_dir": str(out_dir),
            "sample_column": cfg.sample_column,
            "slice_order_column": cfg.slice_order_column,
            "coordinate_columns": {"x": cfg.x_column, "y": cfg.y_column, "z": cfg.z_column},
            "scale_ratio_ok_range": [cfg.scale_ratio_ok_min, cfg.scale_ratio_ok_max],
            **global_summary,
        },
    )

    bbox_overlay_png = out_dir / "xy_bbox_overlay.png"
    _plot_bbox_overlay(global_summary, bbox_overlay_png)

    centroid_delta_png = out_dir / "centroid_delta_by_slice.png"
    _plot_centroid_delta(by_slice, centroid_delta_png)

    contour_path = (
        Path(cfg.contour_csv).expanduser()
        if cfg.contour_csv is not None
        else Path(cfg.stack_root).expanduser() / "aligned_contour_3d_points.csv"
    )
    contour = pd.read_csv(contour_path) if contour_path.exists() else pd.DataFrame()
    overlay_slices = _resolve_overlay_slices(cfg.overlay_slices, by_slice["slice_order"].to_numpy(dtype=int))
    overlay_pngs = tuple(
        _plot_slice_density_overlay(
            histoseg=hs,
            spacemap=sm,
            contour=contour,
            slice_order=slice_order,
            cfg=cfg,
            out_png=out_dir / f"slice_density_vs_contour_{slice_order:03d}.png",
        )
        for slice_order in overlay_slices
    )

    return SpaceMapCoordinateQcResult(
        out_dir=out_dir,
        summary_json=summary_json,
        by_slice_csv=by_slice_csv,
        bbox_overlay_png=bbox_overlay_png,
        centroid_delta_png=centroid_delta_png,
        overlay_pngs=overlay_pngs,
    )


def _coordinate_qc_by_slice(
    histoseg: pd.DataFrame,
    spacemap: pd.DataFrame,
    cfg: SpaceMapCoordinateQcConfig,
) -> pd.DataFrame:
    hs_groups = histoseg.groupby(cfg.slice_order_column, sort=True)
    sm_groups = spacemap.groupby(cfg.slice_order_column, sort=True)
    rows: list[dict] = []
    for slice_order in sorted(set(hs_groups.groups) | set(sm_groups.groups)):
        hs = hs_groups.get_group(slice_order) if slice_order in hs_groups.groups else histoseg.iloc[0:0]
        sm = sm_groups.get_group(slice_order) if slice_order in sm_groups.groups else spacemap.iloc[0:0]
        hs_bbox = _bbox(hs, cfg)
        sm_bbox = _bbox(sm, cfg)
        hs_centroid = _centroid(hs, cfg)
        sm_centroid = _centroid(sm, cfg)
        dx = sm_centroid[0] - hs_centroid[0]
        dy = sm_centroid[1] - hs_centroid[1]
        rows.append(
            {
                "slice_order": int(slice_order),
                "sample_id_histoseg": _sample_id_for_slice(hs, cfg),
                "sample_id_spacemap": _sample_id_for_slice(sm, cfg),
                "histoseg_cell_count": int(len(hs)),
                "spacemap_cell_count": int(len(sm)),
                "cell_count_delta": int(len(sm) - len(hs)),
                "histoseg_x_min_um": hs_bbox["x_min"],
                "histoseg_x_max_um": hs_bbox["x_max"],
                "histoseg_y_min_um": hs_bbox["y_min"],
                "histoseg_y_max_um": hs_bbox["y_max"],
                "spacemap_x_min_um": sm_bbox["x_min"],
                "spacemap_x_max_um": sm_bbox["x_max"],
                "spacemap_y_min_um": sm_bbox["y_min"],
                "spacemap_y_max_um": sm_bbox["y_max"],
                "histoseg_width_um": hs_bbox["width"],
                "histoseg_height_um": hs_bbox["height"],
                "spacemap_width_um": sm_bbox["width"],
                "spacemap_height_um": sm_bbox["height"],
                "scale_ratio_width": _safe_div(sm_bbox["width"], hs_bbox["width"]),
                "scale_ratio_height": _safe_div(sm_bbox["height"], hs_bbox["height"]),
                "histoseg_centroid_x_um": hs_centroid[0],
                "histoseg_centroid_y_um": hs_centroid[1],
                "spacemap_centroid_x_um": sm_centroid[0],
                "spacemap_centroid_y_um": sm_centroid[1],
                "centroid_delta_x_um": dx,
                "centroid_delta_y_um": dy,
                "centroid_delta_norm_um": float(np.hypot(dx, dy)) if np.isfinite(dx) and np.isfinite(dy) else np.nan,
                "spacemap_cells_inside_histoseg_bbox_fraction": _bbox_coverage(sm, hs_bbox, cfg),
            }
        )
    return pd.DataFrame(rows)


def _global_coordinate_summary(
    histoseg: pd.DataFrame,
    spacemap: pd.DataFrame,
    by_slice: pd.DataFrame,
    cfg: SpaceMapCoordinateQcConfig,
) -> dict:
    hs_bbox = _bbox(histoseg, cfg)
    sm_bbox = _bbox(spacemap, cfg)
    width_ratio = _safe_div(sm_bbox["width"], hs_bbox["width"])
    height_ratio = _safe_div(sm_bbox["height"], hs_bbox["height"])
    centroid_dx = float(by_slice["centroid_delta_x_um"].mean())
    centroid_dy = float(by_slice["centroid_delta_y_um"].mean())
    centroid_sd_x = float(by_slice["centroid_delta_x_um"].std(ddof=0))
    centroid_sd_y = float(by_slice["centroid_delta_y_um"].std(ddof=0))
    scale_ok = (
        cfg.scale_ratio_ok_min <= width_ratio <= cfg.scale_ratio_ok_max
        and cfg.scale_ratio_ok_min <= height_ratio <= cfg.scale_ratio_ok_max
    )
    cell_counts_ok = bool((by_slice["cell_count_delta"] == 0).all())
    coverage_median = float(by_slice["spacemap_cells_inside_histoseg_bbox_fraction"].median())
    diagnosis: list[str] = []
    if not cell_counts_ok:
        diagnosis.append("cell_count_mismatch")
    if not scale_ok:
        diagnosis.append("scale_or_unit_mismatch_suspected")
    if coverage_median < 0.8:
        diagnosis.append("low_bbox_coverage")
    if not diagnosis:
        diagnosis.append("coordinate_frame_appears_compatible")
    return {
        "histoseg": {"cell_count": int(len(histoseg)), "bbox_um": hs_bbox},
        "spacemap": {"cell_count": int(len(spacemap)), "bbox_um": sm_bbox},
        "global_scale_ratio": {"width": width_ratio, "height": height_ratio},
        "global_translation_vector_um": {
            "mean_dx": centroid_dx,
            "mean_dy": centroid_dy,
            "mean_norm": float(np.hypot(centroid_dx, centroid_dy)),
            "sd_dx": centroid_sd_x,
            "sd_dy": centroid_sd_y,
        },
        "median_spacemap_cells_inside_histoseg_bbox_fraction": coverage_median,
        "cell_counts_match_by_slice": cell_counts_ok,
        "scale_ratio_within_threshold": bool(scale_ok),
        "diagnosis": diagnosis,
    }


def _bbox(frame: pd.DataFrame, cfg: SpaceMapCoordinateQcConfig) -> dict[str, float]:
    if frame.empty:
        return {key: np.nan for key in ("x_min", "x_max", "y_min", "y_max", "width", "height")}
    x = frame[cfg.x_column].to_numpy(dtype=float)
    y = frame[cfg.y_column].to_numpy(dtype=float)
    finite = np.isfinite(x) & np.isfinite(y)
    if not finite.any():
        return {key: np.nan for key in ("x_min", "x_max", "y_min", "y_max", "width", "height")}
    x, y = x[finite], y[finite]
    x_min, x_max = float(np.min(x)), float(np.max(x))
    y_min, y_max = float(np.min(y)), float(np.max(y))
    return {
        "x_min": x_min,
        "x_max": x_max,
        "y_min": y_min,
        "y_max": y_max,
        "width": x_max - x_min,
        "height": y_max - y_min,
    }


def _centroid(frame: pd.DataFrame, cfg: SpaceMapCoordinateQcConfig) -> tuple[float, float]:
    if frame.empty:
        return np.nan, np.nan
    x = frame[cfg.x_column].to_numpy(dtype=float)
    y = frame[cfg.y_column].to_numpy(dtype=float)
    finite = np.isfinite(x) & np.isfinite(y)
    if not finite.any():
        return np.nan, np.nan
    return float(np.mean(x[finite])), float(np.mean(y[finite]))


def _bbox_coverage(points: pd.DataFrame, bbox: dict[str, float], cfg: SpaceMapCoordinateQcConfig) -> float:
    if points.empty or not np.isfinite(list(bbox.values())).all():
        return np.nan
    x = points[cfg.x_column].to_numpy(dtype=float)
    y = points[cfg.y_column].to_numpy(dtype=float)
    inside = (
        (x >= bbox["x_min"])
        & (x <= bbox["x_max"])
        & (y >= bbox["y_min"])
        & (y <= bbox["y_max"])
        & np.isfinite(x)
        & np.isfinite(y)
    )
    return float(np.mean(inside))


def _sample_id_for_slice(frame: pd.DataFrame, cfg: SpaceMapCoordinateQcConfig) -> str | None:
    if frame.empty or cfg.sample_column not in frame.columns:
        return None
    values = frame[cfg.sample_column].dropna().astype(str).unique()
    return values[0] if len(values) else None


def _resolve_overlay_slices(requested: Sequence[int], available: np.ndarray) -> tuple[int, ...]:
    if requested:
        return tuple(int(value) for value in requested)
    if available.size == 0:
        return ()
    available = np.sort(np.unique(available.astype(int)))
    candidates = [available[0], available[len(available) // 2], available[-1]]
    return tuple(dict.fromkeys(int(value) for value in candidates))


def _plot_bbox_overlay(summary: dict, out_png: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 7))
    _draw_bbox(ax, summary["histoseg"]["bbox_um"], color="#2b6cb0", label="HistoSeg")
    _draw_bbox(ax, summary["spacemap"]["bbox_um"], color="#c53030", label="Space-map")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("X (µm)")
    ax.set_ylabel("Y (µm)")
    ax.set_title("Global XY bounding box audit")
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=220)
    plt.close(fig)


def _draw_bbox(ax, bbox: dict[str, float], *, color: str, label: str) -> None:
    if not np.isfinite(list(bbox.values())).all():
        return
    xs = [bbox["x_min"], bbox["x_max"], bbox["x_max"], bbox["x_min"], bbox["x_min"]]
    ys = [bbox["y_min"], bbox["y_min"], bbox["y_max"], bbox["y_max"], bbox["y_min"]]
    ax.plot(xs, ys, color=color, linewidth=2.2, label=label)


def _plot_centroid_delta(by_slice: pd.DataFrame, out_png: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    sc = axes[0].scatter(
        by_slice["centroid_delta_x_um"],
        by_slice["centroid_delta_y_um"],
        c=by_slice["slice_order"],
        cmap="viridis",
        s=42,
    )
    axes[0].axhline(0, color="0.6", linewidth=1)
    axes[0].axvline(0, color="0.6", linewidth=1)
    axes[0].set_aspect("equal", adjustable="box")
    axes[0].set_xlabel("Centroid dx (Space-map - HistoSeg, um)")
    axes[0].set_ylabel("Centroid dy (Space-map - HistoSeg, um)")
    axes[0].set_title("Centroid delta trajectory")
    fig.colorbar(sc, ax=axes[0], label="Slice order")

    axes[1].plot(by_slice["slice_order"], by_slice["centroid_delta_norm_um"], marker="o", linewidth=1.5)
    axes[1].set_xlabel("Slice order")
    axes[1].set_ylabel("Centroid delta norm (um)")
    axes[1].set_title("Centroid delta magnitude")
    axes[1].grid(alpha=0.25)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=220)
    plt.close(fig)


def _plot_slice_density_overlay(
    *,
    histoseg: pd.DataFrame,
    spacemap: pd.DataFrame,
    contour: pd.DataFrame,
    slice_order: int,
    cfg: SpaceMapCoordinateQcConfig,
    out_png: Path,
) -> Path:
    hs = histoseg.loc[histoseg[cfg.slice_order_column] == slice_order]
    sm = spacemap.loc[spacemap[cfg.slice_order_column] == slice_order]
    if len(hs) > cfg.max_points_per_overlay:
        hs = hs.sample(n=cfg.max_points_per_overlay, random_state=cfg.random_state)
    if len(sm) > cfg.max_points_per_overlay:
        sm = sm.sample(n=cfg.max_points_per_overlay, random_state=cfg.random_state + 1)

    fig, ax = plt.subplots(figsize=(8, 7))
    ax.scatter(hs[cfg.x_column], hs[cfg.y_column], s=0.35, c="#2b6cb0", alpha=0.18, label="HistoSeg cells")
    ax.scatter(sm[cfg.x_column], sm[cfg.y_column], s=0.35, c="#c53030", alpha=0.18, label="Space-map cells")
    if not contour.empty and {"slice_order", "structure", "polyline_id", "point_index", "x_um", "y_um"}.issubset(contour.columns):
        sub = contour.loc[contour["slice_order"].astype(int) == int(slice_order)]
        for structure_index, (_, group) in enumerate(sub.groupby("structure", sort=True)):
            color = f"C{structure_index % 10}"
            label_used = False
            for _, line in group.sort_values(["polyline_id", "point_index"]).groupby("polyline_id", sort=False):
                ax.plot(
                    line["x_um"],
                    line["y_um"],
                    color=color,
                    linewidth=1.0,
                    alpha=0.85,
                    label=str(group["structure"].iloc[0]) if not label_used else None,
                )
                label_used = True
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("X (µm)")
    ax.set_ylabel("Y (µm)")
    ax.set_title(f"Slice {slice_order:03d}: cell density vs HistoSeg contours")
    ax.legend(markerscale=6, fontsize=8, frameon=False, loc="best")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=220)
    plt.close(fig)
    return out_png


def _safe_div(numerator: float, denominator: float) -> float:
    if not np.isfinite(numerator) or not np.isfinite(denominator) or denominator == 0:
        return np.nan
    return float(numerator / denominator)


def _require_columns(df: pd.DataFrame, columns: Sequence[str]) -> None:
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=str), encoding="utf-8")


def _parse_int_list(values: Sequence[str]) -> tuple[int, ...]:
    parsed: list[int] = []
    for value in values:
        parsed.extend(int(part.strip()) for part in str(value).split(",") if part.strip())
    return tuple(parsed)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Audit Space-map coordinates against a HistoSeg 3D coordinate frame.")
    parser.add_argument("--spacemap-aligned-cells-parquet", required=True)
    parser.add_argument("--histoseg-aligned-cells-parquet", required=True)
    parser.add_argument("--stack-root", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--sample-column", default="sample_id")
    parser.add_argument("--slice-order-column", default="slice_order")
    parser.add_argument("--x-column", default="x_3d_um")
    parser.add_argument("--y-column", default="y_3d_um")
    parser.add_argument("--z-column", default="z_um")
    parser.add_argument("--overlay-slices", nargs="*", default=())
    parser.add_argument("--max-points-per-overlay", type=int, default=80000)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--contour-csv", default=None)
    parser.add_argument("--scale-ratio-ok-min", type=float, default=0.95)
    parser.add_argument("--scale-ratio-ok-max", type=float, default=1.05)
    args = parser.parse_args(argv)
    result = run_spacemap_coordinate_frame_qc(
        SpaceMapCoordinateQcConfig(
            spacemap_aligned_cells_parquet=args.spacemap_aligned_cells_parquet,
            histoseg_aligned_cells_parquet=args.histoseg_aligned_cells_parquet,
            stack_root=args.stack_root,
            out_dir=args.out_dir,
            sample_column=args.sample_column,
            slice_order_column=args.slice_order_column,
            x_column=args.x_column,
            y_column=args.y_column,
            z_column=args.z_column,
            overlay_slices=_parse_int_list(args.overlay_slices),
            max_points_per_overlay=args.max_points_per_overlay,
            random_state=args.random_state,
            contour_csv=args.contour_csv,
            scale_ratio_ok_min=args.scale_ratio_ok_min,
            scale_ratio_ok_max=args.scale_ratio_ok_max,
        )
    )
    print(json.dumps(asdict(result), indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
