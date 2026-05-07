from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class SpaceMapHarmonizationConfig:
    spacemap_aligned_cells_parquet: str | Path
    histoseg_aligned_cells_parquet: str | Path
    out_parquet: str | Path
    out_summary_json: str | Path | None = None
    sample_column: str = "sample_id"
    barcode_column: str = "barcode"
    x_column: str = "x_3d_um"
    y_column: str = "y_3d_um"
    max_points: int = 10000
    random_state: int = 0
    allow_reflection: bool = False


@dataclass(frozen=True)
class SpaceMapHarmonizationResult:
    out_parquet: Path
    summary_json: Path
    matched_cell_count: int
    points_used: int
    scale: float
    rotation_degrees: float
    rmse_used_before_um: float
    rmse_used_after_um: float


def harmonize_spacemap_frame(cfg: SpaceMapHarmonizationConfig) -> SpaceMapHarmonizationResult:
    """Apply a global 2D similarity transform to Space-map coordinates in HistoSeg's frame."""

    if cfg.max_points <= 1:
        raise ValueError("max_points must be greater than 1.")
    hs = pd.read_parquet(Path(cfg.histoseg_aligned_cells_parquet).expanduser())
    sm = pd.read_parquet(Path(cfg.spacemap_aligned_cells_parquet).expanduser())
    _require_columns(hs, [cfg.sample_column, cfg.barcode_column, cfg.x_column, cfg.y_column])
    _require_columns(sm, [cfg.sample_column, cfg.barcode_column, cfg.x_column, cfg.y_column])

    key_columns = [cfg.sample_column, cfg.barcode_column]
    if hs.duplicated(key_columns).any() or sm.duplicated(key_columns).any():
        raise ValueError("Both HistoSeg and Space-map cell tables must have unique sample/barcode keys.")

    hs_match = hs[[*key_columns, cfg.x_column, cfg.y_column]].rename(
        columns={cfg.x_column: "_histoseg_x_um", cfg.y_column: "_histoseg_y_um"}
    )
    sm_match = sm[[*key_columns, cfg.x_column, cfg.y_column]].rename(
        columns={cfg.x_column: "_spacemap_x_um", cfg.y_column: "_spacemap_y_um"}
    )
    matched = sm_match.merge(hs_match, on=key_columns, how="inner", validate="one_to_one")
    if len(matched) != len(sm) or len(matched) != len(hs):
        raise ValueError(
            "Harmonization requires exact matched cell identities. "
            f"matched={len(matched)}, spacemap={len(sm)}, histoseg={len(hs)}"
        )

    finite = np.isfinite(
        matched[["_spacemap_x_um", "_spacemap_y_um", "_histoseg_x_um", "_histoseg_y_um"]].to_numpy(dtype=float)
    ).all(axis=1)
    matched = matched.loc[finite].reset_index(drop=True)
    if len(matched) < 2:
        raise ValueError("At least two finite matched cells are required for similarity harmonization.")

    if len(matched) > cfg.max_points:
        used = matched.sample(n=cfg.max_points, random_state=cfg.random_state).reset_index(drop=True)
    else:
        used = matched

    src = used[["_spacemap_x_um", "_spacemap_y_um"]].to_numpy(dtype=float)
    dst = used[["_histoseg_x_um", "_histoseg_y_um"]].to_numpy(dtype=float)
    transform = _estimate_similarity_umeyama(src, dst, allow_reflection=cfg.allow_reflection)

    all_src = sm[[cfg.x_column, cfg.y_column]].to_numpy(dtype=float)
    transformed = _apply_similarity(all_src, transform)
    out = sm.copy()
    out["x_spacemap_pre_harmonized_um"] = out[cfg.x_column].astype(float)
    out["y_spacemap_pre_harmonized_um"] = out[cfg.y_column].astype(float)
    out[cfg.x_column] = transformed[:, 0]
    out[cfg.y_column] = transformed[:, 1]
    out["coordinate_backend"] = "spacemap_global_similarity_harmonized"

    out_parquet = Path(cfg.out_parquet).expanduser()
    out_parquet.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(out_parquet, index=False)

    used_before = _rmse(src, dst)
    used_after = _rmse(_apply_similarity(src, transform), dst)
    matched_src = matched[["_spacemap_x_um", "_spacemap_y_um"]].to_numpy(dtype=float)
    matched_dst = matched[["_histoseg_x_um", "_histoseg_y_um"]].to_numpy(dtype=float)
    matched_before = _rmse(matched_src, matched_dst)
    matched_after = _rmse(_apply_similarity(matched_src, transform), matched_dst)
    rotation_degrees = float(np.degrees(np.arctan2(transform["rotation_matrix"][1, 0], transform["rotation_matrix"][0, 0])))

    summary_json = (
        Path(cfg.out_summary_json).expanduser()
        if cfg.out_summary_json is not None
        else out_parquet.with_suffix(".harmonization.json")
    )
    _write_json(
        summary_json,
        {
            "method": "global_similarity_transform_umeyama_2d",
            "formula": "x_histoseg_frame = scale * rotation_matrix @ x_spacemap + translation_um",
            "local_nonrigid_warp_applied": False,
            "allow_reflection": cfg.allow_reflection,
            "spacemap_aligned_cells_parquet": str(Path(cfg.spacemap_aligned_cells_parquet).expanduser()),
            "histoseg_aligned_cells_parquet": str(Path(cfg.histoseg_aligned_cells_parquet).expanduser()),
            "out_parquet": str(out_parquet),
            "matched_cell_count": int(len(matched)),
            "points_used": int(len(used)),
            "random_state": cfg.random_state,
            "scale": float(transform["scale"]),
            "rotation_degrees": rotation_degrees,
            "rotation_matrix": transform["rotation_matrix"].tolist(),
            "translation_um": transform["translation"].tolist(),
            "rmse_used_before_um": used_before,
            "rmse_used_after_um": used_after,
            "rmse_all_matched_before_um": matched_before,
            "rmse_all_matched_after_um": matched_after,
        },
    )

    return SpaceMapHarmonizationResult(
        out_parquet=out_parquet,
        summary_json=summary_json,
        matched_cell_count=int(len(matched)),
        points_used=int(len(used)),
        scale=float(transform["scale"]),
        rotation_degrees=rotation_degrees,
        rmse_used_before_um=used_before,
        rmse_used_after_um=used_after,
    )


def _estimate_similarity_umeyama(
    src: np.ndarray,
    dst: np.ndarray,
    *,
    allow_reflection: bool,
) -> dict[str, np.ndarray | float]:
    if src.shape != dst.shape or src.ndim != 2 or src.shape[1] != 2:
        raise ValueError("src and dst must both be N x 2 arrays.")
    src_mean = src.mean(axis=0)
    dst_mean = dst.mean(axis=0)
    src_centered = src - src_mean
    dst_centered = dst - dst_mean
    src_var = float(np.mean(np.sum(src_centered * src_centered, axis=1)))
    if src_var <= 0:
        raise ValueError("Source coordinates have zero variance.")
    covariance = (dst_centered.T @ src_centered) / len(src)
    u, singular_values, vt = np.linalg.svd(covariance)
    sign = np.ones(2)
    if not allow_reflection and np.linalg.det(u @ vt) < 0:
        sign[-1] = -1
    correction = np.diag(sign)
    rotation = u @ correction @ vt
    scale = float(np.sum(singular_values * sign) / src_var)
    translation = dst_mean - scale * (rotation @ src_mean)
    return {
        "scale": scale,
        "rotation_matrix": rotation,
        "translation": translation,
    }


def _apply_similarity(points: np.ndarray, transform: dict[str, np.ndarray | float]) -> np.ndarray:
    rotation = np.asarray(transform["rotation_matrix"], dtype=float)
    scale = float(transform["scale"])
    translation = np.asarray(transform["translation"], dtype=float)
    return (scale * (rotation @ points.T)).T + translation


def _rmse(src: np.ndarray, dst: np.ndarray) -> float:
    diff = src - dst
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def _require_columns(df: pd.DataFrame, columns: list[str]) -> None:
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=str), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Harmonize Space-map aligned coordinates to the HistoSeg frame using only a global similarity transform."
    )
    parser.add_argument("--spacemap-aligned-cells-parquet", required=True)
    parser.add_argument("--histoseg-aligned-cells-parquet", required=True)
    parser.add_argument("--out-parquet", required=True)
    parser.add_argument("--out-summary-json", default=None)
    parser.add_argument("--sample-column", default="sample_id")
    parser.add_argument("--barcode-column", default="barcode")
    parser.add_argument("--x-column", default="x_3d_um")
    parser.add_argument("--y-column", default="y_3d_um")
    parser.add_argument("--max-points", type=int, default=10000)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--allow-reflection", action="store_true")
    args = parser.parse_args(argv)
    result = harmonize_spacemap_frame(
        SpaceMapHarmonizationConfig(
            spacemap_aligned_cells_parquet=args.spacemap_aligned_cells_parquet,
            histoseg_aligned_cells_parquet=args.histoseg_aligned_cells_parquet,
            out_parquet=args.out_parquet,
            out_summary_json=args.out_summary_json,
            sample_column=args.sample_column,
            barcode_column=args.barcode_column,
            x_column=args.x_column,
            y_column=args.y_column,
            max_points=args.max_points,
            random_state=args.random_state,
            allow_reflection=args.allow_reflection,
        )
    )
    print(json.dumps(asdict(result), indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
