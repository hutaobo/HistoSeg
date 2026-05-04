from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Literal

import pandas as pd


CoordinateUnit = Literal["um", "pixel"]


@dataclass(frozen=True)
class SpaceMapImportConfig:
    spacemap_aligned_csv: str | Path
    slice_manifest: str | Path
    out_parquet: str | Path
    histoseg_aligned_cells_parquet: str | Path | None = None
    h5ad: str | Path | None = None
    out_summary_json: str | Path | None = None
    coordinate_unit: CoordinateUnit = "um"
    pixel_size_um: float = 0.2125
    x_column: str = "x"
    y_column: str = "y"
    z_column: str = "z"
    layer_column: str = "layer"
    sample_column: str = "sample_id"
    barcode_column: str = "barcode"
    allow_row_order_match: bool = False
    use_spacemap_z: bool = False
    z_spacing_um: float = 5.0


@dataclass(frozen=True)
class SpaceMapImportResult:
    out_parquet: Path
    summary_json: Path
    cell_count: int
    slice_count: int
    match_strategy: str
    coordinate_unit: str
    layer_mapping_mode: str


def import_spacemap_aligned_coordinates(cfg: SpaceMapImportConfig) -> SpaceMapImportResult:
    """Normalize Space-map aligned coordinates to HistoSeg's aligned cell schema."""

    if cfg.coordinate_unit not in {"um", "pixel"}:
        raise ValueError("coordinate_unit must be 'um' or 'pixel'.")
    if cfg.pixel_size_um <= 0:
        raise ValueError("pixel_size_um must be positive.")
    if cfg.z_spacing_um <= 0:
        raise ValueError("z_spacing_um must be positive.")
    if cfg.histoseg_aligned_cells_parquet is None and cfg.h5ad is None:
        raise ValueError("Provide either histoseg_aligned_cells_parquet or h5ad as the reference cell table.")
    if cfg.histoseg_aligned_cells_parquet is not None and cfg.h5ad is not None:
        raise ValueError("Provide only one of histoseg_aligned_cells_parquet or h5ad.")

    manifest = _load_slice_manifest(Path(cfg.slice_manifest), cfg)
    reference = _load_reference_cells(cfg)
    spacemap = pd.read_csv(Path(cfg.spacemap_aligned_csv))
    _require_columns(spacemap, [cfg.x_column, cfg.y_column])
    spacemap_norm, layer_mode = _normalise_spacemap_table(spacemap, manifest, cfg)
    output, match_strategy = _match_to_reference(reference, spacemap_norm, cfg)

    out_parquet = Path(cfg.out_parquet).expanduser()
    out_parquet.parent.mkdir(parents=True, exist_ok=True)
    output.to_parquet(out_parquet, index=False)

    summary_json = (
        Path(cfg.out_summary_json).expanduser()
        if cfg.out_summary_json is not None
        else out_parquet.with_suffix(".summary.json")
    )
    summary = {
        "spacemap_aligned_csv": str(Path(cfg.spacemap_aligned_csv).expanduser()),
        "slice_manifest": str(Path(cfg.slice_manifest).expanduser()),
        "histoseg_aligned_cells_parquet": (
            str(Path(cfg.histoseg_aligned_cells_parquet).expanduser())
            if cfg.histoseg_aligned_cells_parquet is not None
            else None
        ),
        "h5ad": str(Path(cfg.h5ad).expanduser()) if cfg.h5ad is not None else None,
        "out_parquet": str(out_parquet),
        "cell_count": int(len(output)),
        "slice_count": int(output[cfg.sample_column].nunique()),
        "match_strategy": match_strategy,
        "coordinate_unit": cfg.coordinate_unit,
        "pixel_size_um": cfg.pixel_size_um,
        "z_source": "spacemap_csv" if cfg.use_spacemap_z and cfg.z_column in spacemap.columns else "slice_manifest",
        "layer_mapping_mode": layer_mode,
        "allow_row_order_match": cfg.allow_row_order_match,
        "counts_by_slice_order": {
            str(int(order)): int(count)
            for order, count in output.groupby("slice_order", sort=True).size().items()
        },
    }
    _write_json(summary_json, summary)
    return SpaceMapImportResult(
        out_parquet=out_parquet,
        summary_json=summary_json,
        cell_count=int(len(output)),
        slice_count=int(output[cfg.sample_column].nunique()),
        match_strategy=match_strategy,
        coordinate_unit=cfg.coordinate_unit,
        layer_mapping_mode=layer_mode,
    )


def _load_slice_manifest(path: Path, cfg: SpaceMapImportConfig) -> pd.DataFrame:
    manifest = pd.read_csv(path.expanduser())
    order_column = "order" if "order" in manifest.columns else "slice_order"
    _require_columns(manifest, [order_column, cfg.sample_column])
    result = manifest.copy()
    result["slice_order"] = result[order_column].astype(int)
    result[cfg.sample_column] = result[cfg.sample_column].astype(str)
    if "z_um" not in result.columns:
        result["z_um"] = (result["slice_order"] - 1) * cfg.z_spacing_um
    result["z_um"] = result["z_um"].astype(float)
    if result["slice_order"].duplicated().any():
        raise ValueError("Slice manifest contains duplicate slice orders.")
    if result[cfg.sample_column].duplicated().any():
        raise ValueError(f"Slice manifest contains duplicate {cfg.sample_column!r} values.")
    return result[["slice_order", cfg.sample_column, "z_um"]].sort_values("slice_order").reset_index(drop=True)


def _load_reference_cells(cfg: SpaceMapImportConfig) -> pd.DataFrame:
    if cfg.histoseg_aligned_cells_parquet is not None:
        reference = pd.read_parquet(Path(cfg.histoseg_aligned_cells_parquet).expanduser())
    else:
        try:
            import anndata as ad
        except Exception as exc:  # pragma: no cover - dependency availability.
            raise ImportError("anndata is required when --h5ad is used as the reference table.") from exc
        adata = ad.read_h5ad(Path(cfg.h5ad).expanduser(), backed="r")
        try:
            reference = adata.obs.copy()
        finally:
            if getattr(adata, "isbacked", False):
                adata.file.close()
        if cfg.barcode_column not in reference.columns:
            reference[cfg.barcode_column] = reference.index.astype(str)
        reference = reference.reset_index(drop=True)

    _require_columns(reference, [cfg.sample_column])
    if cfg.barcode_column not in reference.columns and not cfg.allow_row_order_match:
        raise ValueError(
            f"Reference table lacks {cfg.barcode_column!r}; use --allow-row-order-match only if Space-map rows "
            "are in the same within-slice order as the reference cells."
        )
    reference = reference.copy()
    reference[cfg.sample_column] = reference[cfg.sample_column].astype(str)
    if cfg.barcode_column in reference.columns:
        reference[cfg.barcode_column] = reference[cfg.barcode_column].astype(str)
    reference["_reference_order"] = range(len(reference))
    return reference


def _normalise_spacemap_table(
    spacemap: pd.DataFrame,
    manifest: pd.DataFrame,
    cfg: SpaceMapImportConfig,
) -> tuple[pd.DataFrame, str]:
    table = spacemap.copy()
    if cfg.layer_column not in table.columns:
        if cfg.sample_column not in table.columns:
            raise ValueError(
                f"Space-map table must contain {cfg.layer_column!r} or {cfg.sample_column!r} for slice mapping."
            )
        table[cfg.layer_column] = table[cfg.sample_column]

    mapped, mode = _map_layers(table[cfg.layer_column], manifest, cfg.sample_column)
    table[cfg.sample_column] = mapped[cfg.sample_column].to_numpy()
    table["slice_order"] = mapped["slice_order"].to_numpy()
    if cfg.use_spacemap_z and cfg.z_column in table.columns:
        table["z_um"] = pd.to_numeric(table[cfg.z_column], errors="raise").astype(float)
    else:
        table["z_um"] = mapped["z_um"].to_numpy(dtype=float)

    scale = cfg.pixel_size_um if cfg.coordinate_unit == "pixel" else 1.0
    x_raw = pd.to_numeric(table[cfg.x_column], errors="raise").astype(float)
    y_raw = pd.to_numeric(table[cfg.y_column], errors="raise").astype(float)
    table["x_spacemap_raw"] = x_raw
    table["y_spacemap_raw"] = y_raw
    table["x_spacemap_um"] = x_raw * scale
    table["y_spacemap_um"] = y_raw * scale
    table["spacemap_layer"] = table[cfg.layer_column].astype(str)
    if cfg.barcode_column in table.columns:
        table[cfg.barcode_column] = table[cfg.barcode_column].astype(str)
    return table, mode


def _map_layers(layer_values: pd.Series, manifest: pd.DataFrame, sample_column: str) -> tuple[pd.DataFrame, str]:
    unique = pd.Series(layer_values.dropna().unique()).astype(str)
    by_sample = manifest.set_index(sample_column)
    if set(unique).issubset(set(by_sample.index.astype(str))):
        mapped = layer_values.astype(str).map(by_sample.to_dict(orient="index"))
        return pd.DataFrame(list(mapped)), "sample_id"

    numeric = unique.map(_integer_like_or_none)
    if numeric.isna().any():
        raise ValueError("Could not map Space-map layer values to manifest sample IDs or numeric slice orders.")
    numeric_values = numeric.astype(int)
    by_order = manifest.set_index("slice_order")
    if set(numeric_values).issubset(set(by_order.index.astype(int))):
        layer_to_order = dict(zip(unique, numeric_values))
        orders = layer_values.astype(str).map(layer_to_order).astype(int)
        return by_order.loc[orders].reset_index(drop=False).reset_index(drop=True), "one_based_order"

    shifted = numeric_values + 1
    if set(shifted).issubset(set(by_order.index.astype(int))):
        layer_to_order = dict(zip(unique, shifted))
        orders = layer_values.astype(str).map(layer_to_order).astype(int)
        return by_order.loc[orders].reset_index(drop=False).reset_index(drop=True), "zero_based_order"

    raise ValueError("Numeric Space-map layer values do not match manifest slice orders.")


def _match_to_reference(
    reference: pd.DataFrame,
    spacemap: pd.DataFrame,
    cfg: SpaceMapImportConfig,
) -> tuple[pd.DataFrame, str]:
    if cfg.barcode_column in reference.columns and cfg.barcode_column in spacemap.columns:
        key_columns = [cfg.sample_column, cfg.barcode_column]
        if spacemap.duplicated(key_columns).any():
            duplicates = spacemap.loc[spacemap.duplicated(key_columns, keep=False), key_columns].head().to_dict("records")
            raise ValueError(f"Space-map table has duplicate sample/barcode keys: {duplicates}")
        match_columns = [
            *key_columns,
            "slice_order",
            "z_um",
            "x_spacemap_raw",
            "y_spacemap_raw",
            "x_spacemap_um",
            "y_spacemap_um",
            "spacemap_layer",
        ]
        space_match = spacemap[match_columns].rename(
            columns={"slice_order": "slice_order_spacemap", "z_um": "z_spacemap_um"}
        )
        merged = reference.merge(space_match, on=key_columns, how="left", validate="one_to_one")
        strategy = "sample_barcode"
    else:
        if not cfg.allow_row_order_match:
            raise ValueError(
                "Space-map and reference tables must both contain barcode for safe matching. "
                "Use --allow-row-order-match only for verified same-order exports."
            )
        ref = reference.copy()
        sm = spacemap.copy()
        ref["_within_slice_index"] = ref.groupby(cfg.sample_column, sort=False).cumcount()
        sm["_within_slice_index"] = sm.groupby(cfg.sample_column, sort=False).cumcount()
        expected = ref.groupby(cfg.sample_column).size().sort_index()
        observed = sm.groupby(cfg.sample_column).size().sort_index()
        if not expected.equals(observed):
            raise ValueError(
                "Row-order matching requires identical per-slice counts. "
                f"Reference counts={expected.to_dict()}, Space-map counts={observed.to_dict()}"
            )
        key_columns = [cfg.sample_column, "_within_slice_index"]
        match_columns = [
            *key_columns,
            "slice_order",
            "z_um",
            "x_spacemap_raw",
            "y_spacemap_raw",
            "x_spacemap_um",
            "y_spacemap_um",
            "spacemap_layer",
        ]
        space_match = sm[match_columns].rename(
            columns={"slice_order": "slice_order_spacemap", "z_um": "z_spacemap_um"}
        )
        merged = ref.merge(space_match, on=key_columns, how="left", validate="one_to_one")
        strategy = "within_slice_row_order"

    if merged["x_spacemap_um"].isna().any() or merged["y_spacemap_um"].isna().any():
        missing = int(merged["x_spacemap_um"].isna().sum())
        raise ValueError(f"Failed to match {missing} reference cells to Space-map coordinates.")

    merged = merged.sort_values("_reference_order").reset_index(drop=True)
    output = merged.copy()
    for column, renamed in (
        ("x_3d_um", "x_histoseg_3d_um"),
        ("y_3d_um", "y_histoseg_3d_um"),
        ("z_um", "z_histoseg_um"),
    ):
        if column in output.columns:
            output[renamed] = output[column]
    output["slice_order"] = output["slice_order_spacemap"].astype(int)
    output["x_3d_um"] = output["x_spacemap_um"].astype(float)
    output["y_3d_um"] = output["y_spacemap_um"].astype(float)
    output["z_um"] = output["z_spacemap_um"].astype(float)
    output["coordinate_backend"] = "spacemap"
    drop_columns = [
        column
        for column in ("_reference_order", "_within_slice_index", "slice_order_spacemap", "z_spacemap_um")
        if column in output.columns
    ]
    output = output.drop(columns=drop_columns)
    return output, strategy


def _integer_like_or_none(value: str) -> int | float:
    try:
        number = float(value)
    except ValueError:
        return math.nan
    rounded = round(number)
    if abs(number - rounded) > 1e-6:
        return math.nan
    return int(rounded)


def _require_columns(df: pd.DataFrame, columns: list[str]) -> None:
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=str), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Normalize Space-map aligned CSV coordinates to HistoSeg cell Parquet.")
    parser.add_argument("--spacemap-aligned-csv", required=True)
    parser.add_argument("--slice-manifest", required=True)
    parser.add_argument("--out-parquet", required=True)
    parser.add_argument("--histoseg-aligned-cells-parquet", default=None)
    parser.add_argument("--h5ad", default=None)
    parser.add_argument("--out-summary-json", default=None)
    parser.add_argument("--coordinate-unit", choices=["um", "pixel"], default="um")
    parser.add_argument("--pixel-size-um", type=float, default=0.2125)
    parser.add_argument("--x-column", default="x")
    parser.add_argument("--y-column", default="y")
    parser.add_argument("--z-column", default="z")
    parser.add_argument("--layer-column", default="layer")
    parser.add_argument("--sample-column", default="sample_id")
    parser.add_argument("--barcode-column", default="barcode")
    parser.add_argument("--allow-row-order-match", action="store_true")
    parser.add_argument("--use-spacemap-z", action="store_true")
    parser.add_argument("--z-spacing-um", type=float, default=5.0)
    args = parser.parse_args(argv)
    result = import_spacemap_aligned_coordinates(
        SpaceMapImportConfig(
            spacemap_aligned_csv=args.spacemap_aligned_csv,
            slice_manifest=args.slice_manifest,
            out_parquet=args.out_parquet,
            histoseg_aligned_cells_parquet=args.histoseg_aligned_cells_parquet,
            h5ad=args.h5ad,
            out_summary_json=args.out_summary_json,
            coordinate_unit=args.coordinate_unit,
            pixel_size_um=args.pixel_size_um,
            x_column=args.x_column,
            y_column=args.y_column,
            z_column=args.z_column,
            layer_column=args.layer_column,
            sample_column=args.sample_column,
            barcode_column=args.barcode_column,
            allow_row_order_match=args.allow_row_order_match,
            use_spacemap_z=args.use_spacemap_z,
            z_spacing_um=args.z_spacing_um,
        )
    )
    print(json.dumps(asdict(result), indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
