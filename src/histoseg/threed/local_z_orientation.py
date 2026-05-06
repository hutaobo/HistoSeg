"""Transcript-level local-z orientation correction for HistoSeg 3D stacks."""

from __future__ import annotations

from dataclasses import dataclass
from html import escape
import math
from pathlib import Path
from typing import Any, Mapping, Sequence, Union

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from .cell_cloud import (
    SliceCellTransform,
    apply_similarity_to_points,
    load_cell_alignment_transforms,
)
from .multislice import discover_xenium_slices


PathLike = Union[str, Path]

LOCAL_Z_ORIENTATION_SCHEMA_VERSION = 1
LOCAL_Z_ORIENTATION_STATES = ("preserve", "reverse")
ALIGNED_TRANSCRIPT_OPTIONAL_COLUMNS = (
    "transcript_id",
    "cell_id",
    "overlaps_nucleus",
    "qv",
    "structure",
)
ALIGNED_TRANSCRIPT_COLUMNS = (
    "sample_id",
    "slice_order",
    "gene",
    "x_raw_um",
    "y_raw_um",
    "z_raw_um",
    "x_3d_um",
    "y_3d_um",
    "z_3d_um",
    "slice_stack_z_um",
    "local_z_corrected_um",
    "local_z_mid_um",
    "local_z_orientation",
    "local_z_flip_applied",
    *ALIGNED_TRANSCRIPT_OPTIONAL_COLUMNS,
)


@dataclass(frozen=True)
class LocalZOrientationConfig:
    """Configuration for transcript-only local-z orientation inference.

    This workflow does not change HistoSeg contour/cell stack reconstruction.
    It applies the existing aligned-stack XY transforms to transcript coordinates,
    infers whether each slice's internal transcript z-axis should be preserved or
    reversed, and writes a transcript-level 3D table.
    """

    xenium_root: PathLike
    stack_root: PathLike
    out_dir: PathLike | None = None
    sample_glob: str = "*"
    transcript_relpath: str = "transcripts.parquet"
    gene_column: str | None = None
    x_column: str | None = None
    y_column: str | None = None
    z_column: str | None = None
    transcript_id_column: str | None = None
    structure_column: str | None = None
    pixel_size_um: float = 0.2125
    z_band_fraction: float = 0.25
    min_band_transcripts: int = 10
    max_signature_genes: int = 256
    low_confidence_margin: float = 0.02
    vertical_qc_backend: str = "ovrlpy"
    apply_local_z_flip: bool = True
    ovrlpy_kde_bandwidth: float = 2.5
    ovrlpy_n_components: int = 20
    ovrlpy_n_workers: int = 1
    ovrlpy_fit_umap: bool = True
    ovrlpy_min_transcripts: float = 10.0
    doublet_min_signal: float = 4.0
    doublet_integrity_sigma: float = 1.0
    doublet_exclusion_radius_um: float = 20.0
    chunk_size: int = 100000


@dataclass
class LocalZOrientationResult:
    out_dir: Path
    manifest_csv: Path
    aligned_transcripts_parquet: Path
    biological_report_html: Path
    marker_gradients_csv: Path
    vertical_qc_dir: Path
    slice_count: int
    transcript_count: int
    best_score: float
    second_score: float
    confidence_margin: float
    low_confidence: bool


@dataclass(frozen=True)
class _TranscriptSource:
    sample_id: str
    sample_dir: Path
    xenium_dir: Path


@dataclass
class _OvrlpyQC:
    status: str
    doublets_csv: Path | None = None
    doublet_count: int = 0
    signal_integrity_npy: Path | None = None
    signal_map_npy: Path | None = None
    signal_integrity_mean: float = math.nan
    signal_integrity_p05: float = math.nan
    signal_integrity_p95: float = math.nan
    signal_mean: float = math.nan
    signal_p05: float = math.nan
    signal_p95: float = math.nan
    error: str | None = None


@dataclass
class _SliceLocalZSummary:
    sample_id: str
    order: int
    stack_z_um: float
    transcript_path: Path
    z_min_um: float
    z_max_um: float
    z_mid_um: float
    transcript_count: int
    signature_transcript_count: int
    excluded_doublet_transcripts: int
    low_band_count: int
    high_band_count: int
    raw_low_signature: dict[str, float]
    raw_high_signature: dict[str, float]
    ovrlpy_qc: _OvrlpyQC


@dataclass(frozen=True)
class _OrientationInference:
    states_by_sample: dict[str, str]
    best_score: float
    second_score: float
    confidence_margin: float
    low_confidence: bool
    edge_scores: dict[tuple[str, str], float]


def run_local_z_orientation_correction(
    cfg: LocalZOrientationConfig,
) -> LocalZOrientationResult:
    """Infer and apply transcript local-z orientation for an aligned HistoSeg stack."""

    _validate_config(cfg)
    stack_root = Path(cfg.stack_root).expanduser().resolve()
    xenium_root = Path(cfg.xenium_root).expanduser().resolve()
    out_dir = Path(cfg.out_dir).expanduser().resolve() if cfg.out_dir else stack_root
    out_dir.mkdir(parents=True, exist_ok=True)
    vertical_qc_dir = out_dir / "vertical_qc"
    vertical_qc_dir.mkdir(parents=True, exist_ok=True)

    transforms = load_cell_alignment_transforms(
        stack_root,
        pixel_size_um=cfg.pixel_size_um,
        chunk_size=cfg.chunk_size,
    )
    sources = _discover_transcript_sources(xenium_root, sample_glob=cfg.sample_glob)

    summaries: list[_SliceLocalZSummary] = []
    for transform in sorted(transforms, key=lambda item: item.order):
        transcript_path = _find_transcript_path(transform.sample_id, sources, xenium_root, cfg)
        summary = _summarize_slice_transcripts(
            transform=transform,
            transcript_path=transcript_path,
            vertical_qc_dir=vertical_qc_dir,
            cfg=cfg,
        )
        summaries.append(summary)

    inference = _infer_orientation_sequence(summaries, cfg)
    manifest = _build_manifest_frame(summaries, inference, cfg)
    manifest_path = out_dir / "local_z_orientation_manifest.csv"
    manifest.to_csv(manifest_path, index=False)

    aligned_path = out_dir / "aligned_transcripts_3d.parquet"
    transcript_count = _write_aligned_transcripts(
        summaries=summaries,
        transforms=transforms,
        inference=inference,
        output_parquet=aligned_path,
        cfg=cfg,
    )

    marker_gradients = _build_marker_gradient_frame(summaries, inference)
    marker_gradients_path = vertical_qc_dir / "top_bottom_marker_gradients.csv"
    marker_gradients.to_csv(marker_gradients_path, index=False)

    report_path = vertical_qc_dir / "biological_z_report.html"
    _write_biological_report(
        report_path,
        manifest=manifest,
        marker_gradients=marker_gradients,
        inference=inference,
        cfg=cfg,
    )

    return LocalZOrientationResult(
        out_dir=out_dir,
        manifest_csv=manifest_path,
        aligned_transcripts_parquet=aligned_path,
        biological_report_html=report_path,
        marker_gradients_csv=marker_gradients_path,
        vertical_qc_dir=vertical_qc_dir,
        slice_count=len(summaries),
        transcript_count=transcript_count,
        best_score=inference.best_score,
        second_score=inference.second_score,
        confidence_margin=inference.confidence_margin,
        low_confidence=inference.low_confidence,
    )


def _validate_config(cfg: LocalZOrientationConfig) -> None:
    if cfg.vertical_qc_backend not in {"none", "ovrlpy"}:
        raise ValueError("vertical_qc_backend must be 'none' or 'ovrlpy'.")
    if cfg.pixel_size_um <= 0:
        raise ValueError("pixel_size_um must be greater than 0.")
    if not (0.0 < cfg.z_band_fraction < 0.5):
        raise ValueError("z_band_fraction must be greater than 0 and less than 0.5.")
    if cfg.min_band_transcripts < 1:
        raise ValueError("min_band_transcripts must be at least 1.")
    if cfg.max_signature_genes < 1:
        raise ValueError("max_signature_genes must be at least 1.")
    if cfg.low_confidence_margin < 0:
        raise ValueError("low_confidence_margin must be non-negative.")
    if cfg.ovrlpy_kde_bandwidth <= 0:
        raise ValueError("ovrlpy_kde_bandwidth must be greater than 0.")
    if cfg.ovrlpy_n_components < 1:
        raise ValueError("ovrlpy_n_components must be at least 1.")
    if cfg.ovrlpy_n_workers < 1:
        raise ValueError("ovrlpy_n_workers must be at least 1.")
    if cfg.ovrlpy_min_transcripts <= 0:
        raise ValueError("ovrlpy_min_transcripts must be greater than 0.")
    if cfg.doublet_exclusion_radius_um < 0:
        raise ValueError("doublet_exclusion_radius_um must be non-negative.")
    if cfg.chunk_size <= 0:
        raise ValueError("chunk_size must be greater than 0.")


def _discover_transcript_sources(
    xenium_root: Path,
    *,
    sample_glob: str,
) -> dict[str, list[_TranscriptSource]]:
    sources: dict[str, list[_TranscriptSource]] = {}
    try:
        for item in discover_xenium_slices(xenium_root, sample_glob=sample_glob):
            sources.setdefault(item.sample_id, []).append(
                _TranscriptSource(
                    sample_id=item.sample_id,
                    sample_dir=item.sample_dir,
                    xenium_dir=item.xenium_dir,
                )
            )
    except Exception:
        pass

    candidates = [xenium_root]
    if xenium_root.exists() and xenium_root.is_dir():
        candidates.extend(child for child in xenium_root.glob(sample_glob) if child.is_dir())
    for candidate in candidates:
        sample_id = _sample_id_from_transcript_source(candidate)
        sources.setdefault(candidate.name, []).append(
            _TranscriptSource(
                sample_id=sample_id,
                sample_dir=candidate,
                xenium_dir=candidate,
            )
        )
        if sample_id != candidate.name:
            sources.setdefault(sample_id, []).append(
                _TranscriptSource(
                    sample_id=sample_id,
                    sample_dir=candidate,
                    xenium_dir=candidate,
                )
            )
    return sources


def _sample_id_from_transcript_source(path: Path) -> str:
    name = path.name
    if name.endswith(".pyxenium.slide.zarr"):
        return name[: -len(".pyxenium.slide.zarr")]
    return name


def _find_transcript_path(
    sample_id: str,
    sources: Mapping[str, Sequence[_TranscriptSource]],
    xenium_root: Path,
    cfg: LocalZOrientationConfig,
) -> Path:
    candidate_dirs: list[Path] = []
    for source in sources.get(sample_id, ()):
        candidate_dirs.extend([source.xenium_dir, source.sample_dir])
    candidate_dirs.extend([xenium_root / sample_id, xenium_root])

    relpaths = [cfg.transcript_relpath] if cfg.transcript_relpath else []
    relpaths.extend(
        [
            "transcripts.parquet",
            "transcripts.csv.gz",
            "transcripts.csv",
            "transcripts.tsv.gz",
            "transcripts.tsv",
        ]
    )
    seen: set[Path] = set()
    for base in candidate_dirs:
        if base.exists() and base.is_dir() and base.name.endswith(".pyxenium.slide.zarr"):
            return base
        for relpath in relpaths:
            candidate = (base / relpath).expanduser().resolve()
            if candidate in seen:
                continue
            seen.add(candidate)
            if candidate.exists():
                return candidate
    raise FileNotFoundError(
        f"Could not find a transcript table for sample {sample_id!r} below {xenium_root}."
    )


def _summarize_slice_transcripts(
    *,
    transform: SliceCellTransform,
    transcript_path: Path,
    vertical_qc_dir: Path,
    cfg: LocalZOrientationConfig,
) -> _SliceLocalZSummary:
    raw = _read_transcript_table(transcript_path)
    transcripts = _normalize_transcripts(raw, cfg)
    if transcripts.empty:
        raise ValueError(f"{transcript_path} does not contain valid transcript x/y/z rows.")

    sample_qc_dir = vertical_qc_dir / str(transform.sample_id)
    sample_qc_dir.mkdir(parents=True, exist_ok=True)
    ovrlpy_qc = _run_vertical_qc(transcripts, sample_qc_dir, cfg)
    filtered, excluded_count = _filter_doublet_regions(transcripts, ovrlpy_qc, cfg)
    signatures = _build_z_band_signatures(filtered, cfg)

    return _SliceLocalZSummary(
        sample_id=str(transform.sample_id),
        order=int(transform.order),
        stack_z_um=float(transform.z_um),
        transcript_path=transcript_path,
        z_min_um=float(transcripts["z_raw_um"].min()),
        z_max_um=float(transcripts["z_raw_um"].max()),
        z_mid_um=float((transcripts["z_raw_um"].min() + transcripts["z_raw_um"].max()) / 2.0),
        transcript_count=int(len(transcripts)),
        signature_transcript_count=int(len(filtered)),
        excluded_doublet_transcripts=int(excluded_count),
        low_band_count=int(signatures["low_count"]),
        high_band_count=int(signatures["high_count"]),
        raw_low_signature=signatures["low_signature"],
        raw_high_signature=signatures["high_signature"],
        ovrlpy_qc=ovrlpy_qc,
    )


def _read_transcript_table(path: Path) -> pd.DataFrame:
    suffixes = "".join(path.suffixes).lower()
    if path.is_dir() and path.name.endswith(".pyxenium.slide.zarr"):
        return _read_pyxenium_slide_transcripts(path)
    if suffixes.endswith(".parquet"):
        return pd.read_parquet(path)
    if suffixes.endswith(".tsv") or suffixes.endswith(".tsv.gz"):
        return pd.read_csv(path, sep="\t")
    if suffixes.endswith(".csv") or suffixes.endswith(".csv.gz"):
        return pd.read_csv(path)
    raise ValueError(f"Unsupported transcript table format: {path}")


def _read_pyxenium_slide_transcripts(path: Path) -> pd.DataFrame:
    try:
        from pyXenium import read_slide
    except Exception as exc:
        raise ImportError(
            "Reading pyXeniumSlide zarr transcript stores requires pyXenium."
        ) from exc

    slide = read_slide(path)
    try:
        transcripts = slide.points["transcripts"]
    except Exception as exc:
        fallback = _read_pyxenium_slide_sidecar_transcripts(path)
        if fallback is not None:
            return fallback
        raise ValueError(f"{path} does not contain a transcripts point source.") from exc

    frame = transcripts.copy()
    if "valid" in frame.columns:
        frame = frame.loc[frame["valid"].astype(bool)].copy()
    rename = {
        "gene_name": "gene",
        "quality_score": "qv",
    }
    present = {src: dst for src, dst in rename.items() if src in frame.columns and dst not in frame.columns}
    if present:
        frame = frame.rename(columns=present)
    return frame


def _read_pyxenium_slide_sidecar_transcripts(path: Path) -> pd.DataFrame | None:
    sample_id = _sample_id_from_transcript_source(path)
    candidate_dirs = [
        path,
        path.parent / sample_id,
    ]
    names = [
        "transcripts.parquet",
        "transcripts.csv.gz",
        "transcripts.csv",
        "transcripts.tsv.gz",
        "transcripts.tsv",
    ]
    for base in candidate_dirs:
        for name in names:
            candidate = base / name
            if candidate.exists() and candidate.is_file():
                return _read_transcript_table(candidate)
    return None


def _normalize_transcripts(frame: pd.DataFrame, cfg: LocalZOrientationConfig) -> pd.DataFrame:
    gene_col = _resolve_column(frame, cfg.gene_column, ["gene", "feature_name", "gene_name"], "gene")
    x_col = _resolve_column(frame, cfg.x_column, ["x", "x_um", "x_location", "x_centroid"], "x")
    y_col = _resolve_column(frame, cfg.y_column, ["y", "y_um", "y_location", "y_centroid"], "y")
    z_col = _resolve_column(frame, cfg.z_column, ["z", "z_um", "z_location"], "z")
    transcript_id_col = _optional_column(
        frame,
        cfg.transcript_id_column,
        ["transcript_id", "id", "molecule_id"],
    )
    structure_col = _optional_column(frame, cfg.structure_column, ["structure"])

    result = pd.DataFrame(
        {
            "gene": frame[gene_col].astype(str),
            "x_raw_um": pd.to_numeric(frame[x_col], errors="coerce"),
            "y_raw_um": pd.to_numeric(frame[y_col], errors="coerce"),
            "z_raw_um": pd.to_numeric(frame[z_col], errors="coerce"),
        }
    )
    if transcript_id_col is not None:
        result["transcript_id"] = frame[transcript_id_col].astype(str)
    if structure_col is not None:
        result["structure"] = frame[structure_col].astype(str)
    optional_columns = {
        "cell_id": _optional_column(frame, None, ["cell_id"]),
        "overlaps_nucleus": _optional_column(frame, None, ["overlaps_nucleus"]),
        "qv": _optional_column(frame, None, ["qv", "quality_score"]),
    }
    for optional, column in optional_columns.items():
        if column is not None:
            result[optional] = frame[column].to_numpy()

    finite = np.isfinite(result[["x_raw_um", "y_raw_um", "z_raw_um"]].to_numpy(dtype=float)).all(axis=1)
    non_empty_gene = result["gene"].astype(str).str.strip() != ""
    return result.loc[finite & non_empty_gene].reset_index(drop=True)


def _resolve_column(
    frame: pd.DataFrame,
    requested: str | None,
    candidates: Sequence[str],
    label: str,
) -> str:
    if requested is not None:
        if requested not in frame.columns:
            raise ValueError(f"Transcript table is missing requested {label} column {requested!r}.")
        return requested
    for candidate in candidates:
        if candidate in frame.columns:
            return candidate
    raise ValueError(
        f"Transcript table is missing a {label} column. Tried: {', '.join(candidates)}."
    )


def _optional_column(
    frame: pd.DataFrame,
    requested: str | None,
    candidates: Sequence[str],
) -> str | None:
    if requested is not None:
        if requested not in frame.columns:
            raise ValueError(f"Transcript table is missing requested column {requested!r}.")
        return requested
    for candidate in candidates:
        if candidate in frame.columns:
            return candidate
    return None


def _run_vertical_qc(
    transcripts: pd.DataFrame,
    sample_qc_dir: Path,
    cfg: LocalZOrientationConfig,
) -> _OvrlpyQC:
    sample_qc_dir.mkdir(parents=True, exist_ok=True)
    if cfg.vertical_qc_backend == "none":
        empty_doublets = sample_qc_dir / "doublets.csv"
        pd.DataFrame(columns=["x", "y"]).to_csv(empty_doublets, index=False)
        return _OvrlpyQC(status="skipped", doublets_csv=empty_doublets)

    cached_qc = _load_cached_vertical_qc(sample_qc_dir)
    if cached_qc is not None:
        return cached_qc

    try:
        import ovrlpy
    except Exception as exc:
        raise ImportError(
            "Local-z orientation was requested with vertical_qc_backend='ovrlpy', "
            "but ovrlpy is not installed. Install the optional vertical QC dependency "
            "or use vertical_qc_backend='none'."
        ) from exc

    coordinate_df = transcripts.loc[:, ["gene", "x_raw_um", "y_raw_um", "z_raw_um"]].rename(
        columns={"x_raw_um": "x", "y_raw_um": "y", "z_raw_um": "z"}
    )
    try:
        dataset = ovrlpy.Ovrlp(
            coordinate_df,
            KDE_bandwidth=float(cfg.ovrlpy_kde_bandwidth),
            n_components=int(cfg.ovrlpy_n_components),
            n_workers=int(cfg.ovrlpy_n_workers),
        )
        dataset.analyse(
            min_transcripts=float(cfg.ovrlpy_min_transcripts),
            fit_umap=bool(cfg.ovrlpy_fit_umap),
        )
    except Exception as exc:
        return _OvrlpyQC(status="analysis_failed", error=str(exc))

    integrity_value = getattr(dataset, "integrity_map", None)
    signal_value = getattr(dataset, "signal_map", None)
    integrity_path = _save_dataset_array(integrity_value, sample_qc_dir / "signal_integrity.npy")
    signal_path = _save_dataset_array(signal_value, sample_qc_dir / "signal_map.npy")
    integrity_stats = _array_stats(integrity_value)
    signal_stats = _array_stats(signal_value)

    try:
        doublets = dataset.detect_doublets(
            min_signal=float(cfg.doublet_min_signal),
            integrity_sigma=float(cfg.doublet_integrity_sigma),
        )
        doublets_df = _as_pandas_frame(doublets)
    except Exception as exc:
        doublets_df = pd.DataFrame([{"error": str(exc)}])
        status = "doublet_detection_failed"
    else:
        status = "ok"

    doublets_path = sample_qc_dir / "doublets.csv"
    doublets_df.to_csv(doublets_path, index=False)
    return _OvrlpyQC(
        status=status,
        doublets_csv=doublets_path,
        doublet_count=int(len(doublets_df)) if {"x", "y"}.issubset(doublets_df.columns) else 0,
        signal_integrity_npy=integrity_path,
        signal_map_npy=signal_path,
        signal_integrity_mean=integrity_stats["mean"],
        signal_integrity_p05=integrity_stats["p05"],
        signal_integrity_p95=integrity_stats["p95"],
        signal_mean=signal_stats["mean"],
        signal_p05=signal_stats["p05"],
        signal_p95=signal_stats["p95"],
    )


def _load_cached_vertical_qc(sample_qc_dir: Path) -> _OvrlpyQC | None:
    integrity_path = sample_qc_dir / "signal_integrity.npy"
    signal_path = sample_qc_dir / "signal_map.npy"
    doublets_path = sample_qc_dir / "doublets.csv"
    if not (integrity_path.exists() and signal_path.exists() and doublets_path.exists()):
        return None
    try:
        integrity = np.load(integrity_path, mmap_mode="r")
        signal = np.load(signal_path, mmap_mode="r")
        doublets = pd.read_csv(doublets_path)
    except Exception:
        return None
    integrity_stats = _array_stats(integrity)
    signal_stats = _array_stats(signal)
    return _OvrlpyQC(
        status="ok",
        doublets_csv=doublets_path,
        doublet_count=int(len(doublets)) if {"x", "y"}.issubset(doublets.columns) else 0,
        signal_integrity_npy=integrity_path,
        signal_map_npy=signal_path,
        signal_integrity_mean=integrity_stats["mean"],
        signal_integrity_p05=integrity_stats["p05"],
        signal_integrity_p95=integrity_stats["p95"],
        signal_mean=signal_stats["mean"],
        signal_p05=signal_stats["p05"],
        signal_p95=signal_stats["p95"],
    )


def _save_dataset_array(value: Any, output_path: Path) -> Path | None:
    if value is None:
        return None
    try:
        np.save(output_path, np.asarray(value))
    except Exception:
        return None
    return output_path


def _array_stats(value: Any) -> dict[str, float]:
    if value is None:
        return {"mean": math.nan, "p05": math.nan, "p95": math.nan}
    try:
        array = np.asarray(value, dtype=float)
    except Exception:
        return {"mean": math.nan, "p05": math.nan, "p95": math.nan}
    finite = array[np.isfinite(array)]
    if finite.size == 0:
        return {"mean": math.nan, "p05": math.nan, "p95": math.nan}
    return {
        "mean": float(np.mean(finite)),
        "p05": float(np.quantile(finite, 0.05)),
        "p95": float(np.quantile(finite, 0.95)),
    }


def _as_pandas_frame(value: Any) -> pd.DataFrame:
    if isinstance(value, pd.DataFrame):
        return value.copy()
    if hasattr(value, "to_pandas"):
        return value.to_pandas()
    return pd.DataFrame(value)


def _filter_doublet_regions(
    transcripts: pd.DataFrame,
    ovrlpy_qc: _OvrlpyQC,
    cfg: LocalZOrientationConfig,
) -> tuple[pd.DataFrame, int]:
    if (
        cfg.doublet_exclusion_radius_um <= 0
        or ovrlpy_qc.doublets_csv is None
        or not ovrlpy_qc.doublets_csv.exists()
    ):
        return transcripts, 0
    doublets = pd.read_csv(ovrlpy_qc.doublets_csv)
    if doublets.empty:
        return transcripts, 0
    x_col = _optional_column(doublets, None, ["x", "x_um", "x_location"])
    y_col = _optional_column(doublets, None, ["y", "y_um", "y_location"])
    if x_col is None or y_col is None:
        return transcripts, 0
    xy = transcripts[["x_raw_um", "y_raw_um"]].to_numpy(dtype=float)
    doublet_xy = doublets[[x_col, y_col]].apply(pd.to_numeric, errors="coerce").dropna().to_numpy(dtype=float)
    if len(doublet_xy) == 0:
        return transcripts, 0
    tree = cKDTree(doublet_xy)
    near_doublet = np.asarray(
        tree.query_ball_point(xy, r=float(cfg.doublet_exclusion_radius_um), return_length=True)
    ) > 0
    kept = transcripts.loc[~near_doublet].reset_index(drop=True)
    if kept.empty:
        return transcripts, 0
    return kept, int(near_doublet.sum())


def _build_z_band_signatures(
    transcripts: pd.DataFrame,
    cfg: LocalZOrientationConfig,
) -> dict[str, Any]:
    z = transcripts["z_raw_um"].to_numpy(dtype=float)
    low_cut = float(np.quantile(z, cfg.z_band_fraction))
    high_cut = float(np.quantile(z, 1.0 - cfg.z_band_fraction))
    low = transcripts.loc[transcripts["z_raw_um"] <= low_cut]
    high = transcripts.loc[transcripts["z_raw_um"] >= high_cut]
    if len(low) < cfg.min_band_transcripts:
        low = transcripts.nsmallest(cfg.min_band_transcripts, "z_raw_um")
    if len(high) < cfg.min_band_transcripts:
        high = transcripts.nlargest(cfg.min_band_transcripts, "z_raw_um")
    return {
        "low_count": int(len(low)),
        "high_count": int(len(high)),
        "low_signature": _gene_signature(low["gene"], cfg.max_signature_genes),
        "high_signature": _gene_signature(high["gene"], cfg.max_signature_genes),
    }


def _gene_signature(genes: pd.Series, max_genes: int) -> dict[str, float]:
    counts = genes.astype(str).value_counts().head(int(max_genes))
    total = float(counts.sum())
    if total <= 0:
        return {}
    return {str(gene): float(count / total) for gene, count in counts.items()}


def _infer_orientation_sequence(
    summaries: Sequence[_SliceLocalZSummary],
    cfg: LocalZOrientationConfig,
) -> _OrientationInference:
    if not summaries:
        raise ValueError("At least one slice is required for local-z orientation inference.")
    if len(summaries) == 1:
        sample_id = summaries[0].sample_id
        return _OrientationInference(
            states_by_sample={sample_id: "preserve"},
            best_score=0.0,
            second_score=0.0,
            confidence_margin=0.0,
            low_confidence=True,
            edge_scores={},
        )

    paths: dict[str, list[tuple[float, tuple[str, ...]]]] = {
        state: [(0.0, (state,))] for state in LOCAL_Z_ORIENTATION_STATES
    }
    for index in range(1, len(summaries)):
        next_paths: dict[str, list[tuple[float, tuple[str, ...]]]] = {}
        for state in LOCAL_Z_ORIENTATION_STATES:
            candidates: list[tuple[float, tuple[str, ...]]] = []
            for prev_state in LOCAL_Z_ORIENTATION_STATES:
                edge = _edge_score(summaries[index - 1], prev_state, summaries[index], state)
                for prev_score, prev_path in paths[prev_state]:
                    candidates.append((float(prev_score + edge), prev_path + (state,)))
            next_paths[state] = _top_unique_paths(candidates, limit=2)
        paths = next_paths

    all_paths = _top_unique_paths([candidate for group in paths.values() for candidate in group], limit=2)
    best_score, best_path = all_paths[0]
    second_score = all_paths[1][0] if len(all_paths) > 1 else best_score
    edge_count = max(len(summaries) - 1, 1)
    average_best = float(best_score) / edge_count
    average_second = float(second_score) / edge_count
    margin = average_best - average_second

    states_by_sample = {
        summary.sample_id: best_path[index] for index, summary in enumerate(summaries)
    }
    edge_scores: dict[tuple[str, str], float] = {}
    for index in range(1, len(summaries)):
        prev = summaries[index - 1]
        cur = summaries[index]
        edge_scores[(prev.sample_id, cur.sample_id)] = _edge_score(
            prev,
            states_by_sample[prev.sample_id],
            cur,
            states_by_sample[cur.sample_id],
        )
    return _OrientationInference(
        states_by_sample=states_by_sample,
        best_score=average_best,
        second_score=average_second,
        confidence_margin=float(margin),
        low_confidence=bool(margin < cfg.low_confidence_margin),
        edge_scores=edge_scores,
    )


def _top_unique_paths(
    candidates: Sequence[tuple[float, tuple[str, ...]]],
    *,
    limit: int,
) -> list[tuple[float, tuple[str, ...]]]:
    unique: dict[tuple[str, ...], float] = {}
    for score, path in candidates:
        previous = unique.get(path)
        if previous is None or score > previous:
            unique[path] = float(score)
    ordered = sorted(
        ((score, path) for path, score in unique.items()),
        key=lambda item: (-item[0], item[1]),
    )
    return ordered[:limit]


def _edge_score(
    prev_summary: _SliceLocalZSummary,
    prev_state: str,
    next_summary: _SliceLocalZSummary,
    next_state: str,
) -> float:
    prev_high = _surface_signature(prev_summary, prev_state, "high")
    next_low = _surface_signature(next_summary, next_state, "low")
    return _cosine_similarity(prev_high, next_low)


def _surface_signature(
    summary: _SliceLocalZSummary,
    state: str,
    side: str,
) -> dict[str, float]:
    if state not in LOCAL_Z_ORIENTATION_STATES:
        raise ValueError(f"Unsupported local-z orientation state: {state!r}.")
    if side == "low":
        return summary.raw_low_signature if state == "preserve" else summary.raw_high_signature
    if side == "high":
        return summary.raw_high_signature if state == "preserve" else summary.raw_low_signature
    raise ValueError("side must be 'low' or 'high'.")


def _cosine_similarity(first: Mapping[str, float], second: Mapping[str, float]) -> float:
    genes = sorted(set(first) | set(second))
    if not genes:
        return 0.0
    a = np.array([float(first.get(gene, 0.0)) for gene in genes], dtype=float)
    b = np.array([float(second.get(gene, 0.0)) for gene in genes], dtype=float)
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denom <= 0:
        return 0.0
    return float(np.dot(a, b) / denom)


def _build_manifest_frame(
    summaries: Sequence[_SliceLocalZSummary],
    inference: _OrientationInference,
    cfg: LocalZOrientationConfig,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for index, summary in enumerate(summaries):
        selected = inference.states_by_sample[summary.sample_id]
        previous_edge = math.nan
        next_edge = math.nan
        if index > 0:
            previous = summaries[index - 1]
            previous_edge = inference.edge_scores.get((previous.sample_id, summary.sample_id), math.nan)
        if index + 1 < len(summaries):
            nxt = summaries[index + 1]
            next_edge = inference.edge_scores.get((summary.sample_id, nxt.sample_id), math.nan)
        rows.append(
            {
                "schema_version": LOCAL_Z_ORIENTATION_SCHEMA_VERSION,
                "sample_id": summary.sample_id,
                "order": summary.order,
                "stack_z_um": summary.stack_z_um,
                "transcript_path": str(summary.transcript_path),
                "selected_orientation": selected,
                "applied": bool(cfg.apply_local_z_flip),
                "z_min_um": summary.z_min_um,
                "z_max_um": summary.z_max_um,
                "z_mid_um": summary.z_mid_um,
                "transcript_count": summary.transcript_count,
                "signature_transcript_count": summary.signature_transcript_count,
                "excluded_doublet_transcripts": summary.excluded_doublet_transcripts,
                "low_band_count": summary.low_band_count,
                "high_band_count": summary.high_band_count,
                "ovrlpy_status": summary.ovrlpy_qc.status,
                "ovrlpy_doublet_count": summary.ovrlpy_qc.doublet_count,
                "signal_integrity_mean": summary.ovrlpy_qc.signal_integrity_mean,
                "signal_integrity_p05": summary.ovrlpy_qc.signal_integrity_p05,
                "signal_integrity_p95": summary.ovrlpy_qc.signal_integrity_p95,
                "signal_mean": summary.ovrlpy_qc.signal_mean,
                "signal_p05": summary.ovrlpy_qc.signal_p05,
                "signal_p95": summary.ovrlpy_qc.signal_p95,
                "best_score": inference.best_score,
                "second_score": inference.second_score,
                "confidence_margin": inference.confidence_margin,
                "low_confidence": inference.low_confidence,
                "edge_score_to_previous": previous_edge,
                "edge_score_to_next": next_edge,
            }
        )
    return pd.DataFrame(rows)


def _write_aligned_transcripts(
    *,
    summaries: Sequence[_SliceLocalZSummary],
    transforms: Sequence[SliceCellTransform],
    inference: _OrientationInference,
    output_parquet: Path,
    cfg: LocalZOrientationConfig,
) -> int:
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except Exception as exc:  # pragma: no cover - pyarrow is a package dependency.
        raise ImportError("Writing aligned transcript Parquet requires pyarrow.") from exc

    transforms_by_sample = {str(transform.sample_id): transform for transform in transforms}
    output_parquet.parent.mkdir(parents=True, exist_ok=True)
    writer: Any | None = None
    total = 0
    try:
        for summary in summaries:
            transform = transforms_by_sample[summary.sample_id]
            frame = _normalize_transcripts(_read_transcript_table(summary.transcript_path), cfg)
            aligned = _aligned_transcript_frame(
                frame,
                summary=summary,
                transform=transform,
                orientation=inference.states_by_sample[summary.sample_id],
                cfg=cfg,
            )
            table = pa.Table.from_pandas(aligned, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(output_parquet, table.schema)
            writer.write_table(table)
            total += int(len(aligned))
    finally:
        if writer is not None:
            writer.close()

    if writer is None:
        pd.DataFrame().to_parquet(output_parquet, index=False)
    return total


def _aligned_transcript_frame(
    frame: pd.DataFrame,
    *,
    summary: _SliceLocalZSummary,
    transform: SliceCellTransform,
    orientation: str,
    cfg: LocalZOrientationConfig,
) -> pd.DataFrame:
    xy_native = frame[["x_raw_um", "y_raw_um"]].to_numpy(dtype=float) / float(cfg.pixel_size_um)
    if transform.hard is not None:
        xy_native = apply_similarity_to_points(xy_native, transform.hard)
    if transform.tps is not None:
        xy_native = transform.tps.warp(xy_native, chunk_size=cfg.chunk_size)
    xy_um = xy_native * float(cfg.pixel_size_um)

    z_raw = frame["z_raw_um"].to_numpy(dtype=float)
    should_reverse = bool(cfg.apply_local_z_flip and orientation == "reverse")
    if should_reverse:
        local_z = summary.z_min_um + summary.z_max_um - z_raw
    else:
        local_z = z_raw.copy()
    z_3d = float(summary.stack_z_um) + (local_z - float(summary.z_mid_um))

    output = pd.DataFrame(
        {
            "sample_id": summary.sample_id,
            "slice_order": int(summary.order),
            "gene": frame["gene"].astype(str).to_numpy(),
            "x_raw_um": frame["x_raw_um"].to_numpy(dtype=float),
            "y_raw_um": frame["y_raw_um"].to_numpy(dtype=float),
            "z_raw_um": z_raw,
            "x_3d_um": xy_um[:, 0],
            "y_3d_um": xy_um[:, 1],
            "z_3d_um": z_3d,
            "slice_stack_z_um": float(summary.stack_z_um),
            "local_z_corrected_um": local_z,
            "local_z_mid_um": float(summary.z_mid_um),
            "local_z_orientation": orientation,
            "local_z_flip_applied": should_reverse,
        }
    )
    row_count = len(output)
    for optional in ALIGNED_TRANSCRIPT_OPTIONAL_COLUMNS:
        if optional in {"transcript_id", "cell_id", "structure"}:
            if optional in frame.columns:
                output[optional] = frame[optional].astype("string").reset_index(drop=True)
            else:
                output[optional] = pd.array([pd.NA] * row_count, dtype="string")
        elif optional == "overlaps_nucleus":
            if optional in frame.columns:
                output[optional] = (
                    pd.to_numeric(frame[optional], errors="coerce")
                    .fillna(0)
                    .astype("int64")
                    .to_numpy()
                )
            else:
                output[optional] = np.zeros(row_count, dtype="int64")
        elif optional == "qv":
            if optional in frame.columns:
                output[optional] = pd.to_numeric(frame[optional], errors="coerce").astype("float64").to_numpy()
            else:
                output[optional] = np.full(row_count, np.nan, dtype="float64")
    return output.loc[:, list(ALIGNED_TRANSCRIPT_COLUMNS)]


def _build_marker_gradient_frame(
    summaries: Sequence[_SliceLocalZSummary],
    inference: _OrientationInference,
    *,
    max_markers_per_slice: int = 12,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for summary in summaries:
        state = inference.states_by_sample[summary.sample_id]
        bottom = _surface_signature(summary, state, "low")
        top = _surface_signature(summary, state, "high")
        genes = sorted(set(bottom) | set(top))
        ranked = sorted(
            genes,
            key=lambda gene: abs(float(top.get(gene, 0.0)) - float(bottom.get(gene, 0.0))),
            reverse=True,
        )[:max_markers_per_slice]
        for gene in ranked:
            rows.append(
                {
                    "sample_id": summary.sample_id,
                    "order": summary.order,
                    "gene": gene,
                    "corrected_top_fraction": float(top.get(gene, 0.0)),
                    "corrected_bottom_fraction": float(bottom.get(gene, 0.0)),
                    "top_minus_bottom": float(top.get(gene, 0.0)) - float(bottom.get(gene, 0.0)),
                }
            )
    return pd.DataFrame(rows)


def _write_biological_report(
    output_path: Path,
    *,
    manifest: pd.DataFrame,
    marker_gradients: pd.DataFrame,
    inference: _OrientationInference,
    cfg: LocalZOrientationConfig,
) -> None:
    manifest_rows = "\n".join(
        "<tr>"
        f"<td>{escape(str(row.sample_id))}</td>"
        f"<td>{int(row.order)}</td>"
        f"<td>{escape(str(row.selected_orientation))}</td>"
        f"<td>{escape(str(row.ovrlpy_status))}</td>"
        f"<td>{int(row.ovrlpy_doublet_count)}</td>"
        f"<td>{float(row.signal_integrity_mean):.4g}</td>"
        f"<td>{float(row.signal_mean):.4g}</td>"
        f"<td>{float(row.edge_score_to_previous):.3f}</td>"
        f"<td>{float(row.edge_score_to_next):.3f}</td>"
        "</tr>"
        for row in manifest.itertuples(index=False)
    )
    marker_rows = "\n".join(
        "<tr>"
        f"<td>{escape(str(row.sample_id))}</td>"
        f"<td>{escape(str(row.gene))}</td>"
        f"<td>{float(row.corrected_top_fraction):.3f}</td>"
        f"<td>{float(row.corrected_bottom_fraction):.3f}</td>"
        f"<td>{float(row.top_minus_bottom):.3f}</td>"
        "</tr>"
        for row in marker_gradients.head(60).itertuples(index=False)
    )
    low_confidence_note = (
        "The orientation sequence is marked low confidence; forced correction still used "
        "the best-scoring preserve/reverse path."
        if inference.low_confidence
        else "The best orientation path exceeded the configured confidence margin."
    )
    z3d_min = (
        manifest["stack_z_um"].astype(float)
        + manifest["z_min_um"].astype(float)
        - manifest["z_mid_um"].astype(float)
    ).min()
    z3d_max = (
        manifest["stack_z_um"].astype(float)
        + manifest["z_max_um"].astype(float)
        - manifest["z_mid_um"].astype(float)
    ).max()
    total_doublets = int(pd.to_numeric(manifest["ovrlpy_doublet_count"], errors="coerce").fillna(0).sum())
    integrity_median = pd.to_numeric(
        manifest["signal_integrity_mean"], errors="coerce"
    ).median()
    html = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>HistoSeg local-z biological QC</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 28px; color: #17202a; }}
    table {{ border-collapse: collapse; margin: 18px 0; width: 100%; }}
    th, td {{ border: 1px solid #d5dbe3; padding: 6px 8px; text-align: left; }}
    th {{ background: #eef3f8; }}
    .note {{ background: #fff8dc; border: 1px solid #ead27a; padding: 10px; }}
    code {{ background: #f2f4f7; padding: 1px 4px; }}
  </style>
</head>
<body>
  <h1>HistoSeg local-z biological QC</h1>
  <p>
    This report summarizes transcript-level local-z orientation correction.
    Contour and cell stack reconstruction are unchanged; only transcript
    <code>local_z_corrected_um</code> and <code>z_3d_um</code> use the selected
    preserve/reverse states.
  </p>
  <p class="note">{escape(low_confidence_note)}</p>
  <ul>
    <li>Vertical QC backend: {escape(cfg.vertical_qc_backend)}</li>
    <li>Best adjacent-surface continuity score: {inference.best_score:.4f}</li>
    <li>Second-best score: {inference.second_score:.4f}</li>
    <li>Confidence margin: {inference.confidence_margin:.4f}</li>
    <li>Local-z flip applied: {str(bool(cfg.apply_local_z_flip)).lower()}</li>
    <li>Total Ovrlpy doublet calls: {total_doublets}</li>
    <li>Median signal integrity mean: {float(integrity_median):.4g}</li>
    <li>Transcript z_3d range: {float(z3d_min):.3f} to {float(z3d_max):.3f} um</li>
  </ul>

  <h2>Per-slice orientation and vertical QC</h2>
  <table>
    <thead>
      <tr>
        <th>Sample</th><th>Order</th><th>Orientation</th><th>Ovrlpy status</th>
        <th>Doublets</th><th>Integrity mean</th><th>Signal mean</th>
        <th>Prev continuity</th><th>Next continuity</th>
      </tr>
    </thead>
    <tbody>{manifest_rows}</tbody>
  </table>

  <h2>Top-vs-bottom marker gradients</h2>
  <table>
    <thead>
      <tr>
        <th>Sample</th><th>Gene</th><th>Corrected top fraction</th>
        <th>Corrected bottom fraction</th><th>Top minus bottom</th>
      </tr>
    </thead>
    <tbody>{marker_rows}</tbody>
  </table>

  <p>
    Interpret these outputs as QC and hypothesis-generation evidence for
    vertical overlaps, tissue folds, doublets, and layer-like molecular signal.
    They are not standalone proof of biological structure.
  </p>
</body>
</html>
"""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html, encoding="utf-8")


__all__ = [
    "LOCAL_Z_ORIENTATION_SCHEMA_VERSION",
    "LOCAL_Z_ORIENTATION_STATES",
    "LocalZOrientationConfig",
    "LocalZOrientationResult",
    "run_local_z_orientation_correction",
]
