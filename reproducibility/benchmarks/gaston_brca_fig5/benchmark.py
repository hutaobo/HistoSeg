"""Benchmark HistoSeg against GASTON Fig. 5 on Xenium breast cancer Rep1.

The module is intentionally self-contained so it can be copied to an A100
project folder and run with ``python -m`` without changing HistoSeg public
package entry points.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time
from typing import Iterable, Mapping, Sequence
from urllib.error import URLError
from urllib.request import urlopen

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import pearsonr, spearmanr
from shapely.geometry import Point, Polygon, shape
from shapely.ops import unary_union


DEFAULT_XENIUM_OUTS = Path(
    r"Y:\long\10X_datasets\Xenium\Xenium_Breast_Cancer"
    r"\Xenium_FFPE_Human_Breast_Cancer_Rep1_outs"
)
DEFAULT_OUT_DIR = Path(r"D:\histoseg_gaston_brca_fig5_benchmark")
REMOTE_OUT_DIR = "/data/taobo.hu/projects/histoseg_gaston_brca_fig5_benchmark"

TUTORIAL_SOURCE_URL = (
    "https://gaston.readthedocs.io/en/latest/notebooks/tutorials/"
    "xenium_brca_tutorial.html"
)
NATURE_FIG5_URL = "https://www.nature.com/articles/s41592-024-02503-3/figures/5"
GASTON_REPO_URL = "https://github.com/raphael-group/GASTON"

CROP_BOUNDS = (1500.0, 3250.0, 2000.0, 4000.0)
GASTON_DOMAIN_NAMES = {
    0: "Invasive tumor",
    1: "Stromal",
    2: "Immune",
    3: "DCIS/myoepithelial/perivascular",
}
GASTON_DOMAIN_COLORS = {
    "Invasive tumor": "#e9967a",
    "Stromal": "#00bfff",
    "Immune": "#32cd32",
    "DCIS/myoepithelial/perivascular": "#ff0000",
}
HISTOSEG_STRUCTURE_SPECS = [
    {
        "structure_id": 1,
        "structure_name": "Invasive tumor",
        "structure_color": "#e9967a",
        "cluster_ids": [1, 3, 7, 9, 12, 16, 17, 18, 19],
    },
    {
        "structure_id": 2,
        "structure_name": "Stromal",
        "structure_color": "#00bfff",
        "cluster_ids": [2, 4, 10, 20],
    },
    {
        "structure_id": 3,
        "structure_name": "Immune",
        "structure_color": "#32cd32",
        "cluster_ids": [5, 6, 14, 15],
    },
    {
        "structure_id": 4,
        "structure_name": "DCIS/myoepithelial/perivascular",
        "structure_color": "#ff0000",
        "cluster_ids": [8, 11, 13],
    },
]
DOMAIN_NAME_TO_ID = {str(spec["structure_name"]): int(spec["structure_id"]) for spec in HISTOSEG_STRUCTURE_SPECS}
GASTON_Q_VALS = [0.1, 0.1, 0.1, 0.05]
CELL_TYPE_PANEL_ORDER = [
    "CD8+_T_Cells",
    "CD4+_T_Cells",
    "Invasive_Tumor",
    "DCIS_2",
    "Myoepi_ACTA2+",
    "Stromal",
    "B_Cells",
]
GENE_GRADIENT_GENES = ["SELL", "TCL1A"]


class BenchmarkInputError(RuntimeError):
    """Raised when a required benchmark input is missing or inconsistent."""


@dataclass(frozen=True)
class CropBounds:
    xmin: float = CROP_BOUNDS[0]
    xmax: float = CROP_BOUNDS[1]
    ymin: float = CROP_BOUNDS[2]
    ymax: float = CROP_BOUNDS[3]

    def contains(self, cells: pd.DataFrame) -> pd.Series:
        return (
            (cells["x_centroid"] > self.xmin)
            & (cells["x_centroid"] < self.xmax)
            & (cells["y_centroid"] > self.ymin)
            & (cells["y_centroid"] < self.ymax)
        )


def parse_seed_spec(value: str | Sequence[int]) -> list[int]:
    if isinstance(value, str):
        seeds: list[int] = []
        for part in value.split(","):
            token = part.strip()
            if not token:
                continue
            if "-" in token:
                start_s, end_s = token.split("-", 1)
                seeds.extend(range(int(start_s), int(end_s) + 1))
            else:
                seeds.append(int(token))
        return sorted(dict.fromkeys(seeds))
    return sorted(dict.fromkeys(int(seed) for seed in value))


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_cells(xenium_outs: Path) -> pd.DataFrame:
    path = xenium_outs / "cells.parquet"
    if not path.exists():
        raise BenchmarkInputError(f"Missing Xenium cells table: {path}")
    cells = pd.read_parquet(path)
    required = {"cell_id", "x_centroid", "y_centroid"}
    missing = required.difference(cells.columns)
    if missing:
        raise BenchmarkInputError(f"cells.parquet missing columns: {sorted(missing)}")
    cells = cells.copy()
    cells["Barcode"] = cells["cell_id"].astype(str)
    return cells


def select_adh_crop(cells: pd.DataFrame, bounds: CropBounds) -> pd.DataFrame:
    crop = cells.loc[bounds.contains(cells)].copy()
    if crop.empty:
        raise BenchmarkInputError(f"ADH crop is empty for bounds {asdict(bounds)}")
    crop["crop_index"] = np.arange(len(crop), dtype=np.int64)
    return crop


def load_graphclust(xenium_outs: Path) -> pd.DataFrame:
    path = xenium_outs / "analysis" / "clustering" / "gene_expression_graphclust" / "clusters.csv"
    if not path.exists():
        raise BenchmarkInputError(f"Missing GraphClust clusters.csv: {path}")
    clusters = pd.read_csv(path)
    if not {"Barcode", "Cluster"}.issubset(clusters.columns):
        raise BenchmarkInputError(f"GraphClust table must contain Barcode and Cluster: {path}")
    clusters = clusters[["Barcode", "Cluster"]].copy()
    clusters["Barcode"] = clusters["Barcode"].astype(str)
    clusters["Cluster"] = clusters["Cluster"].astype(str)
    return clusters


def assign_histoseg_seed_domains(crop_cells: pd.DataFrame, graphclust: pd.DataFrame) -> pd.DataFrame:
    cluster_to_domain: dict[str, tuple[int, str]] = {}
    for spec in HISTOSEG_STRUCTURE_SPECS:
        for cluster_id in spec["cluster_ids"]:
            cluster_to_domain[str(cluster_id)] = (int(spec["structure_id"]), str(spec["structure_name"]))

    graphclust = graphclust.copy()
    graphclust["Barcode"] = graphclust["Barcode"].astype(str)
    graphclust["Cluster"] = graphclust["Cluster"].astype(str)
    merged = crop_cells.merge(graphclust, on="Barcode", how="left", validate="one_to_one")
    if merged["Cluster"].isna().any():
        missing = int(merged["Cluster"].isna().sum())
        raise BenchmarkInputError(f"{missing} ADH crop cells did not align to GraphClust clusters")

    domain_payload = merged["Cluster"].map(cluster_to_domain)
    if domain_payload.isna().any():
        missing_clusters = sorted(merged.loc[domain_payload.isna(), "Cluster"].unique().tolist())
        raise BenchmarkInputError(f"Unmapped GraphClust clusters: {missing_clusters}")

    merged["histoseg_seed_structure_id"] = [int(item[0]) for item in domain_payload]
    merged["histoseg_seed_structure_name"] = [str(item[1]) for item in domain_payload]
    return merged


def download_cell_type_labels(urls: Sequence[str], out_dir: Path) -> Path:
    if not urls:
        raise BenchmarkInputError(
            "No public Cell_Barcode_Type_Matrices.csv URL was provided. "
            "Pass --cell-type-url with an exact public source."
        )
    ensure_dir(out_dir)
    errors: list[str] = []
    for url in urls:
        label_path = out_dir / "Cell_Barcode_Type_Matrices.csv"
        try:
            with urlopen(url, timeout=60) as response:
                payload = response.read()
            label_path.write_bytes(payload)
            return label_path
        except (OSError, URLError) as exc:
            errors.append(f"{url}: {exc}")
    raise BenchmarkInputError("Could not fetch public cell-type labels:\n" + "\n".join(errors))


def validate_cell_type_labels(label_csv: Path, all_cells: pd.DataFrame, crop_cells: pd.DataFrame) -> pd.DataFrame:
    labels = pd.read_csv(label_csv)
    if "Barcode" not in labels.columns or "Cluster" not in labels.columns:
        raise BenchmarkInputError(
            "Cell_Barcode_Type_Matrices.csv must contain Barcode and Cluster columns "
            "for exact barcode alignment."
        )
    labels = labels[["Barcode", "Cluster"]].copy()
    labels["Barcode"] = labels["Barcode"].astype(str)
    labels["cell_type"] = labels["Cluster"].astype(str)
    if labels["Barcode"].duplicated().any():
        duplicates = labels.loc[labels["Barcode"].duplicated(), "Barcode"].head().tolist()
        raise BenchmarkInputError(f"Cell-type label file contains duplicate barcodes: {duplicates}")

    all_barcodes = set(all_cells["Barcode"].astype(str))
    label_barcodes = set(labels["Barcode"].astype(str))
    missing_all = all_barcodes.difference(label_barcodes)
    if missing_all:
        raise BenchmarkInputError(
            f"Cell-type labels do not cover all Xenium cells; missing {len(missing_all)} barcodes"
        )

    merged = crop_cells[["Barcode"]].merge(labels[["Barcode", "cell_type"]], on="Barcode", how="left")
    if merged["cell_type"].isna().any():
        raise BenchmarkInputError("Cell-type labels failed to align to all ADH crop cells")
    return merged


def _read_10x_h5_sparse(path: Path) -> tuple[sparse.csr_matrix, np.ndarray, np.ndarray]:
    import h5py

    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"]
        data = matrix["data"][:]
        indices = matrix["indices"][:]
        indptr = matrix["indptr"][:]
        shape = tuple(int(item) for item in matrix["shape"][:])
        barcodes = np.asarray([item.decode() for item in matrix["barcodes"][:]])
        gene_labels = np.asarray([item.decode() for item in matrix["features/name"][:]])
    # 10x HDF5 stores features x cells in CSC-like form; transpose to cells x genes.
    counts = sparse.csc_matrix((data, indices, indptr), shape=shape).T.tocsr()
    return counts, barcodes, gene_labels


def load_counts_for_crop(
    xenium_outs: Path,
    crop_cells: pd.DataFrame,
) -> tuple[sparse.csr_matrix, np.ndarray, np.ndarray]:
    counts, barcodes, gene_labels = _read_10x_h5_sparse(xenium_outs / "cell_feature_matrix.h5")
    barcode_to_pos = pd.Series(np.arange(len(barcodes)), index=barcodes.astype(str))
    missing = set(crop_cells["Barcode"].astype(str)).difference(barcode_to_pos.index)
    if missing:
        raise BenchmarkInputError(f"{len(missing)} crop barcodes missing from cell_feature_matrix.h5")
    positions = barcode_to_pos.loc[crop_cells["Barcode"].astype(str)].to_numpy()
    return counts[positions].tocsr(), gene_labels, barcodes[positions]


def compute_feature_matrix(
    counts: sparse.csr_matrix,
    *,
    method: str,
    num_dims: int,
    glmpca_iters: int,
    pearson_clip: float,
    random_state: int,
) -> np.ndarray:
    method = method.lower()
    if method == "glmpca":
        try:
            from glmpca import glmpca
        except ImportError as exc:
            raise BenchmarkInputError("glmpca is required for --feature-method glmpca") from exc
        sums = np.asarray(counts.sum(axis=0)).ravel()
        top = np.argsort(sums)[-min(10000, counts.shape[1]) :]
        counts_dense = np.asarray(counts[:, top].todense())
        result = glmpca.glmpca(
            counts_dense.T,
            num_dims,
            fam="poi",
            penalty=1,
            verbose=True,
            ctl={"maxIter": int(glmpca_iters), "eps": 1e-4, "optimizeTheta": True},
        )
        return np.asarray(result["factors"], dtype=np.float32)
    if method == "pearson":
        try:
            from gaston import parse_adata
        except ImportError as exc:
            raise BenchmarkInputError("gaston-spatial is required for --feature-method pearson") from exc
        coords_placeholder = np.column_stack(
            [np.arange(counts.shape[0], dtype=float), np.zeros(counts.shape[0], dtype=float)]
        )
        return np.asarray(
            parse_adata.get_top_pearson_residuals(
                num_dims,
                np.asarray(counts.todense()),
                coords_placeholder,
                clip=float(pearson_clip),
            ),
            dtype=np.float32,
        )
    if method == "log1p-pca":
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler

        x = np.log1p(np.asarray(counts.todense()))
        x = StandardScaler(with_mean=True, with_std=True).fit_transform(x)
        return PCA(n_components=num_dims, random_state=random_state).fit_transform(x).astype(np.float32)
    raise ValueError(f"Unknown feature method: {method}")


def gene_count_table(counts: sparse.csr_matrix, gene_labels: np.ndarray, genes: Sequence[str]) -> pd.DataFrame:
    gene_to_idx = {str(gene): index for index, gene in enumerate(gene_labels)}
    payload: dict[str, np.ndarray] = {}
    for gene in genes:
        if gene not in gene_to_idx:
            raise BenchmarkInputError(f"Required gene {gene!r} is absent from Xenium panel")
        payload[gene] = np.asarray(counts[:, gene_to_idx[gene]].todense()).ravel()
    return pd.DataFrame(payload)


def write_inputs_manifest(out_dir: Path, payload: Mapping[str, object]) -> Path:
    manifest_path = out_dir / "inputs_manifest.json"
    manifest_path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")
    return manifest_path


def prepare_inputs(
    *,
    xenium_outs: Path,
    out_dir: Path,
    bounds: CropBounds,
    feature_method: str,
    num_dims: int,
    glmpca_iters: int,
    pearson_clip: float,
    cell_type_urls: Sequence[str],
    cell_type_csv: Path | None,
    require_cell_types: bool,
    random_state: int = 0,
) -> dict[str, Path]:
    t0 = time.perf_counter()
    data_dir = ensure_dir(out_dir / "xenium_tumor_data")
    raw_dir = ensure_dir(out_dir / "public_sources")
    all_cells = load_cells(xenium_outs)
    crop_cells = select_adh_crop(all_cells, bounds)
    graphclust = load_graphclust(xenium_outs)
    seed_cells = assign_histoseg_seed_domains(crop_cells, graphclust)
    counts, gene_labels, bounded_barcodes = load_counts_for_crop(xenium_outs, crop_cells)
    features = compute_feature_matrix(
        counts,
        method=feature_method,
        num_dims=num_dims,
        glmpca_iters=glmpca_iters,
        pearson_clip=pearson_clip,
        random_state=random_state,
    )
    gene_counts = gene_count_table(counts, gene_labels, GENE_GRADIENT_GENES)

    label_source = None
    cell_types: pd.DataFrame | None = None
    try:
        label_path = cell_type_csv if cell_type_csv is not None else download_cell_type_labels(cell_type_urls, raw_dir)
        cell_types = validate_cell_type_labels(label_path, all_cells, crop_cells)
        label_source = {
            "path": str(label_path),
            "sha256": sha256_file(label_path),
            "url": cell_type_urls[0] if cell_type_csv is None and cell_type_urls else None,
        }
    except BenchmarkInputError as exc:
        missing_path = out_dir / "missing_public_cell_type_labels.json"
        missing_path.write_text(
            json.dumps(
                {
                    "status": "missing",
                    "reason": str(exc),
                    "required_for": [
                        "cell_type_proportion_metrics.csv",
                        "GASTON Fig. 5d/e cell-type-specific panels",
                    ],
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        if require_cell_types:
            raise

    arrays = {
        "coords_mat_bounded": data_dir / "coords_mat_bounded.npy",
        "glmpca_bounded": data_dir / "glmpca_bounded.npy",
        "counts_mat_bounded": data_dir / "counts_mat_bounded.npy",
        "gene_labels": data_dir / "gene_labels.npy",
        "cell_labels_bounded": data_dir / "cell_labels_bounded.npy",
    }
    np.save(arrays["coords_mat_bounded"], crop_cells[["x_centroid", "y_centroid"]].to_numpy(dtype=np.float32))
    np.save(arrays["glmpca_bounded"], features.astype(np.float32))
    np.save(arrays["counts_mat_bounded"], np.asarray(counts.todense(), dtype=np.int32))
    np.save(arrays["gene_labels"], gene_labels.astype(str))
    if cell_types is not None:
        np.save(arrays["cell_labels_bounded"], cell_types["cell_type"].astype(str).to_numpy())

    cells_out = crop_cells.copy()
    cells_out = cells_out.reset_index(drop=True)
    cells_out["bounded_barcode"] = bounded_barcodes.astype(str)
    cells_out = pd.concat([cells_out, gene_counts], axis=1)
    if cell_types is not None:
        cells_out["cell_type"] = cell_types["cell_type"].to_numpy()
    cells_out.to_parquet(data_dir / "crop_cells.parquet", index=False)
    seed_cells.to_parquet(data_dir / "histoseg_seed_cells.parquet", index=False)
    graphclust.to_csv(data_dir / "graphclust_clusters.csv", index=False)

    write_inputs_manifest(
        out_dir,
        {
            "created_at": pd.Timestamp.utcnow().isoformat(),
            "xenium_outs": str(xenium_outs),
            "crop_bounds": asdict(bounds),
            "crop_cell_count": int(len(crop_cells)),
            "feature_method": feature_method,
            "feature_dims": int(num_dims),
            "sources": {
                "nature_fig5": NATURE_FIG5_URL,
                "gaston_tutorial": TUTORIAL_SOURCE_URL,
                "gaston_repo": GASTON_REPO_URL,
                "cell_type_labels": label_source,
            },
            "input_hashes": {
                "cells.parquet": sha256_file(xenium_outs / "cells.parquet"),
                "cell_feature_matrix.h5": sha256_file(xenium_outs / "cell_feature_matrix.h5"),
                "graphclust_clusters.csv": sha256_file(
                    xenium_outs
                    / "analysis"
                    / "clustering"
                    / "gene_expression_graphclust"
                    / "clusters.csv"
                ),
            },
            "runtime": runtime_fingerprint(),
            "elapsed_seconds": round(time.perf_counter() - t0, 3),
        },
    )
    return {**arrays, "crop_cells": data_dir / "crop_cells.parquet"}


def runtime_fingerprint() -> dict[str, object]:
    payload: dict[str, object] = {
        "python": sys.version,
        "platform": platform.platform(),
        "packages": {},
    }
    for package in ["numpy", "pandas", "scipy", "sklearn", "torch", "gaston", "glmpca", "histoseg"]:
        try:
            module = __import__(package)
            payload["packages"][package] = getattr(module, "__version__", "unknown")
        except Exception:
            payload["packages"][package] = None
    try:
        proc = subprocess.run(
            ["nvidia-smi", "--query-gpu=name,memory.total,memory.free", "--format=csv,noheader"],
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )
        payload["nvidia_smi"] = proc.stdout.strip() if proc.returncode == 0 else None
    except Exception:
        payload["nvidia_smi"] = None
    return payload


def run_gaston_training(
    *,
    out_dir: Path,
    seeds: Sequence[int],
    epochs: int,
    checkpoint: int,
    hidden_layers: Sequence[int],
    optimizer: str,
    device: str,
) -> Path:
    try:
        import torch
        from gaston import neural_net
    except ImportError as exc:
        raise BenchmarkInputError("gaston-spatial and torch are required for GASTON training") from exc

    data_dir = out_dir / "xenium_tumor_data"
    coords = np.load(data_dir / "coords_mat_bounded.npy")
    features = np.load(data_dir / "glmpca_bounded.npy")
    s_torch, a_torch = neural_net.load_rescale_input_data(coords, features)
    model_dir = ensure_dir(out_dir / "gaston_models")
    runtime_rows: list[dict[str, object]] = []

    for seed in seeds:
        seed_dir = ensure_dir(model_dir / f"seed{int(seed)}")
        if (seed_dir / "final_model.pt").exists():
            runtime_rows.append({"method": "GASTON", "stage": "train_seed", "seed": seed, "status": "exists"})
            continue
        start = time.perf_counter()
        neural_net.train(
            s_torch,
            a_torch,
            S_hidden_list=list(hidden_layers),
            A_hidden_list=list(hidden_layers),
            epochs=int(epochs),
            checkpoint=int(checkpoint),
            device=device,
            save_dir=str(seed_dir),
            optim=optimizer,
            seed=int(seed),
            save_final=True,
        )
        runtime_rows.append(
            {
                "method": "GASTON",
                "stage": "train_seed",
                "seed": int(seed),
                "status": "ok",
                "elapsed_seconds": round(time.perf_counter() - start, 3),
                "device": device,
                "cuda_available": bool(torch.cuda.is_available()),
            }
        )
        pd.DataFrame(runtime_rows).to_csv(out_dir / "runtime_memory_metrics.csv", index=False)
    return model_dir


def _torch_load(path: Path):
    import torch

    try:
        return torch.load(path, map_location="cpu", weights_only=False)
    except TypeError:
        return torch.load(path, map_location="cpu")


def select_best_gaston_seed(model_dir: Path, seeds: Sequence[int]) -> tuple[int, Path, float]:
    from gaston.neural_net import get_loss

    best_seed = None
    best_model_path = None
    best_loss = np.inf
    for seed in seeds:
        seed_dir = model_dir / f"seed{int(seed)}"
        if not seed_dir.exists():
            continue
        model_path = seed_dir / "final_model.pt"
        if not model_path.exists():
            epoch_models = sorted(seed_dir.glob("model_epoch_*.pt"))
            if not epoch_models:
                continue
            model_path = epoch_models[-1]
        model = _torch_load(model_path).cpu()
        s_torch = _torch_load(seed_dir / "Storch.pt").cpu()
        a_torch = _torch_load(seed_dir / "Atorch.pt").cpu()
        loss = float(get_loss(model, s_torch, a_torch).detach().cpu().item())
        if loss < best_loss:
            best_seed = int(seed)
            best_model_path = model_path
            best_loss = loss
    if best_seed is None or best_model_path is None:
        raise BenchmarkInputError(f"No complete GASTON seed folders found in {model_dir}")
    return best_seed, best_model_path, best_loss


def process_gaston_output(*, out_dir: Path, seeds: Sequence[int], num_domains: int) -> Path:
    try:
        from gaston import dp_related, isodepth_scaling
    except ImportError as exc:
        raise BenchmarkInputError("gaston-spatial is required to process GASTON output") from exc

    model_dir = out_dir / "gaston_models"
    best_seed, best_model_path, best_loss = select_best_gaston_seed(model_dir, seeds)
    model = _torch_load(best_model_path).cpu()
    seed_dir = best_model_path.parent
    features_scaled = _torch_load(seed_dir / "Atorch.pt").cpu().detach().numpy()
    coords_scaled = _torch_load(seed_dir / "Storch.pt").cpu().detach().numpy()
    coords_um = np.load(out_dir / "xenium_tumor_data" / "coords_mat_bounded.npy")
    cells = pd.read_parquet(out_dir / "xenium_tumor_data" / "crop_cells.parquet")

    gaston_isodepth, gaston_labels = dp_related.get_isodepth_labels(
        model,
        features_scaled,
        coords_scaled,
        int(num_domains),
    )
    gaston_isodepth = np.max(gaston_isodepth) - gaston_isodepth
    gaston_labels = (int(num_domains) - 1) - gaston_labels.astype(int)
    gaston_isodepth_um, scaling_factors = isodepth_scaling.adjust_isodepth(
        gaston_isodepth,
        gaston_labels,
        coords_um,
        num_domains=int(num_domains),
        q_vals=GASTON_Q_VALS,
        visualize=False,
        return_scaling_factors=True,
    )

    out = cells.copy()
    out["gaston_label"] = gaston_labels.astype(int)
    out["gaston_domain_name"] = [GASTON_DOMAIN_NAMES[int(label)] for label in gaston_labels]
    out["gaston_isodepth_um"] = gaston_isodepth_um.astype(float)
    out_path = out_dir / "gaston_reference_cells.parquet"
    out.to_parquet(out_path, index=False)
    (out_dir / "gaston_model_selection.json").write_text(
        json.dumps(
            {
                "best_seed": best_seed,
                "best_model_path": str(best_model_path),
                "best_loss": best_loss,
                "num_domains": int(num_domains),
                "q_vals": GASTON_Q_VALS,
                "scaling_factors": [float(item) for item in scaling_factors],
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    return out_path


def run_histoseg_domains(out_dir: Path) -> Path:
    from histoseg.contour.multi_structure import MultiStructureContourConfig, run_multi_structure_contours

    data_dir = out_dir / "xenium_tumor_data"
    h_dir = ensure_dir(out_dir / "histoseg")
    crop_cells = pd.read_parquet(data_dir / "crop_cells.parquet")
    crop_cells_path = h_dir / "crop_cells.parquet"
    crop_cells.to_parquet(crop_cells_path, index=False)
    clusters_path = data_dir / "graphclust_clusters.csv"
    start = time.perf_counter()
    result = run_multi_structure_contours(
        MultiStructureContourConfig(
            clusters_csv=clusters_path,
            cells_parquet=crop_cells_path,
            out_dir=h_dir,
            structures=HISTOSEG_STRUCTURE_SPECS,
            bins_x=900,
            bins_y=700,
            gaussian_sigma=2.25,
            min_cells=500,
            xenium_pixel_size_um=0.2125,
        )
    )
    assigned = pd.read_parquet(result.partition_table)
    assigned = assigned.rename(
        columns={
            "isoline_structure_id": "histoseg_structure_id",
            "isoline_structure_name": "histoseg_structure_name",
        }
    )
    if "histoseg_structure_name" not in assigned.columns:
        assigned["histoseg_structure_name"] = assigned["histoseg_structure_id"].map(
            {value: key for key, value in DOMAIN_NAME_TO_ID.items()}
        )
    unassigned = assigned["histoseg_structure_id"].isna() | (assigned["histoseg_structure_id"].astype(int) <= 0)
    if bool(unassigned.any()):
        raise BenchmarkInputError(
            f"HistoSeg left {int(unassigned.sum())} ADH crop cells unassigned; "
            "full benchmark acceptance requires every crop cell to map to one structure."
        )
    assigned = attach_signed_distances(assigned, Path(result.geojson), xenium_pixel_size_um=0.2125)
    out_path = out_dir / "histoseg_domains.parquet"
    assigned.to_parquet(out_path, index=False)
    append_runtime_row(
        out_dir,
        {
            "method": "HistoSeg",
            "stage": "multi_structure_contours",
            "status": "ok",
            "elapsed_seconds": round(time.perf_counter() - start, 3),
        },
    )
    return out_path


def append_runtime_row(out_dir: Path, row: Mapping[str, object]) -> None:
    path = out_dir / "runtime_memory_metrics.csv"
    existing = pd.read_csv(path) if path.exists() else pd.DataFrame()
    pd.concat([existing, pd.DataFrame([row])], ignore_index=True).to_csv(path, index=False)


def _load_domain_polygons(geojson_path: Path, xenium_pixel_size_um: float) -> dict[str, object]:
    payload = json.loads(geojson_path.read_text(encoding="utf-8"))
    grouped: dict[str, list[Polygon]] = {}
    for feature in payload.get("features", []):
        props = feature.get("properties", {})
        name = str(props.get("assigned_structure") or props.get("name"))
        geom = shape(feature["geometry"])
        scaled = Polygon(
            [(float(x) * xenium_pixel_size_um, float(y) * xenium_pixel_size_um) for x, y in geom.exterior.coords]
        )
        if scaled.is_valid and not scaled.is_empty:
            grouped.setdefault(name, []).append(scaled)
    return {name: unary_union(polys) for name, polys in grouped.items() if polys}


def attach_signed_distances(
    assigned: pd.DataFrame,
    geojson_path: Path,
    *,
    xenium_pixel_size_um: float,
) -> pd.DataFrame:
    polygons = _load_domain_polygons(geojson_path, xenium_pixel_size_um)
    out = assigned.copy()
    signed = []
    for row in out.itertuples(index=False):
        name = getattr(row, "histoseg_structure_name", None)
        if not name or str(name) not in polygons:
            signed.append(np.nan)
            continue
        point = Point(float(getattr(row, "x_centroid")), float(getattr(row, "y_centroid")))
        polygon = polygons[str(name)]
        distance = float(polygon.boundary.distance(point))
        signed.append(distance if polygon.covers(point) else -distance)
    out["histoseg_signed_distance_um"] = signed
    return out


def merge_reference_tables(out_dir: Path) -> pd.DataFrame:
    gaston = pd.read_parquet(out_dir / "gaston_reference_cells.parquet")
    histoseg = pd.read_parquet(out_dir / "histoseg_domains.parquet")
    keep_cols = [
        "Barcode",
        "histoseg_structure_id",
        "histoseg_structure_name",
        "histoseg_signed_distance_um",
    ]
    histoseg = histoseg[[col for col in keep_cols if col in histoseg.columns]].copy()
    merged = gaston.merge(histoseg, on="Barcode", how="inner", validate="one_to_one")
    if len(merged) != len(gaston):
        raise BenchmarkInputError(
            f"GASTON/HistoSeg merge lost cells: gaston={len(gaston)} merged={len(merged)}"
        )
    return merged


def compute_domain_overlap_metrics(cells: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for domain in GASTON_DOMAIN_NAMES.values():
        g = cells["gaston_domain_name"] == domain
        h = cells["histoseg_structure_name"] == domain
        intersection = int((g & h).sum())
        union = int((g | h).sum())
        rows.append(
            {
                "domain": domain,
                "gaston_cells": int(g.sum()),
                "histoseg_cells": int(h.sum()),
                "intersection": intersection,
                "union": union,
                "jaccard": intersection / union if union else np.nan,
                "dice": (2 * intersection) / (int(g.sum()) + int(h.sum())) if int(g.sum()) + int(h.sum()) else np.nan,
                "recall_vs_gaston": intersection / int(g.sum()) if int(g.sum()) else np.nan,
                "precision_vs_gaston": intersection / int(h.sum()) if int(h.sum()) else np.nan,
            }
        )
    rows.append(
        {
            "domain": "overall",
            "gaston_cells": int(len(cells)),
            "histoseg_cells": int(len(cells)),
            "intersection": int((cells["gaston_domain_name"] == cells["histoseg_structure_name"]).sum()),
            "union": int(len(cells)),
            "jaccard": np.nan,
            "dice": np.nan,
            "recall_vs_gaston": np.nan,
            "precision_vs_gaston": np.nan,
            "accuracy": float((cells["gaston_domain_name"] == cells["histoseg_structure_name"]).mean()),
        }
    )
    return pd.DataFrame(rows)


def _safe_corr(x: Iterable[float], y: Iterable[float], method: str) -> float:
    x_arr = np.asarray(list(x), dtype=float)
    y_arr = np.asarray(list(y), dtype=float)
    mask = np.isfinite(x_arr) & np.isfinite(y_arr)
    if mask.sum() < 3 or np.unique(x_arr[mask]).size < 2 or np.unique(y_arr[mask]).size < 2:
        return np.nan
    if method == "spearman":
        return float(spearmanr(x_arr[mask], y_arr[mask]).statistic)
    return float(pearsonr(x_arr[mask], y_arr[mask]).statistic)


def compute_isodepth_sdf_metrics(cells: pd.DataFrame) -> pd.DataFrame:
    rows = [
        {
            "domain": "global",
            "n_cells": int(len(cells)),
            "pearson_isodepth_vs_sdf": _safe_corr(
                cells["gaston_isodepth_um"], cells["histoseg_signed_distance_um"], "pearson"
            ),
            "spearman_isodepth_vs_sdf": _safe_corr(
                cells["gaston_isodepth_um"], cells["histoseg_signed_distance_um"], "spearman"
            ),
            "median_abs_signed_distance_um": float(np.nanmedian(np.abs(cells["histoseg_signed_distance_um"]))),
        }
    ]
    for domain in GASTON_DOMAIN_NAMES.values():
        subset = cells.loc[cells["gaston_domain_name"] == domain]
        rows.append(
            {
                "domain": domain,
                "n_cells": int(len(subset)),
                "pearson_isodepth_vs_sdf": _safe_corr(
                    subset["gaston_isodepth_um"], subset["histoseg_signed_distance_um"], "pearson"
                ),
                "spearman_isodepth_vs_sdf": _safe_corr(
                    subset["gaston_isodepth_um"], subset["histoseg_signed_distance_um"], "spearman"
                ),
                "median_abs_signed_distance_um": float(np.nanmedian(np.abs(subset["histoseg_signed_distance_um"]))),
            }
        )
    return pd.DataFrame(rows)


def _domain_depth_bins(cells: pd.DataFrame, depth_col: str, domain_col: str) -> pd.Series:
    bin_index = pd.Series(index=cells.index, dtype="Int64")
    for label, domain in GASTON_DOMAIN_NAMES.items():
        n_bins = [7, 10, 7, 7][int(label)]
        subset = cells.loc[cells[domain_col] == domain, depth_col]
        if subset.empty:
            continue
        ranks = subset.rank(method="first")
        bins = pd.qcut(ranks, q=min(n_bins, len(subset)), labels=False, duplicates="drop")
        bin_index.loc[subset.index] = bins.astype("Int64")
    return bin_index


def compute_cell_type_proportion_metrics(cells: pd.DataFrame) -> pd.DataFrame:
    if "cell_type" not in cells.columns:
        raise BenchmarkInputError("cell_type column missing; cannot compute Fig. 5 cell-type proportions")
    rows: list[dict[str, object]] = []
    for method, depth_col, domain_col in [
        ("GASTON", "gaston_isodepth_um", "gaston_domain_name"),
        ("HistoSeg", "histoseg_signed_distance_um", "histoseg_structure_name"),
    ]:
        bins = _domain_depth_bins(cells, depth_col, domain_col)
        tmp = cells.copy()
        tmp["depth_bin"] = bins
        for (domain, depth_bin), subset in tmp.dropna(subset=["depth_bin"]).groupby([domain_col, "depth_bin"]):
            total = len(subset)
            proportions = subset["cell_type"].value_counts(normalize=True)
            for cell_type in sorted(set(CELL_TYPE_PANEL_ORDER).intersection(proportions.index)):
                rows.append(
                    {
                        "method": method,
                        "domain": domain,
                        "depth_bin": int(depth_bin),
                        "cell_type": cell_type,
                        "proportion": float(proportions[cell_type]),
                        "n_cells": int(total),
                    }
                )
    return pd.DataFrame(rows)


def compute_gene_gradient_metrics(cells: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for method, depth_col, domain_col in [
        ("GASTON", "gaston_isodepth_um", "gaston_domain_name"),
        ("HistoSeg", "histoseg_signed_distance_um", "histoseg_structure_name"),
    ]:
        for domain in GASTON_DOMAIN_NAMES.values():
            subset = cells.loc[cells[domain_col] == domain]
            for gene in GENE_GRADIENT_GENES:
                expr = np.log1p(pd.to_numeric(subset[gene], errors="coerce").to_numpy(dtype=float))
                depth = pd.to_numeric(subset[depth_col], errors="coerce").to_numpy(dtype=float)
                mask = np.isfinite(expr) & np.isfinite(depth)
                slope = np.nan
                if mask.sum() >= 3 and np.unique(depth[mask]).size > 1:
                    slope = float(np.polyfit(depth[mask], expr[mask], 1)[0])
                rows.append(
                    {
                        "method": method,
                        "domain": domain,
                        "gene": gene,
                        "n_cells": int(mask.sum()),
                        "mean_log1p_expression": float(np.nanmean(expr)) if len(expr) else np.nan,
                        "slope_log1p_per_um": slope,
                        "spearman_depth_expression": _safe_corr(depth, expr, "spearman"),
                    }
                )
    return pd.DataFrame(rows)


def write_metrics(out_dir: Path) -> dict[str, Path]:
    cells = merge_reference_tables(out_dir)
    cells.to_parquet(out_dir / "benchmark_merged_cells.parquet", index=False)
    outputs = {
        "domain_overlap": out_dir / "domain_overlap_metrics.csv",
        "isodepth_sdf": out_dir / "isodepth_sdf_metrics.csv",
        "gene_gradient": out_dir / "gene_gradient_metrics.csv",
    }
    compute_domain_overlap_metrics(cells).to_csv(outputs["domain_overlap"], index=False)
    compute_isodepth_sdf_metrics(cells).to_csv(outputs["isodepth_sdf"], index=False)
    compute_gene_gradient_metrics(cells).to_csv(outputs["gene_gradient"], index=False)
    if "cell_type" in cells.columns:
        outputs["cell_type_proportion"] = out_dir / "cell_type_proportion_metrics.csv"
        compute_cell_type_proportion_metrics(cells).to_csv(outputs["cell_type_proportion"], index=False)
    else:
        missing = out_dir / "cell_type_proportion_metrics.missing.json"
        missing.write_text(
            json.dumps(
                {
                    "status": "missing",
                    "reason": "Exact public cell-type labels were not available/aligned.",
                },
                indent=2,
            ),
            encoding="utf-8",
        )
    return outputs


def render_figures(out_dir: Path) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    cells = merge_reference_tables(out_dir)
    fig_dir = ensure_dir(out_dir / "figures")
    outputs: list[Path] = []

    panels = [
        ("gaston_domains", "gaston_domain_name", None),
        ("histoseg_domains", "histoseg_structure_name", None),
        ("gaston_isodepth", "gaston_isodepth_um", "viridis"),
        ("histoseg_sdf", "histoseg_signed_distance_um", "coolwarm"),
    ]
    for stem, column, cmap in panels:
        fig, ax = plt.subplots(figsize=(6, 6))
        if cmap is None:
            colors = cells[column].map(GASTON_DOMAIN_COLORS).fillna("#999999")
            ax.scatter(-cells["x_centroid"], cells["y_centroid"], c=colors, s=2, linewidths=0)
        else:
            im = ax.scatter(-cells["x_centroid"], cells["y_centroid"], c=cells[column], cmap=cmap, s=2, linewidths=0)
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        ax.set_aspect("equal")
        ax.set_axis_off()
        fig.tight_layout()
        for suffix in ["png", "svg"]:
            path = fig_dir / f"{stem}.{suffix}"
            fig.savefig(path, dpi=240)
            outputs.append(path)
        plt.close(fig)

    overlap = pd.read_csv(out_dir / "domain_overlap_metrics.csv")
    overlap = overlap.loc[overlap["domain"] != "overall"]
    fig, ax = plt.subplots(figsize=(7, 4))
    sns.barplot(data=overlap, x="domain", y="jaccard", ax=ax, color="#4c78a8")
    ax.set_ylim(0, 1)
    ax.set_xlabel("")
    ax.set_ylabel("Cell-level Jaccard")
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    for suffix in ["png", "svg"]:
        path = fig_dir / f"domain_jaccard_summary.{suffix}"
        fig.savefig(path, dpi=240)
        outputs.append(path)
    plt.close(fig)

    gradients = pd.read_csv(out_dir / "gene_gradient_metrics.csv")
    fig, ax = plt.subplots(figsize=(8, 4))
    sns.barplot(data=gradients, x="domain", y="spearman_depth_expression", hue="gene", ax=ax)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xlabel("")
    ax.set_ylabel("Spearman(depth, log1p expression)")
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    for suffix in ["png", "svg"]:
        path = fig_dir / f"gene_gradient_summary.{suffix}"
        fig.savefig(path, dpi=240)
        outputs.append(path)
    plt.close(fig)
    return outputs


def run_smoke(out_dir: Path) -> None:
    ensure_dir(out_dir)
    rng = np.random.default_rng(0)
    n = 120
    x = rng.uniform(1500, 3250, n)
    y = rng.uniform(2000, 4000, n)
    domains = np.asarray(list(GASTON_DOMAIN_NAMES.values()))
    labels = rng.integers(0, len(domains), n)
    cells = pd.DataFrame(
        {
            "Barcode": [str(i + 1) for i in range(n)],
            "cell_id": [str(i + 1) for i in range(n)],
            "x_centroid": x,
            "y_centroid": y,
            "gaston_label": labels,
            "gaston_domain_name": domains[labels],
            "gaston_isodepth_um": x - x.min(),
            "histoseg_structure_id": labels + 1,
            "histoseg_structure_name": domains[labels],
            "histoseg_signed_distance_um": y - np.median(y),
            "cell_type": rng.choice(CELL_TYPE_PANEL_ORDER, n),
            "SELL": rng.poisson(1, n),
            "TCL1A": rng.poisson(1, n),
        }
    )
    gaston_cols = [
        "Barcode",
        "cell_id",
        "x_centroid",
        "y_centroid",
        "gaston_label",
        "gaston_domain_name",
        "gaston_isodepth_um",
        "cell_type",
        "SELL",
        "TCL1A",
    ]
    cells[gaston_cols].to_parquet(out_dir / "gaston_reference_cells.parquet", index=False)
    cells[["Barcode", "histoseg_structure_id", "histoseg_structure_name", "histoseg_signed_distance_um"]].to_parquet(
        out_dir / "histoseg_domains.parquet", index=False
    )
    write_inputs_manifest(
        out_dir,
        {
            "mode": "smoke",
            "crop_bounds": asdict(CropBounds()),
            "runtime": runtime_fingerprint(),
        },
    )
    write_metrics(out_dir)
    render_figures(out_dir)


def run_pipeline(args: argparse.Namespace) -> None:
    out_dir = ensure_dir(Path(args.out_dir).expanduser().resolve())
    bounds = CropBounds(*[float(v) for v in args.crop_bounds])
    seeds = parse_seed_spec(args.gaston_seeds)
    mode = args.mode
    if mode == "smoke":
        run_smoke(out_dir)
        return
    if mode in {"prepare", "all", "histoseg", "gaston"}:
        prepare_inputs(
            xenium_outs=Path(args.xenium_outs).expanduser().resolve(),
            out_dir=out_dir,
            bounds=bounds,
            feature_method=args.feature_method,
            num_dims=args.feature_dims,
            glmpca_iters=args.glmpca_iters,
            pearson_clip=args.pearson_clip,
            cell_type_urls=args.cell_type_url,
            cell_type_csv=Path(args.cell_type_csv).expanduser().resolve() if args.cell_type_csv else None,
            require_cell_types=not args.allow_missing_cell_types,
        )
    if mode in {"gaston", "all"}:
        run_gaston_training(
            out_dir=out_dir,
            seeds=seeds,
            epochs=args.gaston_epochs,
            checkpoint=args.gaston_checkpoint,
            hidden_layers=args.gaston_hidden_layers,
            optimizer=args.gaston_optimizer,
            device=args.gaston_device,
        )
        process_gaston_output(out_dir=out_dir, seeds=seeds, num_domains=args.num_domains)
    if mode in {"histoseg", "all"}:
        run_histoseg_domains(out_dir)
    if mode in {"metrics", "figures", "all"}:
        if mode != "figures":
            write_metrics(out_dir)
        render_figures(out_dir)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the GASTON Fig. 5 / HistoSeg Xenium breast benchmark."
    )
    parser.add_argument(
        "--mode",
        choices=["prepare", "histoseg", "gaston", "metrics", "figures", "all", "smoke"],
        default="all",
    )
    parser.add_argument("--xenium-outs", default=str(DEFAULT_XENIUM_OUTS))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument(
        "--crop-bounds",
        nargs=4,
        type=float,
        default=list(CROP_BOUNDS),
        metavar=("XMIN", "XMAX", "YMIN", "YMAX"),
    )
    parser.add_argument("--cell-type-url", action="append", default=[])
    parser.add_argument("--cell-type-csv", default=None)
    parser.add_argument(
        "--allow-missing-cell-types",
        action="store_true",
        help="Allow non-cell-type parts of the benchmark to run without exact public cell labels.",
    )
    parser.add_argument("--feature-method", choices=["glmpca", "pearson", "log1p-pca"], default="glmpca")
    parser.add_argument("--feature-dims", type=int, default=20)
    parser.add_argument("--glmpca-iters", type=int, default=30)
    parser.add_argument("--pearson-clip", type=float, default=0.01)
    parser.add_argument("--gaston-seeds", default="0-29")
    parser.add_argument("--gaston-epochs", type=int, default=100000)
    parser.add_argument("--gaston-checkpoint", type=int, default=500)
    parser.add_argument("--gaston-hidden-layers", nargs="+", type=int, default=[20, 20])
    parser.add_argument("--gaston-optimizer", default="adam")
    parser.add_argument("--gaston-device", default="cuda")
    parser.add_argument("--num-domains", type=int, default=4)
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    run_pipeline(args)


if __name__ == "__main__":
    main()
