```markdown
# HistoSeg

[![PyPI](https://img.shields.io/pypi/v/histoseg.svg)](https://pypi.org/project/histoseg/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/histoseg.svg)](https://pypi.org/project/histoseg/)
[![Docs](https://readthedocs.org/projects/histoseg/badge/?version=latest)](https://histoseg.readthedocs.io/en/latest/)
[![Publish to PyPI](https://github.com/hutaobo/HistoSeg/actions/workflows/publish.yml/badge.svg)](https://github.com/hutaobo/HistoSeg/actions/workflows/publish.yml)
[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)

HistoSeg is a Python toolkit for **spatial transcriptomics segmentation / geometry extraction**.

The current focus is **Pattern1 isoline (0.5)** contour generation from cell clusters (e.g. 10x Xenium GraphClust output):
- pick a set of “target clusters” (Pattern1)
- fit a KNN regressor to estimate *P(target)* over space
- smooth the probability field
- extract a contour (isoline) at **level = 0.5**
- save contour vertices and a quick preview plot

## Quick links

- Documentation: https://histoseg.readthedocs.io
- Source code: https://github.com/hutaobo/HistoSeg
- Issue tracker: https://github.com/hutaobo/HistoSeg/issues

> ⚠️ **License note**
>
> This repository is licensed under **CC BY-NC 4.0** (non-commercial).
> Commercial use is not permitted without prior permission. See `LICENSE`.

## Installation

### Install from PyPI (recommended)

```bash
pip install -U histoseg
```

### Install from source (for development)

```bash
git clone https://github.com/hutaobo/HistoSeg.git
cd HistoSeg
pip install -U pip
pip install -e .
```

### Dependencies

The Pattern1 isoline workflow uses:
- numpy, pandas
- scipy
- scikit-learn
- matplotlib
- a Parquet engine (pyarrow is recommended)

If you run into missing imports, install them explicitly:

```bash
pip install -U numpy pandas pyarrow scipy scikit-learn matplotlib
```

Optional:
- Hugging Face downloader: `pip install -U huggingface_hub`

## Tutorial: Pattern1 isoline (0.5)

### What you need (inputs)

The isoline workflow expects the following files:

1. `clusters.csv`  
   - Typically from GraphClust: `analysis/clustering/gene_expression_graphclust/clusters.csv`
   - Must contain columns: `Barcode`, `Cluster`

2. `cells.parquet`  
   - A cell-level table with spatial coordinates (x/y-like columns)
   - Must contain at least:
     - coordinate columns (e.g. `x`/`y` or `x_centroid`/`y_centroid`)
     - an id column that can be aligned with `clusters.csv:Barcode` (the code tries several common column names)

3. `tissue_boundary.csv` (optional but recommended if you enable synthetic background)
   - Must contain columns `x,y` **or** `X,Y`

### What you get (outputs)

By default, the pipeline writes into `out_dir`:
- `params.json` — all parameters + inferred join columns
- `pattern1_isoline_<level>_<i>.npy` — contour vertices (Nx2 arrays)
- `pattern1_isoline_<level>.png` — quick preview plot

### One-liner (from a Hugging Face dataset repo)

This follows the example notebook in `examples/contour_generation_pattern1_from_hf.ipynb`.

```python
# pip install -U histoseg
# pip install -U huggingface_hub pandas pyarrow numpy scipy scikit-learn matplotlib

from histoseg import run_pattern1_isoline_from_hf

PATTERN1 = (10, 23, 19, 27, 14, 20, 25, 26)

result = run_pattern1_isoline_from_hf(
    repo_id="hutaobo/output-XETG00082_C105",
    revision="main",  # or a commit hash for strict reproducibility
    out_dir="outputs/pattern1_isoline0p5_from_graphclust",
    pattern1_clusters=PATTERN1,

    # Defaults are intentionally exposed for tuning:
    grid_n=1200,
    knn_k=30,
    smooth_sigma=5.0,
    min_cells_inside=10,
)

print("Outputs folder:", result.out_dir)
print("Preview image:", result.preview_png)
print("Contours:", len(result.contours))
```

### Run on local files

```python
from histoseg import Pattern1IsolineConfig, run_pattern1_isoline

PATTERN1 = (10, 23, 19, 27, 14, 20, 25, 26)

cfg = Pattern1IsolineConfig(
    clusters_csv="/path/to/analysis/clustering/gene_expression_graphclust/clusters.csv",
    cells_parquet="/path/to/cells.parquet",
    tissue_boundary_csv="/path/to/tissue_boundary.csv",
    out_dir="outputs/pattern1_isoline0p5",
    pattern1_clusters=PATTERN1,

    # Optional tuning:
    grid_n=1200,
    knn_k=30,
    smooth_sigma=5.0,
    min_cells_inside=10,
)

result = run_pattern1_isoline(cfg)
print(result)
```

### How it works (workflow overview)

```mermaid
flowchart TD
  A[clusters.csv\nBarcode/Cluster] --> C[Align barcodes\nwith cells.parquet]
  B[cells.parquet\nx/y + id-like column] --> C
  C --> D[Select target clusters\n(Pattern1)]
  D --> E[Sample background points\n(other cells)]
  F[tissue_boundary.csv] --> G[Generate synthetic background\n(optional)]
  G --> E
  D --> H[KNN regression\npredict P(target)]
  E --> H
  H --> I[Predict on mesh grid]
  I --> J[Gaussian smoothing]
  J --> K[Mask by tissue\n(nearest-cell threshold)]
  K --> L[Extract isoline\nlevel = 0.5]
  L --> M[Filter loops\nmin_cells_inside]
  M --> N[Save params.json\n+ contours .npy\n+ preview .png]
```

### Troubleshooting & tuning

If no contour is found, try:
- decrease `min_cells_inside` (e.g. 10 → 3)
- increase `smooth_sigma` (e.g. 5 → 8)
- increase `knn_k` (e.g. 30 → 50)
- reduce `grid_n` to speed up (note: `grid_n=1200` can be heavy)

## API reference (high-level)

### Pattern1 isoline

- `Pattern1IsolineConfig`  
  Dataclass holding all parameters and input paths.

- `run_pattern1_isoline(cfg) -> Pattern1IsolineResult`  
  Runs the full pipeline on local files.

- `run_pattern1_isoline_from_hf(repo_id, revision="main", ...) -> Pattern1IsolineResult`  
  Convenience wrapper that downloads required files from a Hugging Face *dataset repo* and then runs the pipeline.

### Hugging Face I/O helpers

- `download_xenium_outs(repo_id, revision="main", clusters_relpath=..., cache_dir=None)`  
  Downloads `cells.parquet`, `tissue_boundary.csv`, and the specified `clusters.csv` from a dataset repo.

### SFPlot utilities (legacy / optional)

This repository contains a small subset of SFPlot-style utilities and re-exports:
- `compute_cophenetic_distances_from_df(df, ...)`
- `plot_cophenetic_heatmap(matrix, ...)`

## GUI (experimental)

A GUI entry point is configured as:

```bash
histoseg-gui
```

Note:
- The current GUI code path is still in flux and may require extra dependencies (e.g. Pillow) and/or an external `sfplot` installation.
- For production workflows, prefer the Python API shown above.

## Contributing

Issues and pull requests are welcome.

When reporting a bug, please include:
- OS + Python version
- `histoseg` version
- minimal reproducible code (or a small input subset)
- expected vs. actual behavior

## License

CC BY-NC 4.0 (non-commercial). See `LICENSE` for details.
```
