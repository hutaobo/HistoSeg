# Contour Analysis

Contour Analysis is the spatial/cell-coordinate workflow group exposed as
`histoseg.contour`. Use it when the primary input is a table of cells with
spatial coordinates and cluster assignments.

## Workflows

- Pattern1 isoline contours from clustered cell coordinates.
- Multi-structure contour partitioning.
- Xenium Explorer annotation exports.
- Hugging Face dataset helper workflows for Xenium-style public datasets.

## Inputs

Pattern1 and multi-structure contour workflows use:

- `clusters.csv` with barcode and cluster columns.
- `cells.parquet` with cell identifiers and x/y coordinates.
- optional `tissue_boundary.csv` for tissue-aware background generation.
- selected cluster IDs for the target Pattern1 or named structures.

## Outputs

- contour vertices and preview PNGs.
- region partition tables.
- Xenium Explorer-compatible GeoJSON/CSV exports for multi-structure contours.
- `params.json`, summary CSVs, and metrics JSON files.

## Pattern1 Python API

```python
from histoseg.contour import Pattern1IsolineConfig, run_pattern1_isoline

cfg = Pattern1IsolineConfig(
    clusters_csv="/path/to/clusters.csv",
    cells_parquet="/path/to/cells.parquet",
    tissue_boundary_csv="/path/to/tissue_boundary.csv",
    out_dir="outputs/pattern1_isoline0p5",
    pattern1_clusters=(10, 23, 19, 27, 14, 20, 25, 26),
)
result = run_pattern1_isoline(cfg)
```

## Multi-Structure Python API

```python
from histoseg.contour import MultiStructureContourConfig, MultiStructureSpec
from histoseg.contour import run_multi_structure_contours

result = run_multi_structure_contours(
    MultiStructureContourConfig(
        clusters_csv="/path/to/clusters.csv",
        cells_parquet="/path/to/cells.parquet",
        out_dir="outputs/multi_structure",
        structures=[
            MultiStructureSpec("Region A", [1, 2, 3]),
            MultiStructureSpec("Region B", [4, 5, 6]),
        ],
    )
)
```

## CLI

```bash
histoseg-contour pattern1 --clusters-csv clusters.csv --cells-parquet cells.parquet --out-dir outputs/p1 --pattern1-clusters 10,23,19
histoseg-contour multi-structure --clusters-csv clusters.csv --cells-parquet cells.parquet --out-dir outputs/ms --structures-json structures.json
```

## Tutorials

- {doc}`tutorials/pattern1_isoline`
- {doc}`workflows/contour_generation`
