# Pattern1 Contour Generation

This workflow generates a 0.5 isoline contour for selected clusters from
cell-coordinate data. It belongs to the `histoseg.contour` feature group.

## Overview

::::{grid} 1 2 2 2
:::{grid-item-card} Inputs
- `clusters.csv`
- `cells.parquet`
- optional `tissue_boundary.csv`
:::
:::{grid-item-card} Outputs
- contour vertex arrays
- preview PNG
- `params.json`
:::
::::

## Step-by-Step

1. Align cell IDs between `clusters.csv` and `cells.parquet`.
2. Select target Pattern1 clusters.
3. Sample background points from non-target cells and optional synthetic tissue background.
4. Train a KNN regressor on target/background labels.
5. Predict on a mesh grid and smooth the probability field.
6. Extract the 0.5 isoline.
7. Filter small or unsupported loops.
8. Write contours, preview images, and provenance metadata.

## Python API

```python
from histoseg.contour import Pattern1IsolineConfig, run_pattern1_isoline

result = run_pattern1_isoline(
    Pattern1IsolineConfig(
        clusters_csv="clusters.csv",
        cells_parquet="cells.parquet",
        tissue_boundary_csv="tissue_boundary.csv",
        out_dir="outputs/pattern1",
        pattern1_clusters=(10, 23, 19),
    )
)
```

## CLI

```bash
histoseg-contour pattern1 --clusters-csv clusters.csv --cells-parquet cells.parquet --out-dir outputs/pattern1 --pattern1-clusters 10,23,19
```
