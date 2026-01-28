# Contour generation (Pattern1 isoline 0.5)

This page documents the pipeline used to generate a 0.5 isoline contour for a set of clusters (Pattern1).

## Overview

::::{grid} 1 2 2 2
:::{grid-item-card} Inputs
- `clusters.csv`
- `cells.parquet`
- `tissue_boundary.csv`
:::
:::{grid-item-card} Outputs
- `params.json`
- `pattern1_isoline0p5_*.npy`
- `pattern1_isoline0p5.png`
:::
::::

## Steps

1. **Align IDs** between `clusters.csv` (Barcode) and `cells.parquet` (id-like column)
2. **Define Pattern1 clusters** and extract target cells
3. **Sample background points** from other cells + synthetic background (optional)
4. **Train KNN regressor** (target=1, background=0)
5. **Predict on grid** + Gaussian smoothing
6. **Mask by tissue** using nearest-cell distance threshold
7. **Extract contour** at level = 0.5
8. **Filter loops** by minimum number of target cells inside
9. **Save** parameters + contour vertices + diagnostic plot

## Full script (source of truth)

```{literalinclude} ../workflows/contour_generation_pattern1.py
---
language: python
linenos: true
---
