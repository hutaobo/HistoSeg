# Contour generation based on cell clusters

这里写一句话：这个脚本干啥、输入输出是什么。

## Overview

::::{grid} 1 2 2 2
:::{grid-item-card} Inputs
- `clusters.csv`
- `cells.parquet`
- `tissue_boundary.csv`
:::
:::{grid-item-card} Outputs
- `pattern1_isoline0p5_*.npy`
- `pattern1_isoline0p5.png`
- `params.json`
:::
::::

## Step-by-step

1. Align IDs
2. Select Pattern1 clusters
3. Train KNN
4. Predict on grid + smoothing
5. Extract isoline (0.5)
6. Filter loops

## Full script

```{literalinclude} ../../workflows/contour_generation_pattern1.py
---
language: python
linenos: true
---
