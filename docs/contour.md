# Pattern1 Isoline Contours

This page summarizes the Pattern1 contour workflow in the `histoseg.contour`
feature group. For the broader Contour Analysis overview, see
{doc}`contour_analysis`.

## When To Reach For This Workflow

Use Pattern1 isolines when you want a single target structure traced from a
selected set of clusters and exported as contour loops plus quick-review
artifacts.

## Inputs

- `clusters.csv`
- `cells.parquet`
- optional `tissue_boundary.csv`

## Outputs

- `params.json`
- contour vertex arrays
- preview PNG

## Core Steps

1. Align IDs between `clusters.csv` and `cells.parquet`.
2. Select Pattern1 clusters.
3. Sample background points.
4. Train a KNN regressor.
5. Predict on a grid and smooth the field.
6. Extract the 0.5 isoline.
7. Filter loops by minimum target-cell support.
8. Save contours, preview images, and provenance metadata.
