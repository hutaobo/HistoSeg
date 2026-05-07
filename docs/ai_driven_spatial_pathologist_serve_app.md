# AI Driven Spatial Pathologist Serve App

:::{div} hs-section-kicker
Web app
:::

:::{div} hs-page-intro
The [AI Driven Spatial Pathologist Serve app](https://ai-driven-spatial-pathologist.serve.scilifelab.se/)
is a browser-based entry point for dendrogram-guided `histoseg.contour`
multi-structure contour analysis. It helps turn Xenium cluster assignments into
larger, reviewable spatial structures before exporting HistoSeg-style contour
artifacts.
:::

:::{div} hs-metadata
<span>Input: Xenium `cells.parquet` and GraphClust `clusters.csv`</span>
<span>Selection: cophenetic dendrogram branches</span>
<span>Outputs: preview PNGs, contour arrays, partition tables, and Xenium Explorer exports</span>
:::

## What The App Does

The Serve app wraps the multi-structure contour workflow in an interactive
Gradio interface. It first aligns uploaded cell coordinates with cluster labels,
computes Search-and-Find / cophenetic relationships between clusters, cuts the
resulting dendrogram into candidate structures, and lets you choose one or more
branches. The selected cluster groups are then converted into non-overlapping
structure partitions and contour exports.

Use the app when you want an interactive structure-picking step before running
the final contour generation. Use the Python API or CLI from
{doc}`contour_analysis` when you need scripted or batch execution.

## Input Files

The app expects Xenium-derived tables from the same sample.

| File | Required | Expected contents |
| --- | --- | --- |
| `cells.parquet` | yes | Cell identifiers plus x/y coordinates. In a standard Xenium export this is usually `outs/cells.parquet`. |
| `clusters.csv` | yes | GraphClust-style cell-to-cluster assignments. A common Xenium path is `outs/analysis/clustering/gene_expression_graphclust/clusters.csv`. |
| `tissue_boundary.csv` | optional | Upload field exists in the app, but the current multi-structure mode does not use this file. |

`clusters.csv` must include:

- `Barcode`
- `Cluster`

`cells.parquet` must include:

- a joinable cell identifier column, such as `cell_id`, `barcode`, `Barcode`,
  `cell_barcode`, or `id`;
- x/y coordinate columns, such as `x_centroid` and `y_centroid`, or equivalent
  columns like `x`, `y`, `x_location`, `y_location`, `x_centroid_um`, or
  `y_centroid_um`.

Cluster labels may be numeric or text labels. Numeric labels such as `10`,
`10.0`, and `"10"` are normalized to the same cluster ID during matching.

## Workflow

1. Open the Serve app.
2. Upload `cells.parquet` and `clusters.csv`.
3. Choose how many candidate structures the dendrogram should be cut into.
4. Click **Build dendrogram and candidate structures**.
5. Select one or more candidate structures in the checklist, or type one
   structure per line in the cluster-ID box.
6. Click **Run multi-structure contour analysis**.
7. Review the contour preview, run summary, and downloadable output files.

Manual structure lines use comma-separated cluster IDs:

```text
10,23,19
27,14,20
25,26
```

Each line becomes one named structure in the final contour run. A cluster ID can
belong to only one structure in a single run.

## Output Files

A successful run writes a preview, partition table, metrics, contour arrays, and
Xenium Explorer-compatible annotation exports.

| Output | Description |
| --- | --- |
| `multi_structure_contour_preview.png` | Visual preview of assigned cells and generated structure contours. |
| `cophenetic_heatmap_row_coph.png` | Cophenetic heatmap used as the structure-selection reference. |
| `structure_contour_metrics.json` | Selected structures, parameters, isoline metrics, assignment metrics, and cophenetic stats. |
| `cells_with_structure_partition.parquet` or `cells_with_structure_partition.csv` | Cell table with `isoline_structure_id` and `isoline_structure_name` assignments. |
| `structure_<id>_contour_<n>.npy` | NumPy vertex arrays for each generated contour path. |
| `structure_contour_cell_counts.csv` | Per-structure selected cell counts, assigned cell counts, and contour counts. |
| `xenium_explorer_annotations.geojson` | Polygon annotations for Xenium Explorer import. |
| `xenium_explorer_annotations.csv` | Vertex table with `Selection`, `X`, and `Y` columns. |
| `xenium_explorer_annotations_summary.csv` | Per-polygon structure IDs, names, component indices, and vertex counts. |
| `histoseg_outputs.zip` | Optional archive containing the run outputs when the Serve instance has enough disk space. |

The contour model works internally in micron coordinate space. Xenium Explorer
exports are written in pixel space using the app's Xenium pixel-size setting.

## Troubleshooting

- **`clusters.csv` is missing `Barcode` or `Cluster`:** export the GraphClust
  assignment table or rename the columns before upload.
- **Cell and cluster IDs do not match:** confirm that both files come from the
  same Xenium run. The app can handle common `-1` barcode suffix differences,
  but the IDs still need to refer to the same cells.
- **The app reports fewer than two clusters:** the dendrogram step needs at
  least two matched cluster labels after file alignment.
- **A cluster ID is repeated:** remove duplicated cluster IDs across manual
  structure lines or selected structures. One cluster can belong to only one
  structure per run.
- **No selected structure has enough cells:** select more populated structures
  or lower the minimum cells parameter in the advanced controls.
- **Large runs are slow or unstable:** lower the partition grid sizes before
  running the final contour analysis.
