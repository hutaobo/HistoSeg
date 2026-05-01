# Contour Analysis

:::{div} hs-section-kicker
Guide
:::

:::{div} hs-page-intro
Contour Analysis is the spatial or cell-coordinate workflow group exposed as
`histoseg.contour`. Use it when the primary input is a table of cells with
spatial coordinates and cluster assignments, and the goal is to generate
review-ready spatial structures rather than image-first segmentation masks.
:::

:::{div} hs-metadata
<span>Input: cells with x/y coordinates</span>
<span>Outputs: contours, previews, exports</span>
<span>Supports Pattern1 and named structures</span>
:::

## When To Use It

Choose this workflow group when you want to:

- extract Pattern1 isolines from clustered cell coordinates;
- partition cells into multiple named structures with non-overlapping contours;
- export Xenium Explorer-compatible review layers; or
- analyze how generated structures touch, overlap, or enclose one another.

## Workflows

- Pattern1 isoline contours from clustered cell coordinates.
- Multi-structure contour partitioning.
- Xenium Explorer annotation exports.
- Hugging Face dataset helper workflows for Xenium-style public datasets.

## Inputs And Outputs

Pattern1 and multi-structure contour workflows use:

- `clusters.csv` with barcode and cluster columns.
- `cells.parquet` with cell identifiers and x/y coordinates.
- optional `tissue_boundary.csv` for tissue-aware background generation.
- selected cluster IDs for the target Pattern1 or named structures.
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

## Boundary Network Python API

```python
from histoseg.contour import BoundaryNetworkConfig, run_group_boundary_network

result = run_group_boundary_network(
    BoundaryNetworkConfig(
        boundary_csv="/path/to/group_boundary_overlap_filtered.csv",
        out_dir="outputs/boundary_network",
        drop_structures=("1", "6"),
    )
)
```

The boundary network workflow consumes a group-level boundary overlap table,
including either legacy columns such as `total_shared_boundary` and
`mean_shared_boundary` or HistoSeg topology columns such as
`shared_boundary_length_um`.

## Contour Adjacency Python API

```python
from histoseg.contour import ContourAdjacencyConfig, run_contour_adjacency

result = run_contour_adjacency(
    ContourAdjacencyConfig(
        contours="/path/to/contours.csv",
        groupby="structure",
        out_dir="outputs/contour_adjacency",
    )
)
```

The contour adjacency workflow recomputes topology from contour geometries.
Boundary-neighbor pairs contribute shared boundary length, and enclosure pairs
contribute the inner contour boundary length. It writes a long-form edge CSV, a
square adjacency matrix CSV with an empty diagonal, a network plot, and a
heatmap.

## CLI

The CLI surface is aligned with the same workflow concepts used in the Python
API, which makes it straightforward to promote an exploratory run into a batch
pipeline:

```bash
histoseg-contour pattern1 --clusters-csv clusters.csv --cells-parquet cells.parquet --out-dir outputs/p1 --pattern1-clusters 10,23,19
histoseg-contour multi-structure --clusters-csv clusters.csv --cells-parquet cells.parquet --out-dir outputs/ms --structures-json structures.json
histoseg-contour boundary-network --boundary-csv group_boundary_overlap_filtered.csv --out-dir outputs/boundary_network
histoseg-contour adjacency --contours-csv contours.csv --groupby structure --out-dir outputs/contour_adjacency
```

## Related Pages

- {doc}`tutorials/pattern1_isoline`
- {doc}`workflows/contour_generation`
