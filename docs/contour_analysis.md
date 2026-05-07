# Contour Analysis

:::{div} hs-section-kicker
Guide
:::

:::{div} hs-page-intro
Contour Analysis is the spatial or cell-coordinate workflow group exposed as
`histoseg.contour`. Use it when the primary input is a table of cells with
spatial coordinates and cluster assignments, and the goal is to generate
review-ready semantic spatial structures rather than image-first segmentation
masks. Structure groups can be curated directly or selected after sfplot
Search-and-Find / cophenetic StructureMap analysis.
:::

:::{div} hs-metadata
<span>Input: cells with x/y coordinates</span>
<span>Outputs: contours, previews, exports</span>
<span>Supports StructureMap-guided semantic contours, Pattern1, gene isolines, and named structures</span>
:::

## When To Use It

Choose this workflow group when you want to:

- extract Pattern1 isolines from clustered cell coordinates;
- extract gene/transcript-defined isolines from Xenium transcript tables;
- partition cells into multiple named structures with non-overlapping contours;
- convert selected or curated structure groups into continuous semantic
  contours for 3D reconstruction and SDF quantification;
- export Xenium Explorer-compatible review layers; or
- analyze how generated structures touch, overlap, or enclose one another.

## Semantic Contour Chain

HistoSeg's contour workflow is designed to bridge relationship-level structure
analysis and geometry-level reconstruction. sfplot can be used upstream to
compute Search-and-Find relationships and cophenetic StructureMap summaries
between spatial labels. Those results help define or audit structure groups,
while HistoSeg produces the final continuous contour geometry. HistoSeg consumes
selected or curated structure groups and synthesizes continuous semantic
boundaries from the observed cells or transcripts.

The implemented chain is:

1. Use sfplot Search-and-Find / StructureMap analysis, or expert curation, to
   identify a structure group as a set of cluster labels.
2. Use HistoSeg Pattern1, gene-isoline, or multi-structure workflows to turn
   that group into a continuous isoline or partition contour.
3. Export named GeoJSON/CSV contours for review, alignment, 3D reconstruction,
   and SDF-based gene-structure quantification.

The selection step remains explicit: users provide target clusters, transcript
targets, or multi-structure specifications, and HistoSeg builds the geometric
contours from those definitions.

## Workflows

- Pattern1 isoline contours from clustered cell coordinates.
- Gene/transcript isoline contours from Xenium transcript tables.
- Multi-structure contour partitioning.
- Xenium Explorer annotation exports.
- Hugging Face dataset helper workflows for Xenium-style public datasets.

For an interactive dendrogram-guided entry point, use the
{doc}`AI Driven Spatial Pathologist Serve app <ai_driven_spatial_pathologist_serve_app>`.
The web app builds candidate structures from uploaded Xenium cluster tables,
lets you select dendrogram branches, and then runs the multi-structure contour
workflow with downloadable HistoSeg outputs.

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

Gene/transcript isolines use Xenium `transcripts.parquet` plus `cells.parquet`.
Selected transcripts are adapted into the Pattern1 engine as the target class,
with cell centroids providing the background.

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

## Gene/Transcript Isoline Python API

```python
from histoseg.contour import GeneTranscriptIsolineConfig, run_gene_transcript_isoline

result = run_gene_transcript_isoline(
    GeneTranscriptIsolineConfig(
        xenium_root="/path/to/polyp",
        out_dir="outputs/gene_isolines",
        genes=("GREM1", "COL1A1"),
        sample_glob="A079-C-008_*",
        qv_min=20,
    )
)
print(result.run_log_csv)
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
histoseg-contour gene-isoline --xenium-root polyp --sample-glob "A079-C-008_*" --genes GREM1,COL1A1 --out-dir outputs/gene_isolines
histoseg-contour multi-structure --clusters-csv clusters.csv --cells-parquet cells.parquet --out-dir outputs/ms --structures-json structures.json
histoseg-contour boundary-network --boundary-csv group_boundary_overlap_filtered.csv --out-dir outputs/boundary_network
histoseg-contour adjacency --contours-csv contours.csv --groupby structure --out-dir outputs/contour_adjacency
```

## Related Pages

- {doc}`ai_driven_spatial_pathologist_serve_app`
- {doc}`tutorials/pattern1_isoline`
- {doc}`workflows/contour_generation`
