# 3D Reconstruction

3D Reconstruction is the flagship HistoSeg workflow group. It is designed for
same-sample, multi-slice Xenium contour reconstruction: align neighboring 2D
contour annotations, build 3D contour stacks, export per-structure PLY/OBJ
meshes, and inspect QC artifacts before downstream analysis. The Python
namespace is `histoseg.threed` because module names cannot begin with a digit.

The 3D workflow starts from named semantic contours. In a typical HistoSeg run,
sfplot Search-and-Find / StructureMap relationships or manual biological
curation define candidate structure groups. HistoSeg then converts those groups
into continuous 2D isoline or multi-structure partition contours, and the 3D
stack reconstruction aligns these named contours across slices. The aligned
semantic contours are the source for cell-cloud projection, optional mesh
visualization, and SDF-based gene-structure quantification.

The flagship 3D surface currently has these public workflows:

- Pairwise conservative TPS soft alignment for a fixed contour GeoJSON and a
  moving contour GeoJSON that has already been hard-aligned with a similarity
  transform.
- Same-sample Xenium contour stack reconstruction from multiple slice folders,
  using GitHub `pyXenium`/`XeniumSlide` for IO, unified Leiden labels from a
  merged AnnData when available, sequential hard+soft alignment, and 3D
  contour surface mesh export.
- 3D gene spatial module discovery from an aligned cell table and merged
  AnnData: gene enrichment voxels, nested hotspot surfaces, voxel overlap/SDF
  distance metrics, and gene-by-structure clustermaps.
- 3D cell cloud rendering from aligned cell Parquet or merged AnnData into a
  browser-shareable Plotly HTML, with optional contour overlays.
- Per-gland 3D QC atlas rendering from aligned GeoJSON contours, assigning
  cross-slice gland IDs and writing local zoom pages for small-component
  reconstruction review.

## Current Scope

- Treat 3D contour reconstruction as the main HistoSeg workflow, with 2D
  contour generation and H&E analysis supporting upstream preparation and QC.
- Start from named semantic contours generated from selected or curated
  structure groups, rather than arbitrary polygons.
- Align neighboring 2D Xenium contour slices before downstream 3D
  reconstruction.
- Read ordered Xenium slice folders through `pyXenium.io.read_xenium(...,
  as_="slide")`.
- Reconstruct 3D contour stacks as sampled 3D contour points, per-structure
  smoothed Marching Cubes PLY/OBJ meshes, and an interactive Plotly HTML view.
- Preserve GeoJSON feature properties and annotation structure.
- Report union IoU, per-structure IoU, landmark residuals, and geometry
  validity.
- Write QC overlays and a TPS diagnostic report.

The v1 surface mesh is a conservative voxelized contour-stack reconstruction:
aligned 2D polygons are rasterized into a 3D volume, optionally smoothed with a
3D Gaussian filter, and converted to a triangle surface with Marching Cubes.
It is intended for biomedical 3D reconstruction/QC, not CAD-level editing.

## Method Positioning

| Method or resource | Primary input | Primary output | HistoSeg relationship |
| --- | --- | --- | --- |
| HistoSeg 3D | Ordered Xenium slices with named semantic contours | Topology-checked aligned contours, surface meshes, SDF metrics, gland QC tables and HTML atlases | Reconstructs interpretable contour compartments and records per-pair QC for review. |
| CODA | Serial tissue images and masks | Large-tissue 3D image reconstructions | `coda-image` is a CODA-inspired hard seed only: it rasterizes HistoSeg contour unions, estimates Radon rotation plus phase translation, and competes against the contour seed before HistoSeg TPS. It is not a full CODA reimplementation. |
| Space-map | Serial-section spatial transcriptomics atlas data and analysis context | Atlas-level 3D single-cell tissue reconstruction resource | HistoSeg can consume ordered slices from this type of dataset, but focuses on semantic contour reconstruction, QC, SDF analysis and local component tracking rather than recreating the full Space-map atlas pipeline. |

HistoSeg therefore complements image- or atlas-first 3D methods: it preserves
semantic contour labels, exposes conservative hard-seed tournament provenance,
and makes fold, compression, expansion, disconnected fragment, and gland
branch/merge candidates explicit review artifacts.

Limited numbers of 2D sections can support alignment QC, morphology review,
and hypothesis generation about local continuity. They should not be presented
as strong 3D biological claims without denser sectioning, independent samples,
or orthogonal validation.

## Label-Free Cross-Group Alignment

When two adjacent contour slices were clustered independently, their
`assigned_structure` names may not be comparable across files. In that setting,
HistoSeg can run a label-free preflight that treats labels only as within-slice
groups, evaluates cross-group contour constellations, and estimates the moving
slice transform from the best local overlap anchor. The workflow preserves all
original labels and does not perform semantic harmonization.

See the {doc}`validated breast label-free group alignment example
<label_free_alignment>` for the interactive RTD overlay and command used to
align `breastrep1S2.geojson` with `breastrep2S3.geojson`.

The same group-correspondence seed is available inside stack reconstruction.
With the default `--registration-backend auto`, HistoSeg evaluates both the
standard semantic contour seed and the label-free group seed. The label-free
seed is selected only when the anchor transform is accepted, enough anchor
pairs are used, and the median residual is below the configured limit. If the
selected label-free seed matches different fixed and moving group names,
semantic TPS refinement is skipped for that slice pair because cross-slice
biological structure identity has not been harmonized.

## Python API

```python
from histoseg.threed import (
    ThreeDContourReconstructionConfig,
    run_3d_contour_reconstruction,
)

cfg = ThreeDContourReconstructionConfig(
    fixed_geojson="slice_01.geojson",
    moving_hard_aligned_geojson="slice_02_hard_aligned_to_01.geojson",
    out_dir="outputs/3d_soft_alignment",
    group_property="structure",
    diagnostic_structure="Structure 5",
)

result = run_3d_contour_reconstruction(cfg)
print(result.soft_aligned_geojson)
print(result.diagnostic_report_png)
```

For a multi-slice stack:

```python
from histoseg.threed import (
    ThreeDStackReconstructionConfig,
    run_3d_stack_reconstruction,
)

cfg = ThreeDStackReconstructionConfig(
    xenium_root="polyp",
    segmentation_strategy="polyp/contour for alignment/segmentationstrategy.txt",
    merged_h5ad="polyp/pdc_merge_leiden/polyp_32samples_processed_leiden.h5ad",
    merged_cluster_column="leiden_1_0",
    out_dir="outputs/polyp_3d_reconstruction",
    z_spacing_um=5.0,
)

result = run_3d_stack_reconstruction(cfg)
print(result.visualization_html)
print(result.mesh_dir)
```

For 3D gene spatial modules:

```python
from histoseg.threed import (
    SpatialModuleDiscoveryConfig,
    run_spatial_module_discovery,
)

cfg = SpatialModuleDiscoveryConfig(
    h5ad="polyp/pdc_merge_leiden/polyp_32samples_processed_leiden.h5ad",
    aligned_cells_parquet="outputs/polyp_3d_reconstruction/aligned_leiden_3d_cells.parquet",
    stack_root="outputs/polyp_3d_reconstruction",
    genes=("GREM1", "COL1A1", "FAP", "EPCAM", "MUC2", "LGR5"),
)

result = run_spatial_module_discovery(cfg)
print(result.fraction_inside_matrix_csv)
```

For an interactive 3D cell cloud:

```python
from histoseg.threed import CellCloudRenderConfig, render_cell_cloud_html

result = render_cell_cloud_html(
    CellCloudRenderConfig(
        stack_root="outputs/polyp_3d_reconstruction",
        aligned_cells_parquet="outputs/polyp_3d_reconstruction/aligned_leiden_3d_cells.parquet",
        out_html="outputs/polyp_3d_reconstruction/leiden_3d_cells.html",
        label_column="leiden_1_0",
    )
)
print(result.out_html)
```

For a per-gland QC atlas:

```python
from histoseg.threed import GlandQCAtlasConfig, render_gland_qc_atlas

result = render_gland_qc_atlas(
    GlandQCAtlasConfig(
        stack_root="outputs/polyp_3d_reconstruction",
        aligned_cells_parquet="outputs/polyp_3d_reconstruction/aligned_leiden_3d_cells.parquet",
        out_dir="outputs/polyp_3d_reconstruction/gland_qc",
        structures=("Structure 3", "Structure 4"),
        max_gland_pages=250,
    )
)
print(result.atlas_html)
```

For lumen-seeded gland instance segmentation, tracking, and atlas rendering:

```python
from histoseg.threed import (
    GlandInstanceSegmentationConfig,
    GlandInstanceTrackingConfig,
    run_gland_instance_detection,
)

result = run_gland_instance_detection(
    segmentation_config=GlandInstanceSegmentationConfig(
        stack_root="outputs/polyp_3d_reconstruction",
        aligned_cells_parquet="outputs/polyp_3d_reconstruction/aligned_leiden_3d_cells.parquet",
        out_dir="outputs/polyp_3d_reconstruction/gland_instances",
    ),
    tracking_config=GlandInstanceTrackingConfig(
        segmentation_result_dir="outputs/polyp_3d_reconstruction/gland_instances",
        out_dir="outputs/polyp_3d_reconstruction/gland_instances",
    ),
)
print(result.gland_instance_atlas_html)
```

## CLI

```bash
histoseg-3d reconstruct \
  --fixed-geojson slice_01.geojson \
  --moving-hard-aligned-geojson slice_02_hard_aligned_to_01.geojson \
  --out-dir outputs/3d_soft_alignment \
  --group-property structure \
  --diagnostic-structure "Structure 5"
```

The CLI performs only the soft TPS step. The moving sample should already be
hard-aligned to the fixed sample.

For a same-sample stack:

```bash
histoseg-3d reconstruct-stack \
  --xenium-root polyp \
  --segmentation-strategy "polyp/contour for alignment/segmentationstrategy.txt" \
  --merged-h5ad "polyp/pdc_merge_leiden/polyp_32samples_processed_leiden.h5ad" \
  --merged-cluster-column leiden_1_0 \
  --out-dir outputs/polyp_3d_reconstruction \
  --z-spacing-um 5 \
  --mesh-smoothing-sigma-um 40 \
  --mesh-export-formats ply,obj \
  --no-alignment-preview
```

Set `--mesh-smoothing-sigma-um 0` to disable the 3D Gaussian smoothing step.
`--mesh-level` must stay between `0` and `1` for the Marching Cubes surface.

## Registration Backends

`histoseg-3d reconstruct-stack` defaults to `--registration-backend auto`.
The automatic backend runs the established semantic contour seed
(`contour-tps`) and the label-free cross-group seed (`label-free-group`) as hard
alignment candidates. HistoSeg selects `label-free-group` only when its local
anchor transform passes the configured acceptance checks; otherwise it falls
back to `contour-tps`. For stacks where the tissue silhouette gives a better
global initial pose, use `--registration-backend coda-image` explicitly.

`label-free-group` is geometric only. It preserves original contour properties
and can align a moving slice using a cross-named local contour group, but it
does not rename structures or infer biological identity. When the selected
fixed and moving groups have different names, HistoSeg records
`semantic_soft_skipped_reason="cross_named_label_free_group_match"` and uses the
label-free hard result directly for that pair.

### CODA-Inspired Image Registration

The `coda-image` backend is inspired by CODA's image-registration strategy from
Kiemen et al., "CODA: quantitative 3D reconstruction of large tissues at
cellular resolution," *Nature Methods* 19, 1490-1499 (2022),
[doi:10.1038/s41592-022-01650-9](https://doi.org/10.1038/s41592-022-01650-9),
and the [CODA methodology page](https://labs.pathology.jhu.edu/kiemen/coda-3d/).
It is not a full CODA reimplementation. This first HistoSeg backend rasterizes
the contour union as a tissue-mask proxy, estimates a Radon rotation and global
phase-correlation translation, then passes that hard-aligned result into the
same TPS topology guard used by `contour-tps`.

```bash
histoseg-3d reconstruct-stack \
  --xenium-root polyp \
  --segmentation-strategy "polyp/contour for alignment/segmentationstrategy.txt" \
  --out-dir outputs/polyp_3d_reconstruction \
  --registration-backend coda-image \
  --coda-raster-size 512 \
  --coda-angle-step 1.0
```

Alignment summaries record the backend, CODA-inspired credit, DOI, Radon angle,
phase-correlation shift, and hash-relevant preprocessing parameters so AnnData
cell-cloud caches are invalidated when geometry-defining registration state
changes.

`pairwise_alignment_metrics.csv` records the selected hard seed, contour seed
IoU, CODA seed IoU, seed rotation difference, CODA Radon/phase metadata,
soft-topology validity, checked grid cells, folded/compressed/expanded cells,
and min/median/max TPS area ratios. These columns are the first stop for
broken or folded reconstruction triage before inspecting overlay PNGs and HTML
views.

Render aligned cells into a Plotly/WebGL HTML view:

```bash
histoseg-3d render-cell-cloud \
  --stack-root outputs/polyp_3d_reconstruction \
  --aligned-cells-parquet outputs/polyp_3d_reconstruction/aligned_leiden_3d_cells.parquet \
  --out-html outputs/polyp_3d_reconstruction/leiden_3d_cells.html \
  --label-column leiden_1_0 \
  --max-points 300000
```

Render a local QC atlas for gland-like components:

```bash
histoseg-3d render-gland-qc-atlas \
  --stack-root outputs/polyp_3d_reconstruction \
  --aligned-cells-parquet outputs/polyp_3d_reconstruction/aligned_leiden_3d_cells.parquet \
  --out-dir outputs/polyp_3d_reconstruction/gland_qc \
  --structures "Structure 3" "Structure 4" \
  --max-gland-pages 250
```

This command writes `gland_tracks.csv`, `gland_qc_index.csv`,
`gland_qc_atlas.html`, and one `glands/gland_*.html` page per tracked gland.
The tracking uses neighboring-slice, same-structure contour component links.
The original slice-local `component_index` is preserved as metadata but is not
treated as a global gland identity. The CSV tables always cover every tracked
gland; `--max-gland-pages` limits local HTML rendering to the highest-priority
review targets for faster first-pass QC.

Detect lumen-seeded gland/crypt instances, track them across slices, and render
an instance atlas:

```bash
histoseg-3d detect-gland-instances \
  --stack-root outputs/polyp_3d_reconstruction \
  --aligned-cells-parquet outputs/polyp_3d_reconstruction/aligned_leiden_3d_cells.parquet \
  --out-dir outputs/polyp_3d_reconstruction/gland_instances \
  --epithelial-structures "Structure 3" "Structure 4" \
  --markers EPCAM MUC2 LGR5 OLFM4 MKI67
```

The default tracker uses one-to-one Hungarian assignment between neighboring
slices. Add `--allow-many-to-many` only for exploratory mutual-best linking.
Possible two-to-one or one-to-two branch/merge relationships are reported in
`gland_instance_tracks.csv` and `gland_instance_qc_index.csv` instead of being
silently collapsed.

If the aligned Parquet does not exist yet, render directly from AnnData by
supplying `--h5ad` and `--out-parquet` instead of `--aligned-cells-parquet`.
HistoSeg uses its AnnData alignment cache when available and warns when the
render target is large enough to slow typical browser WebGL sessions.

For end-to-end gene spatial module discovery:

```bash
histoseg-3d discover-spatial-modules \
  --h5ad polyp/pdc_merge_leiden/polyp_32samples_processed_leiden.h5ad \
  --aligned-cells-parquet outputs/polyp_3d_reconstruction/aligned_leiden_3d_cells.parquet \
  --stack-root outputs/polyp_3d_reconstruction \
  --genes GREM1 COL1A1 COL1A2 ACTA2 PDGFRA TAGLN FAP EPCAM MUC2 OLFM4 LGR5 MKI67
```

Advanced users can rerun only the downstream steps:

```bash
histoseg-3d quantify-gene-structure \
  --stack-root outputs/polyp_3d_reconstruction \
  --gene-density-dir outputs/polyp_3d_reconstruction/gene_overlays/batch_3d_genes/GREM1_density \
  --gene GREM1

histoseg-3d plot-spatial-modules \
  --batch-dir outputs/polyp_3d_reconstruction/gene_overlays/batch_3d_genes \
  --matrix fraction_inside \
  --hotspot top05
```

## Outputs

- `soft_aligned_contours.geojson`
- `soft_tps_alignment_metrics.csv`
- `soft_tps_landmarks.csv`
- `soft_tps_diagnostic_residuals.csv`
- `soft_tps_alignment_summary.json`
- `soft_tps_overlay_hard_before.png`
- `soft_tps_overlay_soft_after.png`
- `soft_tps_landmarks_qc.png`
- `soft_tps_diagnostic_report.png`

Multi-slice reconstruction writes:

- `xenium_slice_manifest.csv`
- `aligned_slice_manifest.csv`
- `pairwise_alignment_metrics.csv`
- `aligned_contour_3d_points.csv`
- `histoseg_3d_contour_stack.html`
- `meshes/Structure_*.ply`
- `meshes/Structure_*.obj`
- `meshes/mesh_manifest.csv`
- `meshes/mesh_qc_summary.json`
- `leiden_3d_cells.html` from `render-cell-cloud`
- `gland_qc_atlas.html` and `glands/gland_*.html` from
  `render-gland-qc-atlas`
- `gland_instances/slice_gland_instances.geojson`,
  `gland_instance_tracks.csv`, `gland_instance_qc_index.csv`, and
  `gland_instance_atlas.html` from `detect-gland-instances`

Mesh cleanup removes tiny disconnected fragments by default when a minimum
component volume is provided. Use `--min-mesh-component-volume-um3` to set the
volume threshold, or `--no-mesh-cleanup` to retain all disconnected components.
`meshes/mesh_manifest.csv` records connected-component counts before and after
filtering so fractured surfaces can be audited.

Spatial module discovery writes:

- `<GENE>_density/<GENE>_3d_enrichment_voxels.csv`
- `<GENE>_density/isosurfaces/<GENE>_enrichment_top*.ply`
- `<GENE>_density/structure_relationships/<GENE>_structure_3d_overlap_metrics.csv`
- `<GENE>_density/structure_relationships/<GENE>_structure_3d_distance_metrics.csv`
- `batch_gene_status.csv`
- `gene_structure_overlap_fraction_matrix.csv`
- `gene_structure_signed_distance_matrix.csv`
- `gene_structure_fraction_inside_matrix.csv`
- `fraction_inside_top05_spatial_clustermap.png`

See the [3D contour soft alignment tutorial](tutorials/3d_soft_alignment)
for a complete polyp contour example with bundled data and rendered QC figures.
See the [3D contour stack reconstruction tutorial](tutorials/3d_stack_reconstruction)
for the 32-slice pyXenium/Leiden stack workflow.
