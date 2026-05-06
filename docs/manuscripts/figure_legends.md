# HistoSeg Figure Legends Draft

**Draft baseline:** manuscript evidence anchored to `publication-alpha-20260503`
and the 32-slice polyp validation artifacts regenerated after the measured SDF
anisotropy sweep.

**Status:** first Nature Methods-style figure legend draft. Panel composition
should be revisited after final figure assembly, but all quantitative values
below are tied to current manuscript artifacts.

## Figure 1 | HistoSeg reconstructs multi-slice Xenium tissue architecture as an auditable 3D analysis object

**A,** Ordered multi-slice Xenium inputs are converted into per-slice
multi-structure 2D semantic contour annotations. sfplot Search-and-Find /
cophenetic StructureMap relationships can be used to define or audit structure
groups, and HistoSeg converts the selected groups into continuous isoline or
partition contours. Contours are stored as explicit GeoJSON-derived records,
preserving slice identity, structure labels and physical slice order before 3D
reconstruction. **B,** Consecutive slices are registered with topology-aware
alignment safeguards. The hard step is similarity registration, estimating
rotation, translation and uniform scale before optional TPS refinement. The
workflow records hard-aligned and soft-aligned contour states, displacement
diagnostics and quality-control summaries so that accepted transformations can
be reviewed rather than treated as hidden preprocessing. **C,** Accepted aligned contours
are rasterized onto a 3D `(z, y, x)` grid to form binary structure masks.
Distances used for gene-structure quantification are computed from these masks
with `scipy.ndimage.distance_transform_edt(..., sampling=(z_um, y_um, x_um))`,
with positive values outside a structure and negative values inside. Optional
surface meshes are used for visualization and are not the primitive used for
SDF quantification. **D,** The accepted stack links multiple downstream
representations: aligned 3D cell clouds, sampled contour points, optional
surface meshes, and gene-by-structure spatial-module matrices. In the 32-slice
polyp validation, HistoSeg rendered a browser-shareable 3D Leiden cell cloud
from 2,785,128 aligned cells by deterministic sampling to 300,000 Plotly WebGL
points, with 32 Leiden labels and five contour traces.

**Source artifacts:** `aligned_slice_manifest.csv`;
`docs/_static/manuscript/figure1_workflow_schematic.png`;
`docs/_static/manuscript/figure1_workflow_schematic.svg`;
`docs/_static/manuscript/figure1_histology_contour_overlay_template.png`;
`docs/_static/manuscript/figure1_histology_contour_overlay_template.svg`;
`pairwise_alignment_metrics.csv`;
`pairwise_alignments/*/soft_tps/soft_tps_diagnostic_report.png`;
`aligned_contour_3d_points.csv`; `histoseg_3d_contour_stack.html`;
`meshes/mesh_manifest.csv`;
`reproducibility/results/figure1_leiden_3d_cells_300k.html`;
`reproducibility/results_manifest.json`.

## Figure 2 | Physical signed-distance fields provide stable 3D gene-structure metrics under anisotropic sampling

**A,** HistoSeg defines each tissue structure as a binary 3D mask $M_j$ on a
physical `(z, y, x)` voxel grid. The outside transform
$D_{\mathrm{out},j}$ is computed from $1-M_j$, the inside transform
$D_{\mathrm{in},j}$ is computed from $M_j$, and the signed-distance field
$D_j$ is assembled as positive outside the mask and negative inside the mask.
All distances are reported in microns through physical EDT sampling. **B,** A
5x5x5 anisotropic truth-table fixture validates the discrete SDF contract using
`spacing_zyx_um=(5,1,1)`. The fixture contains one hotspot voxel inside the
structure, one voxel one $x$-step outside and one voxel one $z$-step
outside; all expected values passed exactly, including median signed distance
of 1.0 um, maximum signed distance of 5.0 um and
`fraction_inside_structure=1/3`. **C,** Empty-mask contract validation confirms
that undefined distance cases are represented explicitly. Empty hotspot and
empty structure inputs both returned zero inside/touching fractions and NaN
distance summaries, avoiding false zero-distance or infinite-distance
interpretations. **D,** A synthetic anisotropy sweep discretized the same
tilted ellipsoid structure and hotspot at $1\times$, $2\times$,
$5\times$ and $10\times$ $z$-spacing while preserving 1 um in-plane
sampling and evaluating the public `compute_hotspot_sdf_metrics` API. Across
the sweep, `fraction_inside_structure` ranged from 0.9674 to 0.9698, with a
maximum absolute drift of 0.0025, equivalent to 0.248 percentage points,
relative to the $1\times$ baseline. Median signed distance shifted by at most
0.566 um, and the maximum mean absolute drift across the median, 5th percentile
and 95th percentile signed-distance summaries was 1.477 um.

**Source artifacts:** `reproducibility/validation/sdf_anisotropy_sweep.py`;
`reproducibility/validation/results/sdf_truth_table.csv`;
`reproducibility/validation/results/empty_mask_contract.csv`;
`reproducibility/validation/results/anisotropy_sensitivity.csv`;
`reproducibility/validation/results/anisotropy_sensitivity_fraction_inside.csv`;
`reproducibility/validation/results/anisotropy_sensitivity_signed_distance.csv`;
`reproducibility/validation/results/anisotropy_sensitivity.png`;
`reproducibility/validation/results/validation_manifest.json`.

## Figure 3 | A 32-slice polyp reconstruction resolves stromal-immune, epithelial and macrophage-associated 3D compartments

**A,** Scale view of the 32-slice polyp reconstruction. The aligned cell table
contains 2,785,128 cells across 32 slices and five reconstructed contour
structures; the interactive cell-cloud validation samples 300,000 cells for
browser rendering while preserving the aligned 3D coordinate contract. **B,**
Gene-by-structure clustering using `fraction_inside_structure` for the top 5%
hotspot identifies structure-associated marker programs. Structure 5 contains
the dominant embedded signals for stromal matrix, activated fibroblast,
contractile and immune markers, including GREM1, FAP, ACTA2, COL1A1, PTPRC,
CD3D, CD3E, MS4A1 and CD79A. **C,** The signed-distance clustermap provides a
complementary view of physical embedding. Negative median signed distances
support embedded hotspots within a structure, whereas positive values indicate
nearby but externally positioned hotspots. This separates the strongly embedded
Structure 5 stromal-immune niche from more boundary-associated or distributed
signals. **D,** GREM1 provides the clearest nested-hotspot example in the
current 24-gene panel. The GREM1 top 5% hotspot is Structure 5-dominant
(`fraction_inside_structure=0.862`) and has a median signed distance of
-35 um, linking the 3D hotspot surface to the SDF quantification table. **E,**
Compartment-level summary of marker programs. Structure 5 is interpreted as a
mixed stromal-immune niche; Structure 3 is enriched for epithelial
differentiation/proliferation signals including MUC2, OLFM4 and MKI67;
Structure 4 carries an epithelial/stem-like LGR5/EPCAM readout; and Structure 2
is macrophage/perivascular-associated, with C1QA, PDGFRA and CD68 as the
strongest supporting markers. These assignments are marker-panel
interpretations and should not be read as cell-type deconvolution.

**Source artifacts:** `gene_overlays/batch_3d_genes_starter_panel/batch_gene_status.csv`;
`gene_overlays/batch_3d_genes_starter_panel/gene_structure_fraction_inside_matrix.csv`;
`gene_overlays/batch_3d_genes_starter_panel/gene_structure_signed_distance_matrix.csv`;
`docs/_static/threed/polyp/24gene/fraction_inside_top05_spatial_clustermap.png`;
`docs/_static/threed/polyp/24gene/signed_distance_top05_spatial_clustermap.png`;
`docs/_static/threed/polyp/24gene/GREM1_nested_3d_hotspot_surfaces.png`;
`docs/_static/threed/polyp/24gene/GREM1_structure5_focus.png`;
`docs/_static/threed/polyp/24gene/signed_distance_marker_group_violin.png`;
`docs/_static/threed/polyp/24gene/compartment_summary_table.csv`;
`docs/_static/threed/polyp/24gene/compartment_summary_table.png`;
`reproducibility/results/figure1_leiden_3d_cells_300k.html`;
`reproducibility/results_manifest.json`;
`docs/manuscripts/polyp_24_gene_evidence.tsv`.

## Figure 4 | Benchmark matrix for geometry, topology and expression-domain consistency

**A,** Planned method matrix for the Nature Methods Article benchmark. The
matrix includes naive stacking, hard-only HistoSeg, full HistoSeg, PASTE,
PASTE2, GPSA, STAligner and SPACEL/Scube, with adapter status, version policy,
default-parameter policy and seed recording specified in
`reproducibility/benchmarks/method_matrix.json`. **B,** Synthetic
known-transform smoke benchmark for geometry recovery. The CI-safe smoke run
reports known-transform error, union IoU, mean per-structure IoU, centroid
drift and landmark distance, allowing the benchmark machinery to be tested
without private data. Published external methods are explicitly marked as
`synthetic_proxy` in this smoke run until full adapters are executed on the
archived public bundle. **C,** Topology and expression-domain consistency
summary. Fold and compression failures are reported as metrics rather than
filtered away, and expression/domain consistency is retained as a separate
biological plausibility readout. **D,** Runtime and memory summary. Wall-clock
runtime and peak memory are recorded for each method and should be regenerated
for both the public demo subset and the full archived dataset before
submission.

**Source artifacts:** `reproducibility/benchmarks/method_matrix.json`;
`reproducibility/benchmarks/method_benchmark_metrics.csv`;
`reproducibility/benchmarks/benchmark_summary.csv`;
`reproducibility/benchmarks/benchmark_manifest.json`;
`reproducibility/benchmarks/figure4_benchmark_matrix.png`;
`reproducibility/benchmarks/figure4_benchmark_matrix.svg`.
