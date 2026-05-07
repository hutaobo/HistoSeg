---
orphan: true
---

# HistoSeg Main Text Draft

**Frozen baseline:** `publication-alpha-20260503` at commit `2361fbc`

**Draft status:** Nature Methods Article working draft with Discussion,
references, availability statements, benchmark scaffold and figure-generation
contracts in place. Public accession and DOI fields remain pending until final
GEO/Zenodo deposits are assigned.

## Working Title

HistoSeg: StructureMap-guided semantic contours for signed-distance
reconstruction and 3D spatial quantification of multi-slice spatial
transcriptomics tissue architecture

## Abstract

Spatial transcriptomics is typically analyzed as a set of two-dimensional
fields, even when tissue architecture is inherently three-dimensional.
HistoSeg reconstructs ordered multi-slice Xenium data into auditable 3D tissue
objects by combining sfplot-guided semantic contour definition,
topology-aware contour alignment, deterministic cell-cloud projection and
physically sampled signed-distance fields (SDFs). StructureMap relationships
support the selection of biologically meaningful structure groups; HistoSeg then
converts those groups into continuous contours for 3D reconstruction. SDF
metrics are computed directly from 3D structure masks with anisotropic voxel
spacing, enabling stable gene-structure quantification across sparse slice
sampling. In synthetic validation, increasing $z$-spacing tenfold changed
fraction-inside summaries by at most 0.248 percentage points and median signed
distance by at most 0.566 $\mu m$. Applied to a 32-slice polyp stack with
2,785,128 aligned cells, HistoSeg resolved marker-supported stromal-immune,
epithelial and macrophage/perivascular-associated compartments, including a
GREM1-associated Structure 5 niche.

## Introduction

Spatial transcriptomics has made it possible to profile gene expression in
histological context, but most computational workflows still treat each tissue
section as a self-contained two-dimensional plane. This is a practical
simplification rather than a biological truth. Epithelial folds, stromal
interfaces, vascular tracts and immune niches extend through depth, and their
interpretation can change when a marker is viewed as an isolated 2D signal
instead of as part of a 3D anatomical neighborhood. Multi-slice assays therefore
create an opportunity: if serial sections can be registered, quantified and
audited in physical coordinates, spatial transcriptomics can begin to describe
tissue architecture as a volume rather than as a gallery of adjacent images.

This opportunity is technically difficult. Spatial transcriptomics slices are
sparse along the $z$-axis relative to their in-plane resolution, and
section-to-section deformation can be large enough that simple stacking
produces misleading geometry. Registration errors are especially harmful when
they are propagated into downstream gene-structure statistics: a visually
plausible surface can still assign a hotspot to the wrong tissue compartment,
and a missing or empty structure can silently become a false distance. A
methods framework for 3D spatial transcriptomics must therefore do more than
render a volume. It must preserve slice provenance, guard against topological
fold-over, expose its alignment state, and define gene-structure association in
physical units that remain interpretable under anisotropic sampling.

HistoSeg addresses this problem by treating 3D reconstruction as a set of
explicit contracts, beginning with semantic contour definition. The upstream
sfplot layer computes Search-and-Find relationships and cophenetic StructureMap
summaries that help define or audit biologically meaningful structure groups.
HistoSeg then consumes selected or curated groups and converts them into
continuous isoline or multi-structure partition contours. These contours are
semantic anatomical objects: they carry structure labels, can be reviewed as
GeoJSON/CSV annotations, and can be propagated through alignment and 3D
quantification.

Ordered 2D semantic contours are aligned with a conservative hard similarity
registration step followed by optional thin-plate spline soft alignment, each
accepted only when objective geometry and topology checks pass. The hard step
estimates rotation, translation and uniform scale; it is therefore not a
pure rigid-only registration, although rigid registration is its scale-fixed
special case. Accepted contours define a stack manifest that anchors downstream
cell-cloud projection, surface visualization and gene-structure quantification.
For quantification, HistoSeg rasterizes aligned contours into 3D structure masks
and computes signed-distance fields using
`scipy.ndimage.distance_transform_edt` with physical sampling
`(z_um, y_um, x_um)`. Distances are negative inside a structure, positive
outside it, and undefined empty-mask cases are reported as NaN summaries with
zero inside/touching fractions rather than as hidden failures.

Here we describe the HistoSeg 3D workflow and validate its core SDF
quantification contract. We first show how topology-aware alignment and
manifest-based provenance convert a serial Xenium stack into an auditable 3D
analysis object. We then demonstrate that the implemented SDF metrics have
stable behavior in unit-scale truth tables, empty-mask edge cases and a
synthetic anisotropy sweep. Finally, we apply the workflow to a 32-slice polyp
dataset containing 2,785,128 aligned cells and a 24-gene marker panel. The
analysis resolves interpretable 3D compartments, including a mixed
stromal-immune Structure 5 niche with embedded GREM1, fibroblast, extracellular
matrix, T-cell and B-cell marker signals. These results position HistoSeg as a
reproducible foundation for moving spatial transcriptomics analysis from
section-level visualization toward volumetric tissue architecture.

## Results

### HistoSeg converts semantic structure maps into auditable 3D tissue objects

HistoSeg was designed to convert ordered multi-slice Xenium outputs into a
single 3D analytical object whose geometry, cell coordinates, and downstream
gene-structure measurements remain traceable to explicit intermediate files
(Fig. 1). The workflow begins with a semantic structure layer. sfplot
Search-and-Find / cophenetic StructureMap relationships can be used to identify
or audit groups of spatial labels that behave as coherent tissue structures.
HistoSeg consumes those selected or curated groups, converts them into
continuous 2D semantic contours, aligns consecutive slices, and exports an
aligned stack manifest. This manifest defines the slice order, physical
$z$-positions, accepted contour paths, and alignment role for each slice.
Downstream 3D outputs, including contour point clouds, optional surface meshes,
aligned cell-cloud tables, and gene-by-structure matrices, are generated from
this accepted manifest rather than from ad hoc visualization state.

The alignment step combines a conservative hard similarity registration with an
optional soft thin-plate spline (TPS) warp (Fig. 1B). For each non-reference
slice, HistoSeg first aligns the moving slice to the previously accepted slice.
The hard step estimates rotation, translation and uniform scale by optimizing
contour-union IoU from PCA/centroid-derived and multi-start seeds; it is
accepted only if the union IoU between matched structures does not decrease.
Soft alignment is then fit from matched structure-boundary landmarks and
accepted only when topology and geometry quality-control checks pass. With the
default per-structure acceptance rule, a structure keeps its soft-aligned
geometry only when its per-structure IoU does not regress; otherwise the
hard-aligned geometry is retained. This conservative acceptance rule makes the
aligned stack a reviewable object: each accepted or rejected pairwise step has
an associated metric record and, when enabled, diagnostic overlays.

The 32-slice polyp reconstruction provides the first large-scale use case for
this contract. The frozen publication-alpha baseline contains a 32-slice stack
with 5 contour structures and 2,785,128 aligned cells. The cell cloud was
rendered with `histoseg-3d render-cell-cloud` into a browser-viewable Plotly
HTML after deterministic downsampling to 300,000 points (Fig. 1A,D). The
reproduction manifest records the rendered cell count, source cell count, label
count, contour trace count, output hashes, and stack alignment hash. In the
current baseline, the reconstructed stack resolves to
`alignment_hash=5a3517710b1f2666ed7bfb0625bb47a7e16fd715e8df56144291acaff33e1017`,
which binds the cell-cloud and spatial-module outputs to a specific geometry
state.

### Physical signed-distance fields define gene-structure association metrics

HistoSeg quantifies gene localization against 3D structure masks with an
implementation-level signed-distance field (SDF) contract rather than with
surface-rendering geometry (Fig. 2A). Each structure is represented as a binary
voxel mask on a `(z, y, x)` grid. Distances are computed with
`scipy.ndimage.distance_transform_edt` using physical sampling
`(z_um, y_um, x_um)`, so the returned distances are in microns and directly
encode anisotropic voxel spacing. For structure mask $M_j$, HistoSeg computes
an outside transform from $1-M_j$ and an inside transform from $M_j$. The
signed field is positive outside the structure and negative inside the
structure. Boundary voxels are not clamped to zero; the implementation uses the
discrete inside and outside transforms directly.

Gene-specific enrichment fields are computed on the same 3D grid. For each
gene, HistoSeg accumulates per-voxel expression and total cell counts, smooths
both fields, and divides smoothed expression by smoothed cellularity only where
the smoothed cell count is valid. Nested hotspot masks are then defined from
empirical quantiles of the smoothed valid enrichment field. The implemented
levels `top15`, `top10`, and `top05` correspond to the 85th, 90th, and 95th
percentiles, respectively. A hotspot voxel is retained only when it is in the
quantifiable valid mask and its smoothed enrichment is greater than or equal to
the corresponding threshold.

For each gene, hotspot level, and structure, HistoSeg reports complementary
overlap and SDF metrics (Fig. 2B,C). The `fraction_inside_structure` metric is
the fraction of hotspot voxels with signed distance $D_j(v) < 0$. The signed
distance distribution is summarized by the minimum, median, mean, maximum, and
5th and 95th percentiles over hotspot voxels. The `signed_distance` module
matrix uses the median signed distance, whereas the `fraction_inside` module
matrix uses `fraction_inside_structure`. Empty cases are explicit: an empty
hotspot or empty structure returns NaN distance summaries and zero inside or
touching fractions. This prevents missing geometry from being interpreted as an
infinite biological distance.

The current manuscript baseline supports these claims through implementation
tests, frozen methods text, and a synthetic anisotropy sweep. Unit-scale
fixtures assert the anisotropic distance contract on small 3D masks, including
the expected behavior for inside voxels, outside $x$-steps, outside
$z$-steps, and empty-mask cases. In the synthetic sweep, the same tilted
ellipsoid scene was discretized at $1\times$, $2\times$, $5\times$, and
$10\times$ $z$-spacing while preserving $1\,\mu m$ in-plane sampling.
The `fraction_inside_structure` summary varied from 0.9674 to 0.9698, with a
maximum absolute drift of 0.0025 (0.248 percentage points) relative to the
$1\times$ baseline. The median signed distance shifted by at most 0.566
$\mu m$, and the three-summary signed-distance mean absolute drift
(median, 5th percentile, 95th percentile) was at most 1.477 $\mu m$. Empty
hotspot and empty structure cases both returned zero inside/touching fractions
and NaN distance summaries. Fig. 2 should therefore be framed as a measured
validation of the implemented physical-distance contract and its edge-case
behavior, not as an external reconstruction-method benchmark.

### The 24-gene polyp panel resolves biologically interpretable 3D compartments

We next applied the frozen HistoSeg 3D workflow to a 24-gene starter panel in
the 32-slice polyp stack (Fig. 3). All 24 genes completed the spatial-module
batch run. The resulting `top05` matrices revealed reproducible gene-structure
patterns across stromal, epithelial, immune, macrophage-associated, and
vascular marker classes. These interpretations are marker-panel readouts, not
cell-type deconvolution claims, and structure names remain numeric pending
orthogonal histology review.

Structure 5 emerged as the dominant mixed stromal-immune niche (Fig. 3B,C).
Multiple extracellular-matrix, activated-fibroblast, and contractile stromal
markers had their strongest `top05` fraction-inside values in Structure 5 and
negative median signed distances, placing their high-enrichment hotspots inside
the reconstructed contour. Representative values include GREM1
(`fraction_inside_structure=0.862`, median signed distance $-35$ um), FAP
(0.793, $-35$ um), ACTA2 (0.773, $-25$ um), COL1A1 (0.767, $-25$ um),
and COL1A2 (0.757, $-25$ um). The same structure also embedded immune
markers, including PTPRC (0.872, $-35$ um), CD3D (0.873, $-35$ um), CD3E
(0.876, $-35$ um), MS4A1 (0.805, $-25$ um), and CD79A (0.791, $-25$ um).
This combination supports the interpretation of Structure 5 as a
stromal-immune niche rather than a purely stromal compartment.

The epithelial markers separated into at least two spatial programs. Structure
3 carried the strongest MUC2 signal (`fraction_inside_structure=0.603`,
median signed distance $-10$ um), consistent with a secretory epithelial
component. OLFM4 showed a near-boundary pattern with Structure 5 and Structure
3 signals and should be interpreted as a distributed differentiation or
interface marker rather than forced into a single compartment. MKI67 had its
largest epithelial-associated signal in Structure 3 but a positive median
signed distance in the current `top05` summary, suggesting proximity to this
compartment rather than deep embedding. Together these patterns support a
Structure 3 epithelial differentiation/proliferation axis that should be
presented with this boundary nuance.

Structure 4 carried a distinct epithelial/stem-like readout. LGR5 had its
strongest `top05` fraction-inside value in Structure 4 (0.579, $-5$ um), and
EPCAM also peaked in Structure 4 (0.415, $10$ um). The separation between
MUC2-heavy Structure 3 and LGR5/EPCAM-associated Structure 4 suggests that the
32-slice reconstruction preserves epithelial state heterogeneity in 3D rather
than collapsing epithelial markers into a single contiguous compartment.

Structure 2 showed the clearest macrophage/perivascular-associated signal in
this panel. C1QA (0.270, $30$ um), PDGFRA (0.224, $45$ um), and CD68
(0.157, $80$ um) were Structure 2-dominant by `top05` fraction-inside. LYZ
was more distributed and lower-fraction, so it should be treated as supportive
myeloid context rather than a defining Structure 2 marker. The positive signed
distances for several of these markers indicate that this compartment is
associated with, but not necessarily deeply embedded within, the Structure 2
mask at the `top05` threshold. This finding should therefore be framed as a
macrophage/perivascular-associated compartment, not a pure macrophage region.

### GREM1 provides a nested 3D hotspot example for Structure 5

GREM1 provides the clearest visual and quantitative example of HistoSeg's
biological readout (Fig. 3D). In the 24-gene panel, the GREM1 `top05` hotspot
was Structure 5-dominant by fraction-inside and had a negative median signed
distance, indicating that the highest-enrichment voxels were embedded inside
the Structure 5 mask rather than merely adjacent to it. The nested hotspot
surface gives a visual counterpart to the SDF table: GREM1, FAP, ACTA2,
COL1A1, and COL1A2 mark the stromal arm of the niche, while PTPRC, CD3D, CD3E,
MS4A1, and CD79A show that immune markers occupy the same 3D structural
habitat. This combination makes GREM1/Structure 5 the current best candidate
for a hero panel linking topology, signed-distance quantification, and
marker-based interpretation.

The GREM1 example also illustrates the value of reporting both fraction-inside
and signed-distance summaries. Fraction-inside establishes that a large share
of the high-enrichment hotspot falls within Structure 5. Signed distance
establishes that the hotspot is embedded in the discrete SDF sense. The two
metrics therefore answer different questions: how much of the hotspot belongs
to a structure, and how physically embedded the hotspot is within that
structure. This distinction is important for boundary markers such as OLFM4 and
LYZ, where a single fraction or distance metric can otherwise overstate a
compartment label.

### A reproducibility wrapper binds manuscript figures to the frozen baseline

The publication-alpha baseline includes a reproducibility wrapper,
`reproducibility/run_paper_pipeline.py`, that regenerates the current
paper-facing artifacts through public `histoseg-3d` commands. By default, the
wrapper resolves required inputs from
`reproducibility/paper_data_manifest.json` under `--data-root`, using local
paths only when explicitly supplied as overrides. It uses the validated
24-gene starter-panel matrices rather than rerunning the heavier discovery
step, renders the 300,000-point cell-cloud HTML, regenerates the `top05`
fraction-inside and signed-distance clustermaps, and writes
`reproducibility/results_manifest.json` with command records, manifest-linked
input records, output sizes, output SHA256 hashes, and the recomputed alignment
hash. The same wrapper can rerun `discover-spatial-modules` with
`--run-discovery` when full regeneration is required.

This design keeps the manuscript artifacts linked to the exact HistoSeg CLI
surface while avoiding hidden analysis code in the paper directory. Large raw
inputs, including `.h5ad`, Parquet, and interactive HTML artifacts, are
referenced through manifest records, public DOI/accession fields, size and hash
metadata rather than copied into the repository. The resulting manuscript
package therefore separates three layers: the frozen software baseline, the
public data bundle identified by permanent records, and the derived figure
artifacts used in the paper.

## Discussion

HistoSeg addresses a narrow but recurring gap in multi-slice spatial
transcriptomics: how to turn curated two-dimensional tissue structures into a
three-dimensional object whose geometry, topology, cell coordinates and
gene-structure distances are auditable. Most established spatial-transcriptomic
methods operate on a different primary target. PASTE and PASTE2 align spot
clouds across adjacent slices through optimal-transport formulations and are
well suited to transcriptome-driven stack alignment, including partial overlap
in PASTE2. GPSA estimates shared coordinates through Gaussian-process warps.
STAligner and SPACEL use graph or deep-learning representations to integrate
spatial domains across slices and, in SPACEL/Scube, to construct 3D stacks.
Squidpy and BANKSY provide scalable spatial-omics analysis, neighborhood
statistics, cell typing and domain segmentation rather than an explicit
contour-to-SDF physical-distance contract. SEQUOIA illustrates a complementary
image-to-transcriptome direction: it predicts expression patterns from
histology rather than reconstructing measured multi-slice molecular geometry.
HistoSeg is therefore not a replacement for these tools. Its contribution is a
structure-first layer that preserves semantic contours, exposes hard and soft
registration diagnostics, and computes downstream association metrics from
physically sampled 3D masks.

The signed-distance-field representation is deliberately conservative. By
computing distances from voxelized masks with explicit `(z, y, x)` sampling,
HistoSeg makes the unit of every reported distance a property of the data
contract rather than a by-product of a visualization mesh. This avoids several
failure modes of mesh-based proximity measurements, including sensitivity to
surface triangulation, local holes, smoothing choices and inconsistent
inside/outside orientation. The tradeoff is that SDF accuracy is limited by
voxel discretization and by the sparse sampling of physical sections in the z
dimension. The anisotropy sweep quantifies this limit for a smooth synthetic
scene: increasing z spacing tenfold changed fraction-inside by less than 0.25
percentage points and median signed distance by less than 0.6 um. These values
do not prove universal accuracy; they establish that the implemented EDT
contract is numerically stable under a controlled geometry and that future
biological claims should report voxel size, slice spacing and sensitivity
summaries alongside any signed-distance interpretation.

Topology-aware alignment is also a safeguard rather than a guarantee. Hard
similarity registration keeps the transformation easy to inspect, whereas
TPS/RBF soft alignment can absorb local deformation between sections. The same
flexibility can introduce folds, local compression or biologically implausible
warps if accepted uncritically. HistoSeg therefore treats displacement fields,
area changes and topology QC as first-class outputs. This is especially
important for sparse z sampling: adjacent sections are not repeated
measurements of the same physical plane, and aggressive warping can erase real
through-depth variation. The appropriate benchmark is consequently broader
than a single overlap score. The publication benchmark matrix must report
known-transform recovery, union and per-structure IoU, centroid drift,
landmark distance, topology failures, expression/domain consistency, runtime
and memory for naive stacking, HistoSeg ablations and published alignment
methods. The repository now includes the matrix and a deterministic synthetic
smoke benchmark, but those smoke rows are not a substitute for the full
adapter-based benchmark on the public bundle.

The 32-slice polyp case study demonstrates how the geometric contract supports
biological hypothesis generation, but it remains a single-sample discovery
example. Structure 5 is supported by stromal matrix, activated fibroblast and
immune markers, and GREM1 has the strongest embedded hotspot in the current
panel. This is consistent with prior work linking GREM1 and fibroblast-rich
intestinal niches, and with FAP-positive stromal remodeling in colorectal
tumor microenvironments. However, the evidence is a marker-panel
interpretation. It is not cell-type deconvolution, it does not estimate
cell-state proportions, and it does not establish that Structure 5 is a
histologically reviewed compartment. Final claims require H&E or morphology
overlays with structure outlines, orthogonal pathology review, and either a
second polyp cohort or a second tissue system showing that the same
contour-to-SDF workflow behaves predictably outside this specimen.

The software and data-release requirements are therefore part of the method,
not administrative afterthoughts. A Nature Methods Article needs public data
and code records that allow reviewers to rerun the analysis without private
Windows paths. The reproduction wrapper now resolves inputs through a
paper-data manifest with DOI/accession fields and fails by default if required
public records are still pending. This makes the remaining gap explicit:
raw/processed expression data, aligned contours, masks, GeoJSON, figure source
data, benchmark outputs and the frozen software release must be archived with
permanent identifiers before submission. If the publication code snapshot can
be relicensed under an OSI-approved license, the Code Availability statement
should say so and cite the archived release. If rights-holder approval is not
granted, the manuscript should state the restrictive access terms plainly and
seek editorial guidance before full submission.

## Online Methods

The implementation-faithful Online Methods draft is maintained in
`docs/manuscripts/histoseg_online_methods_sdf_alignment.md`.

## Figure Legends

The current draft legends are maintained in
`docs/manuscripts/figure_legends.md`.

## Data availability

The 32-slice polyp Xenium dataset, processed expression object, aligned
contours, GeoJSON annotations, binary masks, meshes, figure source data,
benchmark outputs and reproduction manifests will be deposited before
submission. Expression-style data should be deposited in GEO where appropriate;
derived reconstruction artifacts and figure source data should be deposited in
Zenodo. The repository-level data contract is maintained in
`reproducibility/paper_data_manifest.json`, which records `role`,
`relative_path`, `accession_or_doi`, `sha256`, `size_bytes`,
`license_or_access_terms` and `generated_by` for each required input or output.
Current accession fields are marked pending and must be replaced with final
GEO/Zenodo identifiers before peer review.

## Code availability

The HistoSeg source code is available from the project repository and the
publication version will be frozen as a tagged release with a Zenodo DOI before
submission. The current repository license is PolyForm Noncommercial 1.0.0.
Because Springer Nature encourages licenses approved by the Open Source
Initiative for shared research code, the preferred publication plan is to
release the manuscript code snapshot under Apache-2.0 or BSD-3-Clause if
rights-holder approval is granted. If relicensing is not possible, the final
Code Availability statement must state the non-commercial access terms and the
authors should query the editor before submission.

## References

1. Ståhl, P. L. et al. Visualization and analysis of gene expression in tissue
   sections by spatial transcriptomics. Science 353, 78-82 (2016).
   https://doi.org/10.1126/science.aaf2403
2. Moses, L. & Pachter, L. Museum of spatial transcriptomics. Nat. Methods 19,
   534-546 (2022). https://doi.org/10.1038/s41592-022-01409-2
3. Janesick, A. et al. High resolution mapping of the tumor microenvironment
   using integrated single-cell, spatial and in situ analysis. Nat. Commun. 14,
   8353 (2023). https://doi.org/10.1038/s41467-023-43458-x
4. Palla, G. et al. Squidpy: a scalable framework for spatial omics analysis.
   Nat. Methods 19, 171-178 (2022).
   https://doi.org/10.1038/s41592-021-01358-2
5. Singhal, V. et al. BANKSY unifies cell typing and tissue domain segmentation
   for scalable spatial omics data analysis. Nat. Genet. 56, 431-441 (2024).
   https://doi.org/10.1038/s41588-024-01664-3
6. Zeira, R., Land, M. & Raphael, B. J. Alignment and integration of spatial
   transcriptomics data. Nat. Methods 19, 567-575 (2022).
   https://doi.org/10.1038/s41592-022-01459-6
7. Liu, X., Zeira, R. & Raphael, B. J. Partial alignment of multislice
   spatially resolved transcriptomics data. Genome Res. 33, 1124-1132 (2023).
   https://doi.org/10.1101/gr.277670.123
8. Jones, A. et al. Alignment of spatial genomics data using deep Gaussian
   processes. Nat. Methods 20, 1379-1387 (2023).
   https://doi.org/10.1038/s41592-023-01972-2
9. Zhou, X., Dong, K. & Zhang, S. Integrating spatial transcriptomics data
   across different conditions, technologies and developmental stages. Nat.
   Comput. Sci. 3, 894-906 (2023).
   https://doi.org/10.1038/s43588-023-00528-w
10. Weng, W. et al. SPACEL: deep learning-based characterization of spatial
    transcriptome architectures. Nat. Commun. 14, 7603 (2023).
    https://doi.org/10.1038/s41467-023-43220-3
11. Clifton, K. et al. STalign: alignment of spatial transcriptomics data using
    diffeomorphic metric mapping. Nat. Commun. 14, 8123 (2023).
    https://doi.org/10.1038/s41467-023-43915-7
12. Pizurica, M. et al. Digital profiling of gene expression from histology
    images with linearized attention. Nat. Commun. 15 (2024).
    https://doi.org/10.1038/s41467-024-54182-5
13. Biancalani, T. et al. Deep learning and alignment of spatially resolved
    single-cell transcriptomes with Tangram. Nat. Methods 18, 1352-1362
    (2021). https://doi.org/10.1038/s41592-021-01264-7
14. Hu, J. et al. SpaGCN: integrating gene expression, spatial location and
    histology to identify spatial domains and spatially variable genes by graph
    convolutional network. Nat. Methods 18, 1342-1351 (2021).
    https://doi.org/10.1038/s41592-021-01255-8
15. Zhao, E. et al. Spatial transcriptomics at subspot resolution with
    BayesSpace. Nat. Biotechnol. 39, 1375-1384 (2021).
    https://doi.org/10.1038/s41587-021-00935-2
16. Dries, R. et al. Giotto: a toolbox for integrative analysis and
    visualization of spatial expression data. Genome Biol. 22, 78 (2021).
    https://doi.org/10.1186/s13059-021-02286-2
17. Virtanen, P. et al. SciPy 1.0: fundamental algorithms for scientific
    computing in Python. Nat. Methods 17, 261-272 (2020).
    https://doi.org/10.1038/s41592-019-0686-2
18. Maurer, C. R. Jr, Qi, R. & Raghavan, V. A linear time algorithm for
    computing exact Euclidean distance transforms of binary images in arbitrary
    dimensions. IEEE Trans. Pattern Anal. Mach. Intell. 25, 265-270 (2003).
    https://doi.org/10.1109/TPAMI.2003.1177156
19. Bookstein, F. L. Principal warps: thin-plate splines and the decomposition
    of deformations. IEEE Trans. Pattern Anal. Mach. Intell. 11, 567-585
    (1989). https://doi.org/10.1109/34.24792
20. Duchon, J. Splines minimizing rotation-invariant semi-norms in Sobolev
    spaces. In Constructive Theory of Functions of Several Variables 85-100
    (Springer, 1977). https://doi.org/10.1007/BFb0086566
21. Schindelin, J. et al. Fiji: an open-source platform for biological-image
    analysis. Nat. Methods 9, 676-682 (2012).
    https://doi.org/10.1038/nmeth.2019
22. van der Walt, S. et al. scikit-image: image processing in Python. PeerJ 2,
    e453 (2014). https://doi.org/10.7717/peerj.453
23. Hunter, J. D. Matplotlib: a 2D graphics environment. Comput. Sci. Eng. 9,
    90-95 (2007). https://doi.org/10.1109/MCSE.2007.55
24. Pedregosa, F. et al. Scikit-learn: machine learning in Python. J. Mach.
    Learn. Res. 12, 2825-2830 (2011).
25. Kosinski, C. et al. Gene expression patterns of human colon tops and basal
    crypts and BMP antagonists as intestinal stem cell niche factors. Proc.
    Natl Acad. Sci. USA 104, 15418-15423 (2007).
    https://doi.org/10.1073/pnas.0707210104
26. Davis, H. et al. Aberrant epithelial GREM1 expression initiates colonic
    tumorigenesis from cells outside the stem cell niche. Nat. Med. 21, 62-70
    (2015). https://doi.org/10.1038/nm.3750
27. Karagiannis, G. S. et al. Fibroblast-derived Gremlin1 localises to
    epithelial cells at the base of the intestinal crypt. Oncotarget 10,
    4630-4639 (2019). https://doi.org/10.18632/oncotarget.27048
28. Henry, L. R. et al. Clinical implications of fibroblast activation protein
    in patients with colon cancer. Clin. Cancer Res. 13, 1736-1741 (2007).
    https://doi.org/10.1158/1078-0432.CCR-06-1746
29. Gerling, M. et al. Stromal Hedgehog signalling is downregulated in colon
    cancer and its restoration restrains tumour growth. Nat. Commun. 7, 12321
    (2016). https://doi.org/10.1038/ncomms12321
30. Nature Methods. What makes a Nature Methods paper. Nat. Methods 19,
    1125-1126 (2022). https://doi.org/10.1038/s41592-022-01558-4
31. Springer Nature. Code policy. https://www.springernature.com/gp/open-science/code-policy
32. Nature Portfolio. Data availability statements and data citation policy:
    guidance for authors.
    https://www.nature.com/documents/nr-data-availability-statements-data-citations.pdf
