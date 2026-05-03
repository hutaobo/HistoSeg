# HistoSeg Main Text Draft

**Frozen baseline:** `publication-alpha-20260503` at commit `2361fbc`

**Draft status:** first Results-centered manuscript draft. Abstract,
Introduction, Discussion, references, and final figure legends remain to be
completed after the Results claims are locked.

## Working Title

HistoSeg: signed-distance reconstruction and 3D spatial quantification of
multi-slice spatial transcriptomics tissue architecture

## Abstract

[To be written after the Results section is finalized.]

## Introduction

[To be written after the Results claims and figure panels are finalized.]

## Results

### HistoSeg reconstructs ordered Xenium slices as auditable 3D tissue objects

HistoSeg was designed to convert ordered multi-slice Xenium outputs into a
single 3D analytical object whose geometry, cell coordinates, and downstream
gene-structure measurements remain traceable to explicit intermediate files
(Fig. 1). The workflow starts from per-slice cell coordinates and structure
annotations, extracts 2D tissue contours, aligns consecutive slices, and exports
an aligned stack manifest. This manifest defines the slice order, physical
\(z\)-positions, accepted contour paths, and alignment role for each slice.
Downstream 3D outputs, including contour point clouds, optional surface meshes,
aligned cell-cloud tables, and gene-by-structure matrices, are generated from
this accepted manifest rather than from ad hoc visualization state.

The alignment step combines a conservative hard similarity transform with an
optional soft thin-plate spline (TPS) warp (Fig. 1B). For each non-reference
slice, HistoSeg first aligns the moving slice to the previously accepted slice.
The hard alignment is accepted only if the union intersection-over-union (IoU)
between matched structures does not decrease. Soft alignment is then fit from
matched structure-boundary landmarks and accepted only when three conditions
are simultaneously satisfied: union IoU does not worsen relative to the hard
alignment, the topology quality-control grid detects no fold-over or excessive
cell-area distortion, and the warped geometries remain valid after repair.
This conservative acceptance rule makes the aligned stack a reviewable object:
each accepted or rejected pairwise step has an associated metric record and,
when enabled, diagnostic overlays.

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
encode anisotropic voxel spacing. For structure mask \(M_j\), HistoSeg computes
an outside transform from \(1-M_j\) and an inside transform from \(M_j\). The
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
the fraction of hotspot voxels with signed distance \(D_j(v) < 0\). The signed
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
the expected behavior for inside voxels, outside \(x\)-steps, outside
\(z\)-steps, and empty-mask cases. In the synthetic sweep, the same tilted
ellipsoid scene was discretized at \(1\times\), \(2\times\), \(5\times\), and
\(10\times\) \(z\)-spacing while preserving \(1\,\mu m\) in-plane sampling.
The `fraction_inside_structure` summary varied from 0.9674 to 0.9698, with a
maximum absolute drift of 0.0025 (0.248 percentage points) relative to the
\(1\times\) baseline. The median signed distance shifted by at most 0.566
\(\mu m\), and the three-summary signed-distance mean absolute drift
(median, 5th percentile, 95th percentile) was at most 1.477 \(\mu m\). Empty
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
(`fraction_inside_structure=0.862`, median signed distance \(-35\) um), FAP
(0.793, \(-35\) um), ACTA2 (0.773, \(-25\) um), COL1A1 (0.767, \(-25\) um),
and COL1A2 (0.757, \(-25\) um). The same structure also embedded immune
markers, including PTPRC (0.872, \(-35\) um), CD3D (0.873, \(-35\) um), CD3E
(0.876, \(-35\) um), MS4A1 (0.805, \(-25\) um), and CD79A (0.791, \(-25\) um).
This combination supports the interpretation of Structure 5 as a
stromal-immune niche rather than a purely stromal compartment.

The epithelial markers separated into at least two spatial programs. Structure
3 carried the strongest MUC2 signal (`fraction_inside_structure=0.603`,
median signed distance \(-10\) um), consistent with a secretory epithelial
component. OLFM4 showed a near-boundary pattern with Structure 5 and Structure
3 signals and should be interpreted as a distributed differentiation or
interface marker rather than forced into a single compartment. MKI67 had its
largest epithelial-associated signal in Structure 3 but a positive median
signed distance in the current `top05` summary, suggesting proximity to this
compartment rather than deep embedding. Together these patterns support a
Structure 3 epithelial differentiation/proliferation axis that should be
presented with this boundary nuance.

Structure 4 carried a distinct epithelial/stem-like readout. LGR5 had its
strongest `top05` fraction-inside value in Structure 4 (0.579, \(-5\) um), and
EPCAM also peaked in Structure 4 (0.415, \(10\) um). The separation between
MUC2-heavy Structure 3 and LGR5/EPCAM-associated Structure 4 suggests that the
32-slice reconstruction preserves epithelial state heterogeneity in 3D rather
than collapsing epithelial markers into a single contiguous compartment.

Structure 2 showed the clearest macrophage/perivascular-associated signal in
this panel. C1QA (0.270, \(30\) um), PDGFRA (0.224, \(45\) um), and CD68
(0.157, \(80\) um) were Structure 2-dominant by `top05` fraction-inside. LYZ
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
wrapper uses the validated 24-gene starter-panel matrices rather than rerunning
the heavier discovery step. It renders the 300,000-point cell-cloud HTML,
regenerates the `top05` fraction-inside and signed-distance clustermaps, and
writes `reproducibility/results_manifest.json` with command records, input
paths, output sizes, output SHA256 hashes, and the recomputed alignment hash.
The same wrapper can rerun `discover-spatial-modules` with `--run-discovery`
when full regeneration is required.

This design keeps the manuscript artifacts linked to the exact HistoSeg CLI
surface while avoiding hidden analysis code in the paper directory. Large raw
inputs, including `.h5ad`, Parquet, and interactive HTML artifacts, are
referenced by path and hash rather than copied into the repository. The
resulting manuscript package therefore separates three layers: the frozen
software baseline, the local data state identified by alignment provenance, and
the derived figure artifacts used in the paper.

## Discussion

[To be written after the Results claims, figure panels, and validation gaps are
reviewed. Candidate themes: SDFs as a transparent physical-distance contract;
topology-aware alignment as a conservative safeguard for sparse multi-slice
data; marker-panel interpretation versus cell-type deconvolution; limitations
of single-sample biological discovery; required second-cohort or second-tissue
generalization.]

## Online Methods

The implementation-faithful Online Methods draft is maintained in
`docs/manuscripts/histoseg_online_methods_sdf_alignment.md`.

## Figure Legends

[To be written after final figure assembly.]

## References

[Reference X: Xenium spatial transcriptomics]

[Reference X: signed distance fields and Euclidean distance transforms]

[Reference X: thin-plate spline registration]

[Reference X: colorectal polyp biology and stromal-immune niches]
