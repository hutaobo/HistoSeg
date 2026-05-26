# Breast Partial-Anchor Alignment Methods

This note describes the manuscript-ready workflow used for the fully automatic
breast partial-anchor alignment example. The workflow starts from two raw
Xenium output folders, discovers spatial structures without manual landmarks,
generates HistoSeg contour GeoJSON files, aligns the two contour slices, fits an
anchor-only residual TPS soft deformation from the accepted anchors, and renders
a clean before/after figure.

The example is a pairwise 3D alignment quality-control demonstration. It shows
that two related breast cancer sections can be placed into a shared coordinate
frame when only part of the tissue overlaps. It should not be interpreted as a
complete dense 3D biological reconstruction from two sections alone.

## Input Data

Each Xenium output folder contributes two standard files:

- `cells.parquet`, containing cell identifiers and physical x/y centroids.
- `analysis/clustering/gene_expression_graphclust/clusters.csv`, containing
  GraphClust cell-to-cluster assignments with `Barcode` and `Cluster` columns.

The reproducibility script accepts either a direct Xenium `outs` folder or a
parent folder containing one direct child ending in `_outs`. This supports the
breast Rep2 layout where the user-facing folder contains
`Xenium_FFPE_Human_Breast_Cancer_Rep2_outs`.

## Plain-Language Method

1. **Load cell coordinates and GraphClust labels.** HistoSeg joins
   `cells.parquet` with `clusters.csv` using the cell barcode or cell ID. The
   joined table provides one x/y coordinate and one GraphClust label per cell.

2. **Automatically discover spatial structures.** GraphClust labels are treated
   as the initial fine labels. HistoSeg computes a cluster-to-cluster spatial
   distance matrix using the Search-and-Find relationship implemented through
   `Cell-GPS`. The matrix is converted to a cophenetic relationship matrix and
   hierarchically cut into a small number of coarse structures. When
   `cluster_count=auto`, HistoSeg evaluates several candidate cuts and selects
   the structure count that balances separation, within-structure compactness,
   and minimum structure cell fraction.

3. **Generate continuous contours from the selected structures.** For each
   automatically selected `Structure N`, HistoSeg builds a smoothed density
   field from the cells assigned to the structure. The workflow rasterizes all
   structures on a shared x/y grid, applies Gaussian smoothing, forms
   structure-support masks, resolves overlaps into an exclusive partition, and
   extracts polygon boundaries from the partition masks. The result is a
   non-overlapping multi-structure contour layer written as
   `xenium_explorer_annotations.geojson`.

4. **Search for partial-anchor correspondences.** HistoSeg does not assume that
   all contours in the two sections have counterparts. It evaluates every
   fixed-structure and moving-structure pair as a candidate local overlap. For
   each candidate pair, it compares individual contour descriptors, including
   centroid layout, area, shape, and local neighborhood context.

5. **Estimate the hard alignment from automatic anchors only.** The accepted
   transform is fitted from high-confidence RANSAC inlier anchors. The user does
   not provide landmarks and does not specify which structures should match.
   Non-overlapping contours are carried as passive geometry: they receive the
   final transform but do not pull the transform toward themselves.

6. **Fit anchor-only residual TPS soft alignment.** After the hard transform,
   HistoSeg uses only the final accepted anchor points to fit a thin-plate
   residual field. The TPS corrects local residual deformation in the
   overlapping region while adding zero-residual identity padding anchors around
   the tissue bounding box. This constrains deformation away from the evidence
   region and keeps passive or no-counterpart tissue from attracting the warp.

7. **Write auditable outputs.** The run writes the fixed and moving
   auto-structure summaries, the two contour GeoJSON files, the selected anchor
   table, the group-correspondence matrix, the hard-aligned moving GeoJSON, the
   soft-aligned moving GeoJSON, the hard and soft alignment summaries, and a
   clean before/after manuscript panel.

## Hard Seed Versus Soft Deformation

The breast workflow deliberately separates the alignment into two layers:

| Layer | What it estimates | Evidence used | What moves |
| --- | --- | --- | --- |
| Hard seed | One global similarity transform \(T\) with rotation, translation, and scale | final RANSAC inlier anchors from label-free group correspondence | every moving contour |
| Anchor-only residual TPS | A local residual displacement field \(u(x)\) after the hard seed | the same accepted hard-aligned anchors, plus zero-displacement identity padding anchors | every moving contour, but non-anchor contours remain passive |

This distinction is important for partial-overlap tissue. The hard seed solves
the global placement problem. The soft step only relaxes the remaining local
residuals supported by accepted anchors. It does not introduce a new semantic
boundary matching objective, and it does not let unmatched tissue pull the warp.

In implementation terms, the manuscript run first writes
`alignment/moving_group_overlap_aligned.geojson`. That file is the hard
similarity result. If anchor-only residual TPS passes QC, the final manuscript
figure uses
`anchor_only_soft_tps/anchor_only_soft_aligned_contours.geojson` instead.

## Mathematical Description

Let the fixed-slice contour set be

\[
F = \{F_i\}_{i=1}^{n_f}
\]

and the moving-slice contour set be

\[
M = \{M_j\}_{j=1}^{n_m}.
\]

Each contour carries a within-slice structure label \(g(F_i)\) or \(g(M_j)\),
but these labels are not assumed to be harmonized across slices. HistoSeg first
builds contour descriptors

\[
d(C) = \left[c(C), a(C), s(C), \mathcal{N}(C)\right],
\]

where \(c(C)\) is the centroid, \(a(C)\) is area, \(s(C)\) is a shape summary,
and \(\mathcal{N}(C)\) describes the local neighboring-contour layout.

For each fixed group \(G_f\) and moving group \(G_m\), candidate anchor pairs
are scored from descriptor compatibility. RANSAC then searches for an inlier
set

\[
A = \{(i,j): F_i \in G_f,\ M_j \in G_m\}
\]

that supports a similarity transform

\[
T(x) = R_{\theta}x + t.
\]

The residual for an anchor pair is

\[
r_{ij} = \left\|c(F_i) - T(c(M_j))\right\|_2.
\]

The selected transform is the candidate that maximizes inlier support and
descriptor agreement while keeping robust residual summaries low. The transform
is first applied to every moving contour,

\[
M^{(0)}_j = T(M_j),
\]

but only the inlier set \(A\) contributes force to the fit. Contours without a
counterpart therefore remain visible as transformed passive tissue rather than
being forced to overlap unrelated fixed contours.

The soft step then fits a residual displacement field from the same accepted
anchors. Let

\[
a_k = T(c(M_{j_k}))
\]

be the hard-aligned moving anchor centroid and

\[
b_k = c(F_{i_k})
\]

be its fixed target. HistoSeg fits a thin-plate radial basis interpolator
\(u(x)\) such that

\[
u(a_k) \approx b_k - a_k.
\]

To limit extrapolation, identity padding anchors \(p_l\) are placed around an
expanded fixed/moving bounding box with zero displacement,

\[
u(p_l) = 0.
\]

The final soft transform is therefore

\[
T_{\mathrm{soft}}(x) = T(x) + u(T(x)).
\]

The soft result is accepted only if the post-warp anchor residuals are low, no
invalid geometries are produced, and the sampled Jacobian grid has an acceptable
negative-Jacobian fraction. In this example the residual TPS uses evidence from
the 20 accepted anchors and keeps all other contours passive.

## HistoSeg Code Mapping

| Method step | HistoSeg API or CLI | Main output |
| --- | --- | --- |
| Resolve raw Xenium folders | `resolve_xenium_output_folder()` | resolved `outs` path |
| Auto-discover structures | `discover_auto_structures(AutoStructureDiscoveryConfig(...))` or `histoseg-contour auto-structure` | `structures.json`, cluster assignment CSV, matrix CSVs |
| Generate contours | `run_multi_structure_contours(MultiStructureContourConfig(...))` | `xenium_explorer_annotations.geojson` |
| Align contours | `align_contours_label_free(LabelFreeContourAlignmentConfig(..., partial_correspondence=True, group_correspondence=True))` | aligned GeoJSON, anchors CSV, group matrix, summary JSON |
| Fit soft deformation | `run_anchor_only_residual_tps(AnchorOnlyResidualTPSConfig(...))` | soft-aligned GeoJSON, TPS landmarks CSV, TPS summary JSON |
| Render manuscript panel | `render_label_free_before_after_panel(LabelFreeBeforeAfterFigureConfig(...))` | soft before/after PNG and SVG |
| End-to-end reproduction | `python reproducibility/run_breast_partial_anchor_from_xenium.py` | manifest JSON and figure assets |

## Soft Alignment Output Files

The script writes the soft deformation artifacts under
`anchor_only_soft_tps/`:

| File | Purpose |
| --- | --- |
| `anchor_only_soft_aligned_contours.geojson` | final moving-slice contours after hard seed plus residual TPS |
| `anchor_only_tps_landmarks.csv` | accepted anchor rows and identity padding anchors used to fit \(u(x)\) |
| `anchor_only_tps_summary.json` | acceptance status, residual summaries, identity padding count, invalid geometry count, and Jacobian QC |
| `anchor_only_tps_before.png` | hard-aligned moving contours before residual TPS |
| `anchor_only_tps_after.png` | moving contours after residual TPS |
| `anchor_only_tps_review.html` | interactive before/after review with residual links |

The manuscript panel is written to
`figure/breast_partial_anchor_before_after.png`. When soft TPS is accepted, this
panel uses the soft-aligned GeoJSON. The script also writes
`figure/breast_partial_anchor_before_after_hard.png` so the hard-only result can
be compared directly.

## Reproducible Command

Run from the HistoSeg repository root:

```bash
python reproducibility/run_breast_partial_anchor_from_xenium.py \
  --out-dir reproducibility/results/breast_partial_anchor_from_xenium \
  --cluster-count auto \
  --min-structure-count 3 \
  --max-structure-count 8
```

The script defaults to the local breast Xenium paths used during method
development:

```text
Y:\long\10X_datasets\Xenium\Xenium_Breast_Cancer\Xenium_FFPE_Human_Breast_Cancer_Rep1_outs
Y:\long\10X_datasets\Xenium\Xenium_Breast_Cancer\Xenium_FFPE_Human_Breast_Cancer_Rep2
```

Override them with `--fixed-xenium-output` and `--moving-xenium-output` when
running on another machine. The script applies anchor-only residual TPS by
default; pass `--no-anchor-only-soft-tps` to keep the hard similarity output as
the manuscript figure.

The main output manifest is:

```text
reproducibility/results/breast_partial_anchor_from_xenium/breast_partial_anchor_manifest.json
```

The manuscript figure is:

```text
reproducibility/results/breast_partial_anchor_from_xenium/figure/breast_partial_anchor_before_after.png
```

## Frozen Raw-Folder Run

The full raw-folder workflow was executed from the repository root with:

```bash
python reproducibility/run_breast_partial_anchor_from_xenium.py --overwrite
```

The run regenerated both contour inputs from the raw Xenium folders and then
performed label-free group-correspondence alignment. The frozen output lives in:

```text
reproducibility/results/breast_partial_anchor_from_xenium/
```

Key results:

| Quantity | Value |
| --- | --- |
| Fixed auto-structure count | 3 |
| Moving auto-structure count | 4 |
| Alignment method | `label_free_group_correspondence_ransac_alignment` |
| Selected fixed group | `Structure 3` |
| Selected moving group | `Structure 3` |
| Anchor pairs selected | 27 |
| Anchor pairs used for transform | 20 |
| Median anchor residual | 43.57 coordinate units |
| P90 anchor residual | 135.00 coordinate units |
| Maximum anchor residual | 146.14 coordinate units |
| Candidate pair count | 760 |
| RANSAC trial count | 1015 |
| Soft method | `anchor_only_residual_tps` |
| Soft TPS accepted | true |
| Soft TPS anchor count | 20 |
| Soft TPS identity padding anchors | 16 |
| Soft post-warp median anchor residual | 0.87 coordinate units |
| Soft post-warp P90 anchor residual | 1.96 coordinate units |
| Soft negative-Jacobian ratio | 0.0 |

The soft residuals are much smaller than the hard residuals because the TPS is
fitted directly to the residual field left after the similarity transform. This
does not mean all tissue is forced to overlap. It means the accepted local
anchor constellation is made internally consistent while unmatched tissue is
carried through the same field as passive geometry.

These values are from the fully regenerated raw-folder workflow. They differ
from the earlier RTD showcase numbers because the manuscript workflow now
reselects coarse structures automatically from the two Xenium output folders
instead of reusing the previously generated interactive Serve-app GeoJSON
files.

## Figure Legend Draft

Fully automatic partial-anchor soft alignment of two breast cancer Xenium
sections. HistoSeg started from the two raw Xenium output folders, joined cell
coordinates with GraphClust labels, automatically discovered coarse spatial
structures, generated multi-structure contour GeoJSON files, and searched all
fixed/moving-structure combinations for a local anchor constellation. A hard
similarity transform was first fitted from automatically detected
high-confidence anchors. HistoSeg then used the same accepted anchors to fit an
anchor-only residual TPS field, allowing local slice deformation while keeping
non-overlapping tissue passive through identity padding and Jacobian QC. Blue,
fixed contours; red, moving contours before alignment; green, moving contours
after anchor-only residual TPS soft alignment; purple, automatically selected
anchor links.
