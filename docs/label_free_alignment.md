# Fully Automatic Partial-Anchor Breast Alignment

This page showcases a two-slice breast cancer contour alignment where the
overlapping tissue is only a local part of the two slices. It is a pairwise 3D
alignment QC example: HistoSeg estimates how the moving contour slice should be
placed relative to the fixed slice before downstream stack reconstruction, while
also making non-overlapping tissue explicit instead of forcing it to match.

The important point is that this run is fully automatic. HistoSeg did not use
manual landmarks, did not require the user to specify which structures should
match, and did not harmonize labels between slices. The method found the local
anchor constellation from the contour geometry itself.

## Result

The before/after panels below use clean documentation renders with no internal
legend or axis text. Blue is the fixed breast slice, red is the moving slice
before alignment, green is the moving slice after alignment, and purple marks
the automatically detected anchor links used by the accepted transform.

```{raw} html
<div style="display:grid; grid-template-columns:repeat(auto-fit,minmax(320px,1fr)); gap:18px; align-items:start; margin:1.2rem 0;">
  <figure style="margin:0;">
    <img
      src="_static/threed/breast/label_free_group_correspondence_20260507/breast_partial_anchor_before_clean.png"
      alt="Breast contour slices before fully automatic partial-anchor alignment"
      style="width:100%; border:1px solid #d9dde3; border-radius:6px; background:white;"
    />
    <figcaption style="font-size:0.92rem; color:#536273; margin-top:0.45rem;">
      <strong>Before.</strong> The moving slice is still in its original local
      position. The large non-overlapping context is not treated as a failure.
    </figcaption>
  </figure>
  <figure style="margin:0;">
    <img
      src="_static/threed/breast/label_free_group_correspondence_20260507/breast_partial_anchor_after_clean.png"
      alt="Breast contour slices after fully automatic partial-anchor alignment"
      style="width:100%; border:1px solid #d9dde3; border-radius:6px; background:white;"
    />
    <figcaption style="font-size:0.92rem; color:#536273; margin-top:0.45rem;">
      <strong>After.</strong> The local overlap is aligned by automatically
      discovered anchors; non-overlapping tissue remains passive geometry.
    </figcaption>
  </figure>
</div>
```

The accepted alignment is between `breastrep1S2.geojson` and
`breastrep2S3.geojson`. HistoSeg selected a cross-group local anchor match from
independently clustered contour labels:

| Metric | Value |
| --- | --- |
| Selected fixed group | `Structure 2` |
| Selected moving group | `Structure 3` |
| Anchor pairs selected | `26` |
| Anchor pairs used for transform | `22` |
| Median anchor residual | `51.74` coordinate units |
| P90 anchor residual | `148.83` coordinate units |
| Rotation | `1.94` degrees |
| Translation | `dx=-1410.1`, `dy=8958.6` |

## Interactive Review

Open the standalone overlay when you want to zoom, pan, and inspect the anchor
links:

```{raw} html
<p>
  <a href="_static/threed/breast/label_free_group_correspondence_20260507/group_overlap_alignment_overlay.html"
     target="_blank" rel="noopener">
    Open the standalone interactive overlay
  </a>
</p>
```

```{raw} html
<iframe
  src="_static/threed/breast/label_free_group_correspondence_20260507/group_overlap_alignment_overlay.html"
  title="Fully automatic partial-anchor breast contour alignment overlay"
  style="width:100%; height:760px; border:1px solid #d9dde3; border-radius:6px;"
  loading="lazy">
</iframe>
```

## Why This Is A Partial-Anchor Alignment

The two slices are not treated as if every contour in one file must have a
counterpart in the other file. HistoSeg first treats `assigned_structure` as a
within-slice grouping variable, evaluates all fixed-group and moving-group
constellations, and then estimates the transform from the best local RANSAC
inlier set.

Only the high-confidence `matched_anchor` pairs fit the transform. The
`matched_review` and `no_counterpart` contours are transformed as passive
geometry, so tissue outside the overlapping region can remain non-overlapping.
This is the behavior needed when two slices share a local 3D neighborhood but
also contain tissue that is absent, shifted, torn, or outside the counterpart
field of view.

## What Was Automated

- Landmark discovery was automatic: the transform used `22` accepted anchor
  pairs selected from `26` candidate inliers.
- Group selection was automatic: the best supported local correspondence was
  `Structure 2` in the fixed slice to `Structure 3` in the moving slice.
- Label handling was conservative: original contour labels were preserved, and
  no semantic harmonization was claimed.
- The full moving slice was transformed, but only anchor contours generated
  alignment force; non-counterpart regions were carried as passive geometry.

This makes the example useful as a 3D alignment QC demonstration for sparse or
partial serial-section overlap. It should not be read as standalone evidence
for dense 3D biology from two slices alone.

## Raw-Folder Soft TPS Result

The static showcase above uses the original GeoJSON overlay produced for the
partial-anchor demonstration. The manuscript reproduction run starts one step
earlier, from the two raw Xenium output folders, and the final figure uses an
additional soft alignment stage:

1. HistoSeg first finds a label-free group-correspondence seed and estimates a
   hard similarity transform from the final RANSAC inlier anchors.
2. It then keeps exactly those accepted anchors as evidence points and fits an
   anchor-only residual TPS displacement field. The soft field corrects local
   deformation after the hard transform.
3. Contours that were not accepted as anchors remain passive geometry. They are
   warped by the fitted field but do not pull the fit toward themselves.
4. Zero-displacement identity anchors are added around the tissue bounding box
   to limit extrapolation away from the overlap evidence.

```{raw} html
<figure style="margin:1.2rem 0;">
  <img
    src="_static/threed/breast/raw_xenium_anchor_only_tps_20260519/breast_partial_anchor_soft_before_after.png"
    alt="Raw Xenium breast partial-anchor alignment with anchor-only residual TPS soft deformation"
    style="width:100%; border:1px solid #d9dde3; border-radius:6px; background:white;"
  />
  <figcaption style="font-size:0.92rem; color:#536273; margin-top:0.45rem;">
    Frozen raw-folder manuscript run. The moving slice is first placed by a
    hard similarity transform and then locally deformed by anchor-only residual
    TPS. Only automatically selected anchors drive the deformation.
  </figcaption>
</figure>
```

For the frozen raw-folder run, HistoSeg regenerated both contour inputs from the
two Xenium folders and selected a different automatically discovered structure
constellation from the older Serve-app GeoJSON showcase. These are the soft-TPS
metrics for the frozen manuscript workflow:

| Metric | Value |
| --- | --- |
| Selected fixed group | `Structure 3` |
| Selected moving group | `Structure 3` |
| Anchor pairs selected | `27` |
| Anchor pairs used by hard transform | `20` |
| Hard median anchor residual | `43.57` coordinate units |
| Hard P90 anchor residual | `135.00` coordinate units |
| Soft method | `anchor_only_residual_tps` |
| Soft TPS identity padding anchors | `16` |
| Soft post-warp median anchor residual | `0.87` coordinate units |
| Soft post-warp P90 anchor residual | `1.96` coordinate units |
| Negative-Jacobian ratio | `0.0` |

The important distinction is that this is not semantic boundary attraction.
Global contour overlap is not the acceptance target. The soft stage is accepted
when the anchor residuals improve, output geometries stay valid, and the sampled
Jacobian grid does not show fold-over.

## Reproduce The Run

### From Raw Xenium Output Folders

For the manuscript workflow, HistoSeg can regenerate the contour inputs instead
of starting from pre-existing GeoJSON files. The end-to-end script starts from
two Xenium output folders, automatically discovers coarse `Structure N` groups,
generates multi-structure contour GeoJSON files, runs the same label-free
partial-anchor hard alignment, fits anchor-only residual TPS from the accepted
anchors, and writes a clean soft-aligned before/after figure plus a provenance
manifest.

```bash
python reproducibility/run_breast_partial_anchor_from_xenium.py \
  --out-dir reproducibility/results/breast_partial_anchor_from_xenium \
  --cluster-count auto \
  --min-structure-count 3 \
  --max-structure-count 8
```

The script defaults to the local breast Rep1/Rep2 Xenium folders used for this
example and resolves nested `_outs` folders automatically. Use
`--fixed-xenium-output` and `--moving-xenium-output` to point at another copy of
the data. Pass `--no-anchor-only-soft-tps` to disable the soft deformation and
keep the hard similarity output as the manuscript figure. The detailed method
description, contour-generation and residual-TPS formulae, API mapping, and
figure-legend draft are in
{doc}`manuscripts/breast_partial_anchor_alignment_methods`.

The default command writes both hard and soft outputs:

| Output | Meaning |
| --- | --- |
| `alignment/moving_group_overlap_aligned.geojson` | hard similarity result from label-free anchors |
| `anchor_only_soft_tps/anchor_only_soft_aligned_contours.geojson` | final soft-aligned moving contours |
| `anchor_only_soft_tps/anchor_only_tps_landmarks.csv` | accepted anchors plus identity padding anchors used for residual TPS |
| `anchor_only_soft_tps/anchor_only_tps_summary.json` | residual, Jacobian, validity, and acceptance metrics |
| `anchor_only_soft_tps/anchor_only_tps_review.html` | interactive review of hard-before-soft and soft-after alignment |
| `figure/breast_partial_anchor_before_after.png` | manuscript panel using the accepted soft output |
| `figure/breast_partial_anchor_before_after_hard.png` | hard-only comparison panel |

### From Existing GeoJSON Files

```bash
histoseg-3d align-contours-label-free \
  --fixed-geojson breastrep1S2.geojson \
  --moving-geojson breastrep2S3.geojson \
  --out-dir histoseg_label_free_breast_partial_anchor_20260507 \
  --partial-correspondence \
  --group-correspondence \
  --overwrite
```

The same label-free group-correspondence seed is also available in
`histoseg-3d reconstruct-stack`. With `--registration-backend auto`, HistoSeg
compares the standard semantic contour seed with the label-free group seed and
selects the label-free candidate only when the automatically discovered anchor
transform passes the configured count and residual thresholds. When the selected
hard seed is `label-free-group`, the stack soft-alignment stage uses
anchor-only residual TPS instead of semantic boundary attraction.

## Source Artifacts

- {download}`Alignment summary JSON <_static/threed/breast/label_free_group_correspondence_20260507/group_overlap_alignment_summary.json>`
- {download}`Automatically selected anchor table <_static/threed/breast/label_free_group_correspondence_20260507/group_ransac_anchors.csv>`
- {download}`Group correspondence matrix CSV <_static/threed/breast/label_free_group_correspondence_20260507/group_correspondence_matrix.csv>`
- {download}`Raw-folder soft TPS summary JSON <_static/threed/breast/raw_xenium_anchor_only_tps_20260519/anchor_only_tps_summary.json>`
- {download}`Raw-folder soft TPS landmarks CSV <_static/threed/breast/raw_xenium_anchor_only_tps_20260519/anchor_only_tps_landmarks.csv>`
- {download}`Raw-folder hard-only comparison PNG <_static/threed/breast/raw_xenium_anchor_only_tps_20260519/breast_partial_anchor_hard_before_after.png>`

```{raw} html
<p>
  <a href="_static/threed/breast/label_free_group_correspondence_20260507/group_correspondence_matrix.html"
     target="_blank" rel="noopener">
    Open the group correspondence matrix HTML
  </a>
</p>
<p>
  <a href="https://github.com/hutaobo/HistoSeg/blob/main/reproducibility/results/breast_partial_anchor_from_xenium/anchor_only_soft_tps/anchor_only_tps_review.html"
     target="_blank" rel="noopener">
    Open the raw-folder anchor-only TPS review HTML on GitHub
  </a>
</p>
```
