:orphan:

# Label-Free Contour Group Alignment

This page documents a validated label-free contour alignment example for two
independently clustered Xenium contour slices. The example is intentionally
kept separate from the same-sample stack reconstruction workflow because the
slice-level `assigned_structure` names are not assumed to match across files.

## Validated Breast Example

The interactive overlay below shows the accepted cross-group alignment result
for `breastrep1S2.geojson` and `breastrep2S3.geojson`. The command preserved
all original contour labels and used them only as within-slice groups. It
automatically selected a local cross-group correspondence:

- fixed group: `Structure 2`
- moving group: `Structure 3`
- accepted inlier pairs: 26
- refined anchor pairs: 22
- median anchor residual: 51.7 coordinate units
- rotation: 1.94 degrees
- translation: `dx=-1410.1`, `dy=8958.6`

[Open the standalone interactive overlay](_static/threed/breast/label_free_group_correspondence_20260507/group_overlap_alignment_overlay.html)

```{raw} html
<iframe
  src="_static/threed/breast/label_free_group_correspondence_20260507/group_overlap_alignment_overlay.html"
  title="Label-free breast contour group alignment overlay"
  style="width:100%; height:760px; border:1px solid #d9dde3; border-radius:6px;"
  loading="lazy">
</iframe>
```

## Command

```bash
histoseg-3d align-contours-label-free \
  --fixed-geojson breastrep1S2.geojson \
  --moving-geojson breastrep2S3.geojson \
  --out-dir histoseg_label_free_breast_group_correspondence_cli_20260507 \
  --partial-correspondence \
  --group-correspondence \
  --overwrite
```

## Interpretation

This mode is for cases where two contour slices have a real spatial overlap but
cannot be aligned by forcing every contour in one file to correspond to a
contour in the other file. HistoSeg first treats `assigned_structure` as a
within-slice grouping variable, evaluates all fixed-group and moving-group
combinations, and estimates the final rigid transform from the best local
RANSAC inlier set. Labels are preserved and no semantic harmonization is
claimed.

The output should be read as a geometric preflight and QC result. It identifies
which local contour constellation can act as the alignment anchor, then moves
the full moving slice according to that local transform. Regions outside the
overlapping tissue area are allowed to remain non-overlapping.
