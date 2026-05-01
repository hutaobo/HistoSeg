# 3D Analysis

3D Analysis is the HistoSeg workflow group for same-sample, multi-slice Xenium
contour reconstruction. It is exposed as `histoseg.threed` because Python
module names cannot begin with a digit.

The first released 3D Analysis workflow is a conservative TPS soft-alignment
step. It takes a fixed contour GeoJSON and a moving contour GeoJSON that has
already been hard-aligned with a similarity transform, then applies a local
thin-plate-spline displacement field to the moving contours.

## Current Scope

- Align neighboring 2D Xenium contour slices before downstream 3D
  reconstruction.
- Preserve GeoJSON feature properties and annotation structure.
- Report union IoU, per-structure IoU, landmark residuals, and geometry
  validity.
- Write QC overlays and a TPS diagnostic report.

This v1 workflow does not yet generate a 3D mesh or volume. It provides the
registration/pre-reconstruction layer needed before stacking multiple slices.

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

See the [3D contour soft alignment tutorial](tutorials/3d_soft_alignment)
for a complete polyp contour example with bundled data and rendered QC figures.
