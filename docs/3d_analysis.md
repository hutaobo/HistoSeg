# 3D Analysis

3D Analysis is the HistoSeg workflow group for same-sample, multi-slice Xenium
contour reconstruction. It is exposed as `histoseg.threed` because Python
module names cannot begin with a digit.

3D Analysis currently has two public workflows:

- Pairwise conservative TPS soft alignment for a fixed contour GeoJSON and a
  moving contour GeoJSON that has already been hard-aligned with a similarity
  transform.
- Same-sample Xenium contour stack reconstruction from multiple slice folders,
  using GitHub `pyXenium`/`XeniumSlide` for IO, unified Leiden labels from a
  merged AnnData when available, sequential hard+soft alignment, and 3D
  contour surface mesh export.

## Current Scope

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

See the [3D contour soft alignment tutorial](tutorials/3d_soft_alignment)
for a complete polyp contour example with bundled data and rendered QC figures.
See the [3D contour stack reconstruction tutorial](tutorials/3d_stack_reconstruction)
for the 32-slice pyXenium/Leiden stack workflow.
