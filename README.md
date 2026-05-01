<div align="center">

# HistoSeg

<p align="center">
  <a href="https://pypi.org/project/histoseg/"><img alt="PyPI" src="https://img.shields.io/pypi/v/histoseg.svg"></a>
  <a href="https://histoseg.readthedocs.io/en/latest/"><img alt="Docs" src="https://readthedocs.org/projects/histoseg/badge/?version=latest"></a>
  <a href="https://github.com/hutaobo/HistoSeg/actions/workflows/publish.yml"><img alt="Publish to PyPI" src="https://github.com/hutaobo/HistoSeg/actions/workflows/publish.yml/badge.svg?branch=main"></a>
  <a href="https://polyformproject.org/licenses/noncommercial/1.0.0/"><img alt="License: PolyForm Noncommercial 1.0.0" src="https://img.shields.io/badge/License-PolyForm--Noncommercial%201.0.0-blue.svg"></a>
</p>

**Python toolkit for H&E image analysis, spatial contour analysis, and 3D contour reconstruction.**

</div>

HistoSeg has three public feature groups:

- **HE Analysis** (`histoseg.he`) for image-based H&E tissue segmentation, neutral tissue partitioning, and aligned-image change detection.
- **Contour Analysis** (`histoseg.contour`) for contour extraction from spatial/cell-coordinate data, including Pattern1 isolines and multi-structure Xenium exports.
- **3D Analysis** (`histoseg.threed`) for same-sample, multi-slice Xenium contour alignment and reconstruction workflows.

Full documentation: [histoseg.readthedocs.io](https://histoseg.readthedocs.io)

## When To Use Each Feature Group

Use **HE Analysis** when your input is an H&E image such as PNG, JPG, TIFF, or GeoTIFF and you want masks, overlays, heatmaps, GeoJSON polygons, or region tables.

Use **Contour Analysis** when your input is spatial cell-coordinate data such as Xenium `cells.parquet` plus cluster assignments, and you want geometry extracted from cell neighborhoods or selected cluster groups.

Use **3D Analysis** when you are preparing for multi-slice Xenium contour reconstruction from the same sample. It can soft-align a hard-aligned moving contour GeoJSON to a fixed reference slice, or build a pyXenium-backed multi-slice contour stack with 3D points, PLY meshes, and an interactive HTML view.

## Installation

```bash
pip install -U histoseg
```

For local Hugging Face MedSAM-backed HE segmentation:

```bash
pip install -U "histoseg[he]"
```

For development:

```bash
git clone https://github.com/hutaobo/HistoSeg.git
cd HistoSeg
pip install -U pip
pip install -e ".[he]"
```

For 3D Xenium stack reconstruction, install the GitHub pyXenium extra:

```bash
pip install -e ".[threed]"
```

## HE Analysis Quickstart

```python
from histoseg.he import HESegmentationConfig, run_he_segmentation

result = run_he_segmentation(
    HESegmentationConfig(
        image="/path/to/he.png",
        out_dir="outputs/he_all_elements",
        task="all_elements",
        backend="heuristic",
        n_components=6,
    )
)

print(result.overlay_png)
print(result.geojson)
```

```bash
histoseg-he all-elements \
  --image /path/to/he.png \
  --out-dir outputs/he_all_elements \
  --backend heuristic
```

HE Analysis currently supports:

- `single`: tissue foreground extraction, or user-prompted region extraction from boxes/points
- `all_elements`: neutral tissue component partitioning (`component_1`, `component_2`, ...)
- `change`: aligned before/after H&E change detection

## Contour Analysis Quickstart

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
print(result.preview_png)
print(len(result.contours))
```

```bash
histoseg-contour pattern1 \
  --clusters-csv clusters.csv \
  --cells-parquet cells.parquet \
  --out-dir outputs/pattern1 \
  --pattern1-clusters 10,23,19
```

Contour Analysis currently supports:

- Pattern1 isoline contour generation from clustered cell coordinates
- multi-structure contour partitioning
- Xenium Explorer annotation exports
- Hugging Face dataset helper workflows for Xenium-style inputs

## 3D Analysis Quickstart

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

```bash
histoseg-3d reconstruct \
  --fixed-geojson slice_01.geojson \
  --moving-hard-aligned-geojson slice_02_hard_aligned_to_01.geojson \
  --out-dir outputs/3d_soft_alignment \
  --group-property structure \
  --diagnostic-structure "Structure 5"
```

For a complete same-sample Xenium stack:

```bash
histoseg-3d reconstruct-stack \
  --xenium-root polyp \
  --segmentation-strategy "polyp/contour for alignment/segmentationstrategy.txt" \
  --merged-h5ad "polyp/pdc_merge_leiden/polyp_32samples_processed_leiden.h5ad" \
  --merged-cluster-column leiden_1_0 \
  --out-dir outputs/polyp_3d_reconstruction \
  --z-spacing-um 5
```

3D Analysis writes aligned per-slice GeoJSON files, pairwise alignment metrics,
sampled 3D contour points, per-structure PLY meshes, and an interactive Plotly
HTML visualization.

## Outputs

HistoSeg workflows write reviewable artifacts such as:

- PNG previews and overlays
- label maps and heatmaps
- GeoJSON polygons
- CSV/Parquet region or contour tables
- `params.json` and `metrics.json` provenance files

## Documentation

- [HE Analysis](https://histoseg.readthedocs.io/en/latest/he_analysis.html)
- [Contour Analysis](https://histoseg.readthedocs.io/en/latest/contour_analysis.html)
- [3D Analysis](https://histoseg.readthedocs.io/en/latest/3d_analysis.html)
- [3D contour soft alignment tutorial](https://histoseg.readthedocs.io/en/latest/tutorials/3d_soft_alignment.html)
- [3D contour stack reconstruction tutorial](https://histoseg.readthedocs.io/en/latest/tutorials/3d_stack_reconstruction.html)
- [API Reference](https://histoseg.readthedocs.io/en/latest/api/index.html)

## License

This project is distributed under the **PolyForm Noncommercial 1.0.0** license.
Academic and other noncommercial use is permitted.
Any commercial use requires a separate commercial license from **SPATHO AB**.
See [LICENSE](LICENSE) for details.
