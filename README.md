<div align="center">

# HistoSeg

<p align="center">
  <a href="https://pypi.org/project/histoseg/"><img alt="PyPI" src="https://img.shields.io/pypi/v/histoseg.svg"></a>
  <a href="https://histoseg.readthedocs.io/en/latest/"><img alt="Docs" src="https://readthedocs.org/projects/histoseg/badge/?version=latest"></a>
  <a href="https://github.com/hutaobo/HistoSeg/actions/workflows/publish.yml"><img alt="Publish to PyPI" src="https://github.com/hutaobo/HistoSeg/actions/workflows/publish.yml/badge.svg?branch=main"></a>
  <a href="https://polyformproject.org/licenses/noncommercial/1.0.0/"><img alt="License: PolyForm Noncommercial 1.0.0" src="https://img.shields.io/badge/License-PolyForm--Noncommercial%201.0.0-blue.svg"></a>
</p>

**Python toolkit for H&E image analysis and spatial contour analysis.**

</div>

HistoSeg has two public feature groups:

- **HE Analysis** (`histoseg.he`) for image-based H&E tissue segmentation, neutral tissue partitioning, and aligned-image change detection.
- **Contour Analysis** (`histoseg.contour`) for contour extraction from spatial/cell-coordinate data, including Pattern1 isolines and multi-structure Xenium exports.

Full documentation: [histoseg.readthedocs.io](https://histoseg.readthedocs.io)

## When To Use Each Feature Group

Use **HE Analysis** when your input is an H&E image such as PNG, JPG, TIFF, or GeoTIFF and you want masks, overlays, heatmaps, GeoJSON polygons, or region tables.

Use **Contour Analysis** when your input is spatial cell-coordinate data such as Xenium `cells.parquet` plus cluster assignments, and you want geometry extracted from cell neighborhoods or selected cluster groups.

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
- [API Reference](https://histoseg.readthedocs.io/en/latest/api/index.html)

## License

This project is distributed under the **PolyForm Noncommercial 1.0.0** license.
Academic and other noncommercial use is permitted.
Any commercial use requires a separate commercial license from **SPATHO AB**.
See [LICENSE](LICENSE) for details.
