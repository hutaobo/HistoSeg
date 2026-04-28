# HE Analysis

HE Analysis is the image-based H&E workflow group exposed as `histoseg.he`.
Use it when the primary input is a histology image and the desired output is a
mask, overlay, heatmap, GeoJSON polygon, or region table.

## Workflows

- `single`: extract tissue foreground by default, or extract user-prompted regions from boxes/points.
- `all_elements`: partition tissue into neutral region components such as `component_1`, `component_2`, and `component_3`.
- `change`: compare aligned before/after H&E images and export changed regions.

HistoSeg does not assign diagnostic labels such as tumor, stroma, or necrosis by
default. Region names are neutral unless supplied by the user or a downstream
classification workflow.

## Inputs

- PNG, JPG, TIFF, or GeoTIFF H&E images.
- Optional prompt boxes or points for `single` region extraction.
- Aligned before/after images for `change` detection.

GeoTIFF coordinate metadata is preserved in exported GeoJSON when `rasterio` is
available.

## Outputs

- PNG overlay images.
- Label maps and change heatmaps.
- GeoJSON polygons.
- CSV/Parquet region tables.
- `params.json` and `metrics.json` provenance files.

## Python API

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
```

For local Hugging Face MedSAM inference, install `histoseg[he]` and use
`backend="medsam"`.

```python
from histoseg.he import HEChangeDetectionConfig, run_he_change_detection

result = run_he_change_detection(
    HEChangeDetectionConfig(
        before_image="before.png",
        after_image="after.png",
        out_dir="outputs/he_change",
    )
)
```

## CLI

```bash
histoseg-he single --image /path/to/he.png --out-dir outputs/he_single --backend heuristic
histoseg-he all-elements --image /path/to/he.png --out-dir outputs/he_all --backend heuristic
histoseg-he change --before before.png --after after.png --out-dir outputs/he_change
```

## Synthetic Image Tutorial

The dependency-light `heuristic` backend is useful for tests and examples:

```python
from pathlib import Path

import numpy as np
from skimage import draw, io

from histoseg.he import HESegmentationConfig, run_he_segmentation

out = Path("outputs/he_synthetic")
out.mkdir(parents=True, exist_ok=True)

image = np.full((128, 160, 3), 255, dtype=np.uint8)
rr, cc = draw.ellipse(64, 80, 42, 52)
image[rr, cc] = [232, 185, 205]
rr, cc = draw.disk((54, 62), 18, shape=image.shape[:2])
image[rr, cc] = [118, 72, 145]

image_path = out / "synthetic_he.png"
io.imsave(image_path, image, check_contrast=False)

result = run_he_segmentation(
    HESegmentationConfig(
        image=image_path,
        out_dir=out / "all_elements",
        task="all_elements",
        backend="heuristic",
        n_components=3,
        min_region_area_px=20,
    )
)

print(result.overlay_png)
print(result.geojson)
```
