from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from skimage import draw, io

from histoseg.he import (
    HEChangeDetectionConfig,
    HERegionSpec,
    HESegmentationConfig,
    run_he_change_detection,
    run_he_segmentation,
)
from histoseg.he.cli import main as he_cli_main


def _write_synthetic_he_image(path: Path, *, changed: bool = False, shape=(128, 160)) -> Path:
    image = np.full((shape[0], shape[1], 3), 255, dtype=np.uint8)
    rr, cc = draw.ellipse(shape[0] // 2, shape[1] // 2, shape[0] // 3, shape[1] // 3)
    image[rr, cc] = np.array([232, 185, 205], dtype=np.uint8)

    rr, cc = draw.disk((shape[0] // 2 - 14, shape[1] // 2 - 22), 18, shape=shape)
    image[rr, cc] = np.array([118, 72, 145], dtype=np.uint8)

    rr, cc = draw.disk((shape[0] // 2 + 18, shape[1] // 2 + 28), 16, shape=shape)
    image[rr, cc] = np.array([207, 112, 140], dtype=np.uint8)

    if changed:
        image[72:96, 96:126] = np.array([88, 52, 120], dtype=np.uint8)

    io.imsave(path, image, check_contrast=False)
    return path


def test_single_heuristic_extracts_tissue_foreground(tmp_path):
    image_path = _write_synthetic_he_image(tmp_path / "he.png")
    result = run_he_segmentation(
        HESegmentationConfig(
            image=image_path,
            out_dir=tmp_path / "single",
            task="single",
            backend="heuristic",
            min_region_area_px=20,
        )
    )

    assert result.n_regions == 1
    assert result.label_map_png is not None and result.label_map_png.exists()
    assert result.overlay_png is not None and result.overlay_png.exists()
    assert result.geojson is not None and result.geojson.exists()
    assert result.regions_csv is not None and result.regions_csv.exists()

    table = pd.read_csv(result.regions_csv)
    assert table.loc[0, "region_name"] == "tissue_foreground"
    assert int(table.loc[0, "pixel_count"]) > 1000


class FakePromptBackend:
    def __init__(self) -> None:
        self.calls = []

    def segment(self, *, image, boxes, points, point_labels, tissue_mask, point_radius_px):
        self.calls.append(
            {
                "boxes": boxes,
                "points": points,
                "point_labels": point_labels,
                "point_radius_px": point_radius_px,
            }
        )
        mask = np.zeros(tissue_mask.shape, dtype=bool)
        for x0, y0, x1, y1 in boxes:
            mask[int(y0) : int(y1), int(x0) : int(x1)] = True
        return mask


def test_prompt_backend_receives_box_and_exports_prompted_region(tmp_path):
    image_path = _write_synthetic_he_image(tmp_path / "he.png")
    backend = FakePromptBackend()

    result = run_he_segmentation(
        HESegmentationConfig(
            image=image_path,
            out_dir=tmp_path / "prompt",
            task="single",
            backend=backend,
            region_specs=[
                HERegionSpec(
                    region_name="user_defined_region",
                    boxes=[[54, 42, 94, 84]],
                    points=[[70, 60]],
                )
            ],
            min_region_area_px=20,
        )
    )

    assert result.n_regions == 1
    assert backend.calls and backend.calls[0]["boxes"] == [[54.0, 42.0, 94.0, 84.0]]
    table = pd.read_csv(result.regions_csv)
    assert table.loc[0, "region_name"] == "user_defined_region"


def test_all_elements_creates_neutral_component_outputs(tmp_path):
    image_path = _write_synthetic_he_image(tmp_path / "he.png")
    result = run_he_segmentation(
        HESegmentationConfig(
            image=image_path,
            out_dir=tmp_path / "all_elements",
            task="all_elements",
            backend="heuristic",
            n_components=3,
            slic_segments=80,
            min_region_area_px=20,
        )
    )

    assert result.n_regions >= 2
    table = pd.read_csv(result.regions_csv)
    assert table["region_name"].str.startswith("component_").all()
    payload = json.loads(result.geojson.read_text(encoding="utf-8"))
    assert payload["type"] == "FeatureCollection"
    assert len(payload["features"]) >= 2


def test_he_change_detection_exports_changed_region(tmp_path):
    before = _write_synthetic_he_image(tmp_path / "before.png")
    after = _write_synthetic_he_image(tmp_path / "after.png", changed=True)

    result = run_he_change_detection(
        HEChangeDetectionConfig(
            before_image=before,
            after_image=after,
            out_dir=tmp_path / "change",
            change_quantile=0.80,
            min_change_area_px=20,
        )
    )

    assert result.n_regions >= 1
    assert result.heatmap_png.exists()
    assert result.change_mask_png.exists()
    assert result.geojson is not None and result.geojson.exists()

    table = pd.read_csv(result.changes_csv)
    assert int(table["pixel_count"].sum()) >= 20


def test_he_change_detection_rejects_mismatched_sizes(tmp_path):
    before = _write_synthetic_he_image(tmp_path / "before.png")
    after = _write_synthetic_he_image(tmp_path / "after.png", shape=(96, 120))

    with pytest.raises(ValueError, match="same shape"):
        run_he_change_detection(
            HEChangeDetectionConfig(
                before_image=before,
                after_image=after,
                out_dir=tmp_path / "bad_change",
            )
        )


def test_he_cli_single_smoke(tmp_path, capsys):
    image_path = _write_synthetic_he_image(tmp_path / "he.png")
    out_dir = tmp_path / "cli"

    he_cli_main(
        [
            "single",
            "--image",
            str(image_path),
            "--out-dir",
            str(out_dir),
            "--backend",
            "heuristic",
            "--min-region-area-px",
            "20",
        ]
    )

    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert payload["n_regions"] == 1
    assert (out_dir / "he_single_label_map.png").exists()
