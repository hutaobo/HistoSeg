from __future__ import annotations

import json
from pathlib import Path

import pytest

import histoseg.threed.multislice as multislice


def test_coda_image_tournament_promotes_higher_iou_seed(monkeypatch, tmp_path):
    fixed = tmp_path / "fixed.geojson"
    moving = tmp_path / "moving.geojson"
    fixed.write_text("{}", encoding="utf-8")
    moving.write_text("{}", encoding="utf-8")
    out = tmp_path / "moving_hard_aligned.geojson"
    summary_json = tmp_path / "hard_similarity_alignment.json"

    def fake_contour_align(**kwargs):
        Path(kwargs["output_geojson"]).write_text("contour", encoding="utf-8")
        return _summary(
            backend="contour-tps",
            output_geojson=kwargs["output_geojson"],
            iou_after=0.81,
            rotation_degrees=10.0,
        )

    def fake_coda_align(**kwargs):
        Path(kwargs["output_geojson"]).write_text("coda", encoding="utf-8")
        return _summary(
            backend="coda-image",
            output_geojson=kwargs["output_geojson"],
            iou_after=0.92,
            rotation_degrees=120.0,
            coda_image={
                "radon_rotation_degrees": 120.0,
                "phase_shift_y": -3.0,
                "phase_shift_x": 4.0,
                "orientation_disambiguation": "selected_original_orientation",
            },
        )

    monkeypatch.setattr(multislice, "hard_align_geojson", fake_contour_align)
    monkeypatch.setattr(multislice, "_coda_image_align_geojson", fake_coda_align)

    summary = multislice._run_tournament_hard_alignment(
        fixed_geojson=fixed,
        moving_geojson=moving,
        output_geojson=out,
        summary_json=summary_json,
        overwrite=True,
    )

    assert out.read_text(encoding="utf-8") == "coda"
    assert summary_json.exists()
    assert json.loads(summary_json.read_text(encoding="utf-8")) == summary
    assert summary["registration_backend"] == "coda-image"
    assert summary["selected_hard_seed_backend"] == "coda-image"
    assert summary["union_iou_after_hard"] == pytest.approx(0.92)
    assert summary["hard_alignment_tournament"]["selected_backend"] == "coda-image"
    assert summary["hard_alignment_tournament"]["rotation_difference_degrees"] == pytest.approx(
        110.0
    )
    assert [item["backend"] for item in summary["hard_alignment_candidates"]] == [
        "contour-tps",
        "coda-image",
    ]
    assert summary["hard_alignment_candidates"][0]["union_iou_after_hard"] == pytest.approx(
        0.81
    )
    assert summary["hard_alignment_candidates"][1]["union_iou_after_hard"] == pytest.approx(
        0.92
    )
    assert summary["coda_image"]["radon_rotation_degrees"] == pytest.approx(120.0)


def test_coda_image_tournament_keeps_contour_seed_on_tie(monkeypatch, tmp_path):
    fixed = tmp_path / "fixed.geojson"
    moving = tmp_path / "moving.geojson"
    fixed.write_text("{}", encoding="utf-8")
    moving.write_text("{}", encoding="utf-8")
    out = tmp_path / "moving_hard_aligned.geojson"

    def fake_contour_align(**kwargs):
        Path(kwargs["output_geojson"]).write_text("contour", encoding="utf-8")
        return _summary(
            backend="contour-tps",
            output_geojson=kwargs["output_geojson"],
            iou_after=0.88,
            rotation_degrees=0.0,
        )

    def fake_coda_align(**kwargs):
        Path(kwargs["output_geojson"]).write_text("coda", encoding="utf-8")
        return _summary(
            backend="coda-image",
            output_geojson=kwargs["output_geojson"],
            iou_after=0.88,
            rotation_degrees=180.0,
            coda_image={"radon_rotation_degrees": 180.0},
        )

    monkeypatch.setattr(multislice, "hard_align_geojson", fake_contour_align)
    monkeypatch.setattr(multislice, "_coda_image_align_geojson", fake_coda_align)

    summary = multislice._run_tournament_hard_alignment(
        fixed_geojson=fixed,
        moving_geojson=moving,
        output_geojson=out,
        summary_json=None,
        overwrite=True,
    )

    assert out.read_text(encoding="utf-8") == "contour"
    assert summary["registration_backend"] == "coda-image"
    assert summary["selected_hard_seed_backend"] == "contour-tps"
    assert summary["hard_alignment_tournament"]["selected_backend"] == "contour-tps"
    assert summary["transform"]["rotation_degrees"] == pytest.approx(0.0)


def _summary(
    *,
    backend: str,
    output_geojson: Path,
    iou_after: float,
    rotation_degrees: float,
    coda_image: dict | None = None,
) -> dict:
    payload = {
        "fixed_geojson": "fixed.geojson",
        "moving_geojson": "moving.geojson",
        "output_geojson": str(output_geojson),
        "registration_backend": backend,
        "transform": {
            "origin_x": 0.0,
            "origin_y": 0.0,
            "rotation_degrees": rotation_degrees,
            "scale": 1.0,
            "translate_x": 2.0,
            "translate_y": -1.0,
        },
        "affine_params": None,
        "optimization": {"accepted": True},
        "union_iou_before_hard": 0.25,
        "union_iou_after_hard": iou_after,
        "hard_alignment_accepted": True,
        "per_structure_iou_after_hard": {"Structure 1": iou_after},
    }
    if coda_image is not None:
        payload["coda_image"] = coda_image
    return payload
