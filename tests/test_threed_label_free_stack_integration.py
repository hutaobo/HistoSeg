from __future__ import annotations

import json

from shapely import affinity
from shapely.geometry import Polygon, mapping

from histoseg.threed import ThreeDStackReconstructionConfig
from histoseg.threed.multislice import (
    _SliceInput,
    _label_free_group_candidate_guard,
    _pairwise_row,
    _run_hard_alignment_backend,
    _resolve_soft_alignment_mode,
    _semantic_soft_alignment_policy,
)


def test_auto_backend_selects_cross_named_label_free_group(tmp_path):
    fixed = tmp_path / "fixed.geojson"
    moving = tmp_path / "moving.geojson"
    fixed_group = _constellation_geometries()
    fixed_distractor = [affinity.translate(geom, xoff=2600, yoff=1000) for geom in fixed_group[:4]]
    moving_group = [affinity.translate(geom, xoff=1500, yoff=-850) for geom in fixed_group]
    moving_distractor = [affinity.translate(geom, xoff=-2200, yoff=1800) for geom in fixed_group[:4]]
    _write_geojson(
        fixed,
        [*fixed_group, *fixed_distractor],
        ["Structure 2"] * len(fixed_group) + ["Structure 1"] * len(fixed_distractor),
    )
    _write_geojson(
        moving,
        [*moving_group, *moving_distractor],
        ["Structure 3"] * len(moving_group) + ["Structure 4"] * len(moving_distractor),
    )

    summary = _run_hard_alignment_backend(
        _stack_cfg(
            tmp_path,
            registration_backend="auto",
            label_free_group_min_component_area_um2=0.0,
        ),
        fixed_geojson=fixed,
        moving_geojson=moving,
        output_geojson=tmp_path / "moving_hard_aligned.geojson",
        summary_json=tmp_path / "hard_summary.json",
    )

    assert summary["registration_backend"] == "auto"
    assert summary["selected_hard_seed_backend"] == "label-free-group"
    assert summary["label_free_fixed_group"] == "Structure 2"
    assert summary["label_free_moving_group"] == "Structure 3"
    assert summary["label_free_used_anchor_pair_count"] >= 3
    assert (tmp_path / "moving_hard_aligned.geojson").exists()
    allowed, reason = _semantic_soft_alignment_policy(summary)
    assert allowed is False
    assert reason == "cross_named_label_free_group_match"


def test_auto_backend_falls_back_when_label_free_has_too_few_anchors(tmp_path):
    fixed = tmp_path / "fixed.geojson"
    moving = tmp_path / "moving.geojson"
    geoms = _constellation_geometries()
    _write_geojson(fixed, geoms, ["Structure 2"] * len(geoms))
    _write_geojson(
        moving,
        [affinity.translate(geom, xoff=1500, yoff=-850) for geom in geoms],
        ["Structure 3"] * len(geoms),
    )

    summary = _run_hard_alignment_backend(
        _stack_cfg(
            tmp_path,
            registration_backend="auto",
            label_free_min_anchor_count=99,
            label_free_group_min_component_area_um2=0.0,
        ),
        fixed_geojson=fixed,
        moving_geojson=moving,
        output_geojson=tmp_path / "moving_hard_aligned.geojson",
        summary_json=tmp_path / "hard_summary.json",
    )

    assert summary["registration_backend"] == "auto"
    assert summary["selected_hard_seed_backend"] == "contour-tps"
    assert summary["label_free_fixed_group"] is None
    assert summary["label_free_moving_group"] is None
    assert summary["label_free_guarded_candidate_accepted"] is False


def test_guard_rejects_two_anchor_similarity_rotation(tmp_path):
    ok, guard = _label_free_group_candidate_guard(
        _label_free_summary_for_guard(anchor_count=2, rotation_degrees=32.0),
        _stack_cfg(tmp_path, label_free_min_anchor_count=2),
    )

    assert ok is False
    assert guard["reason"] == "too_few_anchors_for_similarity_transform"


def test_guard_rejects_low_anchor_unconstrained_similarity(tmp_path):
    ok, guard = _label_free_group_candidate_guard(
        _label_free_summary_for_guard(anchor_count=5, rotation_degrees=3.0),
        _stack_cfg(tmp_path, label_free_min_anchor_count=4),
    )

    assert ok is False
    assert guard["reason"] == "too_few_anchors_for_unconstrained_similarity"


def test_guard_rejects_high_rotation_without_strong_evidence(tmp_path):
    ok, guard = _label_free_group_candidate_guard(
        _label_free_summary_for_guard(
            anchor_count=7,
            rotation_degrees=48.0,
            residual_median=20.0,
            coverage=0.6,
            quadrants=4,
        ),
        _stack_cfg(tmp_path, label_free_min_anchor_count=4),
    )

    assert ok is False
    assert guard["reason"] == "rotation_exceeds_prior_without_enough_anchors"


def test_guard_accepts_strong_low_rotation_label_free_candidate(tmp_path):
    ok, guard = _label_free_group_candidate_guard(
        _label_free_summary_for_guard(
            anchor_count=8,
            rotation_degrees=6.0,
            residual_median=35.0,
            coverage=0.3,
            quadrants=3,
        ),
        _stack_cfg(tmp_path, label_free_min_anchor_count=4),
    )

    assert ok is True
    assert guard["accepted"] is True


def test_guard_rejects_near_180_rotation(tmp_path):
    ok, guard = _label_free_group_candidate_guard(
        _label_free_summary_for_guard(
            anchor_count=12,
            rotation_degrees=178.0,
            residual_median=10.0,
            coverage=0.8,
            quadrants=4,
        ),
        _stack_cfg(tmp_path, label_free_min_anchor_count=4),
    )

    assert ok is False
    assert guard["reason"] == "near_180_rotation_rejected"


def test_forced_label_free_backend_exposes_hard_summary_schema(tmp_path):
    fixed = tmp_path / "fixed.geojson"
    moving = tmp_path / "moving.geojson"
    geoms = _constellation_geometries()
    _write_geojson(fixed, geoms, ["Structure 2"] * len(geoms))
    _write_geojson(
        moving,
        [affinity.translate(geom, xoff=1500, yoff=-850) for geom in geoms],
        ["Structure 3"] * len(geoms),
    )

    summary = _run_hard_alignment_backend(
        _stack_cfg(
            tmp_path,
            registration_backend="label-free-group",
            label_free_group_min_component_area_um2=0.0,
        ),
        fixed_geojson=fixed,
        moving_geojson=moving,
        output_geojson=tmp_path / "moving_hard_aligned.geojson",
        summary_json=tmp_path / "hard_summary.json",
    )

    assert summary["registration_backend"] == "label-free-group"
    assert summary["selected_hard_seed_backend"] == "label-free-group"
    assert {"rotation_degrees", "scale", "translate_x", "translate_y"}.issubset(
        summary["transform"]
    )
    assert isinstance(summary["hard_alignment_accepted"], bool)
    assert "union_iou_before_hard" in summary
    assert "union_iou_after_hard" in summary
    assert summary["label_free_anchor_landmarks_csv"]
    assert summary["label_free_partial_matches_csv"]


def test_same_named_label_free_group_uses_anchor_only_soft_tps(tmp_path):
    summary = {
        "registration_backend": "auto",
        "selected_hard_seed_backend": "label-free-group",
        "label_free_fixed_group": "Structure 2",
        "label_free_moving_group": "Structure 2",
    }

    allowed, reason = _semantic_soft_alignment_policy(summary)
    assert allowed is False
    assert reason == "label_free_group_uses_anchor_only_soft_alignment"
    active, skipped = _resolve_soft_alignment_mode(
        _stack_cfg(tmp_path),
        summary,
        semantic_soft_allowed=allowed,
        semantic_soft_skipped_reason=reason,
    )
    assert active == "anchor-only"
    assert skipped is None


def test_pairwise_row_records_label_free_stack_metrics(tmp_path):
    hard_summary = {
        "registration_backend": "auto",
        "selected_hard_seed_backend": "label-free-group",
        "union_iou_before_hard": 0.8,
        "union_iou_after_hard": 0.4,
        "hard_alignment_candidates": [
            {"backend": "contour-tps", "union_iou_after_hard": 0.2},
            {
                "backend": "label-free-group",
                "label_free_residual_median": 51.7,
            },
        ],
        "hard_alignment_tournament": {},
        "transform": {
            "rotation_degrees": 1.0,
            "scale": 1.0,
            "translate_x": -10.0,
            "translate_y": 20.0,
        },
        "hard_alignment_accepted": True,
        "label_free_fixed_group": "Structure 2",
        "label_free_moving_group": "Structure 3",
        "label_free_used_anchor_pair_count": 22,
        "label_free_residual_median": 51.7,
        "semantic_soft_allowed": False,
        "semantic_soft_skipped_reason": "cross_named_label_free_group_match",
        "soft_alignment_mode_requested": "auto",
        "active_soft_alignment_mode": "anchor-only",
        "soft_alignment_skipped_reason": None,
    }

    row = _pairwise_row(
        _SliceInput(
            order=2,
            sample_id="slice_2",
            sample_dir=tmp_path / "sample",
            xenium_dir=tmp_path / "sample" / "xenium",
        ),
        hard_summary,
        None,
        None,
    )

    assert row["selected_hard_seed_backend"] == "label-free-group"
    assert row["label_free_fixed_group"] == "Structure 2"
    assert row["label_free_moving_group"] == "Structure 3"
    assert row["label_free_used_anchor_pair_count"] == 22
    assert row["hard_candidate_label_free_group_residual_median"] == 51.7
    assert row["semantic_soft_allowed"] is False
    assert row["semantic_soft_skipped_reason"] == "cross_named_label_free_group_match"
    assert row["soft_alignment_mode_requested"] == "auto"
    assert row["active_soft_alignment_mode"] == "anchor-only"


def test_pairwise_row_records_anchor_only_residual_tps_metrics(tmp_path):
    hard_summary = {
        "registration_backend": "auto",
        "selected_hard_seed_backend": "label-free-group",
        "union_iou_before_hard": 0.4,
        "union_iou_after_hard": 0.5,
        "hard_alignment_candidates": [],
        "hard_alignment_tournament": {},
        "transform": {
            "rotation_degrees": 0.0,
            "scale": 1.0,
            "translate_x": 0.0,
            "translate_y": 0.0,
        },
        "hard_alignment_accepted": True,
        "soft_alignment_mode_requested": "auto",
        "active_soft_alignment_mode": "anchor-only",
        "soft_alignment_runtime_seconds": 1.25,
    }
    summary_path = tmp_path / "anchor_only_tps_summary.json"
    summary_path.write_text("{}", encoding="utf-8")
    soft_summary = {
        "method": "anchor_only_residual_tps",
        "accepted": True,
        "landmarks": {
            "boundary_landmark_count": 6,
            "anchor_landmark_count": 6,
            "identity_padding_count": 16,
            "input_residual_um": {"median": 21.0, "p90": 42.0},
        },
        "qc": {
            "union_iou_hard_before_soft": 0.5,
            "union_iou_soft_after": 0.48,
            "geometry_status_counts": {"valid": 3},
            "topology_check": {
                "valid": True,
                "checked_cells": 2401,
                "folded_cell_count": 0,
                "negative_jacobian_ratio": 0.0,
                "min_jacobian_ratio": 0.94,
            },
            "jacobian_check": {
                "negative_jacobian_ratio": 0.0,
                "min_jacobian_ratio": 0.94,
            },
            "post_warp_residual_um": {"median": 1.0, "p90": 2.0},
        },
    }
    soft_result = type("SoftResult", (), {"summary_json": summary_path})()

    row = _pairwise_row(
        _SliceInput(
            order=2,
            sample_id="slice_2",
            sample_dir=tmp_path / "sample",
            xenium_dir=tmp_path / "sample" / "xenium",
        ),
        hard_summary,
        soft_summary,
        soft_result,
        soft_accepted=True,
    )

    assert row["active_soft_alignment_mode"] == "anchor-only"
    assert row["soft_accepted"] is True
    assert row["soft_union_iou_after"] == 0.48
    assert row["anchor_only_anchor_count"] == 6
    assert row["anchor_only_identity_padding_count"] == 16
    assert row["anchor_only_input_residual_p90"] == 42.0
    assert row["anchor_only_post_residual_p90"] == 2.0
    assert row["anchor_only_negative_jacobian_ratio"] == 0.0
    assert row["anchor_only_min_jacobian_ratio"] == 0.94


def _stack_cfg(tmp_path, **overrides):
    values = {
        "xenium_root": tmp_path,
        "label_free_knn_neighbors": 3,
        "label_free_min_anchor_count": 3,
        "label_free_group_candidate_count": 4,
        "label_free_group_ransac_trials": 2000,
        "label_free_group_residual_limit_um": 300.0,
        "hard_alignment_maxiter": 0,
        "save_alignment_preview_png": False,
        "overwrite": True,
    }
    values.update(overrides)
    return ThreeDStackReconstructionConfig(**values)


def _label_free_summary_for_guard(
    *,
    anchor_count: int,
    rotation_degrees: float,
    residual_median: float = 25.0,
    coverage: float = 0.25,
    quadrants: int = 3,
):
    return {
        "registration_backend": "label-free-group",
        "hard_alignment_accepted": True,
        "transform": {
            "kind": "similarity",
            "rotation_degrees": rotation_degrees,
            "scale": 1.0,
            "translate_x": 100.0,
            "translate_y": -50.0,
        },
        "label_free_used_anchor_pair_count": anchor_count,
        "label_free_residual_median": residual_median,
        "label_free_anchor_coverage_ratio": coverage,
        "label_free_anchor_occupied_quadrants": quadrants,
    }


def _write_geojson(path, geometries, labels) -> None:
    payload = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "properties": {
                    "assigned_structure": label,
                    "structure": label,
                    "structure_id": index + 1,
                    "name": label,
                },
                "geometry": mapping(geom),
            }
            for index, (geom, label) in enumerate(zip(geometries, labels))
        ],
    }
    path.write_text(json.dumps(payload), encoding="utf-8")


def _constellation_geometries() -> list[Polygon]:
    centers = [(40, 40), (140, 45), (80, 150), (230, 80), (250, 200), (360, 120)]
    sizes = [(30, 24), (36, 30), (28, 34), (40, 26), (32, 36), (30, 30)]
    geoms = []
    for (cx, cy), (w, h) in zip(centers, sizes):
        geoms.append(
            Polygon(
                [
                    (cx - w / 2, cy - h / 2),
                    (cx + w / 2, cy - h / 2),
                    (cx + w / 2, cy + h / 2),
                    (cx - w / 2, cy + h / 2),
                ]
            )
        )
    return geoms
