from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

from shapely import affinity
from shapely.geometry import Polygon, mapping, shape
from shapely.ops import unary_union

from histoseg.threed import (
    LabelFreeContourAlignmentConfig,
    align_contours_label_free,
)
from histoseg.threed.cli import main
from histoseg.threed.label_free_alignment import _topology_similarity


def test_label_free_alignment_ignores_permuted_structure_labels(tmp_path):
    fixed = tmp_path / "fixed.geojson"
    moving = tmp_path / "moving.geojson"
    fixed_geoms = [
        Polygon([(0, 0), (80, 0), (80, 30), (0, 30)]),
        Polygon([(120, 20), (160, 20), (160, 70), (120, 70)]),
    ]
    moving_geoms = [affinity.translate(geom, xoff=25, yoff=-15) for geom in fixed_geoms]
    _write_geojson(fixed, fixed_geoms, ["Structure 1", "Structure 2"])
    _write_geojson(moving, moving_geoms, ["Unrelated A", "Unrelated B"])

    result = align_contours_label_free(
        LabelFreeContourAlignmentConfig(
            fixed_geojson=fixed,
            moving_geojson=moving,
            out_dir=tmp_path / "aligned",
            maxiter=0,
            run_soft_tps=False,
            boundary_sample_count=200,
        )
    )

    payload = json.loads(result.aligned_geojson.read_text(encoding="utf-8"))
    aligned = unary_union([shape(feature["geometry"]) for feature in payload["features"]])
    expected = unary_union(fixed_geoms)
    assert _iou(expected, aligned) > 0.95
    assert payload["features"][0]["properties"]["assigned_structure"] == "Unrelated A"
    assert payload["features"][1]["properties"]["assigned_structure"] == "Unrelated B"
    summary = json.loads(result.summary_json.read_text(encoding="utf-8"))
    assert summary["labels_used_for_matching"] is False
    assert summary["semantic_harmonization_performed"] is False


def test_label_free_alignment_handles_no_shared_structure_names(tmp_path):
    fixed = tmp_path / "fixed.geojson"
    moving = tmp_path / "moving.geojson"
    fixed_geoms = [
        Polygon([(0, 0), (50, 0), (50, 20), (0, 20)]),
        Polygon([(20, 55), (75, 55), (75, 90), (20, 90)]),
        Polygon([(95, 5), (130, 5), (130, 60), (95, 60)]),
    ]
    moving_geoms = [affinity.translate(geom, xoff=-40, yoff=30) for geom in fixed_geoms]
    _write_geojson(fixed, fixed_geoms, ["A", "B", "C"])
    _write_geojson(moving, moving_geoms, ["X", "Y", "Z"])

    result = align_contours_label_free(
        LabelFreeContourAlignmentConfig(
            fixed_geojson=fixed,
            moving_geojson=moving,
            out_dir=tmp_path / "aligned",
            maxiter=0,
            run_soft_tps=False,
            boundary_sample_count=240,
        )
    )
    summary = json.loads(result.summary_json.read_text(encoding="utf-8"))

    assert summary["scores"]["after_final"]["iou"] > 0.95
    assert result.landmarks_csv.exists()
    assert result.overlay_html is not None
    assert result.overlay_html.exists()
    assert "Plotly.newPlot" in result.overlay_html.read_text(encoding="utf-8")
    assert result.component_qc_csv is not None
    assert result.component_qc_csv.exists()
    summary = json.loads(result.summary_json.read_text(encoding="utf-8"))
    assert summary["correspondence_diagnostics"]["warning"] is None


def test_label_free_component_weight_cap_limits_dominant_polygon(tmp_path):
    fixed = tmp_path / "fixed.geojson"
    moving = tmp_path / "moving.geojson"
    fixed_geoms = [
        Polygon([(0, 0), (1000, 0), (1000, 1000), (0, 1000)]),
        Polygon([(100, 100), (130, 100), (130, 130), (100, 130)]),
        Polygon([(700, 700), (730, 700), (730, 730), (700, 730)]),
    ]
    moving_geoms = [affinity.translate(geom, xoff=10, yoff=20) for geom in fixed_geoms]
    _write_geojson(fixed, fixed_geoms, ["large", "small1", "small2"])
    _write_geojson(moving, moving_geoms, ["other_large", "other_small1", "other_small2"])

    result = align_contours_label_free(
        LabelFreeContourAlignmentConfig(
            fixed_geojson=fixed,
            moving_geojson=moving,
            out_dir=tmp_path / "aligned",
            maxiter=0,
            run_soft_tps=False,
            max_component_weight=0.2,
            boundary_sample_count=200,
        )
    )
    qc = json.loads(result.summary_json.read_text(encoding="utf-8"))["component_qc"]

    assert qc["max_component_weight"] <= 0.2


def test_label_free_soft_tps_writes_landmarks_and_preserves_valid_geojson(tmp_path):
    fixed = tmp_path / "fixed.geojson"
    moving = tmp_path / "moving.geojson"
    fixed_geoms = [Polygon([(0, 0), (120, 0), (120, 80), (0, 80)])]
    moving_geoms = [Polygon([(4, 2), (115, 0), (124, 84), (2, 78)])]
    _write_geojson(fixed, fixed_geoms, ["A"])
    _write_geojson(moving, moving_geoms, ["B"])

    result = align_contours_label_free(
        LabelFreeContourAlignmentConfig(
            fixed_geojson=fixed,
            moving_geojson=moving,
            out_dir=tmp_path / "aligned",
            maxiter=0,
            run_soft_tps=True,
            sampling_distance_um=20,
            max_landmark_distance_um=30,
            topology_grid_size=0,
            boundary_sample_count=160,
        )
    )
    landmarks = result.landmarks_csv.read_text(encoding="utf-8")
    payload = json.loads(result.aligned_geojson.read_text(encoding="utf-8"))

    assert "label_free_mutual_kdtree" in landmarks
    assert shape(payload["features"][0]["geometry"]).is_valid


def test_label_free_alignment_cli_smoke(tmp_path):
    fixed = tmp_path / "fixed.geojson"
    moving = tmp_path / "moving.geojson"
    square = Polygon([(0, 0), (40, 0), (40, 40), (0, 40)])
    _write_geojson(fixed, [square], ["A"])
    _write_geojson(moving, [affinity.translate(square, xoff=5, yoff=3)], ["B"])

    main(
        [
            "align-contours-label-free",
            "--fixed-geojson",
            str(fixed),
            "--moving-geojson",
            str(moving),
            "--out-dir",
            str(tmp_path / "aligned"),
            "--maxiter",
            "0",
            "--no-soft-tps",
            "--no-preview",
        ]
    )

    assert (tmp_path / "aligned" / "moving_label_free_aligned.geojson").exists()
    assert (tmp_path / "aligned" / "label_free_alignment_summary.json").exists()
    assert (tmp_path / "aligned" / "label_free_alignment_overlay.html").exists()


def test_label_free_alignment_warns_when_internal_contours_are_not_homologous(tmp_path):
    fixed = tmp_path / "fixed.geojson"
    moving = tmp_path / "moving.geojson"
    envelope = Polygon([(0, 0), (400, 0), (400, 300), (0, 300)])
    fixed_geoms = [
        envelope,
        Polygon([(30, 30), (100, 30), (100, 100), (30, 100)]),
        Polygon([(280, 190), (360, 190), (360, 260), (280, 260)]),
    ]
    moving_geoms = [
        envelope,
        Polygon([(260, 35), (360, 35), (360, 100), (260, 100)]),
        Polygon([(35, 210), (120, 210), (120, 275), (35, 275)]),
    ]
    _write_geojson(fixed, fixed_geoms, ["Envelope", "A", "B"])
    _write_geojson(moving, moving_geoms, ["Envelope other", "X", "Y"])

    result = align_contours_label_free(
        LabelFreeContourAlignmentConfig(
            fixed_geojson=fixed,
            moving_geojson=moving,
            out_dir=tmp_path / "aligned",
            maxiter=0,
            run_soft_tps=False,
            boundary_sample_count=200,
        )
    )
    summary = json.loads(result.summary_json.read_text(encoding="utf-8"))
    html = result.overlay_html.read_text(encoding="utf-8") if result.overlay_html else ""

    assert summary["correspondence_diagnostics"]["warning"]
    assert "Alignment warning" in html


def test_partial_correspondence_permuted_labels_yields_anchors(tmp_path):
    fixed = tmp_path / "fixed.geojson"
    moving = tmp_path / "moving.geojson"
    geoms = _constellation_geometries()
    _write_geojson(fixed, geoms, ["A", "B", "C", "D", "E", "F"])
    _write_geojson(moving, geoms, ["X", "Y", "Z", "U", "V", "W"])

    result = align_contours_label_free(
        LabelFreeContourAlignmentConfig(
            fixed_geojson=fixed,
            moving_geojson=moving,
            out_dir=tmp_path / "partial",
            partial_correspondence=True,
            diagnostic_only=True,
            knn_neighbors=3,
            min_anchor_count=3,
        )
    )
    summary = json.loads(result.summary_json.read_text(encoding="utf-8"))
    nodes = result.partial_nodes_csv.read_text(encoding="utf-8")
    matches = result.partial_matches_csv.read_text(encoding="utf-8")

    assert result.aligned_geojson is None
    assert not (tmp_path / "partial" / "moving_label_free_aligned.geojson").exists()
    assert summary["coordinate_warp_performed"] is False
    assert summary["anchor_count"] >= 4
    assert "matched_anchor" in nodes
    assert "matched_anchor" in matches
    assert result.partial_matrix_html is not None
    assert "Partial correspondence candidate matrix" in result.partial_matrix_html.read_text(encoding="utf-8")


def test_partial_correspondence_missing_and_extra_components_are_no_counterpart(tmp_path):
    fixed = tmp_path / "fixed.geojson"
    moving = tmp_path / "moving.geojson"
    fixed_geoms = _constellation_geometries()
    moving_geoms = fixed_geoms[:4] + [
        Polygon([(700, 700), (740, 700), (740, 740), (700, 740)]),
        Polygon([(820, 40), (860, 40), (860, 80), (820, 80)]),
    ]
    _write_geojson(fixed, fixed_geoms, ["A", "B", "C", "D", "E", "F"])
    _write_geojson(moving, moving_geoms, ["A2", "B2", "C2", "D2", "extra1", "extra2"])

    result = align_contours_label_free(
        LabelFreeContourAlignmentConfig(
            fixed_geojson=fixed,
            moving_geojson=moving,
            out_dir=tmp_path / "partial",
            partial_correspondence=True,
            diagnostic_only=True,
            knn_neighbors=3,
            search_window=120,
            min_anchor_count=3,
        )
    )
    summary = json.loads(result.summary_json.read_text(encoding="utf-8"))
    nodes_text = result.partial_nodes_csv.read_text(encoding="utf-8")

    assert summary["anchor_count"] >= 3
    assert "no_counterpart" in nodes_text


def test_partial_anchor_alignment_uses_only_matched_anchors_for_transform(tmp_path):
    fixed = tmp_path / "fixed.geojson"
    moving = tmp_path / "moving.geojson"
    fixed_geoms = _constellation_geometries()
    translated = [affinity.translate(geom, xoff=32, yoff=-18) for geom in fixed_geoms]
    extra = Polygon([(900, 900), (940, 900), (940, 940), (900, 940)])
    moving_geoms = [*translated, extra]
    _write_geojson(fixed, fixed_geoms, ["A", "B", "C", "D", "E", "F"])
    _write_geojson(moving, moving_geoms, ["X", "Y", "Z", "U", "V", "W", "extra"])

    result = align_contours_label_free(
        LabelFreeContourAlignmentConfig(
            fixed_geojson=fixed,
            moving_geojson=moving,
            out_dir=tmp_path / "partial_anchor",
            partial_correspondence=True,
            diagnostic_only=False,
            knn_neighbors=3,
            min_anchor_count=3,
        )
    )

    assert result.aligned_geojson is not None
    assert result.overlay_html is not None
    payload = json.loads(result.aligned_geojson.read_text(encoding="utf-8"))
    aligned_first_six = unary_union(
        [shape(feature["geometry"]) for feature in payload["features"][:6]]
    )
    assert _iou(unary_union(fixed_geoms), aligned_first_six) > 0.95
    summary = json.loads(result.summary_json.read_text(encoding="utf-8"))
    assert summary["method"] == "label_free_partial_anchor_similarity_alignment"
    assert summary["coordinate_warp_performed"] is True
    assert summary["anchor_transform"]["used_anchor_pair_count"] >= 3
    anchor_rows = result.landmarks_csv.read_text(encoding="utf-8")
    assert "used_for_transform" in anchor_rows


def test_overlap_ransac_recovers_large_partial_overlap_offset(tmp_path):
    fixed = tmp_path / "fixed.geojson"
    moving = tmp_path / "moving.geojson"
    fixed_geoms = _constellation_geometries()
    moving_overlap = [
        affinity.translate(geom, xoff=1800, yoff=-950)
        for geom in fixed_geoms
    ]
    moving_geoms = [
        *moving_overlap,
        Polygon([(2600, 1400), (2650, 1400), (2650, 1450), (2600, 1450)]),
    ]
    _write_geojson(fixed, fixed_geoms, ["A", "B", "C", "D", "E", "F"])
    _write_geojson(moving, moving_geoms, ["X", "Y", "Z", "U", "V", "W", "extra"])

    result = align_contours_label_free(
        LabelFreeContourAlignmentConfig(
            fixed_geojson=fixed,
            moving_geojson=moving,
            out_dir=tmp_path / "overlap",
            partial_correspondence=True,
            diagnostic_only=False,
            overlap_ransac=True,
            overlap_candidate_count=4,
            overlap_ransac_trials=2000,
            knn_neighbors=3,
            search_window=300,
            min_anchor_count=3,
        )
    )

    assert result.aligned_geojson is not None
    payload = json.loads(result.aligned_geojson.read_text(encoding="utf-8"))
    aligned_first_six = unary_union(
        [shape(feature["geometry"]) for feature in payload["features"][:6]]
    )
    assert _iou(unary_union(fixed_geoms), aligned_first_six) > 0.9
    summary = json.loads(result.summary_json.read_text(encoding="utf-8"))
    assert summary["method"] == "label_free_overlap_ransac_similarity_alignment"
    assert summary["anchor_transform"]["used_anchor_pair_count"] >= 3


def test_group_correspondence_matches_cross_named_contour_groups(tmp_path):
    fixed = tmp_path / "fixed.geojson"
    moving = tmp_path / "moving.geojson"
    fixed_group = _constellation_geometries()
    fixed_distractor = [
        affinity.translate(geom, xoff=2600, yoff=1000)
        for geom in _constellation_geometries()[:4]
    ]
    moving_group = [
        affinity.translate(geom, xoff=1500, yoff=-850)
        for geom in fixed_group
    ]
    moving_distractor = [
        affinity.translate(geom, xoff=-2200, yoff=1800)
        for geom in _constellation_geometries()[:4]
    ]
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

    result = align_contours_label_free(
        LabelFreeContourAlignmentConfig(
            fixed_geojson=fixed,
            moving_geojson=moving,
            out_dir=tmp_path / "group",
            partial_correspondence=True,
            diagnostic_only=False,
            group_correspondence=True,
            group_candidate_count=4,
            group_ransac_trials=2000,
            group_min_component_area_um2=0.0,
            knn_neighbors=3,
            search_window=300,
            min_anchor_count=3,
        )
    )

    assert result.aligned_geojson is not None
    assert result.group_matrix_csv is not None
    matrix = result.group_matrix_csv.read_text(encoding="utf-8")
    assert "Structure 2,Structure 3" in matrix
    summary = json.loads(result.summary_json.read_text(encoding="utf-8"))
    assert summary["method"] == "label_free_group_correspondence_ransac_alignment"
    assert summary["anchor_transform"]["fixed_group"] == "Structure 2"
    assert summary["anchor_transform"]["moving_group"] == "Structure 3"
    assert summary["anchor_transform"]["used_anchor_pair_count"] >= 3
    payload = json.loads(result.aligned_geojson.read_text(encoding="utf-8"))
    aligned_group = unary_union(
        [shape(feature["geometry"]) for feature in payload["features"][: len(fixed_group)]]
    )
    assert _iou(unary_union(fixed_group), aligned_group) > 0.9


def test_partial_correspondence_excludes_dominant_envelope(tmp_path):
    fixed = tmp_path / "fixed.geojson"
    moving = tmp_path / "moving.geojson"
    envelope = Polygon([(0, 0), (1000, 0), (1000, 800), (0, 800)])
    geoms = [envelope, *_constellation_geometries()]
    _write_geojson(fixed, geoms, ["envelope", "A", "B", "C", "D", "E", "F"])
    _write_geojson(moving, geoms, ["other envelope", "X", "Y", "Z", "U", "V", "W"])

    result = align_contours_label_free(
        LabelFreeContourAlignmentConfig(
            fixed_geojson=fixed,
            moving_geojson=moving,
            out_dir=tmp_path / "partial",
            partial_correspondence=True,
            diagnostic_only=True,
            knn_neighbors=3,
            min_anchor_count=3,
        )
    )
    nodes = result.partial_nodes_csv.read_text(encoding="utf-8")
    matches = result.partial_matches_csv.read_text(encoding="utf-8")

    assert "envelope_only" in nodes
    assert "envelope" not in matches


def test_partial_correspondence_clustered_anchors_warn(tmp_path):
    fixed = tmp_path / "fixed.geojson"
    moving = tmp_path / "moving.geojson"
    envelope = Polygon([(-50, -50), (500, -50), (500, 500), (-50, 500)])
    geoms = [
        envelope,
        Polygon([(10, 10), (40, 10), (40, 40), (10, 40)]),
        Polygon([(70, 15), (105, 15), (105, 45), (70, 45)]),
        Polygon([(20, 80), (55, 80), (55, 115), (20, 115)]),
        Polygon([(90, 90), (120, 90), (120, 120), (90, 120)]),
    ]
    _write_geojson(fixed, geoms, ["envelope", "A", "B", "C", "D"])
    _write_geojson(moving, geoms, ["envelope", "X", "Y", "Z", "W"])

    result = align_contours_label_free(
        LabelFreeContourAlignmentConfig(
            fixed_geojson=fixed,
            moving_geojson=moving,
            out_dir=tmp_path / "partial",
            partial_correspondence=True,
            diagnostic_only=True,
            knn_neighbors=3,
            min_anchor_count=3,
        )
    )
    summary = json.loads(result.summary_json.read_text(encoding="utf-8"))

    assert "low_anchor_coverage" in summary["warnings"] or "clustered_anchors" in summary["warnings"]


def test_partial_correspondence_topology_similarity_is_circular_shift_stable():
    profile = SimpleNamespace(
        knn_distances=(0.45, 0.8, 1.2, 1.55),
        knn_angle_gaps=(0.15, 0.2, 0.35, 0.3),
        knn_log_area_ratios=(-0.2, 0.0, 0.15, 0.4),
    )
    shifted = SimpleNamespace(
        knn_distances=(1.2, 1.55, 0.45, 0.8),
        knn_angle_gaps=(0.35, 0.3, 0.15, 0.2),
        knn_log_area_ratios=(0.15, 0.4, -0.2, 0.0),
    )
    mismatched = SimpleNamespace(
        knn_distances=(0.1, 2.5, 0.2, 3.0),
        knn_angle_gaps=(0.75, 0.05, 0.15, 0.05),
        knn_log_area_ratios=(1.5, -1.3, 1.2, -1.1),
    )

    assert _topology_similarity(profile, shifted) > 0.99
    assert _topology_similarity(profile, mismatched) < 0.25


def test_partial_correspondence_cli_smoke(tmp_path):
    fixed = tmp_path / "fixed.geojson"
    moving = tmp_path / "moving.geojson"
    geoms = _constellation_geometries()
    _write_geojson(fixed, geoms, ["A", "B", "C", "D", "E", "F"])
    _write_geojson(moving, geoms, ["X", "Y", "Z", "U", "V", "W"])

    main(
        [
            "align-contours-label-free",
            "--fixed-geojson",
            str(fixed),
            "--moving-geojson",
            str(moving),
            "--out-dir",
            str(tmp_path / "partial"),
            "--partial-correspondence",
            "--diagnostic-only",
            "--knn-neighbors",
            "3",
            "--min-anchor-count",
            "3",
        ]
    )

    assert (tmp_path / "partial" / "partial_correspondence_nodes.csv").exists()
    assert (tmp_path / "partial" / "partial_correspondence_matches.csv").exists()
    assert (tmp_path / "partial" / "partial_correspondence_summary.json").exists()
    assert (tmp_path / "partial" / "partial_correspondence_matrix.html").exists()


def _write_geojson(path: Path, geometries: list[Polygon], labels: list[str]) -> None:
    payload = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "properties": {
                    "assigned_structure": label,
                    "structure_id": index + 1,
                    "name": label,
                },
                "geometry": mapping(geom),
            }
            for index, (geom, label) in enumerate(zip(geometries, labels))
        ],
    }
    path.write_text(json.dumps(payload), encoding="utf-8")


def _iou(a, b) -> float:
    return float(a.intersection(b).area / a.union(b).area)


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
