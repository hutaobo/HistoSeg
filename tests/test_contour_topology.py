from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from shapely.geometry import MultiPolygon, Polygon, box

from histoseg.contour import (
    ContourTopologyConfig,
    ContourTopologyResult,
    summarize_contour_topology,
)


def _contour_table(records: list[tuple[str, object, str]]) -> pd.DataFrame:
    return pd.DataFrame.from_records(
        [
            {"contour_id": contour_id, "geometry": geometry, "structure": structure}
            for contour_id, geometry, structure in records
        ]
    )


def test_shared_edge_reports_boundary_overlap_and_group_summary():
    contours = _contour_table(
        [
            ("left", box(0, 0, 10, 10), "A"),
            ("right", box(10, 0, 20, 10), "B"),
        ]
    )

    result = summarize_contour_topology(contours, groupby="structure", boundary_tolerance=0.0)

    assert isinstance(ContourTopologyConfig(), ContourTopologyConfig)
    assert isinstance(result, ContourTopologyResult)
    assert len(result.boundary_overlap) == 1
    row = result.boundary_overlap.iloc[0]
    assert row["contour_id_a"] == "left"
    assert row["contour_id_b"] == "right"
    assert row["exact_shared_boundary_length_um"] == pytest.approx(10.0)
    assert row["shared_boundary_length_um"] == pytest.approx(10.0)
    assert row["overlap_fraction_a"] == pytest.approx(0.25)
    assert row["overlap_fraction_b"] == pytest.approx(0.25)
    assert bool(row["is_boundary_neighbor"]) is True

    summary = result.contour_summary.set_index("contour_id")
    assert summary.loc["left", "n_boundary_neighbors"] == 1
    assert summary.loc["right", "total_shared_boundary_length_um"] == pytest.approx(10.0)
    assert len(result.group_boundary_overlap) == 1
    assert result.group_boundary_overlap.iloc[0]["shared_boundary_length_um"] == pytest.approx(10.0)


def test_boundary_tolerance_captures_nearby_edges_without_exact_touching():
    contours = _contour_table(
        [
            ("left", box(0, 0, 10, 10), "A"),
            ("right", box(10.5, 0, 20.5, 10), "B"),
        ]
    )

    exact = summarize_contour_topology(contours, boundary_tolerance=0.0)
    tolerant = summarize_contour_topology(contours, boundary_tolerance=1.0)

    assert exact.boundary_overlap.empty
    assert len(tolerant.boundary_overlap) == 1
    row = tolerant.boundary_overlap.iloc[0]
    assert row["exact_shared_boundary_length_um"] == pytest.approx(0.0)
    assert row["shared_boundary_length_um"] > 9.0
    assert row["min_distance_um"] == pytest.approx(0.5)


def test_area_overlap_is_reported_without_false_boundary_neighbor():
    contours = _contour_table(
        [
            ("a", box(0, 0, 10, 10), "A"),
            ("b", box(5, 5, 15, 15), "B"),
        ]
    )

    result = summarize_contour_topology(contours, boundary_tolerance=1.0)

    assert len(result.boundary_overlap) == 1
    row = result.boundary_overlap.iloc[0]
    assert row["area_overlap_um2"] == pytest.approx(25.0)
    assert row["shared_boundary_length_um"] == pytest.approx(0.0)
    assert bool(row["has_area_overlap"]) is True
    assert bool(row["is_boundary_neighbor"]) is False


def test_enclosure_reports_outer_to_inner_direction():
    contours = _contour_table(
        [
            ("outer", box(0, 0, 20, 20), "Outer"),
            ("inner", box(5, 5, 10, 10), "Inner"),
        ]
    )

    result = summarize_contour_topology(contours, groupby="structure")

    assert len(result.enclosure) == 1
    row = result.enclosure.iloc[0]
    assert row["outer_contour_id"] == "outer"
    assert row["inner_contour_id"] == "inner"
    assert row["inner_area_covered_fraction"] == pytest.approx(1.0)
    assert bool(row["is_enclosed"]) is True

    summary = result.contour_summary.set_index("contour_id")
    assert summary.loc["outer", "n_contained_contours"] == 1
    assert summary.loc["inner", "n_enclosing_contours"] == 1
    assert len(result.group_enclosure) == 1
    assert result.group_enclosure.iloc[0]["outer_group"] == "Outer"


def test_donut_hole_does_not_enclose_contour_inside_hole():
    donut = Polygon(
        [(0, 0), (20, 0), (20, 20), (0, 20)],
        holes=[[(5, 5), (15, 5), (15, 15), (5, 15)]],
    )
    contours = _contour_table(
        [
            ("donut", donut, "Outer"),
            ("hole_fill", box(7, 7, 13, 13), "Inner"),
        ]
    )

    result = summarize_contour_topology(contours)

    assert result.enclosure.empty
    assert result.boundary_overlap.empty


def test_multipolygon_boundary_overlap_sums_across_parts():
    multi = MultiPolygon([box(0, 0, 10, 10), box(20, 0, 30, 10)])
    contours = _contour_table(
        [
            ("multi", multi, "A"),
            ("middle", box(10, 0, 20, 10), "B"),
        ]
    )

    result = summarize_contour_topology(contours, boundary_tolerance=0.0)

    row = result.boundary_overlap.iloc[0]
    assert row["shared_boundary_length_um"] == pytest.approx(20.0)
    assert row["overlap_fraction_a"] == pytest.approx(0.25)
    assert row["overlap_fraction_b"] == pytest.approx(0.5)


def test_validation_errors_are_clear():
    contours = _contour_table([("a", box(0, 0, 1, 1), "A")])

    with pytest.raises(ValueError, match="boundary_tolerance"):
        summarize_contour_topology(contours, boundary_tolerance=-1.0)
    with pytest.raises(ValueError, match="enclosure_min_fraction"):
        summarize_contour_topology(contours, enclosure_min_fraction=1.5)
    with pytest.raises(ValueError, match="missing required columns"):
        summarize_contour_topology(contours, groupby="missing")
    with pytest.raises(TypeError, match="Shapely"):
        summarize_contour_topology(pd.DataFrame({"contour_id": ["bad"], "geometry": [np.nan]}))
