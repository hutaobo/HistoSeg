from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from shapely.geometry import Polygon, mapping

from histoseg.threed import (
    GlandInstanceTrackingConfig,
    track_gland_instances,
)


def test_three_slice_stable_gland_id(tmp_path: Path):
    seg_dir = _write_tracking_fixture(
        tmp_path,
        [
            _instance("s1_a", 1, 0, 0),
            _instance("s2_a", 2, 4, 1),
            _instance("s3_a", 3, 8, 2),
        ],
    )

    result = track_gland_instances(
        GlandInstanceTrackingConfig(segmentation_result_dir=seg_dir, out_dir=tmp_path / "out")
    )

    tracks = pd.read_csv(result.gland_instance_tracks_csv)
    assert tracks["gland_id"].nunique() == 1
    assert tracks["slice_order"].tolist() == [1, 2, 3]


def test_one_to_one_assignment(tmp_path: Path):
    seg_dir = _write_tracking_fixture(
        tmp_path,
        [
            _instance("s1_a", 1, 0, 0),
            _instance("s1_b", 1, 100, 0),
            _instance("s2_a", 2, 3, 0),
            _instance("s2_b", 2, 103, 0),
        ],
    )

    result = track_gland_instances(
        GlandInstanceTrackingConfig(segmentation_result_dir=seg_dir, out_dir=tmp_path / "out")
    )

    tracks = pd.read_csv(result.gland_instance_tracks_csv)
    assert tracks["gland_id"].nunique() == 2
    assert tracks.groupby("gland_id")["slice_order"].nunique().tolist() == [2, 2]


def test_branch_merge_flagged_not_hidden(tmp_path: Path):
    seg_dir = _write_tracking_fixture(
        tmp_path,
        [
            _instance("s1_a", 1, 0, 0),
            _instance("s1_b", 1, 32, 0),
            _instance("s2_ab", 2, 16, 0, width=44),
        ],
    )

    result = track_gland_instances(
        GlandInstanceTrackingConfig(
            segmentation_result_dir=seg_dir,
            out_dir=tmp_path / "out",
            max_centroid_distance_um=80,
            min_overlap_ratio=0.0,
        )
    )

    tracks = pd.read_csv(result.gland_instance_tracks_csv)
    assert tracks["gland_id"].nunique() == 2
    assert tracks["branch_merge_candidates"].fillna("").astype(str).map(bool).any()


def _write_tracking_fixture(tmp_path: Path, instances: list[dict]) -> Path:
    seg_dir = tmp_path / "seg"
    seg_dir.mkdir()
    features = []
    rows = []
    for payload in instances:
        polygon = payload["polygon"]
        props = {
            "slice_instance_id": payload["slice_instance_id"],
            "slice_order": payload["slice_order"],
            "sample_id": f"slice_{payload['slice_order']}",
            "z_um": float((payload["slice_order"] - 1) * 5),
            "semantic_structure": "Structure 3",
            "area_um2": float(polygon.area),
            "centroid_x_um": float(polygon.centroid.x),
            "centroid_y_um": float(polygon.centroid.y),
            "lumen_area_um2": 25.0,
            "lumen_centroid_x_um": float(polygon.centroid.x),
            "lumen_centroid_y_um": float(polygon.centroid.y),
            "ring_support_score": 0.8,
            "epithelial_marker_score": 2.0,
            "stromal_immune_contamination_score": 0.0,
            "qc_flags": "",
            "flag_no_lumen_seed": False,
            "flag_weak_ring": False,
            "flag_merged_candidate": False,
            "flag_small_fragment": False,
            "marker_profile_json": json.dumps({"EPCAM": 10.0}),
        }
        features.append(
            {
                "type": "Feature",
                "geometry": mapping(polygon),
                "properties": props,
            }
        )
        rows.append(props)
    (seg_dir / "slice_gland_instances.geojson").write_text(
        json.dumps({"type": "FeatureCollection", "features": features}),
        encoding="utf-8",
    )
    pd.DataFrame(rows).to_csv(seg_dir / "slice_gland_instances_qc.csv", index=False)
    return seg_dir


def _instance(
    slice_instance_id: str,
    slice_order: int,
    x: float,
    y: float,
    *,
    width: float = 30,
    height: float = 30,
) -> dict:
    polygon = Polygon(
        [
            (x, y),
            (x + width, y),
            (x + width, y + height),
            (x, y + height),
            (x, y),
        ]
    )
    return {
        "slice_instance_id": slice_instance_id,
        "slice_order": slice_order,
        "polygon": polygon,
    }
