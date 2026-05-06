from __future__ import annotations

from pathlib import Path

import pandas as pd

from histoseg.threed import (
    GlandPositionAtlasConfig,
    render_gland_position_atlas,
)


def test_linear_three_slice_gland_yields_linear_family(tmp_path: Path):
    gland_dir = _write_position_fixture(
        tmp_path,
        [
            _row("s1_a", "gland_0001", 1, 0, 0, next_id="s2_a", next_score=0.9),
            _row("s2_a", "gland_0001", 2, 4, 0, prev_id="s1_a", next_id="s3_a", prev_score=0.9, next_score=0.9),
            _row("s3_a", "gland_0001", 3, 8, 0, prev_id="s2_a", prev_score=0.9),
        ],
    )

    result = render_gland_position_atlas(
        GlandPositionAtlasConfig(gland_instance_dir=gland_dir, out_dir=tmp_path / "atlas")
    )

    families = pd.read_csv(result.gland_family_index_csv)
    assert result.gland_count == 1
    assert result.family_count == 1
    assert families.loc[0, "motif_type"] == "linear"
    assert families.loc[0, "member_gland_ids"] == "gland_0001"


def test_candidate_split_family_preserves_original_gland_ids(tmp_path: Path):
    gland_dir = _write_position_fixture(
        tmp_path,
        [
            _row(
                "s1_a",
                "gland_0001",
                1,
                0,
                0,
                next_id="s2_a",
                next_score=0.9,
                branch="s2_b:0.61",
            ),
            _row("s2_a", "gland_0001", 2, 4, 0, prev_id="s1_a", prev_score=0.9),
            _row("s2_b", "gland_0002", 2, 8, 0),
        ],
        biology=True,
    )

    result = render_gland_position_atlas(
        GlandPositionAtlasConfig(gland_instance_dir=gland_dir, out_dir=tmp_path / "atlas")
    )

    families = pd.read_csv(result.gland_family_index_csv)
    links = pd.read_csv(result.gland_family_links_csv)
    assert set(families.loc[0, "member_gland_ids"].split(";")) == {"gland_0001", "gland_0002"}
    assert families.loc[0, "motif_type"] == "candidate_split"
    assert set(links["link_type"]) == {"accepted", "candidate_branch_merge"}
    assert links["source_gland_id"].nunique() == 1
    assert result.gland_count == 2


def test_candidate_merge_family_is_recorded(tmp_path: Path):
    gland_dir = _write_position_fixture(
        tmp_path,
        [
            _row("s1_a", "gland_0001", 1, 0, 0, next_id="s2_ab", next_score=0.9),
            _row("s1_b", "gland_0002", 1, 12, 0, branch="s2_ab:0.72"),
            _row("s2_ab", "gland_0003", 2, 6, 0, prev_id="s1_a", prev_score=0.9),
        ],
    )

    result = render_gland_position_atlas(
        GlandPositionAtlasConfig(gland_instance_dir=gland_dir, out_dir=tmp_path / "atlas")
    )

    families = pd.read_csv(result.gland_family_index_csv)
    events = pd.read_csv(result.fission_events_csv)
    assert families.loc[0, "motif_type"] == "candidate_merge"
    assert events.loc[0, "event_type"] == "candidate_merge"


def test_low_score_candidate_is_excluded(tmp_path: Path):
    gland_dir = _write_position_fixture(
        tmp_path,
        [
            _row("s1_a", "gland_0001", 1, 0, 0, next_id="s2_a", next_score=0.9, branch="s2_b:0.40"),
            _row("s2_a", "gland_0001", 2, 4, 0, prev_id="s1_a", prev_score=0.9),
            _row("s2_b", "gland_0002", 2, 9, 0),
        ],
    )

    result = render_gland_position_atlas(
        GlandPositionAtlasConfig(gland_instance_dir=gland_dir, out_dir=tmp_path / "atlas")
    )

    families = pd.read_csv(result.gland_family_index_csv)
    links = pd.read_csv(result.gland_family_links_csv)
    assert result.family_count == 2
    assert "candidate_branch_merge" not in set(links["link_type"])
    assert set(families["member_gland_ids"]) == {"gland_0001", "gland_0002"}


def test_position_atlas_html_smoke(tmp_path: Path):
    gland_dir = _write_position_fixture(
        tmp_path,
        [
            _row("s1_a", "gland_0001", 1, 0, 0, next_id="s2_a", next_score=0.9, branch="s2_b:0.7"),
            _row("s2_a", "gland_0001", 2, 5, 0, prev_id="s1_a", prev_score=0.9),
            _row("s2_b", "gland_0002", 2, 10, 0),
        ],
        biology=True,
    )

    result = render_gland_position_atlas(
        GlandPositionAtlasConfig(gland_instance_dir=gland_dir, out_dir=tmp_path / "atlas")
    )

    html = result.gland_position_atlas_html.read_text(encoding="utf-8")
    assert "Plotly.newPlot" in html
    assert "gland_0001" in html
    assert "family_0001" in html
    assert "candidate_branch_merge" in html
    assert (tmp_path / "atlas" / "fission_events" / "event_0001.html").exists()


def _write_position_fixture(tmp_path: Path, rows: list[dict], *, biology: bool = False) -> Path:
    gland_dir = tmp_path / "glands"
    gland_dir.mkdir()
    tracks = pd.DataFrame(rows)
    tracks.to_csv(gland_dir / "gland_instance_tracks.csv", index=False)
    index_rows = []
    for gland_id, group in tracks.groupby("gland_id", sort=True):
        index_rows.append(
            {
                "gland_id": gland_id,
                "semantic_structure": group["semantic_structure"].mode().iloc[0],
                "slice_count": int(group["slice_order"].nunique()),
                "component_count": int(len(group)),
                "first_slice_order": int(group["slice_order"].min()),
                "last_slice_order": int(group["slice_order"].max()),
                "z_min_um": float(group["z_um"].min()),
                "z_max_um": float(group["z_um"].max()),
                "z_span_um": float(group["z_um"].max() - group["z_um"].min()),
                "area_median_um2": float(group["area_um2"].median()),
                "area_cv": 0.0,
                "median_ring_support_score": 1.0,
                "median_epithelial_marker_score": 2.0,
                "max_stromal_immune_contamination_score": 0.0,
                "branch_merge_candidate_count": int(group["branch_merge_candidates"].astype(str).map(bool).sum()),
                "qc_flags": "",
                "page": f"glands/{gland_id}.html",
                "page_rendered": True,
            }
        )
    pd.DataFrame(index_rows).to_csv(gland_dir / "gland_instance_qc_index.csv", index=False)
    if biology:
        bio_dir = gland_dir / "biology_mining_20260505"
        bio_dir.mkdir()
        pd.DataFrame(
            [
                {
                    "gland_id": gland_id,
                    "state_entropy": float(i + 1),
                    "fission_candidate_score": float(i + 1) / 10,
                    "lumen_collapse_score": float(i) / 10,
                    "signature_stromal_niche": float(i + 2),
                    "signature_stem_wnt": float(i + 3),
                    "signature_proliferation": float(i + 4),
                }
                for i, gland_id in enumerate(sorted(tracks["gland_id"].unique()))
            ]
        ).to_csv(bio_dir / "gland_biology_feature_matrix.csv", index=False)
        pd.DataFrame(
            [
                {
                    "gland_id": gland_id,
                    "ligand": "WNT2B",
                    "receptor": "FZD7",
                    "interaction": "stromal WNT to epithelial FZD",
                    "mean_boundary_interaction_score": 0.1,
                    "max_boundary_interaction_score": 0.2,
                    "total_outer_ring_cell_count": 10,
                    "total_inside_cell_count": 20,
                }
                for gland_id in sorted(tracks["gland_id"].unique())
            ]
        ).to_csv(bio_dir / "gland_boundary_communication.csv", index=False)
    return gland_dir


def _row(
    slice_instance_id: str,
    gland_id: str,
    slice_order: int,
    x: float,
    y: float,
    *,
    prev_id: str = "",
    next_id: str = "",
    prev_score: float | str = "",
    next_score: float | str = "",
    branch: str = "",
) -> dict:
    return {
        "slice_instance_id": slice_instance_id,
        "slice_order": slice_order,
        "sample_id": f"slice_{slice_order}",
        "z_um": float((slice_order - 1) * 5),
        "semantic_structure": "Structure 3",
        "area_um2": 100.0,
        "centroid_x_um": float(x),
        "centroid_y_um": float(y),
        "lumen_area_um2": 20.0,
        "lumen_centroid_x_um": float(x),
        "lumen_centroid_y_um": float(y),
        "ring_support_score": 1.0,
        "epithelial_marker_score": 2.0,
        "stromal_immune_contamination_score": 0.0,
        "cell_count": 12,
        "qc_flags": "",
        "flag_no_lumen_seed": False,
        "flag_weak_ring": False,
        "flag_merged_candidate": bool(branch),
        "flag_small_fragment": False,
        "marker_profile_json": "{}",
        "gland_id": gland_id,
        "prev_slice_instance_id": prev_id,
        "next_slice_instance_id": next_id,
        "prev_link_score": prev_score,
        "next_link_score": next_score,
        "link_score": next_score or prev_score or "",
        "branch_merge_candidates": branch,
    }
