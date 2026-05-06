from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from shapely.geometry import Polygon, mapping

from histoseg.threed import GlandBiologyMiningConfig, run_gland_biology_mining
from histoseg.threed.cli import main as threed_cli_main


def test_gland_biology_outputs_mechanistic_tables(tmp_path: Path):
    gland_dir, stack_root, cells_path = _write_biology_fixture(tmp_path)
    out_dir = tmp_path / "biology"

    result = run_gland_biology_mining(
        GlandBiologyMiningConfig(
            gland_instance_dir=gland_dir,
            stack_root=stack_root,
            aligned_cells_parquet=cells_path,
            out_dir=out_dir,
        )
    )

    assert result.gland_count == 2
    assert result.assigned_cell_count >= 5
    feature = pd.read_csv(result.gland_biology_feature_matrix_csv)
    fission = pd.read_csv(result.gland_fission_microenvironment_csv)
    entropy = pd.read_csv(result.gland_entropy_and_plasticity_csv)
    axis = pd.read_csv(result.gland_axis_distortion_csv)
    communication = pd.read_csv(result.gland_boundary_communication_csv)
    latent = pd.read_csv(result.gland_morphomolecular_latent_space_csv)
    assignments = pd.read_csv(result.gland_cell_assignments_csv)

    assert {"signature_stem_wnt", "signature_proliferation", "state_entropy"}.issubset(feature.columns)
    assert len(fission) == 1
    assert fission.loc[0, "gland_id"] == "gland_0001"
    entropy_by_gland = entropy.set_index("gland_id")["state_entropy"]
    assert entropy_by_gland["gland_0001"] > entropy_by_gland["gland_0002"]
    assert axis.set_index("gland_id").loc["gland_0001", "proliferation_longitudinal_mean"] > 0
    pair = communication.loc[
        (communication["ligand"] == "WNT2B")
        & (communication["receptor"] == "FZD7")
        & (communication["gland_id"] == "gland_0001")
    ]
    assert not pair.empty
    assert float(pair["max_boundary_interaction_score"].iloc[0]) > 0
    assert {"morphomolecular_pc1", "morphomolecular_pc2"}.issubset(latent.columns)
    g1_assignments = assignments.loc[assignments["gland_id"] == "gland_0001"]
    assert g1_assignments["radial_coord"].max() > g1_assignments["radial_coord"].min()
    assert g1_assignments["longitudinal_coord"].max() == 1.0


def test_analyze_gland_biology_cli_smoke(tmp_path: Path):
    gland_dir, stack_root, cells_path = _write_biology_fixture(tmp_path)
    out_dir = tmp_path / "cli_biology"

    threed_cli_main(
        [
            "analyze-gland-biology",
            "--gland-instance-dir",
            str(gland_dir),
            "--stack-root",
            str(stack_root),
            "--aligned-cells-parquet",
            str(cells_path),
            "--out-dir",
            str(out_dir),
        ]
    )

    assert (out_dir / "gland_biology_feature_matrix.csv").exists()
    assert (out_dir / "gland_boundary_communication.csv").exists()
    assert (out_dir / "gland_morphomolecular_latent_space.csv").exists()


def _write_biology_fixture(tmp_path: Path) -> tuple[Path, Path, Path]:
    stack_root = tmp_path / "stack"
    gland_dir = tmp_path / "glands"
    stack_root.mkdir()
    gland_dir.mkdir()
    xenium_root = tmp_path / "xenium"
    genes = [
        "LGR5",
        "FZD7",
        "MKI67",
        "MUC2",
        "WNT2B",
        "RSPO3",
        "GREM1",
        "FAP",
        "ACTA2",
        "PTPRC",
    ]
    contexts = []
    cells = []
    transcript_rows = {}
    for order, sample in [(1, "slice_1"), (2, "slice_2")]:
        xenium_dir = xenium_root / sample
        (xenium_dir / "cell_feature_matrix").mkdir(parents=True)
        pd.DataFrame(
            [[f"ENSG_{gene}", gene, "Gene Expression"] for gene in genes]
        ).to_csv(
            xenium_dir / "cell_feature_matrix" / "features.tsv.gz",
            sep="\t",
            header=False,
            index=False,
            compression="gzip",
        )
        contexts.append(
            {
                "order": order,
                "sample_id": sample,
                "sample_dir": str(xenium_dir.parent),
                "xenium_dir": str(xenium_dir),
            }
        )
        transcript_rows[sample] = []

    def add_cell(sample: str, order: int, barcode: str, x: float, y: float, genes_for_cell: list[str]):
        cells.append(
            {
                "sample_id": sample,
                "barcode": barcode,
                "slice_order": order,
                "x_3d_um": x,
                "y_3d_um": y,
                "z_um": float((order - 1) * 5),
            }
        )
        for gene in genes_for_cell:
            transcript_rows[sample].append({"cell_id": barcode, "feature_name": gene})

    add_cell("slice_1", 1, "g1_stem", 10, 10, ["LGR5", "FZD7"])
    add_cell("slice_1", 1, "g1_prolif", 15, 10, ["MKI67"])
    add_cell("slice_1", 1, "g1_goblet", 25, 15, ["MUC2"])
    add_cell("slice_1", 1, "g1_ring", 55, 15, ["WNT2B", "RSPO3", "GREM1"])
    add_cell("slice_1", 1, "g2_a", 110, 10, ["MUC2"])
    add_cell("slice_1", 1, "g2_b", 120, 15, ["MUC2"])
    add_cell("slice_2", 2, "g1_s2", 10, 10, ["MKI67", "FZD7"])
    for sample, rows in transcript_rows.items():
        pd.DataFrame(rows).to_parquet(xenium_root / sample / "transcripts.parquet", index=False)
    pd.DataFrame(contexts).to_csv(stack_root / "xenium_slice_manifest.csv", index=False)
    cells_path = tmp_path / "aligned_cells.parquet"
    pd.DataFrame(cells).to_parquet(cells_path, index=False)

    instances = [
        ("001_g1", "gland_0001", 1, "slice_1", _box(0, 0, 35, 30), "002_g1:0.9"),
        ("002_g1", "gland_0001", 2, "slice_2", _box(0, 0, 35, 30), ""),
        ("001_g2", "gland_0002", 1, "slice_1", _box(95, 0, 135, 30), ""),
    ]
    features = []
    rows = []
    for instance_id, gland_id, order, sample, polygon, branch in instances:
        props = {
            "slice_instance_id": instance_id,
            "slice_order": order,
            "sample_id": sample,
            "z_um": float((order - 1) * 5),
            "semantic_structure": "Structure 3",
            "area_um2": float(polygon.area),
            "centroid_x_um": float(polygon.centroid.x),
            "centroid_y_um": float(polygon.centroid.y),
            "lumen_area_um2": 25.0,
            "lumen_centroid_x_um": 10.0 if gland_id == "gland_0001" else 110.0,
            "lumen_centroid_y_um": 10.0,
            "ring_support_score": 1.0,
            "epithelial_marker_score": 1.0,
            "stromal_immune_contamination_score": 0.0,
            "cell_count": 1,
            "qc_flags": "",
            "flag_no_lumen_seed": False,
            "flag_weak_ring": False,
            "flag_merged_candidate": bool(branch),
            "flag_small_fragment": False,
            "marker_profile_json": "{}",
            "gland_id": gland_id,
            "prev_slice_instance_id": "",
            "next_slice_instance_id": "",
            "prev_link_score": "",
            "next_link_score": "",
            "link_score": "",
            "branch_merge_candidates": branch,
        }
        rows.append(props)
        features.append({"type": "Feature", "geometry": mapping(polygon), "properties": props})
    (gland_dir / "slice_gland_instances.geojson").write_text(
        json.dumps({"type": "FeatureCollection", "features": features}),
        encoding="utf-8",
    )
    pd.DataFrame(rows).to_csv(gland_dir / "gland_instance_tracks.csv", index=False)
    pd.DataFrame(
        [
            {
                "gland_id": gland_id,
                "semantic_structure": "Structure 3",
                "slice_count": len([row for row in rows if row["gland_id"] == gland_id]),
                "component_count": len([row for row in rows if row["gland_id"] == gland_id]),
                "first_slice_order": min(row["slice_order"] for row in rows if row["gland_id"] == gland_id),
                "last_slice_order": max(row["slice_order"] for row in rows if row["gland_id"] == gland_id),
                "z_min_um": min(row["z_um"] for row in rows if row["gland_id"] == gland_id),
                "z_max_um": max(row["z_um"] for row in rows if row["gland_id"] == gland_id),
                "z_span_um": max(row["z_um"] for row in rows if row["gland_id"] == gland_id)
                - min(row["z_um"] for row in rows if row["gland_id"] == gland_id),
                "area_median_um2": 1000.0,
                "area_cv": 0.1,
                "median_ring_support_score": 1.0,
                "median_epithelial_marker_score": 1.0,
                "max_stromal_immune_contamination_score": 0.0,
                "branch_merge_candidate_count": 1 if gland_id == "gland_0001" else 0,
                "qc_flags": "merged_candidate" if gland_id == "gland_0001" else "",
                "page": "",
                "page_rendered": False,
            }
            for gland_id in ["gland_0001", "gland_0002"]
        ]
    ).to_csv(gland_dir / "gland_instance_qc_index.csv", index=False)
    return gland_dir, stack_root, cells_path


def _box(x0: float, y0: float, x1: float, y1: float) -> Polygon:
    return Polygon([(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)])
