from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from reproducibility.benchmarks.gaston_brca_fig5 import benchmark as brca


def test_parse_seed_spec_ranges_and_deduplicates():
    assert brca.parse_seed_spec("0-2,2,5") == [0, 1, 2, 5]


def test_select_adh_crop_uses_strict_gaston_tutorial_bounds():
    cells = pd.DataFrame(
        {
            "cell_id": ["a", "b", "c", "d"],
            "Barcode": ["a", "b", "c", "d"],
            "x_centroid": [1500, 1500.1, 3249.9, 3250],
            "y_centroid": [2500, 2000, 3999.9, 2500],
        }
    )
    crop = brca.select_adh_crop(cells, brca.CropBounds())
    assert crop["cell_id"].tolist() == ["c"]
    assert crop["crop_index"].tolist() == [0]


def test_assign_histoseg_seed_domains_covers_fixed_graphclust_mapping():
    crop = pd.DataFrame(
        {
            "cell_id": ["1", "2", "3", "4"],
            "Barcode": ["1", "2", "3", "4"],
            "x_centroid": [1.0, 2.0, 3.0, 4.0],
            "y_centroid": [1.0, 2.0, 3.0, 4.0],
        }
    )
    clusters = pd.DataFrame(
        {
            "Barcode": ["1", "2", "3", "4"],
            "Cluster": [1, 2, 5, 8],
        }
    )
    out = brca.assign_histoseg_seed_domains(crop, clusters)
    assert out["histoseg_seed_structure_name"].tolist() == [
        "Invasive tumor",
        "Stromal",
        "Immune",
        "DCIS/myoepithelial/perivascular",
    ]


def test_assign_histoseg_seed_domains_rejects_unmapped_cluster():
    crop = pd.DataFrame(
        {
            "cell_id": ["1"],
            "Barcode": ["1"],
            "x_centroid": [1.0],
            "y_centroid": [1.0],
        }
    )
    clusters = pd.DataFrame({"Barcode": ["1"], "Cluster": [99]})
    with pytest.raises(brca.BenchmarkInputError, match="Unmapped GraphClust"):
        brca.assign_histoseg_seed_domains(crop, clusters)


def test_validate_cell_type_labels_requires_barcode_alignment(tmp_path: Path):
    labels = tmp_path / "Cell_Barcode_Type_Matrices.csv"
    labels.write_text("Cluster\nB_Cells\n", encoding="utf-8")
    cells = pd.DataFrame({"Barcode": ["1"], "cell_id": ["1"], "x_centroid": [1.0], "y_centroid": [1.0]})
    with pytest.raises(brca.BenchmarkInputError, match="Barcode and Cluster"):
        brca.validate_cell_type_labels(labels, cells, cells)


def test_validate_cell_type_labels_accepts_exact_barcodes(tmp_path: Path):
    labels = tmp_path / "Cell_Barcode_Type_Matrices.csv"
    labels.write_text("Barcode,Cluster\n1,B_Cells\n2,Stromal\n", encoding="utf-8")
    cells = pd.DataFrame(
        {
            "Barcode": ["1", "2"],
            "cell_id": ["1", "2"],
            "x_centroid": [1.0, 2.0],
            "y_centroid": [1.0, 2.0],
        }
    )
    out = brca.validate_cell_type_labels(labels, cells, cells.iloc[[0]])
    assert out["cell_type"].tolist() == ["B_Cells"]


def test_domain_overlap_metrics_on_toy_labels():
    cells = pd.DataFrame(
        {
            "gaston_domain_name": ["Invasive tumor", "Invasive tumor", "Stromal"],
            "histoseg_structure_name": ["Invasive tumor", "Stromal", "Stromal"],
        }
    )
    metrics = brca.compute_domain_overlap_metrics(cells)
    invasive = metrics.loc[metrics["domain"] == "Invasive tumor"].iloc[0]
    overall = metrics.loc[metrics["domain"] == "overall"].iloc[0]
    assert invasive["intersection"] == 1
    assert invasive["union"] == 2
    assert invasive["jaccard"] == pytest.approx(0.5)
    assert overall["accuracy"] == pytest.approx(2 / 3)


def test_smoke_mode_writes_expected_outputs(tmp_path: Path):
    brca.run_smoke(tmp_path)
    expected = [
        "inputs_manifest.json",
        "gaston_reference_cells.parquet",
        "histoseg_domains.parquet",
        "domain_overlap_metrics.csv",
        "isodepth_sdf_metrics.csv",
        "cell_type_proportion_metrics.csv",
        "gene_gradient_metrics.csv",
        "figures/gaston_domains.png",
        "figures/histoseg_domains.svg",
    ]
    for rel in expected:
        assert (tmp_path / rel).exists(), rel
    overlap = pd.read_csv(tmp_path / "domain_overlap_metrics.csv")
    assert np.isfinite(overlap.loc[overlap["domain"] == "overall", "accuracy"]).all()

