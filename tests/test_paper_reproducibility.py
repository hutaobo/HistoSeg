from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_PIPELINE = REPO_ROOT / "reproducibility" / "run_paper_pipeline.py"


def _load_pipeline_module():
    spec = importlib.util.spec_from_file_location("run_paper_pipeline", RUN_PIPELINE)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_run_paper_pipeline_dry_run_uses_public_manifest(tmp_path, capsys):
    module = _load_pipeline_module()
    data_root = tmp_path / "public_data"
    stack_root = data_root / "polyp" / "histoseg_3d_reconstruction"
    stack_root.mkdir(parents=True)
    (stack_root / "aligned_slice_manifest.csv").write_text("order,sample_id,z_um\n", encoding="utf-8")
    manifest = tmp_path / "paper_data_manifest.json"
    _write_manifest(manifest, accession="10.5281/zenodo.1234567")

    assert module.main(
        [
            "--dry-run",
            "--data-root",
            str(data_root),
            "--data-manifest",
            str(manifest),
            "--skip-cell-cloud",
            "--skip-clustermaps",
            "--skip-h5ad-provenance",
        ]
    ) == 0
    out = capsys.readouterr().out
    assert "histoseg_3d_reconstruction" in out
    assert "10.5281/zenodo.1234567" in out
    assert "Y:/" not in out
    assert "Y:\\" not in out
    assert "spatialpathologist" not in out


def test_run_paper_pipeline_rejects_pending_accessions(tmp_path):
    module = _load_pipeline_module()
    data_root = tmp_path / "public_data"
    stack_root = data_root / "polyp" / "histoseg_3d_reconstruction"
    stack_root.mkdir(parents=True)
    (stack_root / "aligned_slice_manifest.csv").write_text("order,sample_id,z_um\n", encoding="utf-8")
    manifest = tmp_path / "paper_data_manifest.json"
    _write_manifest(manifest, accession="PENDING_ZENODO_DOI")

    with pytest.raises(ValueError, match="DOI/accession"):
        module.main(
            [
                "--dry-run",
                "--data-root",
                str(data_root),
                "--data-manifest",
                str(manifest),
                "--skip-cell-cloud",
                "--skip-clustermaps",
                "--skip-h5ad-provenance",
            ]
        )


def test_public_reproducibility_files_do_not_encode_local_windows_paths():
    paths = [
        REPO_ROOT / "reproducibility" / "run_paper_pipeline.py",
        REPO_ROOT / "reproducibility" / "paper_data_manifest.json",
        REPO_ROOT / "reproducibility" / "results_manifest.json",
    ]
    for path in paths:
        text = path.read_text(encoding="utf-8")
        assert "Y:/" not in text
        assert "Y:\\" not in text
        assert "spatialpathologist" not in text


def _write_manifest(path: Path, *, accession: str) -> None:
    records = [
        ("stack_root", "polyp/histoseg_3d_reconstruction", None),
        ("aligned_slice_manifest", "polyp/histoseg_3d_reconstruction/aligned_slice_manifest.csv", 21),
        ("aligned_cells_parquet", "polyp/histoseg_3d_reconstruction/aligned_leiden_3d_cells.parquet", None),
        ("h5ad", "polyp/expression/polyp_32slice_processed_leiden.h5ad", None),
        (
            "spatial_module_batch_dir",
            "polyp/histoseg_3d_reconstruction/gene_overlays/batch_3d_genes_starter_panel",
            None,
        ),
        (
            "batch_gene_status",
            "polyp/histoseg_3d_reconstruction/gene_overlays/batch_3d_genes_starter_panel/batch_gene_status.csv",
            None,
        ),
        (
            "fraction_inside_matrix",
            "polyp/histoseg_3d_reconstruction/gene_overlays/batch_3d_genes_starter_panel/gene_structure_fraction_inside_matrix.csv",
            None,
        ),
        (
            "signed_distance_matrix",
            "polyp/histoseg_3d_reconstruction/gene_overlays/batch_3d_genes_starter_panel/gene_structure_signed_distance_matrix.csv",
            None,
        ),
    ]
    payload = {
        "schema_version": 1,
        "dataset": "test_public_bundle",
        "archive_status": "test",
        "files": [
            {
                "role": role,
                "relative_path": relative_path,
                "accession_or_doi": accession,
                "sha256": None,
                "size_bytes": size_bytes,
                "license_or_access_terms": "test",
                "generated_by": "pytest",
            }
            for role, relative_path, size_bytes in records
        ],
    }
    path.write_text(json.dumps(payload), encoding="utf-8")
