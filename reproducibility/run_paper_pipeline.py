"""Run the HistoSeg paper-figure reproduction wrapper.

This script intentionally remains a thin orchestrator around the public
``histoseg-3d`` CLI. It does not implement reconstruction, quantification, or
rendering logic itself; it checks expected inputs, calls the CLI commands used
for the paper figure artifacts, and writes a JSON manifest with output hashes
and alignment provenance.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any, Iterable, Sequence


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_ROOT = REPO_ROOT / "reproducibility" / "public_data"
DEFAULT_DATA_MANIFEST = REPO_ROOT / "reproducibility" / "paper_data_manifest.json"
DEFAULT_RESULTS_DIR = REPO_ROOT / "reproducibility" / "results"
DEFAULT_MANIFEST_JSON = REPO_ROOT / "reproducibility" / "results_manifest.json"
DEFAULT_HOTSPOT = "top05"
DEFAULT_GENES = (
    "GREM1",
    "COL1A1",
    "COL1A2",
    "ACTA2",
    "PDGFRA",
    "TAGLN",
    "FAP",
    "EPCAM",
    "MUC2",
    "OLFM4",
    "LGR5",
    "MKI67",
    "PTPRC",
    "CD3D",
    "CD3E",
    "CD8A",
    "CD4",
    "MS4A1",
    "CD79A",
    "LYZ",
    "CD68",
    "C1QA",
    "PECAM1",
    "RGS5",
)
DEFAULT_MATRICES = ("fraction_inside", "signed_distance")
SMALL_INPUT_HASH_LIMIT_BYTES = 256 * 1024 * 1024


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    data_root = Path(args.data_root).expanduser()
    data_manifest_path = Path(args.data_manifest).expanduser()
    data_manifest = _load_paper_data_manifest(data_manifest_path)
    data_records = _manifest_records_by_role(data_manifest)

    stack_root = _resolve_path_arg(args.stack_root, data_records, "stack_root", data_root)
    h5ad = _resolve_path_arg(args.h5ad, data_records, "h5ad", data_root)
    aligned_cells_parquet = _resolve_path_arg(
        args.aligned_cells_parquet,
        data_records,
        "aligned_cells_parquet",
        data_root,
    )
    batch_dir = _resolve_path_arg(args.batch_dir, data_records, "spatial_module_batch_dir", data_root)
    results_dir = Path(args.results_dir).expanduser()
    manifest_json = Path(args.manifest_json).expanduser()
    genes = tuple(_parse_tokens(args.genes))
    matrices = tuple(args.matrices)
    required_roles = _required_input_roles(
        matrices,
        run_discovery=args.run_discovery,
        skip_cell_cloud=args.skip_cell_cloud,
        skip_clustermaps=args.skip_clustermaps,
    )

    outputs = _planned_outputs(results_dir, args.hotspot, matrices)
    inputs = _input_records(
        stack_root=stack_root,
        h5ad=h5ad,
        aligned_cells_parquet=aligned_cells_parquet,
        batch_dir=batch_dir,
        matrices=matrices,
        hash_large_inputs=args.hash_large_inputs,
        data_records=data_records,
    )

    _validate_manifest_accessions(
        data_records,
        required_roles=required_roles,
        allow_pending=bool(args.allow_pending_accessions),
    )
    _validate_inputs(inputs, required_roles=required_roles)
    _validate_cli(args.histoseg_cli, dry_run=args.dry_run)

    manifest: dict[str, Any] = {
        "schema_version": 1,
        "pipeline": "histoseg_polyp_24_gene_paper_figures",
        "created_at": _utc_now(),
        "repo_root": str(REPO_ROOT),
        "dry_run": bool(args.dry_run),
        "parameters": {
            "stack_root": str(stack_root),
            "h5ad": str(h5ad),
            "aligned_cells_parquet": str(aligned_cells_parquet),
            "batch_dir": str(batch_dir),
            "data_root": str(data_root),
            "data_manifest": str(data_manifest_path),
            "allow_pending_accessions": bool(args.allow_pending_accessions),
            "results_dir": str(results_dir),
            "hotspot": args.hotspot,
            "genes": list(genes),
            "matrices": list(matrices),
            "max_points": int(args.max_points),
            "random_state": int(args.random_state),
            "label_column": args.label_column,
            "run_discovery": bool(args.run_discovery),
            "skip_cell_cloud": bool(args.skip_cell_cloud),
            "skip_clustermaps": bool(args.skip_clustermaps),
        },
        "paper_data_manifest": {
            "path": str(data_manifest_path),
            "schema_version": data_manifest.get("schema_version"),
            "dataset": data_manifest.get("dataset"),
            "archive_status": data_manifest.get("archive_status"),
        },
        "inputs": inputs,
        "commands": [],
        "outputs": [],
        "alignment": _alignment_provenance(
            stack_root=stack_root,
            h5ad=h5ad,
            pixel_size_um=float(args.pixel_size_um),
            skip_h5ad_provenance=bool(args.skip_h5ad_provenance),
        ),
    }

    if not args.dry_run:
        results_dir.mkdir(parents=True, exist_ok=True)
        manifest_json.parent.mkdir(parents=True, exist_ok=True)

    if args.run_discovery:
        command = [
            args.histoseg_cli,
            "discover-spatial-modules",
            "--h5ad",
            str(h5ad),
            "--aligned-cells-parquet",
            str(aligned_cells_parquet),
            "--stack-root",
            str(stack_root),
            "--out-dir",
            str(batch_dir),
            "--genes",
            *genes,
            "--mesh-export-formats",
            args.mesh_export_formats,
            "--pixel-size-um",
            str(args.pixel_size_um),
        ]
        if args.force_rebuild_masks:
            command.append("--force-rebuild-masks")
        manifest["commands"].append(_run_command("discover_spatial_modules", command, args.dry_run))

    if not args.skip_cell_cloud:
        command = [
            args.histoseg_cli,
            "render-cell-cloud",
            "--stack-root",
            str(stack_root),
            "--aligned-cells-parquet",
            str(aligned_cells_parquet),
            "--out-html",
            str(outputs["cell_cloud_html"]),
            "--label-column",
            args.label_column,
            "--max-points",
            str(args.max_points),
            "--random-state",
            str(args.random_state),
            "--title",
            args.cell_cloud_title,
        ]
        manifest["commands"].append(_run_command("render_cell_cloud", command, args.dry_run))

    if not args.skip_clustermaps:
        for matrix in matrices:
            out_png = outputs["clustermaps"][matrix]
            command = [
                args.histoseg_cli,
                "plot-spatial-modules",
                "--batch-dir",
                str(batch_dir),
                "--matrix",
                matrix,
                "--hotspot",
                args.hotspot,
                "--out-png",
                str(out_png),
            ]
            manifest["commands"].append(_run_command(f"plot_{matrix}", command, args.dry_run))

    manifest["outputs"] = [
        _file_record(path, role=role, sha256=not args.dry_run, skip_large=False)
        for role, path in _iter_output_paths(outputs, args.skip_cell_cloud, args.skip_clustermaps)
    ]
    _write_manifest(manifest_json, manifest, dry_run=args.dry_run)

    if args.dry_run:
        print(json.dumps(manifest, indent=2, ensure_ascii=True, default=str))
    else:
        print(f"Wrote reproduction manifest: {manifest_json}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate the HistoSeg paper reproduction artifacts by calling "
            "existing histoseg-3d CLI commands."
        )
    )
    parser.add_argument("--data-root", default=str(DEFAULT_DATA_ROOT))
    parser.add_argument("--data-manifest", default=str(DEFAULT_DATA_MANIFEST))
    parser.add_argument(
        "--allow-pending-accessions",
        action="store_true",
        help=(
            "Permit manifest records whose accession_or_doi is still marked as pending. "
            "Use only before the public archive DOI/accession is assigned."
        ),
    )
    parser.add_argument(
        "--stack-root",
        default=None,
        help="Optional local override for the reconstructed stack root; defaults to the data manifest.",
    )
    parser.add_argument(
        "--h5ad",
        default=None,
        help="Optional local override for the processed expression object; defaults to the data manifest.",
    )
    parser.add_argument(
        "--aligned-cells-parquet",
        default=None,
        help="Optional local override for aligned cells Parquet; defaults to the data manifest.",
    )
    parser.add_argument(
        "--batch-dir",
        default=None,
        help="Optional local override for spatial-module batch outputs; defaults to the data manifest.",
    )
    parser.add_argument("--results-dir", default=str(DEFAULT_RESULTS_DIR))
    parser.add_argument("--manifest-json", default=str(DEFAULT_MANIFEST_JSON))
    parser.add_argument("--histoseg-cli", default="histoseg-3d")
    parser.add_argument("--genes", nargs="*", default=list(DEFAULT_GENES))
    parser.add_argument("--hotspot", choices=("top15", "top10", "top05"), default=DEFAULT_HOTSPOT)
    parser.add_argument(
        "--matrices",
        nargs="*",
        choices=("fraction_inside", "overlap_fraction", "signed_distance"),
        default=list(DEFAULT_MATRICES),
    )
    parser.add_argument("--label-column", default="leiden_1_0")
    parser.add_argument("--max-points", type=int, default=300000)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--pixel-size-um", type=float, default=0.2125)
    parser.add_argument("--mesh-export-formats", default="ply")
    parser.add_argument("--cell-cloud-title", default="HistoSeg 32-slice polyp Leiden cell cloud")
    parser.add_argument(
        "--run-discovery",
        action="store_true",
        help=(
            "Rerun discover-spatial-modules before plotting. By default, the "
            "script uses the validated starter-panel batch matrices."
        ),
    )
    parser.add_argument("--force-rebuild-masks", action="store_true")
    parser.add_argument("--skip-cell-cloud", action="store_true")
    parser.add_argument("--skip-clustermaps", action="store_true")
    parser.add_argument("--skip-h5ad-provenance", action="store_true")
    parser.add_argument(
        "--hash-large-inputs",
        action="store_true",
        help="Also SHA256-hash large .h5ad/Parquet inputs. This can take time.",
    )
    parser.add_argument("--dry-run", action="store_true")
    return parser


def _planned_outputs(results_dir: Path, hotspot: str, matrices: Sequence[str]) -> dict[str, Any]:
    return {
        "cell_cloud_html": results_dir / "figure1_leiden_3d_cells_300k.html",
        "clustermaps": {
            matrix: results_dir / f"figure3_{matrix}_{hotspot}_spatial_clustermap.png"
            for matrix in matrices
        },
    }


def _load_paper_data_manifest(path: Path) -> dict[str, Any]:
    expanded = Path(path).expanduser()
    if not expanded.exists():
        raise FileNotFoundError(
            f"Paper data manifest not found: {expanded}. "
            "Pass --data-manifest or create the public data manifest first."
        )
    payload = json.loads(expanded.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Paper data manifest must be a JSON object: {expanded}")
    files = payload.get("files")
    if not isinstance(files, list) or not files:
        raise ValueError(f"Paper data manifest requires a non-empty 'files' list: {expanded}")
    required_fields = {
        "role",
        "relative_path",
        "accession_or_doi",
        "sha256",
        "size_bytes",
        "license_or_access_terms",
        "generated_by",
    }
    for index, record in enumerate(files):
        if not isinstance(record, dict):
            raise ValueError(f"Manifest file entry {index} must be an object.")
        missing = sorted(required_fields.difference(record))
        if missing:
            raise ValueError(
                f"Manifest file entry {index} is missing required field(s): {', '.join(missing)}"
            )
    return payload


def _manifest_records_by_role(manifest: dict[str, Any]) -> dict[str, dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}
    for record in manifest.get("files", []):
        role = str(record["role"])
        if role in records:
            raise ValueError(f"Duplicate role in paper data manifest: {role}")
        records[role] = dict(record)
    return records


def _resolve_path_arg(
    value: str | None,
    data_records: dict[str, dict[str, Any]],
    role: str,
    data_root: Path,
) -> Path:
    if value:
        return Path(value).expanduser()
    record = data_records.get(role)
    if record is None:
        raise KeyError(f"Paper data manifest is missing required role: {role}")
    relative_path = Path(str(record["relative_path"]))
    if relative_path.is_absolute():
        return relative_path.expanduser()
    return (data_root / relative_path).expanduser()


def _validate_manifest_accessions(
    data_records: dict[str, dict[str, Any]],
    *,
    required_roles: set[str],
    allow_pending: bool,
) -> None:
    missing_roles = sorted(role for role in required_roles if role not in data_records)
    if missing_roles:
        raise KeyError(
            "Paper data manifest is missing required role(s): " + ", ".join(missing_roles)
        )
    if allow_pending:
        return
    pending_prefixes = ("PENDING", "TBD", "TODO")
    pending = []
    for role in sorted(required_roles):
        accession = str(data_records[role].get("accession_or_doi", "")).strip().upper()
        if not accession or accession.startswith(pending_prefixes):
            pending.append(f"{role}: {data_records[role].get('accession_or_doi')!r}")
    if pending:
        details = "\n  - ".join(pending)
        raise ValueError(
            "Required public-data manifest records must have DOI/accession-backed "
            "accession_or_doi values before full reproduction:\n  - "
            f"{details}\nPass --allow-pending-accessions only for presubmission local checks."
        )


def _input_records(
    *,
    stack_root: Path,
    h5ad: Path,
    aligned_cells_parquet: Path,
    batch_dir: Path,
    matrices: Sequence[str],
    hash_large_inputs: bool,
    data_records: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    records = [
        _file_record(stack_root, role="stack_root", sha256=False),
        _file_record(stack_root / "aligned_slice_manifest.csv", role="aligned_slice_manifest"),
        _file_record(
            aligned_cells_parquet,
            role="aligned_cells_parquet",
            sha256=hash_large_inputs,
            skip_large=not hash_large_inputs,
        ),
        _file_record(h5ad, role="h5ad", sha256=hash_large_inputs, skip_large=not hash_large_inputs),
        _file_record(batch_dir, role="spatial_module_batch_dir", sha256=False),
        _file_record(batch_dir / "batch_gene_status.csv", role="batch_gene_status"),
    ]
    matrix_filenames = {
        "fraction_inside": "gene_structure_fraction_inside_matrix.csv",
        "overlap_fraction": "gene_structure_overlap_fraction_matrix.csv",
        "signed_distance": "gene_structure_signed_distance_matrix.csv",
    }
    for matrix in matrices:
        records.append(_file_record(batch_dir / matrix_filenames[matrix], role=f"{matrix}_matrix"))
    return [_with_manifest_record(record, data_records.get(record["role"])) for record in records]


def _required_input_roles(
    matrices: Sequence[str],
    *,
    run_discovery: bool,
    skip_cell_cloud: bool,
    skip_clustermaps: bool,
) -> set[str]:
    required_roles = {"stack_root", "aligned_slice_manifest"}
    if run_discovery:
        required_roles.update({"h5ad", "aligned_cells_parquet"})
    if not skip_cell_cloud:
        required_roles.add("aligned_cells_parquet")
    if not skip_clustermaps:
        required_roles.update({"spatial_module_batch_dir", "batch_gene_status"})
        required_roles.update(f"{matrix}_matrix" for matrix in matrices)
    return required_roles


def _validate_inputs(inputs: Sequence[dict[str, Any]], *, required_roles: set[str]) -> None:
    missing = [
        f"{record['role']}: {record['path']}"
        for record in inputs
        if record["role"] in required_roles and not record["exists"]
    ]
    if missing:
        details = "\n  - ".join(missing)
        raise FileNotFoundError(f"Missing required reproduction input(s):\n  - {details}")


def _validate_cli(histoseg_cli: str, *, dry_run: bool) -> None:
    if dry_run:
        return
    command_path = Path(histoseg_cli)
    if command_path.parent != Path(".") and command_path.exists():
        return
    if shutil.which(histoseg_cli) is None:
        raise FileNotFoundError(
            f"Could not find {histoseg_cli!r} on PATH. Activate the HistoSeg environment first."
        )


def _run_command(name: str, command: Sequence[str], dry_run: bool) -> dict[str, Any]:
    record: dict[str, Any] = {
        "name": name,
        "argv": list(command),
        "started_at": _utc_now(),
        "dry_run": bool(dry_run),
    }
    if dry_run:
        record.update({"returncode": None, "stdout": "", "stderr": "", "finished_at": _utc_now()})
        return record

    completed = subprocess.run(
        command,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    record.update(
        {
            "returncode": int(completed.returncode),
            "stdout": _truncate(completed.stdout),
            "stderr": _truncate(completed.stderr),
            "finished_at": _utc_now(),
        }
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"{name} failed with exit code {completed.returncode}.\n"
            f"Command: {' '.join(command)}\n"
            f"stderr:\n{completed.stderr}"
        )
    return record


def _alignment_provenance(
    *,
    stack_root: Path,
    h5ad: Path,
    pixel_size_um: float,
    skip_h5ad_provenance: bool,
) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "computed": {
            "status": "not_available",
            "alignment_hash": None,
            "manifest": None,
        },
        "h5ad_cache": {
            "status": "skipped" if skip_h5ad_provenance else "not_available",
            "alignment_hash": None,
            "provenance": None,
        },
    }
    try:
        from histoseg.threed import build_alignment_manifest, hash_alignment_manifest

        manifest = build_alignment_manifest(stack_root, pixel_size_um=pixel_size_um)
        payload["computed"] = {
            "status": "ok",
            "alignment_hash": hash_alignment_manifest(manifest),
            "manifest": manifest,
        }
    except Exception as exc:  # pragma: no cover - depends on local data state.
        payload["computed"] = {
            "status": "error",
            "error": f"{type(exc).__name__}: {exc}",
            "alignment_hash": None,
            "manifest": None,
        }

    if skip_h5ad_provenance:
        return payload
    if not h5ad.exists():
        payload["h5ad_cache"] = {"status": "missing_h5ad", "alignment_hash": None, "provenance": None}
        return payload
    try:
        import anndata as ad

        adata = ad.read_h5ad(h5ad, backed="r")
        try:
            provenance = dict(getattr(adata, "uns", {}).get("histoseg_3d_alignment", {}))
        finally:
            if getattr(adata, "isbacked", False):
                adata.file.close()
        payload["h5ad_cache"] = {
            "status": "ok" if provenance else "missing_cache",
            "alignment_hash": provenance.get("alignment_hash"),
            "provenance": provenance or None,
        }
    except Exception as exc:  # pragma: no cover - optional local dependency/data state.
        payload["h5ad_cache"] = {
            "status": "error",
            "error": f"{type(exc).__name__}: {exc}",
            "alignment_hash": None,
            "provenance": None,
        }
    return payload


def _iter_output_paths(
    outputs: dict[str, Any],
    skip_cell_cloud: bool,
    skip_clustermaps: bool,
) -> Iterable[tuple[str, Path]]:
    if not skip_cell_cloud:
        yield "cell_cloud_html", outputs["cell_cloud_html"]
    if not skip_clustermaps:
        for matrix, path in outputs["clustermaps"].items():
            yield f"{matrix}_clustermap_png", path


def _file_record(
    path: Path,
    *,
    role: str,
    sha256: bool = True,
    skip_large: bool = True,
) -> dict[str, Any]:
    expanded = Path(path).expanduser()
    exists = expanded.exists()
    is_file = expanded.is_file() if exists else False
    is_dir = expanded.is_dir() if exists else False
    size_bytes = expanded.stat().st_size if is_file else None
    digest = None
    hash_status = "not_requested"
    if sha256 and is_file:
        if skip_large and size_bytes is not None and size_bytes > SMALL_INPUT_HASH_LIMIT_BYTES:
            hash_status = "skipped_large_input"
        else:
            digest = _sha256_path(expanded)
            hash_status = "ok"
    elif sha256 and not is_file:
        hash_status = "not_a_file"
    return {
        "role": role,
        "path": str(expanded),
        "exists": bool(exists),
        "is_file": bool(is_file),
        "is_dir": bool(is_dir),
        "size_bytes": size_bytes,
        "sha256": digest,
        "hash_status": hash_status,
    }


def _with_manifest_record(
    record: dict[str, Any],
    manifest_record: dict[str, Any] | None,
) -> dict[str, Any]:
    if manifest_record is None:
        record["manifest"] = None
        return record
    record["manifest"] = {
        "relative_path": manifest_record.get("relative_path"),
        "accession_or_doi": manifest_record.get("accession_or_doi"),
        "license_or_access_terms": manifest_record.get("license_or_access_terms"),
        "generated_by": manifest_record.get("generated_by"),
        "sha256": manifest_record.get("sha256"),
        "size_bytes": manifest_record.get("size_bytes"),
    }
    return record


def _sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _parse_tokens(values: Sequence[str]) -> list[str]:
    tokens: list[str] = []
    for value in values:
        tokens.extend(part.strip() for part in str(value).split(",") if part.strip())
    return tokens


def _write_manifest(path: Path, manifest: dict[str, Any], *, dry_run: bool) -> None:
    if dry_run:
        return
    path.write_text(json.dumps(manifest, indent=2, ensure_ascii=True, default=str) + "\n", encoding="utf-8")


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _truncate(value: str, limit: int = 20000) -> str:
    if len(value) <= limit:
        return value
    return value[:limit] + f"\n... <truncated {len(value) - limit} chars>"


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
