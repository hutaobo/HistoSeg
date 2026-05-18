from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from importlib import metadata
from pathlib import Path
from typing import Any

from histoseg.contour import (
    AutoStructureDiscoveryConfig,
    MultiStructureContourConfig,
    discover_auto_structures,
    resolve_xenium_output_folder,
    run_multi_structure_contours,
)
from histoseg.contour.auto_structure import XENIUM_CELLS_RELPATH, XENIUM_CLUSTERS_RELPATH
from histoseg.threed import (
    LabelFreeBeforeAfterFigureConfig,
    LabelFreeContourAlignmentConfig,
    align_contours_label_free,
    render_label_free_before_after_panel,
)

DEFAULT_FIXED_XENIUM_OUTPUT = (
    r"Y:\long\10X_datasets\Xenium\Xenium_Breast_Cancer"
    r"\Xenium_FFPE_Human_Breast_Cancer_Rep1_outs"
)
DEFAULT_MOVING_XENIUM_OUTPUT = (
    r"Y:\long\10X_datasets\Xenium\Xenium_Breast_Cancer"
    r"\Xenium_FFPE_Human_Breast_Cancer_Rep2"
)


def main() -> None:
    args = _parse_args()
    out_dir = Path(args.out_dir).expanduser().resolve()

    fixed_outs = resolve_xenium_output_folder(args.fixed_xenium_output)
    moving_outs = resolve_xenium_output_folder(args.moving_xenium_output)
    if args.dry_run:
        plan = {
            "dry_run": True,
            "fixed_xenium_output_resolved": str(fixed_outs),
            "moving_xenium_output_resolved": str(moving_outs),
            "out_dir": str(out_dir),
            "steps": [
                "discover fixed auto structures",
                "discover moving auto structures",
                "generate fixed/moving multi-structure GeoJSON contours",
                "run label-free group-correspondence partial-anchor alignment",
                "render manuscript before/after panel",
                "write provenance manifest",
            ],
        }
        print(json.dumps(plan, indent=2, ensure_ascii=False))
        return

    out_dir.mkdir(parents=True, exist_ok=True)
    fixed_auto = _discover_structures(
        xenium_output=fixed_outs,
        out_dir=out_dir / "fixed_auto_structure",
        cluster_count=args.fixed_cluster_count or args.cluster_count,
        args=args,
    )
    moving_auto = _discover_structures(
        xenium_output=moving_outs,
        out_dir=out_dir / "moving_auto_structure",
        cluster_count=args.moving_cluster_count or args.cluster_count,
        args=args,
    )

    fixed_contours = run_multi_structure_contours(
        MultiStructureContourConfig(
            clusters_csv=fixed_outs / XENIUM_CLUSTERS_RELPATH,
            cells_parquet=fixed_outs / XENIUM_CELLS_RELPATH,
            out_dir=out_dir / "fixed_contours",
            structures=fixed_auto.selected_structures,
            bins_x=args.bins_x,
            bins_y=args.bins_y,
            gaussian_sigma=args.gaussian_sigma,
            min_cells=args.min_cells,
        )
    )
    moving_contours = run_multi_structure_contours(
        MultiStructureContourConfig(
            clusters_csv=moving_outs / XENIUM_CLUSTERS_RELPATH,
            cells_parquet=moving_outs / XENIUM_CELLS_RELPATH,
            out_dir=out_dir / "moving_contours",
            structures=moving_auto.selected_structures,
            bins_x=args.bins_x,
            bins_y=args.bins_y,
            gaussian_sigma=args.gaussian_sigma,
            min_cells=args.min_cells,
        )
    )

    alignment = align_contours_label_free(
        LabelFreeContourAlignmentConfig(
            fixed_geojson=fixed_contours.geojson,
            moving_geojson=moving_contours.geojson,
            out_dir=out_dir / "alignment",
            partial_correspondence=True,
            group_correspondence=True,
            group_candidate_count=args.group_candidate_count,
            group_ransac_trials=args.group_ransac_trials,
            group_min_descriptor_score=args.group_min_descriptor_score,
            group_residual_limit_um=args.group_residual_limit_um,
            group_min_component_area_um2=args.group_min_component_area_um2,
            search_window=args.search_window,
            knn_neighbors=args.knn_neighbors,
            min_anchor_score=args.min_anchor_score,
            min_review_score=args.min_review_score,
            min_anchor_count=args.min_anchor_count,
            dpi=args.dpi,
            overwrite=args.overwrite,
        )
    )
    if alignment.aligned_geojson is None:
        raise RuntimeError(f"Label-free alignment did not produce an aligned GeoJSON: {alignment.summary_json}")

    figure = render_label_free_before_after_panel(
        LabelFreeBeforeAfterFigureConfig(
            fixed_geojson=fixed_contours.geojson,
            moving_geojson=moving_contours.geojson,
            aligned_geojson=alignment.aligned_geojson,
            anchors_csv=alignment.landmarks_csv,
            out_png=out_dir / "figure" / "breast_partial_anchor_before_after.png",
            out_svg=out_dir / "figure" / "breast_partial_anchor_before_after.svg",
            dpi=args.figure_dpi,
        )
    )

    manifest = _build_manifest(
        args=args,
        out_dir=out_dir,
        fixed_outs=fixed_outs,
        moving_outs=moving_outs,
        fixed_auto=fixed_auto,
        moving_auto=moving_auto,
        fixed_contours=fixed_contours,
        moving_contours=moving_contours,
        alignment=alignment,
        figure=figure,
    )
    manifest_path = out_dir / "breast_partial_anchor_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    print(json.dumps({"manifest": str(manifest_path), "figure_png": str(figure.out_png)}, indent=2))


def _discover_structures(
    *,
    xenium_output: Path,
    out_dir: Path,
    cluster_count: str,
    args: argparse.Namespace,
):
    return discover_auto_structures(
        AutoStructureDiscoveryConfig(
            xenium_output=xenium_output,
            out_dir=out_dir,
            cluster_count=_parse_cluster_count(cluster_count),
            min_structure_count=args.min_structure_count,
            max_structure_count=args.max_structure_count,
            min_structure_cell_fraction=args.min_structure_cell_fraction,
            linkage_method=args.linkage_method,
            use_cophenetic=not args.no_cophenetic,
            overwrite=args.overwrite,
        )
    )


def _build_manifest(**kwargs: Any) -> dict[str, Any]:
    args = kwargs["args"]
    fixed_auto = kwargs["fixed_auto"]
    moving_auto = kwargs["moving_auto"]
    fixed_contours = kwargs["fixed_contours"]
    moving_contours = kwargs["moving_contours"]
    alignment = kwargs["alignment"]
    figure = kwargs["figure"]
    try:
        histoseg_version = metadata.version("histoseg")
    except Exception:
        histoseg_version = "unknown"
    alignment_summary = json.loads(Path(alignment.summary_json).read_text(encoding="utf-8"))
    return {
        "workflow": "breast_partial_anchor_from_xenium",
        "histoseg_version": histoseg_version,
        "inputs": {
            "fixed_xenium_output": str(kwargs["fixed_outs"]),
            "moving_xenium_output": str(kwargs["moving_outs"]),
        },
        "parameters": vars(args),
        "auto_structure": {
            "fixed": {
                "structures_json": str(fixed_auto.structures_json),
                "summary_json": str(fixed_auto.summary_json),
                "selected_structure_count": fixed_auto.selected_structure_count,
            },
            "moving": {
                "structures_json": str(moving_auto.structures_json),
                "summary_json": str(moving_auto.summary_json),
                "selected_structure_count": moving_auto.selected_structure_count,
            },
        },
        "contours": {
            "fixed": _paths_from_dataclass(fixed_contours),
            "moving": _paths_from_dataclass(moving_contours),
        },
        "alignment": {
            "summary_json": str(alignment.summary_json),
            "aligned_geojson": str(alignment.aligned_geojson),
            "landmarks_csv": str(alignment.landmarks_csv) if alignment.landmarks_csv else None,
            "group_matrix_csv": str(alignment.group_matrix_csv) if alignment.group_matrix_csv else None,
            "method": alignment_summary.get("method"),
            "anchor_transform": alignment_summary.get("anchor_transform"),
            "transform": alignment_summary.get("transform"),
        },
        "figure": _paths_from_dataclass(figure),
    }


def _paths_from_dataclass(obj: Any) -> dict[str, Any]:
    payload = asdict(obj)
    return {key: str(value) if isinstance(value, Path) else value for key, value in payload.items()}


def _parse_cluster_count(value: str | int) -> str | int:
    text = str(value).strip().lower()
    if text == "auto":
        return "auto"
    return int(text)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Reproduce the fully automatic breast partial-anchor alignment figure "
            "starting from two raw Xenium output folders."
        )
    )
    parser.add_argument("--fixed-xenium-output", default=DEFAULT_FIXED_XENIUM_OUTPUT)
    parser.add_argument("--moving-xenium-output", default=DEFAULT_MOVING_XENIUM_OUTPUT)
    parser.add_argument("--out-dir", default="reproducibility/results/breast_partial_anchor_from_xenium")
    parser.add_argument("--cluster-count", default="auto", help="Auto-structure count for both slices, or 'auto'.")
    parser.add_argument("--fixed-cluster-count", default=None, help="Optional fixed-slice override.")
    parser.add_argument("--moving-cluster-count", default=None, help="Optional moving-slice override.")
    parser.add_argument("--min-structure-count", type=int, default=3)
    parser.add_argument("--max-structure-count", type=int, default=8)
    parser.add_argument("--min-structure-cell-fraction", type=float, default=0.01)
    parser.add_argument("--linkage-method", default="average")
    parser.add_argument("--no-cophenetic", action="store_true")
    parser.add_argument("--bins-x", type=int, default=900)
    parser.add_argument("--bins-y", type=int, default=700)
    parser.add_argument("--gaussian-sigma", type=float, default=2.25)
    parser.add_argument("--min-cells", type=int, default=500)
    parser.add_argument("--search-window", type=float, default=800.0)
    parser.add_argument("--knn-neighbors", type=int, default=6)
    parser.add_argument("--min-anchor-score", type=float, default=0.72)
    parser.add_argument("--min-review-score", type=float, default=0.55)
    parser.add_argument("--min-anchor-count", type=int, default=8)
    parser.add_argument("--group-candidate-count", type=int, default=12)
    parser.add_argument("--group-ransac-trials", type=int, default=15000)
    parser.add_argument("--group-min-descriptor-score", type=float, default=0.35)
    parser.add_argument("--group-residual-limit-um", type=float, default=900.0)
    parser.add_argument("--group-min-component-area-um2", type=float, default=5000.0)
    parser.add_argument("--dpi", type=int, default=180)
    parser.add_argument("--figure-dpi", type=int, default=300)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    main()
