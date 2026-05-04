from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from typing import Sequence

from . import (
    CODA_METHOD_REFERENCE_DOI,
    CellCloudProjectionConfig,
    CellCloudRenderConfig,
    GeneStructureQuantificationConfig,
    SpatialModuleDiscoveryConfig,
    SpatialModulePlotConfig,
    ThreeDContourReconstructionConfig,
    ThreeDStackReconstructionConfig,
    plot_spatial_module_clustermap,
    quantify_gene_structure_relationships,
    render_cell_cloud_html,
    run_3d_contour_reconstruction,
    run_3d_stack_reconstruction,
    run_cell_cloud_projection,
    run_spatial_module_discovery,
)


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Run HistoSeg 3D Analysis workflows.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    reconstruct = subparsers.add_parser(
        "reconstruct",
        help="Soft-align a hard-aligned moving contour GeoJSON to a fixed slice.",
    )
    reconstruct.add_argument(
        "--fixed-geojson",
        required=True,
        help="Reference contour GeoJSON path.",
    )
    reconstruct.add_argument(
        "--moving-hard-aligned-geojson",
        required=True,
        help="Moving contour GeoJSON after hard/similarity alignment to the reference.",
    )
    reconstruct.add_argument(
        "--out-dir",
        required=True,
        help="Output directory for soft-alignment artifacts.",
    )
    reconstruct.add_argument(
        "--group-property",
        default="structure",
        help="Feature property used to match contour structures.",
    )
    reconstruct.add_argument(
        "--sampling-distance-um",
        type=float,
        default=50.0,
        help="Boundary landmark sampling interval in micrometers.",
    )
    reconstruct.add_argument(
        "--max-landmark-distance-um",
        type=float,
        default=180.0,
        help="Reject source-target landmark pairs farther apart than this distance.",
    )
    reconstruct.add_argument(
        "--landmarks-per-structure",
        type=int,
        default=260,
        help="Conservative cap for boundary landmarks per structure. Use 0 for no cap.",
    )
    reconstruct.add_argument(
        "--diagnostic-structure-landmarks",
        type=int,
        default=620,
        help="Landmark cap for the diagnostic structure. Use 0 to use the regular cap.",
    )
    reconstruct.add_argument(
        "--landmark-candidate-count",
        type=int,
        default=8,
        help="K nearest fixed-boundary candidates for normal-aware landmark matching.",
    )
    reconstruct.add_argument(
        "--landmark-candidate-spacing-um",
        type=float,
        default=None,
        help="Fixed-boundary candidate sampling interval. Defaults to a conservative auto value.",
    )
    reconstruct.add_argument(
        "--landmark-normal-weight-um",
        type=float,
        default=0.0,
        help="Normal-alignment penalty weight. Use 0 for legacy nearest-projection matching.",
    )
    reconstruct.add_argument(
        "--landmark-normal-step-um",
        type=float,
        default=None,
        help="Arc-length step used to estimate boundary normals. Defaults to an auto value.",
    )
    reconstruct.add_argument(
        "--rbf-kernel",
        default="thin_plate_spline",
        help="SciPy RBFInterpolator kernel.",
    )
    reconstruct.add_argument(
        "--rbf-neighbors",
        type=int,
        default=96,
        help="Local RBF neighbor count. Use 0 for a global interpolator.",
    )
    reconstruct.add_argument(
        "--rbf-smoothing",
        default="1e-4",
        help="RBF smoothing parameter, or 'auto' for k-fold CV selection.",
    )
    reconstruct.add_argument(
        "--topology-grid-size",
        type=int,
        default=24,
        help="Sparse grid size for TPS topology checks. Use 0 to disable.",
    )
    reconstruct.add_argument(
        "--topology-min-area-ratio",
        type=float,
        default=0.5,
        help="Reject soft TPS when any checked grid cell shrinks below this area ratio.",
    )
    reconstruct.add_argument(
        "--topology-max-area-ratio",
        type=float,
        default=2.0,
        help="Reject soft TPS when any checked grid cell expands above this area ratio.",
    )
    reconstruct.add_argument(
        "--diagnostic-structure",
        default="Structure 5",
        help="Structure label for the local diagnostic zoom. Use '' to show all structures.",
    )
    reconstruct.add_argument("--dpi", type=int, default=200, help="Preview PNG resolution.")
    reconstruct.add_argument(
        "--no-preview",
        action="store_true",
        help="Skip writing overlay and diagnostic PNGs.",
    )
    reconstruct.add_argument(
        "--curvature-landmark-weight",
        type=float,
        default=0.5,
        help="Fraction of boundary length to oversample at high-curvature regions (0–1).",
    )
    reconstruct.add_argument(
        "--no-mutual-nn-check",
        action="store_true",
        help="Disable mutual nearest-neighbour consistency filter for landmarks.",
    )
    reconstruct.add_argument(
        "--icp-iterations",
        type=int,
        default=2,
        help="Number of ICP-style iterative TPS refinement passes.",
    )
    reconstruct.add_argument(
        "--zero-anchor-count",
        type=int,
        default=16,
        help="Number of zero-displacement anchors placed on the convex-hull perimeter.",
    )
    reconstruct.add_argument(
        "--landmark-outlier-mad-threshold",
        type=float,
        default=3.5,
        help="MAD threshold for landmark outlier rejection. Use 0 to disable.",
    )
    reconstruct.add_argument(
        "--no-per-structure-soft-acceptance",
        action="store_true",
        help="Disable per-structure soft/hard acceptance mixing.",
    )

    stack = subparsers.add_parser(
        "reconstruct-stack",
        help="Build and align a same-sample multi-slice Xenium contour stack.",
    )
    stack.add_argument(
        "--xenium-root",
        required=True,
        help="Folder containing ordered Xenium slice folders.",
    )
    stack.add_argument(
        "--out-dir",
        required=True,
        help="Output directory for the 3D contour stack artifacts.",
    )
    stack.add_argument(
        "--segmentation-strategy",
        default=None,
        help="Text file with one comma-separated cluster list per Structure N.",
    )
    stack.add_argument(
        "--sample-glob",
        default="*",
        help="Glob used to discover slice folders below --xenium-root.",
    )
    stack.add_argument(
        "--z-spacing-um",
        type=float,
        default=5.0,
        help="Physical spacing between consecutive slices.",
    )
    stack.add_argument(
        "--xenium-pixel-size-um",
        type=float,
        default=0.2125,
        help="Microns per Xenium Explorer annotation pixel.",
    )
    stack.add_argument(
        "--group-property",
        default="structure",
        help="GeoJSON property used for structure labels.",
    )
    stack.add_argument(
        "--cluster-column-name",
        default="cluster",
        help="Cluster column name requested from pyXenium read_xenium.",
    )
    stack.add_argument(
        "--clusters-relpath",
        default="analysis/clustering/gene_expression_graphclust/clusters.csv",
        help="Cluster CSV path relative to each Xenium output directory.",
    )
    stack.add_argument(
        "--merged-h5ad",
        default=None,
        help="Optional merged AnnData with unified per-cell clusters for all slices.",
    )
    stack.add_argument(
        "--merged-cluster-column",
        default=None,
        help="Cluster column in merged AnnData obs. Defaults to leiden_1_0 or first leiden* column.",
    )
    stack.add_argument("--merged-sample-column", default="sample_id")
    stack.add_argument("--merged-barcode-column", default="barcode")
    stack.add_argument("--bins-x", type=int, default=900, help="Contour raster bins in X.")
    stack.add_argument("--bins-y", type=int, default=700, help="Contour raster bins in Y.")
    stack.add_argument("--gaussian-sigma", type=float, default=2.25)
    stack.add_argument("--min-cells", type=int, default=500)
    stack.add_argument("--min-component-pixels", type=int, default=180)
    stack.add_argument(
        "--registration-backend",
        default="contour-tps",
        choices=["contour-tps", "coda-image"],
        help=(
            "Hard-alignment seed backend. 'contour-tps' uses the existing contour "
            "similarity seed before TPS. 'coda-image' uses a CODA-inspired tissue-mask "
            f"Radon rotation plus phase-correlation translation seed (DOI: {CODA_METHOD_REFERENCE_DOI}) "
            "before the same topology-safe contour TPS refinement."
        ),
    )
    stack.add_argument("--hard-alignment-maxiter", type=int, default=80)
    stack.add_argument(
        "--coda-raster-size",
        type=int,
        default=512,
        help="Raster size for the CODA-inspired contour tissue-mask proxy.",
    )
    stack.add_argument(
        "--coda-angle-step",
        type=float,
        default=1.0,
        help="Radon angle grid step, in degrees, for --registration-backend coda-image.",
    )
    stack.add_argument(
        "--coda-phase-upsample-factor",
        type=int,
        default=1,
        help="Phase-correlation upsample factor for --registration-backend coda-image.",
    )
    stack.add_argument(
        "--coda-mask-padding-fraction",
        type=float,
        default=0.05,
        help="Fractional square-bounds padding for CODA-inspired mask rasterization.",
    )
    stack.add_argument(
        "--sampling-distance-um",
        type=float,
        default=50.0,
        help="TPS boundary sampling interval in GeoJSON coordinate units.",
    )
    stack.add_argument(
        "--max-landmark-distance-um",
        type=float,
        default=180.0,
        help="TPS nearest-boundary rejection distance in GeoJSON coordinate units.",
    )
    stack.add_argument(
        "--landmarks-per-structure",
        type=int,
        default=260,
        help="Conservative cap for boundary landmarks per structure. Use 0 for no cap.",
    )
    stack.add_argument(
        "--diagnostic-structure-landmarks",
        type=int,
        default=620,
        help="Landmark cap for the diagnostic structure. Use 0 to use the regular cap.",
    )
    stack.add_argument("--landmark-candidate-count", type=int, default=8)
    stack.add_argument("--landmark-candidate-spacing-um", type=float, default=None)
    stack.add_argument("--landmark-normal-weight-um", type=float, default=0.0)
    stack.add_argument("--landmark-normal-step-um", type=float, default=None)
    stack.add_argument("--rbf-neighbors", type=int, default=96)
    stack.add_argument("--rbf-smoothing", default="1e-4", help="RBF smoothing or 'auto'.")
    stack.add_argument("--topology-grid-size", type=int, default=24)
    stack.add_argument("--topology-min-area-ratio", type=float, default=0.5)
    stack.add_argument("--topology-max-area-ratio", type=float, default=2.0)
    stack.add_argument(
        "--diagnostic-structure",
        default="Structure 5",
        help="Structure label for the local diagnostic zoom. Use '' to show all structures.",
    )
    stack.add_argument("--point-sample-distance-um", type=float, default=80.0)
    stack.add_argument("--voxel-size-um", type=float, default=80.0)
    stack.add_argument("--mesh-method", default="marching_cubes")
    stack.add_argument(
        "--mesh-smoothing-sigma-um",
        type=float,
        default=40.0,
        help="3D Gaussian smoothing sigma in microns. Use 0 to disable smoothing.",
    )
    stack.add_argument("--mesh-level", type=float, default=0.5)
    stack.add_argument(
        "--mesh-export-formats",
        default="ply,obj",
        help="Comma-separated mesh export formats. Supported: ply,obj.",
    )
    stack.add_argument(
        "--no-mesh-cleanup",
        action="store_true",
        help="Skip trimesh post-processing before writing PLY/OBJ.",
    )
    stack.add_argument(
        "--min-mesh-component-volume-um3",
        type=float,
        default=None,
        help="Drop disconnected volume components smaller than this volume.",
    )
    stack.add_argument("--mesh-max-faces-for-html", type=int, default=25000)
    stack.add_argument("--dpi", type=int, default=180)
    stack.add_argument(
        "--no-soft",
        action="store_true",
        help="Use only hard similarity alignment between slices.",
    )
    stack.add_argument(
        "--save-slice-preview",
        action="store_true",
        help="Write the per-slice 2D contour preview PNGs.",
    )
    stack.add_argument(
        "--no-alignment-preview",
        action="store_true",
        help="Skip per-pair soft-alignment preview PNGs.",
    )
    stack.add_argument("--overwrite", action="store_true", help="Recompute existing artifacts.")
    stack.add_argument(
        "--no-multistart",
        action="store_true",
        help="Disable multi-start similarity alignment (single PCA-seeded run only).",
    )
    stack.add_argument(
        "--affine-fallback-iou-threshold",
        type=float,
        default=0.0,
        help="Try 6-DOF affine fallback when similarity IoU is below this value. 0 = disabled.",
    )
    stack.add_argument(
        "--global-drift-correction",
        action="store_true",
        help="Apply linear centroid drift correction after the full alignment chain.",
    )
    stack.add_argument(
        "--curvature-landmark-weight",
        type=float,
        default=0.5,
        help="Fraction of boundary to oversample at high-curvature regions (0–1).",
    )
    stack.add_argument(
        "--no-mutual-nn-check",
        action="store_true",
        help="Disable mutual nearest-neighbour landmark consistency filter.",
    )
    stack.add_argument("--icp-iterations", type=int, default=2)
    stack.add_argument("--zero-anchor-count", type=int, default=16)
    stack.add_argument("--landmark-outlier-mad-threshold", type=float, default=3.5)
    stack.add_argument(
        "--no-per-structure-soft-acceptance",
        action="store_true",
        help="Use global soft/hard acceptance instead of per-structure acceptance.",
    )

    project_cells = subparsers.add_parser(
        "project-cells",
        help="Project merged AnnData cells into a HistoSeg 3D reconstruction.",
    )
    project_cells.add_argument("--h5ad", required=True, help="Merged AnnData .h5ad path.")
    project_cells.add_argument("--stack-root", required=True, help="HistoSeg 3D reconstruction output directory.")
    project_cells.add_argument("--out-parquet", required=True, help="Output aligned 3D cell Parquet path.")
    project_cells.add_argument("--sample-column", default="sample_id")
    project_cells.add_argument("--barcode-column", default="barcode")
    project_cells.add_argument("--x-column", default="x_centroid")
    project_cells.add_argument("--y-column", default="y_centroid")
    project_cells.add_argument(
        "--label-columns",
        nargs="*",
        default=(),
        help="Optional obs columns to carry into the output table. Comma-separated chunks are accepted.",
    )
    project_cells.add_argument("--pixel-size-um", type=float, default=0.2125)
    project_cells.add_argument("--chunk-size", type=int, default=100000)
    project_cells.add_argument(
        "--ignore-cache",
        action="store_true",
        help="Ignore AnnData HistoSeg coordinate cache and recompute coordinates.",
    )
    project_cells.add_argument(
        "--fail-on-stale-cache",
        action="store_true",
        help="Raise an error instead of recomputing when cached coordinates have a stale alignment hash.",
    )
    project_cells.add_argument(
        "--write-cache",
        action="store_true",
        help="Write HistoSeg 3D coordinates and provenance back to AnnData.",
    )
    project_cells.add_argument(
        "--cache-h5ad",
        default=None,
        help="Optional output .h5ad for --write-cache. Defaults to overwriting --h5ad.",
    )
    project_cells.add_argument(
        "--write-scanpy-spatial",
        action="store_true",
        help="With --write-cache, also write aligned XY to .obsm['spatial'].",
    )

    render_cells = subparsers.add_parser(
        "render-cell-cloud",
        help="Render an aligned 3D cell cloud as an interactive Plotly HTML.",
    )
    render_cells.add_argument("--stack-root", required=True, help="HistoSeg 3D reconstruction output directory.")
    render_cells.add_argument("--out-html", required=True, help="Output interactive Plotly HTML path.")
    render_cells.add_argument(
        "--aligned-cells-parquet",
        default=None,
        help="Existing aligned 3D cell Parquet. If omitted, --h5ad and --out-parquet are required.",
    )
    render_cells.add_argument("--h5ad", default=None, help="Merged AnnData .h5ad path to project before rendering.")
    render_cells.add_argument(
        "--out-parquet",
        default=None,
        help="Aligned 3D cell Parquet to write when rendering from --h5ad.",
    )
    render_cells.add_argument("--sample-column", default="sample_id")
    render_cells.add_argument("--barcode-column", default="barcode")
    render_cells.add_argument("--x-column", default="x_centroid")
    render_cells.add_argument("--y-column", default="y_centroid")
    render_cells.add_argument("--label-column", default="leiden_1_0")
    render_cells.add_argument(
        "--label-columns",
        nargs="*",
        default=(),
        help="Extra AnnData obs columns to carry when projecting from --h5ad.",
    )
    render_cells.add_argument("--pixel-size-um", type=float, default=0.2125)
    render_cells.add_argument("--chunk-size", type=int, default=100000)
    render_cells.add_argument("--ignore-cache", action="store_true")
    render_cells.add_argument("--fail-on-stale-cache", action="store_true")
    render_cells.add_argument("--write-cache", action="store_true")
    render_cells.add_argument("--cache-h5ad", default=None)
    render_cells.add_argument("--write-scanpy-spatial", action="store_true")
    render_cells.add_argument("--max-points", type=int, default=300000)
    render_cells.add_argument("--random-state", type=int, default=0)
    render_cells.add_argument("--z-visual-scale", type=float, default=8.0)
    render_cells.add_argument("--marker-size", type=float, default=1.4)
    render_cells.add_argument("--opacity", type=float, default=0.48)
    render_cells.add_argument("--no-contours", action="store_true")
    render_cells.add_argument("--title", default="HistoSeg 3D cell cloud")
    render_cells.add_argument("--performance-warning-threshold", type=int, default=500000)

    discover = subparsers.add_parser(
        "discover-spatial-modules",
        help="Batch-map genes into 3D enrichment fields, surfaces, SDF metrics, and spatial module matrices.",
    )
    discover.add_argument("--h5ad", required=True, help="Merged AnnData .h5ad with gene expression.")
    discover.add_argument(
        "--aligned-cells-parquet",
        required=True,
        help="Aligned 3D cell table, usually aligned_leiden_3d_cells.parquet.",
    )
    discover.add_argument("--stack-root", required=True, help="HistoSeg 3D reconstruction output directory.")
    discover.add_argument("--out-dir", default=None, help="Batch output directory.")
    discover.add_argument("--genes", nargs="*", default=(), help="Gene names or comma-separated gene chunks.")
    discover.add_argument("--gene-file", default=None, help="Optional newline-delimited gene list.")
    discover.add_argument("--sample-column", default="sample_id")
    discover.add_argument("--barcode-column", default="barcode")
    discover.add_argument("--skip-order-check", action="store_true")
    discover.add_argument("--template-density-summary", default=None)
    discover.add_argument("--xy-voxel-size-um", type=float, default=80.0)
    discover.add_argument("--z-voxel-size-um", type=float, default=5.0)
    discover.add_argument("--smoothing-sigma-xy-um", type=float, default=120.0)
    discover.add_argument("--smoothing-sigma-z-um", type=float, default=10.0)
    discover.add_argument("--surface-smoothing-sigma-voxels", default="1.0,0.9,0.9")
    discover.add_argument("--valid-min-cell-count", type=float, default=8.0)
    discover.add_argument("--min-positive-cells", type=int, default=200)
    discover.add_argument("--min-expression-sum", type=float, default=0.0)
    discover.add_argument("--min-surface-voxels", type=int, default=20)
    discover.add_argument("--mesh-export-formats", default="ply,obj")
    discover.add_argument("--structures", nargs="*", default=None)
    discover.add_argument("--group-property", default="structure")
    discover.add_argument("--pixel-size-um", type=float, default=0.2125)
    discover.add_argument("--force-rebuild-masks", action="store_true")

    quantify = subparsers.add_parser(
        "quantify-gene-structure",
        help="Quantify an existing gene density output against 3D structure masks.",
    )
    quantify.add_argument("--stack-root", required=True)
    quantify.add_argument("--gene-density-dir", required=True)
    quantify.add_argument("--gene", required=True)
    quantify.add_argument("--structures", nargs="*", default=None)
    quantify.add_argument("--group-property", default="structure")
    quantify.add_argument("--pixel-size-um", type=float, default=0.2125)
    quantify.add_argument("--force-rebuild-masks", action="store_true")

    plot = subparsers.add_parser(
        "plot-spatial-modules",
        help="Plot clustered gene-by-structure heatmaps from spatial module matrices.",
    )
    plot.add_argument("--batch-dir", required=True)
    plot.add_argument("--hotspot", default="top05", choices=["top15", "top10", "top05"])
    plot.add_argument(
        "--matrix",
        default="fraction_inside",
        choices=["fraction_inside", "overlap_fraction", "signed_distance"],
    )
    plot.add_argument("--out-png", default=None)
    plot.add_argument("--raw-values", action="store_true", help="Plot raw values instead of row z-scores.")
    plot.add_argument("--cmap", default=None)

    args = parser.parse_args(argv)
    if args.command == "reconstruct":
        rbf_smoothing_reconstruct = (
            args.rbf_smoothing if args.rbf_smoothing == "auto" else float(args.rbf_smoothing)
        )
        result = run_3d_contour_reconstruction(
            ThreeDContourReconstructionConfig(
                fixed_geojson=args.fixed_geojson,
                moving_hard_aligned_geojson=args.moving_hard_aligned_geojson,
                out_dir=args.out_dir,
                group_property=args.group_property,
                sampling_distance_um=args.sampling_distance_um,
                max_landmark_distance_um=args.max_landmark_distance_um,
                landmarks_per_structure=_zero_as_none(args.landmarks_per_structure),
                diagnostic_structure_landmarks=_zero_as_none(args.diagnostic_structure_landmarks),
                landmark_candidate_count=args.landmark_candidate_count,
                landmark_candidate_spacing_um=args.landmark_candidate_spacing_um,
                landmark_normal_weight_um=args.landmark_normal_weight_um,
                landmark_normal_step_um=args.landmark_normal_step_um,
                rbf_kernel=args.rbf_kernel,
                rbf_neighbors=args.rbf_neighbors or None,
                rbf_smoothing=rbf_smoothing_reconstruct,
                topology_grid_size=args.topology_grid_size,
                topology_min_area_ratio=args.topology_min_area_ratio,
                topology_max_area_ratio=args.topology_max_area_ratio,
                diagnostic_structure=args.diagnostic_structure or None,
                dpi=args.dpi,
                save_preview_png=not args.no_preview,
                curvature_landmark_weight=args.curvature_landmark_weight,
                mutual_nn_check=not args.no_mutual_nn_check,
                icp_iterations=args.icp_iterations,
                zero_anchor_count=args.zero_anchor_count,
                landmark_outlier_mad_threshold=(
                    args.landmark_outlier_mad_threshold
                    if args.landmark_outlier_mad_threshold > 0
                    else None
                ),
                per_structure_soft_acceptance=not args.no_per_structure_soft_acceptance,
            )
        )
        print(json.dumps(asdict(result), indent=2, ensure_ascii=False, default=str))
    elif args.command == "project-cells":
        result = run_cell_cloud_projection(
            CellCloudProjectionConfig(
                h5ad=args.h5ad,
                stack_root=args.stack_root,
                out_parquet=args.out_parquet,
                sample_column=args.sample_column,
                barcode_column=args.barcode_column,
                x_column=args.x_column,
                y_column=args.y_column,
                label_columns=tuple(_parse_cli_chunks(args.label_columns)),
                pixel_size_um=args.pixel_size_um,
                chunk_size=args.chunk_size,
                ignore_cache=args.ignore_cache,
                fail_on_stale_cache=args.fail_on_stale_cache,
                write_cache=args.write_cache,
                cache_h5ad=args.cache_h5ad,
                write_scanpy_spatial=args.write_scanpy_spatial,
            )
        )
        print(json.dumps(asdict(result), indent=2, ensure_ascii=False, default=str))
    elif args.command == "render-cell-cloud":
        aligned_cells_parquet = args.aligned_cells_parquet
        if aligned_cells_parquet:
            if args.h5ad or args.out_parquet:
                parser.error("--aligned-cells-parquet cannot be combined with --h5ad or --out-parquet.")
        else:
            if not args.h5ad or not args.out_parquet:
                parser.error("render-cell-cloud requires either --aligned-cells-parquet or both --h5ad and --out-parquet.")
            label_columns = tuple(dict.fromkeys([args.label_column, *_parse_cli_chunks(args.label_columns)]))
            projection = run_cell_cloud_projection(
                CellCloudProjectionConfig(
                    h5ad=args.h5ad,
                    stack_root=args.stack_root,
                    out_parquet=args.out_parquet,
                    sample_column=args.sample_column,
                    barcode_column=args.barcode_column,
                    x_column=args.x_column,
                    y_column=args.y_column,
                    label_columns=label_columns,
                    pixel_size_um=args.pixel_size_um,
                    chunk_size=args.chunk_size,
                    ignore_cache=args.ignore_cache,
                    fail_on_stale_cache=args.fail_on_stale_cache,
                    write_cache=args.write_cache,
                    cache_h5ad=args.cache_h5ad,
                    write_scanpy_spatial=args.write_scanpy_spatial,
                )
            )
            aligned_cells_parquet = projection.out_parquet

        result = render_cell_cloud_html(
            CellCloudRenderConfig(
                aligned_cells_parquet=aligned_cells_parquet,
                stack_root=args.stack_root,
                out_html=args.out_html,
                label_column=args.label_column,
                max_points=args.max_points,
                random_state=args.random_state,
                z_visual_scale=args.z_visual_scale,
                marker_size=args.marker_size,
                opacity=args.opacity,
                include_contours=not args.no_contours,
                title=args.title,
                performance_warning_threshold=args.performance_warning_threshold,
            )
        )
        print(json.dumps(asdict(result), indent=2, ensure_ascii=False, default=str))
    elif args.command == "discover-spatial-modules":
        surface_sigma = tuple(
            float(part.strip())
            for part in args.surface_smoothing_sigma_voxels.split(",")
            if part.strip()
        )
        mesh_export_formats = tuple(
            fmt.strip().lower()
            for fmt in args.mesh_export_formats.split(",")
            if fmt.strip()
        )
        structures = tuple(args.structures) if args.structures else (
            "Structure 1",
            "Structure 2",
            "Structure 3",
            "Structure 4",
            "Structure 5",
        )
        result = run_spatial_module_discovery(
            SpatialModuleDiscoveryConfig(
                h5ad=args.h5ad,
                aligned_cells_parquet=args.aligned_cells_parquet,
                stack_root=args.stack_root,
                out_dir=args.out_dir,
                genes=tuple(args.genes),
                gene_file=args.gene_file,
                sample_column=args.sample_column,
                barcode_column=args.barcode_column,
                skip_order_check=args.skip_order_check,
                template_density_summary=args.template_density_summary,
                xy_voxel_size_um=args.xy_voxel_size_um,
                z_voxel_size_um=args.z_voxel_size_um,
                smoothing_sigma_xy_um=args.smoothing_sigma_xy_um,
                smoothing_sigma_z_um=args.smoothing_sigma_z_um,
                surface_smoothing_sigma_voxels_zyx=surface_sigma,
                valid_min_cell_count=args.valid_min_cell_count,
                min_positive_cells=args.min_positive_cells,
                min_expression_sum=args.min_expression_sum,
                min_surface_voxels=args.min_surface_voxels,
                mesh_export_formats=mesh_export_formats,
                structures=structures,
                group_property=args.group_property,
                pixel_size_um=args.pixel_size_um,
                force_rebuild_masks=args.force_rebuild_masks,
            )
        )
        print(json.dumps(asdict(result), indent=2, ensure_ascii=False, default=str))
    elif args.command == "quantify-gene-structure":
        structures = tuple(args.structures) if args.structures else (
            "Structure 1",
            "Structure 2",
            "Structure 3",
            "Structure 4",
            "Structure 5",
        )
        result = quantify_gene_structure_relationships(
            GeneStructureQuantificationConfig(
                stack_root=args.stack_root,
                gene_density_dir=args.gene_density_dir,
                gene=args.gene,
                structures=structures,
                group_property=args.group_property,
                pixel_size_um=args.pixel_size_um,
                force_rebuild_masks=args.force_rebuild_masks,
            )
        )
        print(json.dumps(asdict(result), indent=2, ensure_ascii=False, default=str))
    elif args.command == "plot-spatial-modules":
        result = plot_spatial_module_clustermap(
            SpatialModulePlotConfig(
                batch_dir=args.batch_dir,
                hotspot=args.hotspot,
                matrix=args.matrix,
                out_png=args.out_png,
                row_zscore=not args.raw_values,
                cmap=args.cmap,
            )
        )
        print(json.dumps(asdict(result), indent=2, ensure_ascii=False, default=str))
    elif args.command == "reconstruct-stack":
        mesh_smoothing_sigma_um = args.mesh_smoothing_sigma_um
        if mesh_smoothing_sigma_um < 0:
            parser.error("--mesh-smoothing-sigma-um must be >= 0. Use 0 to disable smoothing.")
        if mesh_smoothing_sigma_um == 0:
            mesh_smoothing_sigma_um = None
        mesh_export_formats = tuple(
            fmt.strip().lower()
            for fmt in args.mesh_export_formats.split(",")
            if fmt.strip()
        )
        rbf_smoothing_stack = (
            args.rbf_smoothing if args.rbf_smoothing == "auto" else float(args.rbf_smoothing)
        )
        result = run_3d_stack_reconstruction(
            ThreeDStackReconstructionConfig(
                xenium_root=args.xenium_root,
                out_dir=args.out_dir,
                segmentation_strategy=args.segmentation_strategy,
                sample_glob=args.sample_glob,
                z_spacing_um=args.z_spacing_um,
                xenium_pixel_size_um=args.xenium_pixel_size_um,
                group_property=args.group_property,
                cluster_column_name=args.cluster_column_name,
                clusters_relpath=args.clusters_relpath,
                merged_h5ad=args.merged_h5ad,
                merged_cluster_column=args.merged_cluster_column,
                merged_sample_column=args.merged_sample_column,
                merged_barcode_column=args.merged_barcode_column,
                bins_x=args.bins_x,
                bins_y=args.bins_y,
                gaussian_sigma=args.gaussian_sigma,
                min_cells=args.min_cells,
                min_component_pixels=args.min_component_pixels,
                save_slice_preview_png=args.save_slice_preview,
                registration_backend=args.registration_backend,
                hard_alignment_maxiter=args.hard_alignment_maxiter,
                hard_alignment_multistart=not args.no_multistart,
                affine_fallback_iou_threshold=args.affine_fallback_iou_threshold,
                coda_raster_size=args.coda_raster_size,
                coda_angle_step=args.coda_angle_step,
                coda_phase_upsample_factor=args.coda_phase_upsample_factor,
                coda_mask_padding_fraction=args.coda_mask_padding_fraction,
                global_drift_correction=args.global_drift_correction,
                run_soft_alignment=not args.no_soft,
                sampling_distance_um=args.sampling_distance_um,
                max_landmark_distance_um=args.max_landmark_distance_um,
                landmarks_per_structure=_zero_as_none(args.landmarks_per_structure),
                diagnostic_structure_landmarks=_zero_as_none(args.diagnostic_structure_landmarks),
                landmark_candidate_count=args.landmark_candidate_count,
                landmark_candidate_spacing_um=args.landmark_candidate_spacing_um,
                landmark_normal_weight_um=args.landmark_normal_weight_um,
                landmark_normal_step_um=args.landmark_normal_step_um,
                rbf_neighbors=args.rbf_neighbors or None,
                rbf_smoothing=rbf_smoothing_stack,
                topology_grid_size=args.topology_grid_size,
                topology_min_area_ratio=args.topology_min_area_ratio,
                topology_max_area_ratio=args.topology_max_area_ratio,
                diagnostic_structure=args.diagnostic_structure or None,
                save_alignment_preview_png=not args.no_alignment_preview,
                curvature_landmark_weight=args.curvature_landmark_weight,
                mutual_nn_check=not args.no_mutual_nn_check,
                icp_iterations=args.icp_iterations,
                zero_anchor_count=args.zero_anchor_count,
                landmark_outlier_mad_threshold=(
                    args.landmark_outlier_mad_threshold
                    if args.landmark_outlier_mad_threshold > 0
                    else None
                ),
                per_structure_soft_acceptance=not args.no_per_structure_soft_acceptance,
                point_sample_distance_um=args.point_sample_distance_um,
                voxel_size_um=args.voxel_size_um,
                mesh_method=args.mesh_method,
                mesh_smoothing_sigma_um=mesh_smoothing_sigma_um,
                mesh_level=args.mesh_level,
                mesh_export_formats=mesh_export_formats,
                mesh_cleanup=not args.no_mesh_cleanup,
                min_mesh_component_volume_um3=args.min_mesh_component_volume_um3,
                mesh_max_faces_for_html=args.mesh_max_faces_for_html,
                overwrite=args.overwrite,
                dpi=args.dpi,
            )
        )
        print(json.dumps(asdict(result), indent=2, ensure_ascii=False, default=str))


def _parse_cli_chunks(values: Sequence[str]) -> list[str]:
    result: list[str] = []
    for value in values:
        result.extend(part.strip() for part in str(value).split(",") if part.strip())
    return result


def _zero_as_none(value: int | None) -> int | None:
    if value == 0:
        return None
    return value


if __name__ == "__main__":
    main()
