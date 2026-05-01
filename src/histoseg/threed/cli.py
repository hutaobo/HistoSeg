from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from typing import Sequence

from . import (
    ThreeDContourReconstructionConfig,
    ThreeDStackReconstructionConfig,
    run_3d_contour_reconstruction,
    run_3d_stack_reconstruction,
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
        type=float,
        default=1e-4,
        help="RBF smoothing parameter.",
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
    stack.add_argument("--hard-alignment-maxiter", type=int, default=80)
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
    stack.add_argument("--landmarks-per-structure", type=int, default=260)
    stack.add_argument("--diagnostic-structure-landmarks", type=int, default=620)
    stack.add_argument("--rbf-neighbors", type=int, default=96)
    stack.add_argument("--rbf-smoothing", type=float, default=1e-4)
    stack.add_argument("--diagnostic-structure", default="Structure 5")
    stack.add_argument("--point-sample-distance-um", type=float, default=80.0)
    stack.add_argument("--voxel-size-um", type=float, default=80.0)
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

    args = parser.parse_args(argv)
    if args.command == "reconstruct":
        result = run_3d_contour_reconstruction(
            ThreeDContourReconstructionConfig(
                fixed_geojson=args.fixed_geojson,
                moving_hard_aligned_geojson=args.moving_hard_aligned_geojson,
                out_dir=args.out_dir,
                group_property=args.group_property,
                sampling_distance_um=args.sampling_distance_um,
                max_landmark_distance_um=args.max_landmark_distance_um,
                landmarks_per_structure=args.landmarks_per_structure or None,
                diagnostic_structure_landmarks=args.diagnostic_structure_landmarks or None,
                rbf_kernel=args.rbf_kernel,
                rbf_neighbors=args.rbf_neighbors or None,
                rbf_smoothing=args.rbf_smoothing,
                diagnostic_structure=args.diagnostic_structure or None,
                dpi=args.dpi,
                save_preview_png=not args.no_preview,
            )
        )
        print(json.dumps(asdict(result), indent=2, ensure_ascii=False, default=str))
    elif args.command == "reconstruct-stack":
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
                hard_alignment_maxiter=args.hard_alignment_maxiter,
                run_soft_alignment=not args.no_soft,
                sampling_distance_um=args.sampling_distance_um,
                max_landmark_distance_um=args.max_landmark_distance_um,
                landmarks_per_structure=args.landmarks_per_structure or None,
                diagnostic_structure_landmarks=args.diagnostic_structure_landmarks or None,
                rbf_neighbors=args.rbf_neighbors or None,
                rbf_smoothing=args.rbf_smoothing,
                diagnostic_structure=args.diagnostic_structure or None,
                save_alignment_preview_png=not args.no_alignment_preview,
                point_sample_distance_um=args.point_sample_distance_um,
                voxel_size_um=args.voxel_size_um,
                mesh_max_faces_for_html=args.mesh_max_faces_for_html,
                overwrite=args.overwrite,
                dpi=args.dpi,
            )
        )
        print(json.dumps(asdict(result), indent=2, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()
