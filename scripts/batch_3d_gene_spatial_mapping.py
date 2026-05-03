from __future__ import annotations

import argparse
import json
from dataclasses import asdict

from histoseg.threed import SpatialModuleDiscoveryConfig, run_spatial_module_discovery


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Batch-map genes onto a HistoSeg 3D reconstruction.")
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--aligned-cells-parquet", required=True)
    parser.add_argument("--stack-root", required=True)
    parser.add_argument("--out-dir", default=None)
    parser.add_argument("--genes", nargs="*", default=())
    parser.add_argument("--gene-file", default=None)
    parser.add_argument("--sample-column", default="sample_id")
    parser.add_argument("--barcode-column", default="barcode")
    parser.add_argument("--skip-order-check", action="store_true")
    parser.add_argument("--template-density-summary", default=None)
    parser.add_argument("--xy-voxel-size-um", type=float, default=80.0)
    parser.add_argument("--z-voxel-size-um", type=float, default=5.0)
    parser.add_argument("--smoothing-sigma-xy-um", type=float, default=120.0)
    parser.add_argument("--smoothing-sigma-z-um", type=float, default=10.0)
    parser.add_argument("--surface-smoothing-sigma-voxels", default="1.0,0.9,0.9")
    parser.add_argument("--valid-min-cell-count", type=float, default=8.0)
    parser.add_argument("--min-positive-cells", type=int, default=200)
    parser.add_argument("--min-expression-sum", type=float, default=0.0)
    parser.add_argument("--min-surface-voxels", type=int, default=20)
    parser.add_argument("--mesh-export-formats", default="ply,obj")
    parser.add_argument("--structures", nargs="*", default=None)
    parser.add_argument("--group-property", default="structure")
    parser.add_argument("--pixel-size-um", type=float, default=0.2125)
    parser.add_argument("--force-rebuild-masks", action="store_true")
    args = parser.parse_args(argv)

    surface_sigma = tuple(float(part.strip()) for part in args.surface_smoothing_sigma_voxels.split(",") if part.strip())
    export_formats = tuple(fmt.strip().lower() for fmt in args.mesh_export_formats.split(",") if fmt.strip())
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
            mesh_export_formats=export_formats,
            structures=structures,
            group_property=args.group_property,
            pixel_size_um=args.pixel_size_um,
            force_rebuild_masks=args.force_rebuild_masks,
        )
    )
    print(json.dumps(asdict(result), indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
