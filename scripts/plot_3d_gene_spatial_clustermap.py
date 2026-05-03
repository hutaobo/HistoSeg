from __future__ import annotations

import argparse
import json
from dataclasses import asdict

from histoseg.threed import SpatialModulePlotConfig, plot_spatial_module_clustermap


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Plot clustered 3D gene-structure spatial matrices.")
    parser.add_argument("--batch-dir", required=True)
    parser.add_argument("--hotspot", default="top05", choices=["top15", "top10", "top05"])
    parser.add_argument(
        "--matrix",
        default="fraction_inside",
        choices=["fraction_inside", "overlap_fraction", "signed_distance"],
    )
    parser.add_argument("--out-png", default=None)
    parser.add_argument("--raw-values", action="store_true")
    parser.add_argument("--cmap", default=None)
    args = parser.parse_args(argv)
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
    print(json.dumps(asdict(result), indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
