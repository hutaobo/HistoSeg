from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path
from typing import Sequence

from .multi_structure import MultiStructureContourConfig, run_multi_structure_contours
from .pattern1_isoline import Pattern1IsolineConfig, run_pattern1_isoline


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Run HistoSeg contour analysis workflows.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    pattern1 = subparsers.add_parser("pattern1", help="Generate Pattern1 isoline contours.")
    pattern1.add_argument("--clusters-csv", required=True, help="GraphClust clusters.csv path.")
    pattern1.add_argument("--cells-parquet", required=True, help="Cell table parquet path.")
    pattern1.add_argument("--out-dir", required=True, help="Output directory.")
    pattern1.add_argument(
        "--pattern1-clusters",
        required=True,
        help="Comma-separated Pattern1 cluster IDs, for example 10,23,19.",
    )
    pattern1.add_argument("--tissue-boundary-csv", default=None, help="Optional tissue boundary CSV.")
    pattern1.add_argument("--grid-n", type=int, default=1200, help="Square grid size.")
    pattern1.add_argument("--knn-k", type=int, default=30, help="KNN neighbors.")
    pattern1.add_argument("--smooth-sigma", type=float, default=5.0, help="Gaussian smoothing sigma.")
    pattern1.add_argument("--min-cells-inside", type=int, default=10, help="Minimum target cells per loop.")

    multi = subparsers.add_parser("multi-structure", help="Generate multi-structure contours.")
    multi.add_argument("--clusters-csv", required=True, help="GraphClust clusters.csv path.")
    multi.add_argument("--cells-parquet", required=True, help="Cell table parquet path.")
    multi.add_argument("--out-dir", required=True, help="Output directory.")
    multi.add_argument(
        "--structures-json",
        required=True,
        help="JSON file containing a list of structure specs with structure_name and cluster_ids.",
    )
    multi.add_argument("--bins-x", type=int, default=900, help="Rasterization bins along x.")
    multi.add_argument("--bins-y", type=int, default=700, help="Rasterization bins along y.")
    multi.add_argument("--gaussian-sigma", type=float, default=2.25, help="Density smoothing sigma.")
    multi.add_argument("--min-cells", type=int, default=500, help="Minimum assigned cells per contour.")

    args = parser.parse_args(argv)
    if args.command == "pattern1":
        result = run_pattern1_isoline(
            Pattern1IsolineConfig(
                clusters_csv=args.clusters_csv,
                cells_parquet=args.cells_parquet,
                tissue_boundary_csv=args.tissue_boundary_csv,
                out_dir=args.out_dir,
                pattern1_clusters=_parse_cluster_list(args.pattern1_clusters),
                grid_n=args.grid_n,
                knn_k=args.knn_k,
                smooth_sigma=args.smooth_sigma,
                min_cells_inside=args.min_cells_inside,
            )
        )
    else:
        structures = json.loads(Path(args.structures_json).read_text(encoding="utf-8"))
        result = run_multi_structure_contours(
            MultiStructureContourConfig(
                clusters_csv=args.clusters_csv,
                cells_parquet=args.cells_parquet,
                out_dir=args.out_dir,
                structures=structures,
                bins_x=args.bins_x,
                bins_y=args.bins_y,
                gaussian_sigma=args.gaussian_sigma,
                min_cells=args.min_cells,
            )
        )

    print(json.dumps(asdict(result), indent=2, ensure_ascii=False, default=str))


def _parse_cluster_list(value: str) -> list[str]:
    clusters = [part.strip() for part in value.split(",") if part.strip()]
    if not clusters:
        raise SystemExit("--pattern1-clusters must contain at least one cluster ID.")
    return clusters


if __name__ == "__main__":
    main()
