from __future__ import annotations

import argparse
import json
from dataclasses import asdict

from histoseg.threed import (
    GeneStructureQuantificationConfig,
    quantify_gene_structure_relationships,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Quantify 3D gene hotspots against HistoSeg structure masks.")
    parser.add_argument("--stack-root", required=True)
    parser.add_argument("--gene-density-dir", required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--structures", nargs="*", default=None)
    parser.add_argument("--group-property", default="structure")
    parser.add_argument("--pixel-size-um", type=float, default=0.2125)
    parser.add_argument("--force-rebuild-masks", action="store_true")
    args = parser.parse_args(argv)
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
    print(json.dumps(asdict(result), indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
