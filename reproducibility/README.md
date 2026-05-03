# HistoSeg Reproducibility Guide

This directory documents the reproducibility entrypoint for the HistoSeg methods
paper. It defines the environment, provenance checks, and thin CLI wrapper used
to recreate the tutorial-derived figure artifacts.

## Environment Setup

For local paper-figure rendering, use the visualization environment:

```bash
conda env create -f environment-viz.yml
conda activate histoseg-viz
```

The viz environment installs the 3D, visualization, and documentation extras:

```bash
pip install -e ".[threed,viz,docs]"
```

For containerized CPU/headless reproduction, build the Docker `viz` target:

```bash
docker build --target viz -t histoseg:viz .
```

You can also use the compose service:

```bash
docker compose run --rm histoseg-viz --help
```

The `viz` target includes Mesa/Xvfb support for static PyVista-based rendering.
It is intended for CPU/headless workflows and does not require a desktop OpenGL
session.

## Reproduce The 24-Gene Polyp Tutorial

The 24-gene polyp tutorial is the source for the manuscript Figure 3 story:

- Documentation: `docs/tutorials/polyp_24_gene_3d_spatial_modules.md`
- Online methods: `docs/manuscripts/histoseg_online_methods_sdf_alignment.md`
- Figure plan: `docs/manuscripts/figure_plan.md`

The tutorial outputs used for Figure 3 are generated from the 32-slice polyp
stack and the starter-panel batch directory:

```text
gene_overlays/batch_3d_genes_starter_panel/
```

The key figure artifacts are:

```text
batch_gene_status.csv
gene_structure_fraction_inside_matrix.csv
gene_structure_signed_distance_matrix.csv
fraction_inside_top05_spatial_clustermap.png
signed_distance_top05_spatial_clustermap.png
```

The post-merge validation cell-cloud artifact is:

```text
post_merge_validation_20260503/leiden_3d_cells_300k_main.html
```

Raw `.h5ad`, Parquet, and HTML files are not copied into the repository. The
paper package records their paths and derived figure artifacts.

To regenerate local paper artifacts from the validated 32-slice polyp paths,
run:

```bash
python reproducibility/run_paper_pipeline.py
```

By default this writes derived outputs under `reproducibility/results/` and a
manifest at `reproducibility/results_manifest.json`. It renders the 300k-point
cell-cloud HTML and regenerates the `top05` fraction-inside and signed-distance
clustermaps from the existing starter-panel matrices. It does not rerun the
heavier `discover-spatial-modules` step unless explicitly requested:

```bash
python reproducibility/run_paper_pipeline.py --run-discovery
```

Use `--dry-run` to print the exact command plan and provenance payload without
creating files.

## Provenance Verification

HistoSeg stores cell-cloud alignment provenance in AnnData
`.uns["histoseg_3d_alignment"]` when projection is run with `--write-cache`.
The provenance payload includes:

- `alignment_hash`
- `alignment_manifest_schema_version`
- `stack_root`
- `pixel_size_um`
- `manifest`

The `alignment_hash` is the SHA-256 hash of a canonical JSON manifest built
from `aligned_slice_manifest.csv`, hard-alignment summaries, accepted soft TPS
landmarks, RBF settings, and topology settings. To verify that a cached cell
cloud matches the paper data state, compare:

1. The current stack's `aligned_slice_manifest.csv`.
2. The canonical `manifest` stored in `.uns["histoseg_3d_alignment"]`.
3. The stored `alignment_hash`.
4. The regenerated hash from `histoseg.threed.hash_alignment_manifest`.

If the hash differs, regenerate the aligned cell cloud instead of reusing cached
coordinates. The CLI supports strict behavior with `--fail-on-stale-cache` and
cache refresh with `--write-cache`.

## Current Reproduction Boundary

This guide intentionally documents the environment and provenance contract
before adding automation. A later `run_paper_pipeline.py` should generate only
the artifacts named in `docs/manuscripts/figure_plan.md` and should avoid
rewriting validated tutorial source artifacts unless explicitly requested.
