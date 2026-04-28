# AI-Driven Spatial Pathologist: Repository Architecture

## Purpose

This document defines an implementation-ready repository architecture for an
AI-driven spatial pathology system that performs:

1. AI-assisted cell type annotation from spatial transcriptomics data
2. cell type relationship modeling through cophenetic distance embedding (CDE)
3. tissue structure discovery and contour extraction
4. pathology structure assignment against a curated knowledge base
5. automated pathology discovery for low-confidence or unmatched structures
6. pathologist-facing overlays, reports, and review workflows

The design is grounded in the current code assets in:

- `src/histoseg/`
- `D:/GitHub/segmentation_methods/src/tissue_structure_pipeline/`

The goal is not to replace the current code immediately. The goal is to define
how the existing prototypes should evolve into a maintainable product-grade
stack.

## Scope

This architecture covers:

- Python package boundaries
- data contracts between stages
- repository layout
- execution entrypoints
- migration plan from current prototype code
- validation and review outputs

This architecture does not prescribe:

- a specific cloud vendor
- a specific frontend framework
- one fixed LLM provider
- a single organ or disease ontology

## Design Principles

1. Geometry and semantics must stay separate.
   Contours answer "where". Pathology assignment answers "what".

2. Every stage must emit confidence and provenance.
   The downstream AI copilot must explain itself from structured evidence.

3. Discovery is a first-class output.
   Low-confidence or unmatched structures are valid scientific findings, not
   pipeline failures.

4. Stable library code and exploratory disease logic move at different speeds.
   Generic algorithms belong in `histoseg`. Dataset-specific experiments can
   remain outside until they stabilize.

5. Reviewability is a product requirement.
   Every case should be exportable to plots, tables, polygons, and a structured
   summary for a pathologist.

## Current Asset Map

The current codebase already contains the core building blocks.

| Capability | Current asset | Role in target architecture |
| --- | --- | --- |
| KNN + Gaussian isoline contouring | `src/histoseg/contour/pattern1_isoline.py` | Generic geometry engine |
| Tissue boundary extraction | `src/histoseg/geometry/tissue_boundary.py` | Boundary and polygon utilities |
| Searcher-Findee distance and cophenetic matrices | `histoseg.contour` utility exports | CDE foundation |
| Structure discovery and semantic assignment | `segmentation_methods/src/tissue_structure_pipeline/structure_assignment.py` | Prototype pathology engine |
| Hierarchical pathology scoring | `segmentation_methods/src/tissue_structure_pipeline/hierarchical_structure_assignment.py` | Reference-based semantic scorer |
| H&E overlay and Xenium export | `segmentation_methods/src/tissue_structure_pipeline/he_overlay.py` | Review and validation layer |
| Reference universe parsing | `segmentation_methods/src/tissue_structure_pipeline/reference_parser.py` | Ontology and knowledge layer |

## Recommended Repository Strategy

Use a two-speed repository strategy.

1. `HistoSeg` remains the reusable library repo.
   It owns stable algorithms, schemas, generic CLIs, and documentation.

2. `segmentation_methods` remains the exploratory methods repo.
   It owns disease-specific scoring logic, reference drafting, and project
   validation until those modules become stable enough to upstream.

3. A future thin application repo can be added later if needed.
   That repo should orchestrate web serving, report rendering, and deployment.
   It should depend on `histoseg` instead of duplicating algorithm code.

This means the shortest implementation path is:

1. expand `HistoSeg` to cover the generic pipeline interfaces
2. keep lung-specific pathology logic in `segmentation_methods`
3. migrate stable pathology modules upstream in phases

## Target `HistoSeg` Package Layout

The recommended package layout inside this repository is:

```text
src/histoseg/
  __init__.py
  schemas/
    __init__.py
    cells.py
    structures.py
    pathology.py
    reports.py
  annotation/
    __init__.py
    marker_rules.py
    llm_assisted.py
    confidence.py
  cde/
    __init__.py
    distance.py
    cophenetic.py
    structure_discovery.py
  segmentation/
    __init__.py
    isoline.py
    partition.py
    overlay.py
    polygon_export.py
  pathology/
    __init__.py
    reference.py
    assignment.py
    discovery.py
    evidence.py
  reports/
    __init__.py
    case_summary.py
    reviewer_exports.py
  io/
    __init__.py
    xenium.py
    huggingface.py
    parquet.py
  cli/
    __init__.py
    main.py
```

## Mapping From Current Modules to Target Modules

| Current module | Target module | Migration note |
| --- | --- | --- |
| `contour/pattern1_isoline.py` | `segmentation/isoline.py` | Keep current implementation, widen from Pattern1 to generic structures |
| `geometry/tissue_boundary.py` | `segmentation/partition.py` or `geometry/` | Keep as shared geometry utility |
| Searcher-Findee utilities | `cde/distance.py` and `cde/cophenetic.py` | Split plotting from computation |
| `segmentation_methods/.../structure_assignment.py` | `pathology/assignment.py` | Migrate only stable generic pieces |
| `segmentation_methods/.../hierarchical_structure_assignment.py` | `pathology/reference.py` and `pathology/assignment.py` | Keep scorer logic, simplify organ-specific heuristics |
| `segmentation_methods/.../he_overlay.py` | `segmentation/overlay.py` and `reports/reviewer_exports.py` | Split geometry, rendering, and export |
| `segmentation_methods/.../reference_parser.py` | `pathology/reference.py` | Upstream the reference universe model |

## Data Contracts

The pipeline should standardize on five core tables.

### 1. `cells.parquet`

One row per spatial cell.

Required fields:

- `cell_id`
- `x_centroid`
- `y_centroid`
- `sample_id`
- `dataset_id`

Optional fields:

- `z_centroid`
- `transcript_count`
- `he_x_px`
- `he_y_px`

### 2. `cell_annotations.parquet`

One row per cell annotation result.

Required fields:

- `cell_id`
- `annotation_label`
- `annotation_level`
- `annotation_confidence`
- `annotation_source`

Recommended evidence fields:

- `marker_support_json`
- `llm_rationale`
- `retrieval_hits_json`
- `model_name`
- `model_version`

### 3. `cluster_relationships.parquet`

Outputs from CDE and structure discovery.

Required fields:

- `cluster_id`
- `cluster_label`
- `top_group`
- `structure_id`
- `distance_profile_json`

Optional fields:

- `row_cophenetic_score`
- `col_cophenetic_score`
- `structure_cut_version`

### 4. `structures.parquet`

One row per discovered structure.

Required fields:

- `structure_id`
- `sample_id`
- `n_cells`
- `cell_type_composition_json`
- `structure_confidence`
- `contour_path`
- `partition_mask_path`

Optional fields:

- `overlay_png_path`
- `xenium_geojson_path`
- `he_overlay_png_path`

### 5. `pathology_assignments.parquet`

One row per structure-level pathology interpretation.

Required fields:

- `structure_id`
- `assigned_label`
- `assigned_reference_id`
- `behavior`
- `assignment_confidence`
- `top_candidates_json`
- `evidence_json`

Optional fields:

- `fallback_applied`
- `fallback_reason`
- `discovery_flag`
- `review_status`

## End-to-End Pipeline

The target pipeline is:

```mermaid
flowchart LR
  A["Spatial transcriptomics inputs"] --> B["Cell annotation"]
  B --> C["Cluster and cell-type normalization"]
  C --> D["Searcher-Findee distance matrix"]
  D --> E["Cophenetic distance embedding"]
  E --> F["Structure discovery"]
  F --> G["Contour and partition extraction"]
  G --> H["Pathology knowledge-base assignment"]
  H --> I["Discovery engine for low-confidence structures"]
  I --> J["Pathologist report and review exports"]
```

## Stage Responsibilities

### Stage 1. Cell annotation

Inputs:

- raw cell expression or pre-clustered labels
- marker references
- optional LLM prompt templates and retrieval context

Outputs:

- one best label per cell
- confidence score
- evidence bundle

Implementation notes:

- start with deterministic marker and atlas matching
- add LLM assistance only for ambiguous or low-confidence cells
- never drop the deterministic evidence when the LLM is used

### Stage 2. CDE and structure discovery

Inputs:

- spatial cell coordinates
- cell labels or cluster labels

Outputs:

- cluster-to-cluster directed distance matrix
- row and column cophenetic matrices
- cluster-to-structure assignment

Implementation notes:

- upstream the fast implementation from `tissue_structure_pipeline/distance_utils.py`
- keep compute and plotting separate
- treat CDE outputs as reusable artifacts, not temporary notebook state

### Stage 3. Geometry and contour extraction

Inputs:

- cells assigned to structures
- optional tissue boundary
- optional H&E alignment

Outputs:

- contour vertices
- partition masks
- Xenium-compatible polygons
- review overlays

Implementation notes:

- preserve the current KNN + Gaussian contouring path
- generalize from a single Pattern1 target to multiple structures
- ensure every structure can be exported as `.npy`, `.geojson`, and image preview

### Stage 4. Pathology assignment

Inputs:

- structure composition
- hierarchical pathology reference

Outputs:

- best pathology label
- top-k alternatives
- score breakdown
- semantic confidence

Implementation notes:

- upstream the reference universe model
- separate generic scoring from organ-specific harmonization rules
- keep lung-specific overrides in project config until generalized

### Stage 5. Discovery

Inputs:

- low-confidence assignments
- structures with high novelty or poor reference fit

Outputs:

- discovery candidate list
- anomaly score
- nearest known structures
- supporting evidence for expert review

Implementation notes:

- discovery is triggered when confidence is below a configurable threshold
- discovery outputs should never be silently merged into known labels
- every discovery candidate should carry its nearest known pathology labels

### Stage 6. Reporting and review

Inputs:

- structure assignments
- overlays
- evidence

Outputs:

- case summary markdown
- reviewer CSV and Parquet bundles
- Xenium Explorer exports
- web-ready JSON payloads

Implementation notes:

- pathologist-facing text should be generated from structured evidence
- reports must separate facts, model suggestions, and open questions

## CLI and Workflow Entry Points

Add a stable CLI layer so every stage can be run independently.

Recommended commands:

```bash
histoseg annotate-cells --config configs/case.json
histoseg compute-cde --config configs/case.json
histoseg discover-structures --config configs/case.json
histoseg build-contours --config configs/case.json
histoseg assign-pathology --config configs/case.json
histoseg run-case --config configs/case.json
histoseg export-review --config configs/case.json
```

Recommended orchestration behavior:

1. each command reads one config file
2. each stage writes durable outputs to disk
3. downstream stages only consume declared artifacts
4. re-running a stage should be idempotent when inputs have not changed

## Configuration Layout

Recommended case config shape:

```json
{
  "dataset_name": "lung_graphclust",
  "inputs": {
    "cells_parquet": "data/cells.parquet",
    "cluster_csv": "data/clusters.csv",
    "cluster_annotation_csv": "data/cluster_annotations.csv",
    "tissue_boundary_csv": "data/tissue_boundary.csv",
    "he_alignment_csv": "data/he_alignment.csv",
    "he_image_tif": "data/he_image.ome.tif"
  },
  "references": {
    "hierarchical_reference_json": "references/lung_reference_hierarchical.json",
    "celltype_harmonization_json": "references/lung_celltype_harmonization.json"
  },
  "outputs": {
    "root_dir": "outputs/lung_graphclust"
  },
  "annotation": {
    "mode": "atlas_plus_llm"
  },
  "cde": {
    "linkage_method": "average"
  },
  "segmentation": {
    "grid_n": 1200,
    "knn_k": 30,
    "smooth_sigma": 5.0
  },
  "pathology": {
    "top_k": 5,
    "discovery_threshold": 0.55
  }
}
```

## Output Layout

Each case should write to one deterministic root.

```text
outputs/<case_id>/
  annotation/
    cell_annotations.parquet
    annotation_summary.json
  cde/
    cluster_distance_matrix.csv
    row_cophenetic_matrix.csv
    col_cophenetic_matrix.csv
    cluster_structure_lookup.csv
  segmentation/
    contour/
    masks/
    previews/
    xenium_explorer_annotations.geojson
  pathology/
    structure_assignments.csv
    structure_assignment_details.json
    discovery_candidates.csv
  reports/
    case_summary.md
    reviewer_bundle.json
```

## Validation Strategy

Validation must exist at three levels.

### Algorithm validation

- contour stability across parameter sweeps
- structure stability across linkage and cut settings
- assignment stability across harmonization rule changes

### Biological validation

- marker agreement for annotated cells
- expected compartment enrichment
- agreement with known pathology structures in benchmark cases

### Review validation

- H&E overlay quality
- Xenium polygon import success
- pathologist agreement and disagreement capture

## Testing Strategy

Add tests in phases.

1. unit tests for schemas and reference parsing
2. unit tests for CDE matrix generation on tiny synthetic datasets
3. golden tests for contour extraction outputs
4. integration tests for a small end-to-end case
5. regression tests for reference-based pathology assignment

## Migration Plan

### Phase 1. Stabilize contracts inside `HistoSeg`

Deliverables:

- add schemas for cells, structures, and pathology assignments
- add CLI entrypoints
- split Searcher-Findee compute from plot utilities

### Phase 2. Generalize contouring

Deliverables:

- convert Pattern1-only contouring into generic structure contouring
- standardize export formats
- add mask and polygon exporters

### Phase 3. Upstream generic pathology components

Deliverables:

- migrate reference universe parsing
- migrate generic hierarchical scorer
- keep organ-specific harmonization in config-driven extension files

### Phase 4. Add discovery and report layers

Deliverables:

- low-confidence structure detector
- discovery candidate ranking
- structured case summary generator

### Phase 5. Add application layer

Deliverables:

- API or web app wrapper
- review status tracking
- user-facing case browser

## Recommended First MVP

The fastest credible MVP is a single-organ lung workflow that:

1. reads Xenium `cells.parquet` and `clusters.csv`
2. imports or computes cell type annotations
3. computes CDE and discovers 4 to 8 structures
4. builds non-overlapping contours and overlays
5. assigns pathology structure labels with confidence and evidence
6. exports a pathologist review bundle

This MVP is aligned with the current `lung_graphclust` project and requires
the least speculative engineering.

## Definition of Done

The repository architecture is considered implemented when:

1. `histoseg run-case --config <case.json>` executes the full workflow
2. every stage writes durable structured artifacts
3. every structure has geometry, semantics, confidence, and provenance
4. unmatched structures are emitted as discovery candidates
5. a reviewer can open overlays and exported polygons without notebook code

## Immediate Next Steps

1. create schema modules and normalize current case outputs to those schemas
2. add `histoseg` CLI entrypoints for CDE, contouring, and report export
3. extract the generic reference universe and scorer from `segmentation_methods`
4. keep lung-specific harmonization and disease rules in external JSON configs
5. add one small benchmark case for end-to-end regression testing
