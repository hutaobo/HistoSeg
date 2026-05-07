---
orphan: true
---

# Polyp 24-Gene 3D Harvest Memo

## Scope

This memo distills the first biological readout from the 32-slice polyp reconstruction after the post-merge HistoSeg 3D validation. It is a draft research artifact for manuscript shaping, not a final causal claim.

## Data Inputs

- Reconstruction scale: 2,785,128 aligned cells across 32 slices and 5 contour structures.
- Cell-cloud validation: `histoseg-3d render-cell-cloud` generated a 300k-point Plotly HTML from `aligned_leiden_3d_cells.parquet` in `Y:\long\spatialpathologist\3D aligment\polyp\histoseg_3d_reconstruction\post_merge_validation_20260503`.
- Spatial module matrices: fraction-inside and signed-distance matrices from `Y:\long\spatialpathologist\3D aligment\polyp\histoseg_3d_reconstruction\gene_overlays\batch_3d_genes_starter_panel`.
- Evidence table: [polyp_24_gene_evidence.tsv](polyp_24_gene_evidence.tsv).

## Working Model

The 24-gene panel resolves four biologically interpretable 3D compartments:

1. Structure 5 is a mixed stromal-immune niche rather than a purely stromal region.
2. Structure 3 is an epithelial differentiation/proliferation compartment.
3. Structure 4 is an epithelial/stem-like compartment.
4. Structure 2 is a macrophage/perivascular-associated compartment.

## Claim 1: Structure 5 Is A Stromal-Immune Niche

Structure 5 contains the strongest embedded signals for stromal matrix and contractile programs alongside immune lineage markers. GREM1, FAP, ACTA2, COL1A1, and COL1A2 are all Structure 5-dominant by top05 fraction-inside and have negative signed distance, placing their high-density voxels inside the contour rather than merely adjacent to it. The same structure also embeds PTPRC, CD3D, CD3E, MS4A1, and CD79A, which makes the niche immunologically mixed.

Key Structure 5 values: GREM1 0.862 / -35 um, FAP 0.793 / -35 um, ACTA2 0.773 / -25 um, COL1A1 0.767 / -25 um, PTPRC 0.872 / -35 um, CD3D 0.873 / -35 um, MS4A1 0.805 / -25 um. Fractions are top05 fraction-inside; distances are signed distance in um.

Interpretation: Structure 5 should be described as a stromal-immune niche with activated fibroblast, ECM, contractile, T-cell, and B-cell components. It is not cleanly separable into a fibroblast-only compartment by this panel.

## Claim 2: Structure 3 Captures Epithelial Differentiation And Proliferation

Structure 3 carries the strongest secretory epithelial signal in MUC2 and a clear proliferation signal in MKI67. OLFM4 is split between Structure 3 and Structure 5, so it is best treated as a differentiation/interface marker rather than forced into one compartment.

Key Structure 3 values: MUC2 0.603 / -10 um, OLFM4 0.550 / -5 um, MKI67 0.331 / 50 um.

Interpretation: Structure 3 is the best candidate for a MUC2/OLFM4/MKI67 epithelial differentiation-proliferation compartment, with OLFM4 marking a boundary state that also touches the Structure 5 niche.

## Claim 3: Structure 4 Captures An Epithelial/Stem-Like Compartment

Structure 4 is dominated by LGR5 and EPCAM. LGR5 has its maximum top05 fraction-inside in Structure 4, and EPCAM also peaks there, supporting a compact epithelial/stem-like readout that is distinct from the MUC2-heavy Structure 3 state.

Key Structure 4 values: LGR5 0.579 / -5 um, EPCAM 0.415 / 10 um.

Interpretation: Structure 4 is a strong candidate for an epithelial/stem-like compartment. In a manuscript figure, it should be paired with Structure 3 to show epithelial state separation rather than treated as an isolated finding.

## Claim 4: Structure 2 Captures A Macrophage/Perivascular-Associated Compartment

Structure 2 has the highest top05 fractions for C1QA, PDGFRA, and CD68. LYZ is lower and more distributed, but it supports the myeloid direction when read alongside CD68 and C1QA. PDGFRA gives the structure a stromal/perivascular axis rather than a pure macrophage label.

Key Structure 2 values: C1QA 0.270 / 30 um, PDGFRA 0.224 / 45 um, CD68 0.157 / 80 um, LYZ 0.148 / 80 um.

Interpretation: Structure 2 should be framed as macrophage/perivascular-associated. The evidence is strongest for C1QA, PDGFRA, and CD68; LYZ is supportive but not dominant.

## Evidence Discipline

External evidence is used as marker-context support, not as causal proof. The TSV evidence notes favor UniProt reviewed entries, Human Protein Atlas gene pages, Reactome/PubMed follow-up routes, and HistoSeg matrix values. For this memo, the decisive spatial claim comes from HistoSeg's top05 fraction-inside and signed-distance matrices, while external resources anchor marker interpretation.

## Caveats

- These are marker-panel interpretations, not cell-type deconvolution outputs.
- Negative signed distance supports embedding inside a structure, but it does not by itself prove direct cell-cell interaction.
- OLFM4 and LYZ should be reported as boundary/distributed markers in the first manuscript draft.
- Structure names remain numeric until paired with histology review or orthogonal annotation.

## Immediate Figure Direction

Use Figure 1 to establish reconstruction scale and validation, Figure 2 to show matrix-level module separation, and Figure 3 to deep-dive into GREM1/Structure 5 as the best nested-hotspot example. See [polyp_24_gene_figure_plan.md](polyp_24_gene_figure_plan.md).
