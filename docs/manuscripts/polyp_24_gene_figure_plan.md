# Polyp 24-Gene Figure Plan

## Figure 1: 32-Slice 3D Reconstruction And Cell-Cloud Validation

**Main claim:** HistoSeg reconstructs a full 32-slice polyp stack into a shareable 3D analytical object.

Panels:

- **Fig 1A:** 3D contour surface overview using `docs/_static/threed/polyp/24gene/polyp_3d_structure_surfaces.png`.
- **Fig 1B:** 300k-point Leiden cell-cloud validation from `Y:\long\spatialpathologist\3D aligment\polyp\histoseg_3d_reconstruction\post_merge_validation_20260503\leiden_3d_cells_300k_main.html`.
- **Fig 1C:** Scale annotation: 2,785,128 aligned cells, 32 slices, 5 contour structures, sampled to 300,000 Plotly WebGL points.
- **Fig 1D:** Provenance note: cell coordinates are aligned through the HistoSeg stack-root/cached projection contract; raw Parquet and HTML are referenced, not copied into the repo.

## Figure 2: Spatial Module Quantification

**Main claim:** The 24-gene panel separates structures by both fraction-inside and signed-distance behavior.

Panels:

- **Fig 2A:** Fraction-inside clustermap from `Y:\long\spatialpathologist\3D aligment\polyp\histoseg_3d_reconstruction\post_merge_validation_20260503\fraction_inside_top05_main.png` or tracked tutorial equivalent `docs/_static/threed/polyp/24gene/fraction_inside_top05_spatial_clustermap.png`.
- **Fig 2B:** Signed-distance clustermap from `Y:\long\spatialpathologist\3D aligment\polyp\histoseg_3d_reconstruction\post_merge_validation_20260503\signed_distance_top05_main.png` or tracked tutorial equivalent `docs/_static/threed/polyp/24gene/signed_distance_top05_spatial_clustermap.png`.
- **Fig 2C:** Compartment summary table: Structure 5 stromal-immune, Structure 3 epithelial differentiation/proliferation, Structure 4 epithelial/stem-like, Structure 2 macrophage/perivascular-associated.
- **Fig 2D:** Boundary-marker footnote: OLFM4 and LYZ are distributed signals and should be shown as nuance, not overfit into a single category.

## Figure 3: GREM1 / Structure 5 Nested Hotspot Deep Dive

**Main claim:** GREM1 provides the clearest visual and quantitative example of a nested stromal hotspot inside Structure 5.

Panels:

- **Fig 3A:** Nested 3D hotspot overview using `docs/_static/threed/polyp/24gene/GREM1_nested_3d_hotspot_surfaces.png`.
- **Fig 3B:** Structure 5 focus using `docs/_static/threed/polyp/24gene/GREM1_structure5_focus.png`.
- **Fig 3C:** Quantitative callout: GREM1 top05 fraction-inside for Structure 5 = 0.862; signed distance = -35 um.
- **Fig 3D:** Marker context: GREM1/FAP/ACTA2/COL1A1 share the same Structure 5 niche with PTPRC/CD3D/CD3E/MS4A1/CD79A, motivating the stromal-immune interpretation.

## Suggested Manuscript Sequence

1. Show geometry first so reviewers trust the 3D object.
2. Show matrices second so the compartment claims are visibly quantitative.
3. Use GREM1 third as the hero biological example because it links visual topology, SDF embedding, and marker interpretation.

## Assets To Avoid Copying

Do not copy raw `.h5ad`, Parquet, or HTML files into the repository. Reference large validation outputs by absolute analysis path, and embed only selected raster figures in manuscript/deck artifacts.
