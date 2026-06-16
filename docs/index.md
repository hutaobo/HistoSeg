# HistoSeg

```{toctree}
:hidden:
:maxdepth: 2

installation
usage
3d_analysis
label_free_alignment
tear_aware_alignment
contour_analysis
he_analysis
workflows/index
tutorials/index
api/index
manuscripts/breast_partial_anchor_alignment_methods
history
```

<div class="hs-hero">
  <p class="hs-tagline">StructureMap-guided semantic contours and 3D Xenium tissue reconstruction.</p>
  <p class="hs-lead">
    HistoSeg is a Python toolkit centered on same-sample, multi-slice Xenium
    contour reconstruction. It turns Cell-GPS Search-and-Find/StructureMap-guided
    or curated structure groups into continuous semantic contours, then aligns
    neighboring contour slices, builds 3D contour stacks, exports PLY/OBJ
    meshes, and writes QC artifacts. <code>histoseg.contour</code> and
    <code>histoseg.he</code> support the 2D spatial and image workflows that
    prepare those structures for reconstruction.
  </p>
  <div class="hs-link-row">
    <a href="3d_analysis.html">3D reconstruction</a>
    <a href="usage.html">Quickstart</a>
    <a href="installation.html">Installation</a>
    <a href="contour_analysis.html">Contour guide</a>
    <a href="ai_driven_spatial_pathologist_serve_app.html">Serve app</a>
    <a href="he_analysis.html">H&amp;E guide</a>
    <a href="https://github.com/hutaobo/HistoSeg">GitHub</a>
    <a href="https://pypi.org/project/histoseg/">PyPI</a>
  </div>
</div>

:::{div} hs-metadata
<span>Same-sample 3D Xenium reconstruction</span>
<span>TPS soft alignment for neighboring contour slices</span>
<span>PLY/OBJ mesh export and QC reports</span>
<span>StructureMap-guided semantic 2D contours</span>
<span>Xenium Explorer-ready exports</span>
:::

```{note}
Start with **3D Reconstruction** when the goal is a multi-slice Xenium contour
stack or 3D mesh. Use **Contour Analysis** to generate 2D structures from
cell-coordinate data, and **H&E Analysis** when the primary input is an image.
```

## Start Here

::::{grid} 1 2 2 4
:gutter: 2
:class-container: link-grid

:::{grid-item-card} 3D Reconstruction
:link: 3d_analysis
:link-type: doc

Align neighboring Xenium contour slices, reconstruct 3D contour stacks, and
export meshes plus QC artifacts.

+++
Flagship workflow
:::

:::{grid-item-card} Installation
:link: installation
:link-type: doc

Set up the package from PyPI or prepare an editable development environment.

+++
Package setup and extras
:::

:::{grid-item-card} Quickstart
:link: usage
:link-type: doc

Run the first 3D reconstruction, contour, or H&E workflow from Python or the
command line.

+++
First commands and code snippets
:::

:::{grid-item-card} Tutorials
:link: tutorials/index
:link-type: doc

Follow the 3D stack reconstruction and soft-alignment walkthroughs end to end.

+++
Hands-on examples
:::

:::{grid-item-card} API Reference
:link: api/index
:link-type: doc

Browse the public interfaces for `histoseg.threed`, `histoseg.contour`, and
`histoseg.he`.

+++
Configs, result types, and runners
:::
::::

## Explore By Workflow

::::{grid} 1 1 2 2
:gutter: 2
:class-container: link-grid

:::{grid-item-card} 3D Reconstruction
:link: 3d_analysis
:link-type: doc

Soft-align neighboring Xenium contour slices and build multi-slice 3D contour
meshes.
:::

:::{grid-item-card} Automatic Partial-Overlap Breast Alignment
:link: label_free_alignment
:link-type: doc

Inspect a fully automatic two-slice breast contour alignment where local
landmarks are found without manual annotation.
:::

:::{grid-item-card} Tear/Missing-Aware Breast Alignment
:link: tear_aware_alignment
:link-type: doc

Follow a B1_2/B1_3 breast pair from raw cell points through Cell-GPS-guided
contours to contour-guided tear closure.
:::

:::{grid-item-card} Contour Analysis
:link: contour_analysis
:link-type: doc

Extract Pattern1 isolines, build StructureMap-guided named structure partitions,
and export Xenium-ready review layers.
:::

:::{grid-item-card} Spatial Pathologist Serve App
:link: ai_driven_spatial_pathologist_serve_app
:link-type: doc

Run the Gradio web workflow for dendrogram-guided multi-structure contour
selection from Xenium cell-coordinate and cluster tables.
:::

:::{grid-item-card} H&E Analysis
:link: he_analysis
:link-type: doc

Segment tissue, partition neutral components, and detect changes between
aligned H&E images.
:::

:::{grid-item-card} Pattern1 Isoline Guide
:link: contour
:link-type: doc

See the contour-generation steps, expected inputs, and artifact set for the
single-structure Pattern1 pipeline.
:::

:::{grid-item-card} Workflow Guides
:link: workflows/index
:link-type: doc

Use workflow-oriented guides when you need a more operational view of a run and
its outputs.
:::
::::

## Project And Internals

::::{grid} 1 1 3 3
:gutter: 2
:class-container: link-grid

:::{grid-item-card} Architecture
:link: ai_driven_spatial_pathologist_architecture
:link-type: doc

Understand how the contour and structure-analysis building blocks fit together.
:::

:::{grid-item-card} Contributing
:link: contributing
:link-type: doc

See the expectations for issues, pull requests, and reproducible bug reports.
:::

:::{grid-item-card} History
:link: history
:link-type: doc

Track the package direction and the public feature groups exposed today.
:::
::::

```{toctree}
:hidden:
:maxdepth: 1
:caption: Additional Guides

contour
ai_driven_spatial_pathologist_serve_app
```

```{toctree}
:hidden:
:maxdepth: 1
:caption: Reference Details

authors
contributing
```

```{toctree}
:hidden:
:maxdepth: 1
:caption: Architecture

ai_driven_spatial_pathologist_architecture
```
