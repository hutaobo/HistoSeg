# HistoSeg

```{toctree}
:hidden:
:maxdepth: 2

installation
usage
he_analysis
contour_analysis
3d_analysis
workflows/index
tutorials/index
api/index
history
```

<div class="hs-hero">
  <p class="hs-tagline">H&E imaging, 2D spatial contours, and 3D contour alignment in one documentation family.</p>
  <p class="hs-lead">
    HistoSeg is a Python toolkit for histology images, cell-coordinate tables,
    and review-ready spatial transcriptomics exports. This site is organized
    around three primary surfaces, <code>histoseg.he</code>,
    <code>histoseg.contour</code>, and <code>histoseg.threed</code>, then connects them to quickstarts,
    workflow guides, tutorials, and API reference pages.
  </p>
  <div class="hs-link-row">
    <a href="usage.html">Quickstart</a>
    <a href="installation.html">Installation</a>
    <a href="he_analysis.html">H&amp;E guide</a>
    <a href="contour_analysis.html">Contour guide</a>
    <a href="3d_analysis.html">3D guide</a>
    <a href="https://github.com/hutaobo/HistoSeg">GitHub</a>
    <a href="https://pypi.org/project/histoseg/">PyPI</a>
  </div>
</div>

:::{div} hs-metadata
<span>Image-first H&amp;E segmentation</span>
<span>Pattern1 and multi-structure contours</span>
<span>TPS soft alignment for 3D contour workflows</span>
<span>Xenium Explorer-ready exports</span>
<span>Python API and CLI parity</span>
:::

```{note}
Start with **H&E Analysis** when the primary input is an image. Start with
**Contour Analysis** when the primary input is a table of cells with
coordinates and cluster assignments. Start with **3D Analysis** when aligning
neighboring Xenium contour slices from the same sample.
```

## Start Here

::::{grid} 1 2 2 4
:gutter: 2
:class-container: link-grid

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

Run the first H&E workflow or contour workflow from Python or the command line.

+++
First commands and code snippets
:::

:::{grid-item-card} Tutorials
:link: tutorials/index
:link-type: doc

Follow guided walkthroughs when you want to see a full workflow end to end.

+++
Hands-on examples
:::

:::{grid-item-card} API Reference
:link: api/index
:link-type: doc

Browse the public interfaces for `histoseg.he`, `histoseg.contour`, and
`histoseg.threed`.

+++
Configs, result types, and runners
:::
::::

## Explore By Workflow

::::{grid} 1 1 2 2
:gutter: 2
:class-container: link-grid

:::{grid-item-card} H&E Analysis
:link: he_analysis
:link-type: doc

Segment tissue, partition neutral components, and detect changes between
aligned H&E images.
:::

:::{grid-item-card} Contour Analysis
:link: contour_analysis
:link-type: doc

Extract Pattern1 isolines, build named structure partitions, and export
Xenium-ready review layers.
:::

:::{grid-item-card} 3D Analysis
:link: 3d_analysis
:link-type: doc

Soft-align neighboring Xenium contour slices before multi-slice 3D
reconstruction.
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
