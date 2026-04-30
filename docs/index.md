# HistoSeg

::::{div} hero-shell hero-layout
:::{div} hero-copy

### H&E images in. Spatial structure out.

HistoSeg is a Python toolkit for teams working across **histology images**,
**cell-coordinate tables**, and **spatial transcriptomics review workflows**.
It combines image-first segmentation, contour generation, Xenium-ready exports,
and reproducible CLI pipelines in a single package.

- Use `histoseg.he` for H&E foreground extraction, component partitioning, and
  aligned-image change detection.
- Use `histoseg.contour` for Pattern1 isolines, multi-structure contours,
  boundary topology analysis, and review-ready exports.
- Move between Python APIs and CLI commands without changing the underlying
  artifact model.

```{button-ref} installation
:color: primary
:shadow:
Install HistoSeg
```

```{button-ref} usage
:color: secondary
:shadow:
Open the quickstart
```
:::

:::{div} hero-panel

#### What the docs are built for

- Getting from installation to a first successful run quickly.
- Choosing the correct workflow based on the modality of your input data.
- Understanding output artifacts before integrating them downstream.
- Jumping into API details only when you need exact config and result types.

```{note}
Start with **H&E Analysis** when the primary input is an image. Start with
**Contour Analysis** when the primary input is a table of cells with
coordinates and cluster assignments.
```
::::

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

Browse the public interfaces for `histoseg.he` and `histoseg.contour`.

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
:maxdepth: 2
:caption: Getting Started

installation
usage
```

```{toctree}
:hidden:
:maxdepth: 2
:caption: Guides

he_analysis
contour_analysis
contour
workflows/index
```

```{toctree}
:hidden:
:maxdepth: 2
:caption: Tutorials

tutorials/index
```

```{toctree}
:hidden:
:maxdepth: 2
:caption: Reference

api/index
authors
history
contributing
```

```{toctree}
:hidden:
:maxdepth: 1
:caption: Architecture

ai_driven_spatial_pathologist_architecture
```
