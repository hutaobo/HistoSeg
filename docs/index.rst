HistoSeg
========

HistoSeg is a Python toolkit for H&E image analysis and spatial contour
analysis. Its public API is organized into two feature groups:

* **HE Analysis** (``histoseg.he``): image-based H&E tissue segmentation,
  neutral tissue component partitioning, and aligned-image change detection.
* **Contour Analysis** (``histoseg.contour``): contour extraction from
  spatial/cell-coordinate data, including Pattern1 isolines and
  multi-structure Xenium exports.

.. toctree::
   :maxdepth: 2
   :caption: HE Analysis

   he_analysis

.. toctree::
   :maxdepth: 2
   :caption: Contour Analysis

   contour_analysis
   contour
   tutorials/pattern1_isoline
   workflows/index

.. toctree::
   :maxdepth: 2
   :caption: Reference

   installation
   usage
   api/index
   authors
   contributing
   history

.. toctree::
   :maxdepth: 1
   :caption: Architecture

   ai_driven_spatial_pathologist_architecture
