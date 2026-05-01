============
Installation
============

.. raw:: html

   <div class="hs-section-kicker">Getting Started</div>
   <div class="hs-page-intro">HistoSeg ships as a Python package for 3D Xenium contour reconstruction, with optional extras for pyXenium-backed stack reconstruction and local MedSAM-backed H&amp;E workflows. Install the 3D extra when you want the flagship reconstruction pipeline.</div>

Quick Install
-------------

Install HistoSeg from PyPI:

.. code-block:: console

    pip install -U histoseg

Flagship 3D Extra
-----------------

3D Reconstruction uses GitHub ``pyXenium`` for Xenium slide IO and ``trimesh``
for mesh export. Install the extra when you want the
``histoseg-3d reconstruct-stack`` workflow:

.. code-block:: console

    pip install -U "histoseg[threed]"

Optional H&E Extra
------------------

HE Analysis can optionally use local Hugging Face MedSAM inference. Install the
extra dependencies when you need the ``backend="medsam"`` workflow:

.. code-block:: console

    pip install -U "histoseg[he]"

Development Setup
-----------------

Clone the repository and install it in editable mode:

.. code-block:: console

    git clone https://github.com/hutaobo/HistoSeg.git
    cd HistoSeg
    pip install -U pip
    pip install -e ".[threed,he]"

.. note::

   Install ``[threed]`` for the main 3D Xenium contour reconstruction workflow.
   Add ``[he]`` when you need image-first H&E analysis with the MedSAM backend
   available locally.

Feature Groups
--------------

``histoseg.threed``
    Flagship 3D reconstruction workflows for Xenium contour stack alignment,
    sampled 3D points, PLY/OBJ mesh export, and QC visualization.

``histoseg.contour``
    Spatial contour workflows for Pattern1 isolines and multi-structure
    contour exports from cell-coordinate inputs.

``histoseg.he``
    H&E image workflows for tissue foreground extraction, neutral component
    partitioning, and aligned-image change detection.
