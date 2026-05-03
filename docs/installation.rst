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

3D Reconstruction uses GitHub ``pyXenium`` for Xenium slide IO, ``anndata``
for merged cell tables, and ``trimesh`` for mesh export. Install the extra
when you want ``histoseg-3d reconstruct-stack``, ``project-cells``,
``render-cell-cloud``, or ``discover-spatial-modules``:

.. code-block:: console

    pip install -U "histoseg[threed]"

Visualization And Docs Extras
-----------------------------

``render-cell-cloud`` writes Plotly HTML without requiring PyVista. PyVista is
only needed for local static tutorial or publication figure rendering, such as
the Polyp 24-gene hero images:

.. code-block:: console

    pip install -U "histoseg[threed,viz,docs]"

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
   Add ``[viz,docs]`` when you need static PyVista figure rendering or local
   documentation builds. Add ``[he]`` when you need image-first H&E analysis
   with the MedSAM backend available locally.

Conda Environments
------------------

The repository includes two environment files:

.. code-block:: console

    conda env create -f environment.yml
    conda env create -f environment-viz.yml

``environment.yml`` is the CPU/headless 3D analysis environment.
``environment-viz.yml`` adds PyVista, VTK, and documentation dependencies for
reproducing static tutorial figures.

Docker Targets
--------------

The Dockerfile exposes separate production and visualization targets:

.. code-block:: console

    docker build --target core -t histoseg:core .
    docker build --target viz -t histoseg:viz .

``core`` installs ``histoseg[threed]`` for headless 3D analysis. ``viz`` adds
Mesa/Xvfb system libraries and ``histoseg[viz,docs]`` so PyVista can render
static figures without a desktop display. The Docker context ignores local
large data products such as ``data/``, ``output/``, ``tmp/``, ``*.h5ad``, and
``*.parquet``.

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
