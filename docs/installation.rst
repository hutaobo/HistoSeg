============
Installation
============

.. raw:: html

   <div class="hs-section-kicker">Getting Started</div>
   <div class="hs-page-intro">HistoSeg ships as a Python package with a lightweight base install and an optional H&amp;E extra for local MedSAM-backed workflows. Install the base package for contour-first work, then add the H&amp;E extra when you need local image-model inference.</div>

Quick Install
-------------

Install HistoSeg from PyPI:

.. code-block:: console

    pip install -U histoseg

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
    pip install -e ".[he]"

.. note::

   Install the base package when your work is centered on contour workflows and
   cell-coordinate data. Add ``[he]`` when you need image-first H&E analysis
   with the MedSAM backend available locally.

Feature Groups
--------------

``histoseg.he``
    H&E image workflows for tissue foreground extraction, neutral component
    partitioning, and aligned-image change detection.

``histoseg.contour``
    Spatial contour workflows for Pattern1 isolines and multi-structure
    contour exports from cell-coordinate inputs.
