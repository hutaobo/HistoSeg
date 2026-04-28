============
Installation
============

Stable Release
--------------

Install HistoSeg from PyPI:

.. code-block:: console

    pip install -U histoseg

HE Analysis can optionally use local Hugging Face MedSAM inference. Install the
extra dependencies when you need the ``backend="medsam"`` workflow:

.. code-block:: console

    pip install -U "histoseg[he]"

Development Install
-------------------

Clone the repository and install it in editable mode:

.. code-block:: console

    git clone https://github.com/hutaobo/HistoSeg.git
    cd HistoSeg
    pip install -U pip
    pip install -e ".[he]"

Feature Groups
--------------

``histoseg.he``
    H&E image workflows for tissue foreground extraction, neutral component
    partitioning, and aligned-image change detection.

``histoseg.contour``
    Spatial contour workflows for Pattern1 isolines and multi-structure
    contour exports from cell-coordinate inputs.
