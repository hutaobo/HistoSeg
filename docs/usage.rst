=====
Usage
=====

HistoSeg exposes two public feature groups.

HE Analysis
-----------

Use ``histoseg.he`` when your input is an H&E image:

.. code-block:: python

    from histoseg.he import HESegmentationConfig, run_he_segmentation

    result = run_he_segmentation(
        HESegmentationConfig(
            image="/path/to/he.png",
            out_dir="outputs/he_all_elements",
            task="all_elements",
            backend="heuristic",
            n_components=6,
        )
    )

    print(result.overlay_png)

The matching CLI entrypoint is:

.. code-block:: console

    histoseg-he all-elements --image he.png --out-dir outputs/he_all --backend heuristic

Contour Analysis
----------------

Use ``histoseg.contour`` when your input is cell-coordinate data:

.. code-block:: python

    from histoseg.contour import Pattern1IsolineConfig, run_pattern1_isoline

    result = run_pattern1_isoline(
        Pattern1IsolineConfig(
            clusters_csv="clusters.csv",
            cells_parquet="cells.parquet",
            out_dir="outputs/pattern1",
            pattern1_clusters=(10, 23, 19),
        )
    )

    print(result.preview_png)

The matching CLI entrypoint is:

.. code-block:: console

    histoseg-contour pattern1 --clusters-csv clusters.csv --cells-parquet cells.parquet --out-dir outputs/pattern1 --pattern1-clusters 10,23,19
