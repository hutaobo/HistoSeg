==========
Quickstart
==========

.. raw:: html

   <div class="hs-section-kicker">Quickstart</div>
   <div class="hs-page-intro">HistoSeg exposes three public workflow surfaces. Choose the one that matches your primary input modality, then move between Python and CLI runs without changing the conceptual workflow or output artifact model.</div>

HE Analysis
-----------

Use ``histoseg.he`` when your starting point is an H&E image and your desired
output is a tissue mask, neutral component partition, overlay, or change map:

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

    histoseg-he all-elements --image he.png --out-dir outputs/he_all --backend heuristic --n-components 6

Contour Analysis
----------------

Use ``histoseg.contour`` when your starting point is cell-coordinate data plus
cluster assignments and you want Pattern1 isolines or named structure contours:

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

3D Analysis
-----------

Use ``histoseg.threed`` when you have neighboring Xenium contour slices from
the same sample and the moving slice has already been hard-aligned to the fixed
slice:

.. code-block:: python

    from histoseg.threed import ThreeDContourReconstructionConfig, run_3d_contour_reconstruction

    result = run_3d_contour_reconstruction(
        ThreeDContourReconstructionConfig(
            fixed_geojson="slice_01.geojson",
            moving_hard_aligned_geojson="slice_02_hard_aligned_to_01.geojson",
            out_dir="outputs/3d_soft_alignment",
            group_property="structure",
        )
    )

    print(result.soft_aligned_geojson)

The matching CLI entrypoint is:

.. code-block:: console

    histoseg-3d reconstruct --fixed-geojson slice_01.geojson --moving-hard-aligned-geojson slice_02_hard_aligned_to_01.geojson --out-dir outputs/3d_soft_alignment --group-property structure

Next Steps
----------

- Read :doc:`he_analysis` for the image-first workflow family.
- Read :doc:`contour_analysis` for contour generation, exports, and topology workflows.
- Read :doc:`3d_analysis` for the contour soft-alignment workflow.
- Open :doc:`tutorials/index` if you want a guided run instead of a reference page.
