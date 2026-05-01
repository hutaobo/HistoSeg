==========
Quickstart
==========

.. raw:: html

   <div class="hs-section-kicker">Quickstart</div>
   <div class="hs-page-intro">HistoSeg's main workflow is 3D Xenium contour reconstruction. Start with <code>histoseg.threed</code> when your goal is a multi-slice contour stack or mesh, then use <code>histoseg.contour</code> and <code>histoseg.he</code> for the 2D spatial and image inputs that support reconstruction.</div>

3D Reconstruction
-----------------

Use ``histoseg.threed`` when you have neighboring Xenium contour slices from
the same sample and want soft alignment, a 3D contour stack, exported meshes,
and QC artifacts:

.. code-block:: console

    histoseg-3d reconstruct-stack --xenium-root polyp --segmentation-strategy segmentationstrategy.txt --out-dir outputs/polyp_3d

For pairwise soft alignment after hard alignment:

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

.. code-block:: console

    histoseg-3d reconstruct --fixed-geojson slice_01.geojson --moving-hard-aligned-geojson slice_02_hard_aligned_to_01.geojson --out-dir outputs/3d_soft_alignment --group-property structure

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

Next Steps
----------

- Read :doc:`3d_analysis` for the flagship contour reconstruction workflow.
- Read :doc:`contour_analysis` for contour generation, exports, and topology workflows.
- Read :doc:`he_analysis` for the image-first workflow family.
- Open :doc:`tutorials/index` if you want a guided run instead of a reference page.
