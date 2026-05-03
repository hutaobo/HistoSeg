from __future__ import annotations

import importlib


def test_feature_group_import_contracts():
    he = importlib.import_module("histoseg.he")
    contour = importlib.import_module("histoseg.contour")
    threed = importlib.import_module("histoseg.threed")

    assert callable(he.run_he_segmentation)
    assert callable(he.run_he_change_detection)
    assert callable(contour.run_pattern1_isoline)
    assert callable(contour.run_gene_transcript_isoline)
    assert callable(contour.run_multi_structure_contours)
    assert callable(contour.run_group_boundary_network)
    assert callable(contour.run_contour_adjacency)
    assert callable(threed.run_3d_contour_reconstruction)
    assert callable(threed.run_3d_stack_reconstruction)
    assert callable(threed.run_cell_cloud_projection)
    assert callable(threed.run_spatial_module_discovery)
    assert callable(threed.quantify_gene_structure_relationships)
    assert callable(threed.plot_spatial_module_clustermap)
    assert callable(threed.discover_xenium_slices)


def test_top_level_exports_feature_group_namespaces():
    import histoseg

    assert histoseg.he.__name__ == "histoseg.he"
    assert histoseg.contour.__name__ == "histoseg.contour"
    assert histoseg.threed.__name__ == "histoseg.threed"
