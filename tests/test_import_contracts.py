from __future__ import annotations

import importlib


def test_feature_group_import_contracts():
    he = importlib.import_module("histoseg.he")
    contour = importlib.import_module("histoseg.contour")

    assert callable(he.run_he_segmentation)
    assert callable(he.run_he_change_detection)
    assert callable(contour.run_pattern1_isoline)
    assert callable(contour.run_multi_structure_contours)


def test_top_level_exports_feature_group_namespaces():
    import histoseg

    assert histoseg.he.__name__ == "histoseg.he"
    assert histoseg.contour.__name__ == "histoseg.contour"
