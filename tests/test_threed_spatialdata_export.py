from __future__ import annotations

import json
import builtins

import pandas as pd
import pytest
from shapely.geometry import box, mapping

from histoseg.threed.spatialdata_export import (
    SpatialData3DExportConfig,
    _collect_aligned_contour_shape_rows,
    export_stack_to_spatialdata_3d,
)


def test_collect_aligned_contour_shapes_preserves_z_and_scales_xy(tmp_path):
    geojson = tmp_path / "aligned.geojson"
    geojson.write_text(
        json.dumps(
            {
                "type": "FeatureCollection",
                "features": [
                    {
                        "type": "Feature",
                        "properties": {"structure": "Local A"},
                        "geometry": mapping(box(0, 0, 10, 20)),
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    manifest = pd.DataFrame(
        [
            {
                "order": 2,
                "sample_id": "slice_2",
                "z_um": 7.5,
                "aligned_geojson": str(geojson),
            }
        ]
    )

    rows = _collect_aligned_contour_shape_rows(
        manifest,
        group_property="structure",
        xenium_pixel_size_um=0.5,
    )

    assert len(rows) == 1
    assert rows[0]["z"] == 7.5
    assert rows[0]["slice_order"] == 2
    assert rows[0]["structure"] == "Local A"
    assert rows[0]["geometry"].bounds == pytest.approx((0.0, 0.0, 5.0, 10.0))


def test_spatialdata_export_reports_optional_dependency_message(tmp_path, monkeypatch):
    stack = tmp_path / "stack"
    stack.mkdir()
    (stack / "aligned_slice_manifest.csv").write_text(
        "order,sample_id,z_um,aligned_geojson\n",
        encoding="utf-8",
    )

    def blocked_import(name, *args, **kwargs):
        if name == "spatialdata":
            raise ModuleNotFoundError("No module named 'spatialdata'")
        return real_import(name, *args, **kwargs)

    real_import = builtins.__import__
    monkeypatch.setattr(builtins, "__import__", blocked_import)

    with pytest.raises(ImportError, match="Exporting to SpatialData requires optional packages"):
        export_stack_to_spatialdata_3d(SpatialData3DExportConfig(stack_root=stack))
