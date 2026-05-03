from __future__ import annotations

import json

import pandas as pd
from shapely.geometry import box, mapping

from histoseg.threed import (
    GeneStructureQuantificationConfig,
    SpatialModulePlotConfig,
    plot_spatial_module_clustermap,
    quantify_gene_structure_relationships,
)


def test_quantify_gene_structure_relationships_outputs_sdf_metrics(tmp_path):
    stack_root = tmp_path / "stack"
    gene_dir = stack_root / "gene_overlays" / "TEST_density"
    iso_dir = gene_dir / "isosurfaces"
    aligned_dir = stack_root / "aligned"
    iso_dir.mkdir(parents=True)
    aligned_dir.mkdir(parents=True)

    for index in (1, 2):
        geojson = aligned_dir / f"slice{index}.geojson"
        _write_geojson(geojson, [_feature(box(0, 0, 20, 20), "Structure 1")])

    pd.DataFrame(
        [
            {"order": 1, "sample_id": "s1", "z_um": 0.0, "aligned_geojson": str(aligned_dir / "slice1.geojson")},
            {"order": 2, "sample_id": "s2", "z_um": 5.0, "aligned_geojson": str(aligned_dir / "slice2.geojson")},
        ]
    ).to_csv(stack_root / "aligned_slice_manifest.csv", index=False)

    _write_json(
        gene_dir / "TEST_3d_enrichment_summary.json",
        {
            "gene": "TEST",
            "voxelization": {
                "xy_voxel_size_um": 10.0,
                "z_voxel_size_um": 5.0,
                "grid_shape_zyx": [2, 4, 4],
                "x_min_um": 0.0,
                "x_max_um": 40.0,
                "y_min_um": 0.0,
                "y_max_um": 40.0,
                "z_min_um": 0.0,
                "z_max_um": 5.0,
            },
        },
    )
    _write_json(
        iso_dir / "TEST_3d_enrichment_isosurfaces_summary.json",
        {
            "gene": "TEST",
            "surface_smoothing_sigma_voxels_zyx": [0.0, 0.0, 0.0],
            "thresholds": {"top15": 0.5, "top10": 0.5, "top05": 0.5},
            "mesh_manifest": [],
        },
    )
    rows = []
    for z in (0.0, 5.0):
        for y in (5.0, 15.0):
            for x in (5.0, 15.0):
                rows.append(
                    {
                        "x_um": x,
                        "y_um": y,
                        "z_um": z,
                        "cell_count": 10.0,
                        "test_positive_cell_count": 5.0,
                        "test_expression_sum": 10.0,
                        "test_enrichment": 1.0,
                        "displayed": True,
                    }
                )
    pd.DataFrame(rows).to_csv(gene_dir / "TEST_3d_enrichment_voxels.csv", index=False)

    result = quantify_gene_structure_relationships(
        GeneStructureQuantificationConfig(
            stack_root=stack_root,
            gene_density_dir=gene_dir,
            gene="TEST",
            structures=("Structure 1",),
            pixel_size_um=1.0,
        )
    )

    assert result.overlap_metrics_csv.exists()
    assert result.distance_metrics_csv.exists()
    assert result.overlap_heatmap_png.exists()
    overlap = pd.read_csv(result.overlap_metrics_csv)
    distance = pd.read_csv(result.distance_metrics_csv)
    assert float(overlap.loc[overlap["hotspot"] == "top05", "fraction_of_gene_in_structure"].iloc[0]) == 1.0
    assert float(distance.loc[distance["hotspot"] == "top05", "median_signed_distance_um"].iloc[0]) < 0


def test_plot_spatial_module_clustermap_writes_png(tmp_path):
    batch_dir = tmp_path / "batch"
    batch_dir.mkdir()
    pd.DataFrame(
        {
            "top05|Structure 1": [0.9, 0.1, 0.2],
            "top05|Structure 2": [0.1, 0.8, 0.2],
            "top05|Structure 3": [0.2, 0.2, 0.9],
        },
        index=["GENE1", "GENE2", "GENE3"],
    ).to_csv(batch_dir / "gene_structure_fraction_inside_matrix.csv")

    result = plot_spatial_module_clustermap(
        SpatialModulePlotConfig(batch_dir=batch_dir, hotspot="top05", matrix="fraction_inside")
    )

    assert result.out_png.exists()
    assert result.out_png.stat().st_size > 0


def _feature(geom, structure: str):
    return {
        "type": "Feature",
        "properties": {"structure": structure},
        "geometry": mapping(geom),
    }


def _write_geojson(path, features):
    _write_json(path, {"type": "FeatureCollection", "features": features})


def _write_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload), encoding="utf-8")
