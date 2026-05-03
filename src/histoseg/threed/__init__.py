"""3D Analysis public API.

This namespace contains same-sample, multi-slice Xenium contour workflows.
It includes conservative TPS soft alignment for a hard-aligned contour pair
and pyXenium-backed multi-slice contour stack reconstruction.
"""

from __future__ import annotations

from .soft_alignment import (
    ThreeDContourReconstructionConfig,
    ThreeDContourReconstructionResult,
    ThreeDFeatureUnavailableError,
    run_3d_contour_reconstruction,
)
from .multislice import (
    ThreeDStackReconstructionConfig,
    ThreeDStackReconstructionResult,
    discover_xenium_slices,
    hard_align_geojson,
    reconstruct_3d_contour_meshes,
    run_3d_stack_reconstruction,
    write_3d_contour_points,
    write_3d_visualization_html,
)
from .cell_cloud import (
    CELL_CLOUD_ALIGNED_XY_OBSM_KEY,
    CELL_CLOUD_OBS_SLICE_KEY,
    CELL_CLOUD_OBSM_KEY,
    CELL_CLOUD_UNS_KEY,
    CellCloudProjectionConfig,
    CellCloudProjectionResult,
    build_alignment_manifest,
    cell_cloud_cache_status,
    cell_cloud_dataframe_from_coordinates,
    hash_alignment_manifest,
    load_cell_alignment_transforms,
    project_cell_coordinates,
    run_cell_cloud_projection,
    write_cell_cloud_cache,
)
from .spatial_modules import (
    GeneStructureQuantificationConfig,
    GeneStructureQuantificationResult,
    SpatialModuleDiscoveryConfig,
    SpatialModuleDiscoveryResult,
    SpatialModulePlotConfig,
    SpatialModulePlotResult,
    compute_hotspot_overlap_metrics,
    compute_hotspot_sdf_metrics,
    plot_spatial_module_clustermap,
    quantify_gene_structure_relationships,
    run_spatial_module_discovery,
)

__all__ = [
    "CELL_CLOUD_ALIGNED_XY_OBSM_KEY",
    "CELL_CLOUD_OBS_SLICE_KEY",
    "CELL_CLOUD_OBSM_KEY",
    "CELL_CLOUD_UNS_KEY",
    "CellCloudProjectionConfig",
    "CellCloudProjectionResult",
    "GeneStructureQuantificationConfig",
    "GeneStructureQuantificationResult",
    "SpatialModuleDiscoveryConfig",
    "SpatialModuleDiscoveryResult",
    "SpatialModulePlotConfig",
    "SpatialModulePlotResult",
    "ThreeDContourReconstructionConfig",
    "ThreeDContourReconstructionResult",
    "ThreeDFeatureUnavailableError",
    "ThreeDStackReconstructionConfig",
    "ThreeDStackReconstructionResult",
    "compute_hotspot_overlap_metrics",
    "compute_hotspot_sdf_metrics",
    "build_alignment_manifest",
    "cell_cloud_cache_status",
    "cell_cloud_dataframe_from_coordinates",
    "discover_xenium_slices",
    "hard_align_geojson",
    "hash_alignment_manifest",
    "load_cell_alignment_transforms",
    "plot_spatial_module_clustermap",
    "project_cell_coordinates",
    "quantify_gene_structure_relationships",
    "reconstruct_3d_contour_meshes",
    "run_3d_contour_reconstruction",
    "run_cell_cloud_projection",
    "run_spatial_module_discovery",
    "run_3d_stack_reconstruction",
    "write_3d_contour_points",
    "write_3d_visualization_html",
    "write_cell_cloud_cache",
]
