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
from .spatial_modules import (
    GeneStructureQuantificationConfig,
    GeneStructureQuantificationResult,
    SpatialModuleDiscoveryConfig,
    SpatialModuleDiscoveryResult,
    SpatialModulePlotConfig,
    SpatialModulePlotResult,
    plot_spatial_module_clustermap,
    quantify_gene_structure_relationships,
    run_spatial_module_discovery,
)

__all__ = [
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
    "discover_xenium_slices",
    "hard_align_geojson",
    "plot_spatial_module_clustermap",
    "quantify_gene_structure_relationships",
    "reconstruct_3d_contour_meshes",
    "run_3d_contour_reconstruction",
    "run_spatial_module_discovery",
    "run_3d_stack_reconstruction",
    "write_3d_contour_points",
    "write_3d_visualization_html",
]
