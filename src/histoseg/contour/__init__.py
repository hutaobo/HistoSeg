"""Contour Analysis public API."""

from histoseg.sfplot.Searcher_Findee_Score import (
    compute_cophenetic_distances_from_df,
    compute_cophenetic_from_distance_matrix,
    compute_searcher_findee_distance_matrix_from_df,
    plot_cophenetic_heatmap,
)

from .multi_structure import (
    MultiStructureContourConfig,
    MultiStructureContourResult,
    MultiStructureSpec,
    run_multi_structure_contours,
)
from .pattern1_isoline import (
    Pattern1IsolineConfig,
    Pattern1IsolineResult,
    SegmentationConfidenceResult,
    compute_segmentation_confidence_score,
    compute_segmentation_confidence_score_from_merged,
    run_pattern1_isoline,
    run_pattern1_isoline_from_hf,
)
from .topology import (
    ContourTopologyConfig,
    ContourTopologyResult,
    summarize_contour_topology,
)

__all__ = [
    "ContourTopologyConfig",
    "ContourTopologyResult",
    "MultiStructureContourConfig",
    "MultiStructureContourResult",
    "MultiStructureSpec",
    "Pattern1IsolineConfig",
    "Pattern1IsolineResult",
    "SegmentationConfidenceResult",
    "compute_cophenetic_distances_from_df",
    "compute_cophenetic_from_distance_matrix",
    "compute_searcher_findee_distance_matrix_from_df",
    "compute_segmentation_confidence_score",
    "compute_segmentation_confidence_score_from_merged",
    "plot_cophenetic_heatmap",
    "run_multi_structure_contours",
    "run_pattern1_isoline",
    "run_pattern1_isoline_from_hf",
    "summarize_contour_topology",
]
