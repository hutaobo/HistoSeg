"""HistoSeg: utilities for spatial transcriptomics segmentation / geometry extraction."""

from .contours.pattern1_isoline import (
    Pattern1IsolineConfig,
    Pattern1IsolineResult,
    run_pattern1_isoline,
    run_pattern1_isoline_from_hf,
)
from .contours.multi_structure import (
    MultiStructureContourConfig,
    MultiStructureContourResult,
    MultiStructureSpec,
    run_multi_structure_contours,
)

__all__ = [
    "MultiStructureContourConfig",
    "MultiStructureContourResult",
    "MultiStructureSpec",
    "Pattern1IsolineConfig",
    "Pattern1IsolineResult",
    "run_multi_structure_contours",
    "run_pattern1_isoline",
    "run_pattern1_isoline_from_hf",
]

from .sfplot.Searcher_Findee_Score import (
    compute_cophenetic_distances_from_df,
    plot_cophenetic_heatmap,
)

__all__ += [
    "compute_cophenetic_distances_from_df",
    "plot_cophenetic_heatmap",
]

from .spatial_pathologist import (
    FullAutoSpatialPathologistConfig,
    SpatialPathologistConfig,
    load_full_auto_spatial_pathologist_config,
    load_spatial_pathologist_config,
    run_full_auto_spatial_pathologist,
    run_spatial_pathologist,
)

__all__ += [
    "FullAutoSpatialPathologistConfig",
    "SpatialPathologistConfig",
    "load_full_auto_spatial_pathologist_config",
    "load_spatial_pathologist_config",
    "run_full_auto_spatial_pathologist",
    "run_spatial_pathologist",
]

# Optional: Hugging Face helpers (requires `huggingface_hub`)
try:
    from .io.huggingface import (
        XeniumOutsPaths,
        download_files,
        download_xenium_outs,
    )
    __all__ += [
        "XeniumOutsPaths",
        "download_files",
        "download_xenium_outs",
    ]
except Exception:
    # Keep base package usable without huggingface_hub installed.
    pass
