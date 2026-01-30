"""HistoSeg: utilities for spatial transcriptomics segmentation / geometry extraction."""

from .contours.pattern1_isoline import (
    Pattern1IsolineConfig,
    Pattern1IsolineResult,
    run_pattern1_isoline,
    run_pattern1_isoline_from_hf,
)

__all__ = [
    "Pattern1IsolineConfig",
    "Pattern1IsolineResult",
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
