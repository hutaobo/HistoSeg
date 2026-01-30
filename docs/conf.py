import sys
from pathlib import Path

# 确保可以 import histoseg
ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if SRC.exists():
    sys.path.insert(0, str(SRC))

project = "HistoSeg"
author = "HistoSeg authors"

extensions = [
    "myst_nb",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
]

html_theme = "sphinx_rtd_theme"

# 明确声明所有文档类型（避免 RTD / Sphinx 猜错）
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
    ".ipynb": "myst-nb",
}

# ★ 核心：绝不在 RTD 上执行 notebook
nb_execution_mode = "off"

# notebook 中的 markdown 用 myst 渲染
nb_render_markdown_format = "myst"

exclude_patterns = [
    "_build",
    "**/.ipynb_checkpoints",
]
