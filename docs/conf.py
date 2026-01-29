# Configuration file for the Sphinx documentation builder.
from __future__ import annotations

import os
import sys
from datetime import datetime
from pathlib import Path

# -- Path setup --------------------------------------------------------------
DOCS_DIR = Path(__file__).resolve().parent
ROOT = DOCS_DIR.parent
SRC = ROOT / "src"

# 兼容 src-layout；即使 RTD 没装 -e .，autodoc 也能 import
if SRC.exists():
    sys.path.insert(0, str(SRC))
sys.path.insert(0, str(ROOT))

# -- Project information -----------------------------------------------------
project = os.environ.get("SPHINX_PROJECT", "HistoSeg")
author = os.environ.get("SPHINX_AUTHOR", "Taobo Hu")
copyright = f"{datetime.now():%Y}, {author}"

# 版本信息：优先读环境变量；读不到就留空（避免因为包名不匹配导致构建失败）
release = os.environ.get("SPHINX_RELEASE", "")
version = release

# -- General configuration ---------------------------------------------------
extensions = [
    # Sphinx core
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",

    # Markdown + notebooks
    # ⚠️ 重点：不要再加 myst_parser。myst_nb 会自动加载 myst_parser，手动加会冲突。
    "myst_nb",

    # UI components
    "sphinx_design",
    "sphinx_copybutton",
]

# 这些目录如果不存在就置空，避免 warnings（有些项目会把 warning 当 error）
templates_path = ["_templates"] if (DOCS_DIR / "_templates").exists() else []
html_static_path = ["_static"] if (DOCS_DIR / "_static").exists() else []

exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "**/.ipynb_checkpoints",
]

# Sphinx >= 4 用 root_doc
root_doc = "index"

# autosummary / autodoc
autosummary_generate = True
autodoc_default_options = {
    "members": True,
    "show-inheritance": True,
}

# 如果你的包 import 会拉一堆重依赖、在 RTD 上不一定装全，
# 可以把它们加到这里 mock 掉（按需填）
autodoc_mock_imports = [
    # "scanpy", "anndata", "torch", "tensorflow", ...
]

# -- MyST / MyST-NB configuration -------------------------------------------
# 这些扩展对 Markdown 写作很实用，且依赖最少（稳）
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "tasklist",
]
myst_heading_anchors = 3

# Notebook：默认不执行（RTD 最稳），需要执行再用环境变量切换
# 可选值：off / auto / force / cache
nb_execution_mode = os.environ.get("NB_EXECUTION_MODE", "off")
nb_execution_timeout = 300
nb_merge_streams = True

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
    ".ipynb": "myst-nb",
}

# -- HTML output -------------------------------------------------------------
html_theme = os.environ.get("SPHINX_HTML_THEME", "pydata_sphinx_theme")
html_theme_options = {
    "show_toc_level": 2,
    "navigation_with_keys": True,
}

# -- Intersphinx -------------------------------------------------------------
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", {}),
}

# -- Copybutton --------------------------------------------------------------
copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True
