import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if SRC.exists():
    sys.path.insert(0, str(SRC))

project = "HistoSeg"
author = "HistoSeg authors"

extensions = [
    "myst_nb",
    "sphinx_design",
    "sphinx_copybutton",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
]

autosummary_generate = True
autodoc_member_order = "bysource"
autodoc_typehints = "description"

html_theme = "sphinx_book_theme"
html_title = "HistoSeg"
templates_path = ["_templates"]
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "myst-nb",
    ".ipynb": "myst-nb",
}

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
    "html_admonition",
    "html_image",
    "linkify",
]
myst_heading_anchors = 3

nb_execution_mode = "off"
nb_render_markdown_format = "myst"

exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "3d_analysis.md",
    "api/threed.rst",
    "**/.ipynb_checkpoints",
]

suppress_warnings = [
    "misc.highlighting_failure",
]

html_show_sourcelink = False

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}

html_context = {
    "display_github": True,
    "github_user": "hutaobo",
    "github_repo": "HistoSeg",
    "github_version": "main",
    "conf_py_path": "/docs/",
}

html_theme_options = {
    "repository_url": "https://github.com/hutaobo/HistoSeg",
    "repository_branch": "main",
    "path_to_docs": "docs",
    "use_repository_button": True,
    "use_issues_button": True,
    "use_edit_page_button": True,
    "home_page_in_toc": True,
    "show_toc_level": 2,
    "navigation_with_keys": True,
    "announcement": "HistoSeg unifies H&E image analysis and spatial contour workflows in one Python toolkit.",
}
