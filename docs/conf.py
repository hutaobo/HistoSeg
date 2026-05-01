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

html_theme = "pydata_sphinx_theme"
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
    "github_user": "hutaobo",
    "github_repo": "HistoSeg",
    "github_version": "main",
    "doc_path": "docs",
}

html_theme_options = {
    "navbar_start": ["navbar-logo"],
    "navbar_center": ["navbar-nav"],
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "navbar_align": "left",
    "header_links_before_dropdown": 6,
    "navigation_with_keys": True,
    "show_prev_next": False,
    "secondary_sidebar_items": ["page-toc"],
    "use_edit_page_button": True,
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/hutaobo/HistoSeg",
            "icon": "fa-brands fa-github",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/histoseg/",
            "icon": "fa-solid fa-cube",
        },
        {
            "name": "Read the Docs",
            "url": "https://histoseg.readthedocs.io/en/latest/",
            "icon": "fa-solid fa-book-open",
        },
    ],
    "announcement": "HistoSeg unifies H&E image analysis, spatial contour workflows, and 3D contour alignment in one Python toolkit.",
}
