import sys
from pathlib import Path

# Ensure we can import histoseg
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

html_theme = "furo"

# Key: let myst-nb parse both .md and .ipynb
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "myst-nb",
    ".ipynb": "myst-nb",
}

# Core: never execute notebooks on RTD
nb_execution_mode = "off"

# Render notebook markdown via MyST
nb_render_markdown_format = "myst"

exclude_patterns = [
    "_build",
    "**/.ipynb_checkpoints",
]

# Remove the default "View page source" link (/_sources/...)
html_show_sourcelink = False

# Enable GitHub links (e.g., "Edit on GitHub" / "View on GitHub")
# This should point to:
# https://github.com/hutaobo/HistoSeg/blob/master/docs/<pagename><suffix>
html_context = {
    "display_github": True,
    "github_user": "hutaobo",
    "github_repo": "HistoSeg",
    "github_version": "master",
    "conf_py_path": "/docs/",
}

# Ensure GitHub links open file view ("blob") instead of edit/raw
html_theme_options = {
    "source_repository": "https://github.com/hutaobo/HistoSeg/",
    "source_branch": "master",
    "source_directory": "docs/",
    "top_of_page_buttons": ["view", "edit"],
}
