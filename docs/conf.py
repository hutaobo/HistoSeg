# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'HistoSeg'
copyright = '2026, Mengping Long; Taobo Hu; Mats Nilsson'
author = 'Mengping Long; Taobo Hu; Mats Nilsson'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    "myst_nb",
    "sphinx_design",
]

myst_enable_extensions = [
    "colon_fence",
]

nb_execution_mode = "off"

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"

html_theme_options = {
    "navbar_end": ["navbar-icon-links"],
    "navigation_depth": 4,
}

# 让你能用 “卡片/按钮/徽章” 这种漂亮组件

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
