# Configuration file for the Sphinx documentation builder.
#
# Minimal + robust RTD configuration for Markdown + notebooks via MyST-NB.

project = "HistoSeg"
copyright = "2026, Mengping Long; Taobo Hu; Mats Nilsson"
author = "Mengping Long; Taobo Hu; Mats Nilsson"
release = "0.1"

# 核心修复：
# - 只启用 myst_nb
# - 不要在 extensions 里同时启用 myst_parser（myst_nb 会自动导入它）
extensions = [
    "myst_nb",
    "sphinx_design",
]

# MyST 语法扩展（myst_nb 会把 myst_parser 的配置选项一起处理）
myst_enable_extensions = [
    "colon_fence",
]

# 如果你不希望 RTD 构建时执行 notebook，保持 off 是最稳妥的
nb_execution_mode = "off"

# 明确告诉 Sphinx 主入口文档名（默认就是 index，但写出来更清晰）
# Sphinx 官方文档：root_doc 默认就是 'index'
root_doc = "index"

# 模板与排除项
templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**/.ipynb_checkpoints"]

# HTML theme
html_theme = "pydata_sphinx_theme"
html_theme_options = {
    "navbar_end": ["navbar-icon-links"],
    "navigation_depth": 4,
}
