# Configuration file for the Sphinx documentation builder.

import os
import subprocess


project = 'SeQuant'
copyright = '2024, Valeev Research Group'
author = 'Valeev Research Group'

# General Configuration
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', "/external/*"]
# templates_path = ['_templates']

html_static_path = ['_static']
html_theme = 'furo' # "sphinx_rtd_theme", "furo", "pydata_sphinx_theme", "sphinx_book_theme"
html_title = "SeQuant Documentation"

# Extensions 
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    "sphinx_copybutton",
    'breathe',
    'exhale',
    "myst_parser"
]

# Exhale Configuration
# Warning: This generates files in the specified folder
exhale_args = {
    "containmentFolder": "./api",
    "rootFileName": "library_root.rst",
    "doxygenStripFromPath": "..",
    "rootFileTitle": "API Reference",
    "createTreeView": False,
    "contentsDirectives": False, # Not needed for themes like furo
}

# Custom options for sphinx_rtd_theme
html_theme_options = {
    "display_version": True,
    "vcs_pageview_mode": "",
    # Toc options
    "collapse_navigation": False,
    "sticky_navigation": True,
    "navigation_depth": 4,
}

# To include GitHub repo info
html_context = {
    "display_github": True,
    "github_user": "ajay-mk",
    "github_repo": "SeQuant",
    "github_version": "master",
    "conf_py_path": "/doc/sphinx/"
}

# Primary Code Language
primary_domain = "cpp"
# Use xeLatex for support for unicode characters
latex_engine = 'xelatex'

# Breathe Configuration
# We use breathe to generate documentation from Doxygen XML files
breath_projects = {"SeQuant": "./../doxygen/xml"}
breathe_default_project = "SeQuant"

# Current default in Sphinx is MathJax v3, but we are using v4
# see https://stackoverflow.com/questions/78813688/mathjax-rendering-issue-in-sphinx-with-symbol
mathjax_path = 'https://cdn.jsdelivr.net/npm/mathjax@4.0.0-beta.3/tex-mml-chtml.js'


# MyST Enables parsing of markdown files
myst_enable_extensions = [
    "amsmath",
    "attrs_inline",
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_admonition",
    "html_image",
    "replacements",
    "smartquotes",
    "strikethrough",
    "substitution",
    "tasklist",
]
