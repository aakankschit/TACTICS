# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys

sys.path.insert(0, os.path.abspath("../.."))

project = "TACTICS"
copyright = "2024, Aakankschit Nandkeolyar"
author = "Aakankschit Nandkeolyar"
version = "0.0"
release = "0.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.githubpages",
    "sphinx.ext.graphviz",
    "sphinx_design",
    "sphinx_copybutton",
    "sphinx_togglebutton",
    "myst_nb",
]

# Graphviz configuration
graphviz_output_format = "svg"

# MyST-NB configuration
nb_execution_mode = "off"  # Don't execute notebooks during build
myst_enable_extensions = [
    "colon_fence",
    "deflist",
]

# Mock imports for packages that may not be available during doc build
autodoc_mock_imports = [
    "polars",
    "altair",
    "rdkit",
    "numpy",
    "pandas",
    "matplotlib",
    "seaborn",
    "IPython",
    "openeye",
    "tqdm",
    "dill",
    "useful_rdkit_utils",
    "scipy",
    "sklearn",
]

# Autodoc configuration - cleaner output
autodoc_typehints = "description"
autodoc_typehints_format = "short"
autodoc_class_signature = "separated"
autodoc_member_order = "bysource"
autodoc_default_options = {
    "members": True,
    "show-inheritance": True,
}

# Napoleon settings for Google/NumPy docstrings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_type_aliases = None

# Copybutton configuration
copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True

templates_path = ["_templates"]
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]

html_theme_options = {
    "show_toc_level": 2,
    "navigation_depth": 3,
    "show_nav_level": 2,
    "navbar_align": "left",
    "navbar_start": ["navbar-logo"],
    "navbar_center": ["navbar-nav"],
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "secondary_sidebar_items": ["page-toc", "edit-this-page"],
    "footer_start": ["copyright"],
    "footer_end": ["sphinx-version"],
    "pygments_light_style": "default",
    "pygments_dark_style": "monokai",
    "logo": {
        "image_light": "_static/images/TACTICS_logo.png",
        "image_dark": "_static/images/TACTICS_logo.png",
        "text": "TACTICS",
        "alt_text": "TACTICS - Home",
    },
}

html_context = {
    "default_mode": "auto",
}

# Intersphinx mapping
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "polars": ("https://pola-rs.github.io/polars/py-polars/html/", None),
}
