# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

project = 'TACTICS'
copyright = '2024, Aakankschit Nandkeolyar'
author = 'Aakankschit Nandkeolyar'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
]

# Mock imports for packages that may not be available during doc build
autodoc_mock_imports = [
    'polars',
    'altair',
    'rdkit',
    'rdkit.Chem',
    'rdkit.Chem.Draw',
    'numpy',
    'pandas',
    'matplotlib',
    'matplotlib.pyplot',
    'seaborn',
    'IPython',
    'IPython.display',
    'openeye',
    'openeye.oechem',
    'tqdm',
    'dill',
    'useful_rdkit_utils',
    'collections',
    'typing',
    'pathlib',
    'dataclasses',
    'functools',
    'multiprocessing',
    'json',
    'pickle',
    'subprocess',
    'sys',
    'os',
    're',
    'warnings',
    'time',
    'random',
    'itertools'
]

# Additional autodoc configuration
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
    'special-members': '__init__',
}

# Don't skip __init__ methods
def skip(app, what, name, obj, would_skip, options):
    if name == "__init__":
        return False
    return would_skip

def setup(app):
    app.connect("autodoc-skip-member", skip)

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
    'polars': ('https://pola-rs.github.io/polars/py-polars/html/', None),
    'rdkit': ('https://www.rdkit.org/docs/', None),
}
