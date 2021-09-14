# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import shutil
import sys
sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------

project = 'homelette'
copyright = '2021, Philipp Junk, Christina Kiel'
author = 'Philipp Junk, Christina Kiel'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'nbsphinx',
    'sphinx_rtd_theme',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
html_logo = 'logo.png'
html_theme_options = {
        'logo_only': False,
        'style_nav_header_background': '#000000',
        }

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']

# -- Options for LaTex output ------------------------------------------------

latex_elements = {
        'preamble': r'''
\setcounter{tocdepth}{1}

\renewcommand{\hyperref}[2][]{#2}
'''
    }


# -- Extension configuration: autodoc ----------------------------------------
autodoc_default_options = {
    'member-order': 'bysource',
}
autoclass_content = 'class'
autodoc_mock_imports = ['altmod', 'modeller', 'ost', 'promod3', 'qmean',
                        'pandas']

# -- Extension configuration: napoleon ---------------------------------------
napoleon_use_ivar = True

# -- Copy notebooks to include in the documentation --------------------------
notebooks = [
        '../examples/Tutorial1_Basics.ipynb',
        '../examples/Tutorial2_Modelling.ipynb',
        '../examples/Tutorial3_Evaluation.ipynb',
        '../examples/Tutorial4_ExtendingHomelette.ipynb',
        '../examples/Tutorial5_Parallelization.ipynb',
        '../examples/Tutorial6_ComplexModelling.ipynb',
        '../examples/Tutorial7_AssemblingPipelines.ipynb',
]

for notebook in notebooks:
    if os.path.exists(notebook):
        shutil.copy(notebook, '.')


# -- Copy logo ---------------------------------------------------------------
if os.path.exists('../logo/logo.png'):
    shutil.copy('../logo/logo.png', '.')
