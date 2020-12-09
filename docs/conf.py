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
import sys
sys.path.insert(0, os.path.abspath('../'))
#sys.path.insert(0, os.path.abspath('../victoriaepi/'))
#sys.path.insert(0, os.path.abspath('../victoria/configparsing/'))


# -- Project information -----------------------------------------------------

project = 'Victoria'
copyright = '2020, Marcos Capistran, Antonio Capella, Andres Christen'
author = 'Marcos Capistran, Antonio Capella, Andres Christen'

# The full version, including alpha/beta/rc tags
release = '0.0.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.coverage',
              'sphinx.ext.napoleon',
              'autodocsumm',
              'sphinxcontrib.spelling',
              'sphinx.ext.viewcode',
              #'sphinxcontrib.programoutput',
              ]

# extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon']
napoleon_use_param = False
# source_suffix = ['.rst', '.md']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store','Analy']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
#html_theme = 'alabaster'
html_theme = 'default'

html_theme_options = {
"stickysidebar": "true"
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

spelling_show_suggestions = True
spelling_ignore_pypi_package_names=True
spelling_ignore_wiki_words=True
spelling_ignore_acronyms=True
spelling_ignore_python_builtins=True
spelling_ignore_importable_modules=True
spelling_ignore_contributor_names=True




#clear;clear; make clean && make html -b spelling
#
