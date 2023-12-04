import os
import sys
sys.path.insert(0, os.path.abspath('./'))

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Computational Physics Project'
copyright = '2023, Alice Eriksson, Emelie Eriksson, Jakov Krnic, Issa Nseir, Markus Wallin, Gustav Wassbäck'
author = 'Alice Eriksson, Emelie Eriksson, Jakov Krnic, Issa Nseir, Markus Wallin, Gustav Wassbäck'
release = '0.2'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    ]


templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

#html_theme = 'classic'
html_theme = "nature"
html_static_path = ['_static']
