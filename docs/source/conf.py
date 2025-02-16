# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Power iterations in parallel using MPI'
copyright = '2025, Mukund Yedunuthala'
author = 'Mukund Yedunuthala'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration
from sphinx.builders.html import StandaloneHTMLBuilder
import subprocess, os
subprocess.call('doxygen Doxyfile.in', shell=True)

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
    'sphinx.ext.inheritance_diagram',
    'sphinx_sitemap',
    'breathe'
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
highlight_language = 'c++'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'shibuya'
html_baseurl = 'https://mukund-yedunuthala.de/docs/mpi-eigen-value/'
html_static_path = ['_static']
breathe_projects = {
	"Power iterations in parallel using MPI": "build/xml/"
}
breathe_default_project = "Power iterations in parallel using MPI"
breathe_default_members = ('members', 'undoc-members')