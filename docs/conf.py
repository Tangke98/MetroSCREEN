# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'MetroSCREEN'
copyright = '2024, Wang Lab at Tongji'
author = 'Ke Tang'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

extensions = [
    'nbsphinx',
    'sphinx_gallery.load_style'
]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
html_theme = 'press'
html_static_path = ['_static']
html_css_files = [
    "custom.css",
]
# html_logo = "_static/img/MetroSCREEN_logo_blue_px.png"

nbsphinx_thumbnails = {
    "tutorials/Stereo-seq": "_static/img/thumbnail/Mouse_OB.png",
    "tutorials/seqFISH": "_static/img/thumbnail/seqFISH.png"
}
