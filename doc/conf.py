import os
import sys
sys.path.insert(0, os.path.abspath('../src/'))

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'sith'
copyright = '2023, Daniel Sucerquia, Mikaela Farrugia'
author = 'Daniel Sucerquia, Mikaela Farrugia'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.coverage',
              'sphinx.ext.napoleon',
              'sphinx_rtd_theme',
              'sphinxcontrib.mermaid']

html_theme = 'sphinx_rtd_theme'
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

#remote_theme: rundocs/jekyll-rtd-theme
html_static_path = ['_static']
html_favicon = '_static/favicon.ico'
html_title = "sith"
html_logo = '_static/favicon.ico'

html_theme_options = {
    "style_nav_header_background": "#3434348D",
    'logo_only': True,
    'display_version' : True,
    'style_external_links' : True
}

# Add custom CSS file
html_css_files = [
    'custom.css',  # Ensure the file path is correct
]

from docutils import nodes
from docutils.parsers.rst import roles

def bashscript_role(name, rawtext, text, lineno, inliner, options={}, content=[]):
    """
    Custom role :bashscript:`display <target-label>` to render a boxed ref link.
    """
    env = inliner.document.settings.env

    parts = text.split('.')
    target = 'modules/' + '.'.join(parts[:-1])
    ref = '#' + text.replace('.', '-')
    ref = ref.replace('_', '-')
    display_text = text.replace('.', '/') + '.sh'


    # Generate relative URI correctly using positional args
    refuri = env.app.builder.get_relative_uri(env.docname, target) + ref

    # Create reference node
    refnode = nodes.reference(rawtext, display_text, refuri=refuri)
    refnode['classes'].append('bashscript')

    return [refnode], []

# Register the role
roles.register_local_role('bashscript', bashscript_role)