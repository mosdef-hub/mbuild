# -*- coding: utf-8 -*-
#
# mbuild documentation build configuration file, created by
# sphinx-quickstart on Wed Oct  1 08:59:12 2014.
#
# This file is execfile()d with the current directory set to its
# containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

import os
import sys
import pathlib

import mock

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

MOCK_MODULES = ['numpy',
                'mdtraj',
                'mdtraj.core.element',
                'nglview',
                'oset',
                'parmed',
                'parmed.utils',
                'parmed.utils.io',
                'parmed.parameters',
                'parmed.periodic_table',
                'scipy',
                'scipy.spatial',
                'scipy.constants',
                'numpy.linalg',
                'gsd.fl',]
for mod_name in MOCK_MODULES:
    sys.modules[mod_name] = mock.Mock()


sys.path.insert(0, os.path.abspath('..'))
# sys.path.insert(0, os.path.abspath('sphinxext'))

base_path = pathlib.Path(__file__).parent
os.system('python {} --name'.format((base_path / '../setup.py').resolve()))

# -- General configuration ------git ------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'nbsphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.autodoc',
    'sphinxcontrib.inkscapeconverter',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'sphinx.ext.autosummary',
    'numpydoc',
]

# nbsphinx specific configuration
nbsphinx_execute = 'never'
nbsphinx_input_prompt = '%s'
nbsphinx_prompt_width = 1.1
html_scaled_image_link = False
nbsphinx_allow_errors = False



# Prolog: Displayed on top of the notebook

nbsphinx_prolog = r"""
{% set docname = 'docs/' + env.doc2path(env.docname, base=None) %}
.. raw:: html

  <div class="admonition note">
      <p>This page was generated from
        <a class="reference external" href="https://github.com/mosdef-hub/mbuild/blob/{{ env.config.release|e }}/{{ docname|e }}">{{ docname|e }}</a>.
        <br>Interactive online version:
        <a href="https://mybinder.org/v2/gh/mosdef-hub/mbuild_tutorials/master?filepath={{ docname|e }}"><img alt="Binder badge" src="https://mybinder.org/badge_logo.svg" style="vertical-align:text-bottom"></a>.
      </p>
    </div>

.. raw:: latex

    \nbsphinxstartnotebook{\scriptsize\noindent\strut
    \textcolor{gray}{The following section was generated from
    \sphinxcode{\sphinxupquote{\strut {{ docname | escape_latex }}}} \dotfill}}

"""

# Displayed on bottom on notebook
nbsphinx_epilog = r"""
{% set docname = 'docs/' + env.doc2path(env.docname, base=None) %}

.. raw:: latex

    \nbsphinxstopnotebook{\scriptsize\noindent\strut
    \textcolor{gray}{\dotfill\ \sphinxcode{\sphinxupquote{\strut
    {{ docname | escape_latex }}}} ends here.}}
"""

autosummary_generate = True
autodoc_default_flags = ['members', ]
numpydoc_class_members_toctree = False

# stackoverflow.com/questions/12206334
numpydoc_show_class_members = False
numpydoc_show_inherited_class_members = False

_python_doc_base = 'http://docs.python.org/3.7'

intersphinx_mapping = {
    _python_doc_base: None,
    'http://docs.scipy.org/doc/numpy': None,
    'http://docs.scipy.org/doc/scipy/reference': None,
    'http://scikit-learn.org/stable': None
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The encoding of source files.
# source_encoding = 'utf-8-sig'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = 'mbuild'
author = 'Mosdef Team'
copyright = u'2014-2019, Vanderbilt University'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#

version = '0.10.13'
release = '0.10.13'


# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
# language = None

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
# today = ''
# Else, today_fmt is used as the format for a strftime call.
# today_fmt = '%B %d, %Y'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build', '**.ipynb_checkpoints']

# The reST default role (used for this markup: `text`) to use for all
# documents.
# default_role = None

# If true, '()' will be appended to :func: etc. cross-reference text.
# add_function_parentheses = True

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
# add_module_names = True

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
# show_authors = False

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# A list of ignored prefixes for module index sorting.
# modindex_common_prefix = []

# If true, keep warnings as "system message" paragraphs in the built documents.
# keep_warnings = False


# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
# html_theme = 'default'
import sphinx_rtd_theme

html_theme = 'sphinx_rtd_theme'
hhtml_theme_path = [
    sphinx_rtd_theme.get_html_theme_path()
]

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
# html_theme_options = {}

# Add any paths that contain custom themes here, relative to this directory.
# html_theme_path = []

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
# html_title = None

# A shorter title for the navigation bar.  Default is the same as html_title.
# html_short_title = None

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
# html_logo = None

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
# html_favicon = None

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Add any extra paths that contain custom files (such as robots.txt or
# .htaccess) here, relative to this directory. These files are copied
# directly to the root of the documentation.
# html_extra_path = []

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
# html_last_updated_fmt = '%b %d, %Y'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
# html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
html_sidebars = {
    '**': ['globaltoc.html', 'sourcelink.html', 'searchbox.html', 'relations.html', 'sourcelink.html'],
}

# Additional templates that should be rendered to pages, maps page names to
# template names.
# html_additional_pages = {}

# If false, no module index is generated.
# html_domain_indices = True

# If false, no index is generated.
# html_use_index = True

# If true, the index is split into individual pages for each letter.
html_split_index = False

# If true, links to the reST sources are added to the pages.
# html_show_sourcelink = True

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
# html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
# html_show_copyright = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
# html_use_opensearch = ''

# This is the file name suffix for HTML files (e.g. ".xhtml").
# html_file_suffix = None

# Output file base name for HTML help builder.
htmlhelp_basename = 'mbuilddoc'

# -- Options for LaTeX output ---------------------------------------------

# See https://www.sphinx-doc.org/en/master/latex.html
latex_elements = {
    'papersize': 'a4paper',
    'printindex': '',
    'sphinxsetup': r"""
        %verbatimwithframe=false,
        %verbatimwrapslines=false,
        %verbatimhintsturnover=false,
        VerbatimColor={HTML}{F5F5F5},
        VerbatimBorderColor={HTML}{E0E0E0},
        noteBorderColor={HTML}{E0E0E0},
        noteborder=1.5pt,
        warningBorderColor={HTML}{E0E0E0},
        warningborder=1.5pt,
        warningBgColor={HTML}{FBFBFB},
    """,
    'preamble': r"""
\usepackage[sc,osf]{mathpazo}
\linespread{1.05}  % see http://www.tug.dk/FontCatalogue/urwpalladio/
\renewcommand{\sfdefault}{pplj}  % Palatino instead of sans serif
\IfFileExists{zlmtt.sty}{
    \usepackage[light,scaled=1.05]{zlmtt}  % light typewriter font from lmodern
}{
    \renewcommand{\ttdefault}{lmtt}  % typewriter font from lmodern
}
\usepackage{booktabs}  % for Pandas dataframes
""",
}

latex_master_doc = 'indexlatex'

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (latex_master_doc, 'mbuild.tex', project,
     author, 'howto'),
]
latex_show_urls = 'footnote'
latex_show_pagerefs = False

# The name of an image file (relative to this directory) to place at the top of
# the title page.
# latex_logo = None

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
# latex_use_parts = False

# If true, show page references after internal links.
# latex_show_pagerefs = False

# If true, show URL addresses after external links.
# latex_show_urls = False

# Documents to append as an appendix to all manuals.
# latex_appendices = []

# If false, no module index is generated.
# latex_domain_indices = True
