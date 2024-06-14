# -*- coding: utf-8 -*-
#
# asterism documentation build configuration file, created by
# sphinx-quickstart on Thu Apr 21 10:33:01 2016.
#
# This file is execfile()d with the current directory set to its
# containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

import sys
import os
import json
import mock
#import jetset
# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
sys.path.insert(0,os.path.abspath('../'))
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'


autodoc_mock_imports=[]
autodoc_mock_imports = ["jetkernel"]
autodoc_mock_imports.append('_jetkernel')
autodoc_mock_imports.append('gammapy')
autodoc_mock_imports.append('sherpa')

for mod_name in autodoc_mock_imports:
    sys.modules[mod_name] = mock.Mock()


import sphinx_bootstrap_theme
if on_rtd==False:  # only import and set the  if we're building docs locally

    

    #theme='bootstrap'
    theme = 'sphinx_book_theme'
    #theme='sphinx_rtd_theme'
    #
else:
    #theme = 'bootstrap'
    theme = 'sphinx_book_theme'




extensions = [
    #'autoapi.extension'
    'sphinx.ext.autodoc',
    'sphinx_automodapi.automodapi',
    'sphinx_automodapi.smart_resolver',
    'sphinx_gallery.load_style',
    'sphinx.ext.doctest',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
    'sphinx.ext.autosummary',
    'sphinx.ext.graphviz',
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.napoleon',
    'sphinx.ext.inheritance_diagram',
    'sphinx.ext.autosummary',
    'nbsphinx',
#    'sphinxcontrib.bibtex',
    'sphinx.ext.mathjax',
]

#bibtex_bibfiles = ['refs.bib']
#bibtex_bibfiles = ['references.bib']
exclude_patterns = ['_build', 
                    '**.ipynb_checkpoints',
                    '../jetkernel/*',
                    '../jetset/jetkernel/*',
                    'documentation_notebooks',
                    'example_notebooks',
                    'slides']

#autoapi_dirs = ['../jetset']

#templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'
with open('../jetset/pkg_info.json') as fp:
    _info = json.load(fp)
__version__ = _info['version']
#import jetset
#__version__ = jetset.__version__

# General information about the project.
# The short X.Y version.
version = __version__
# The full version, including alpha/beta/rc tags.
project = u'jetset'
copyright = u'2019, andrea tramacere'





add_module_names = False
pygments_style = 'sphinx'





#html_static_path = ['_static']
#html_logo = "_static/logo_small_color_transparent.png"


#def setup(app):
#    #app.add_stylesheet("my_theme.css") # also can be a full URL
#    app.add_css_file("css/my_theme.css")

if theme=='bootstrap':
    html_sidebars = {'**': ['localtoc.html','my_side_bar.html']}

    html_logo = "_static/logo_small_color_neg_transparent.png"
    html_theme = 'bootstrap'
    html_static_path = ['_static']
    if not on_rtd:
        html_theme_path = sphinx_bootstrap_theme.get_html_theme_path()
    #else:
    #    html_static_path = ['_static']


    html_theme_options = {
        # Navigation bar title. (Default: ``project`` value)
        'navbar_title': "JetSeT doc",

        # Tab name for entire site. (Default: "Site")
        'navbar_site_name': "JetSeT",

        # A list of tuples containing pages or urls to link to.
        # Valid tuples should be in the following forms:
        #    (name, page)                 # a link to a page
        #    (name, "/aa/bb", 1)          # a link to an arbitrary relative url
        #    (name, "http://example.com", True) # arbitrary absolute url
        # Note the "1" or "True" value above as the third argument to indicate
        # an arbitrary url.
        #'navbar_links': [
        #    ("Examples", "examples"),
        #    ("Link", "http://example.com", True),
        #],

        # Render the next and previous page links in navbar. (Default: true)
        'navbar_sidebarrel': True,

        # Render the current pages TOC in the navbar. (Default: true)
        'navbar_pagenav': True,

        # Tab name for the current pages TOC. (Default: "Page")
        'navbar_pagenav_name': "Page",

        # Global TOC depth for "site" navbar tab. (Default: 1)
        # Switching to -1 shows all levels.
        'globaltoc_depth': 2,

        # Include hidden TOCs in Site navbar?
        #
        # Note: If this is "false", you cannot have mixed ``:hidden:`` and
        # non-hidden ``toctree`` directives in the same page, or else the build
        # will break.
        #
        # Values: "true" (default) or "false"
        'globaltoc_includehidden': True,

        # HTML navbar class (Default: "navbar") to attach to <div> element.
        # For black navbar, do "navbar navbar-inverse"
        'navbar_class': "navbar navbar-inverse",


        # Fix navigation bar to top of page?
        # Values: "true" (default) or "false"
        'navbar_fixed_top': True,

        # Location of link to source.
        # Options are "nav" (default), "footer" or anything else to exclude.
        'source_link_position': "nav",
        # Bootswatch (http://bootswatch.com/) theme.
        #
        # Options are nothing (default) or the name of a valid theme
        # such as "cosmo" or "sandstone".
        #
        # The set of valid themes depend on the version of Bootstrap
        # that's used (the next config option).
        #
        # Currently, the supported themes are:
        # - Bootstrap 2: https://bootswatch.com/2
        # - Bootstrap 3: https://bootswatch.com/3
        'bootswatch_theme': "spacelab",

        #'nosidebar': True,

        # Choose Bootstrap version.
        # Values: "3" (default) or "2" (in quotes)
        'bootstrap_version': "3",
    }


if theme=='sphinx_book_theme':
    html_theme = "sphinx_book_theme"
    html_static_path = ["_static/css/sphinx_book_theme"]
    html_css_files = ["custom.css"]

    html_theme_options = {
        "icon_links": [
        {
            # Label for this link
            "name": "GitHub",
            # URL where the link will redirect
            "url": "https://github.com/andreatramacere/jetset",  # required
            # Icon class (if "type": "fontawesome"), or path to local image (if "type": "local")
            "icon": "fa-brands fa-square-github",
            # The type of image to be used (see below for details)
            "type": "fontawesome",
            }
         ],
        
        "logo": {
        # In a left-to-right context, screen readers will read the alt text
        # first, then the text, so this example will be read as "P-G-G-P-Y
        # (short pause) Home A pretty good geometry package"
        "alt_text": "JetSeT ",
        "text": "documentation",
        "image_light": "_static/logo_large_no_border.png",
        "image_dark": "_static/logo_large_no_border.png",
        },
        "show_toc_level": 3,

    }

    pass



htmlhelp_basename = 'jetsetdoc'


# -- Options for LaTeX output ---------------------------------------------
latex_engine = 'pdflatex'
latex_elements = {
# The paper size ('letterpaper' or 'a4paper').
#'papersize': 'letterpaper',

# The font size ('10pt', '11pt' or '12pt').
#'pointsize': '10pt',

# Additional stuff for the LaTeX preamble.
#'preamble': '',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
  ('index', 'jetset.tex', u'jetset Documentation',
   u'andrea tramacere', 'manual'),
]


man_pages = [
    ('index', 'jetset', u'jetset Documentation',
     [u'andrea tramacere'], 1)
]

# If true, show URL addresses after external links.
#man_show_urls = False


# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
  ('index', 'jetset', u'jetset Documentation',
   u'andrea tramacere', 'jetset', 'One line description of project.',
   'Miscellaneous'),
]



# Bibliographic Dublin Core info.
epub_title = u'jetset'
epub_author = u'andrea tramacere'
epub_publisher = u'andrea tramacere'
epub_copyright = u'2016, andrea tramacere'


epub_exclude_files = ['search.html']


