import sys, datetime, sphinx_rtd_theme

try:
    from jLM import buildData
    releaseVersion = "{}/{}/{}".format(buildData['version'], buildData['gitHash'][:6], buildData['buildTime'].split()[0])
except:
    print("Install LM before generating documentation")
    sys.exit(-1)

extensions =           ['sphinx.ext.autosummary',
                        'sphinx.ext.autodoc',
                        'sphinx.ext.coverage',
                        'sphinx.ext.mathjax',
                        'sphinx.ext.napoleon',
                        'sphinx.ext.todo',
                        'sphinx.ext.viewcode']
templates_path =        ['_templates']
source_suffix =         '.rst'
master_doc =            'index'
project =               'Lattice Microbes'
copyright =             '{}, Luthey-Schulten Lab'.format(datetime.datetime.now().year)
author =                'Tyler M. Earnest, Michael J. Hallock, Andrew Magis, Joseph R. Peterson, Elijah Roberts'
version =               buildData['version']
release =               releaseVersion
language =              None
exclude_patterns =      []
pygments_style =        'sphinx'
todo_include_todos =    True
html_theme =            'sphinx_rtd_theme'
html_theme_path =       [sphinx_rtd_theme.get_html_theme_path()]
html_static_path =      ['_static']
intersphinx_mapping =   {'python': ('https://docs.python.org/3.5', None),
                         'numpy': ('http://docs.scipy.org/doc/numpy/', None)}
autodoc_default_flags = ['members', 'show-inheritance', 'inherited-members', 'undoc-members']
autosummary_generate =  True
html_show_sourcelink =  True
