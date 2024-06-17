rtd_req=['sphinx-bootstrap-theme',
        'sphinx_book_theme',
        'sphinx-automodapi',
        'sphinx_rtd_theme>=0.3.1',
        'sphinxcontrib-bibtex',
        'sphinx-gallery',
        'sphinx-nbexamples',
        'sphinx>=7,<=8',
        'nbsphinx',
        'numpydoc',
        'graphviz',
        'mock']

_skip_list=['pyqt','jupyter','ipython']
f = open("../requirements.txt",'r')
req=f.readlines()
f.close()

req=[n.strip() for n in req]
for r in req[:]:
    for s in _skip_list[:]:

        if s in r:
            req.remove(r)

f = open("./requirements.txt",'w')
for r in req:
    print(r,file=f)
for r in rtd_req:
    print(r,file=f)
f.close()