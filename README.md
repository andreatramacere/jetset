JetSeT (Jet Sed Tool)
=====================

### A multiwavelenght modeler and fitting tool for the SED of relativistic  astrophysical jets

Installation
------------
requires:
- SWIG (http://www.swig.org/)
    - on linux Ubuntu:
        - sudo apt-get install python-dev
        - sudo apt-get install swig
     - on linux Debian:
        - sudo aptitude install python-dev
        - sudo aptitude install swig
     - on linux Fedora:
        - sudo yum install python-dev
        - sudo yum install swig
     - on mac:
        - ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)" < /dev/null 2> /dev/null
        - brew install swig
     - if you use anaconda:
        - conda install -c anaconda swig
        
 - Python:
     - Anaconda (I strongly encourage to use anaconda)
     - otherwise:
        - python 2.7
        - scipy
        - numpy
        - astropy
       
      - iminuit
        
        https://github.com/iminuit/iminuit
        
 - C compiler:
    - gcc compiler

python setup.py install

Jupyter lab
--------------------
if you ecounter problem with jupyter not running interactive plot
"jupyter labextension install @jupyter-widgets/jupyterlab-manager jupyter-matplotlib"

Build documentation
-------------------
 requires: 
    
 - sphinx
 - pylint
 - sphinx-pyreverse: "https://github.com/alendit/sphinx-pyreverse"
 - nbsphinx: "conda install -c conda-forge nbsphinx"
 - sphinx_rtd_theme: conda install -c anaconda sphinx_rtd_theme 
 - sphinx-better-theme: 'https://github.com/irskep/sphinx-better-theme' 
 - sphinx-bootstrap-theme: 'https://github.com/ryan-roemer/sphinx-bootstrap-theme'
 - sphinx automod: 'https://github.com/astropy/sphinx-automodapi'    