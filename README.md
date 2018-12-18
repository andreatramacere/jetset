![img](./logo/logo_large.png)


A multiwavelenght modeler and fitting tool for the SED of relativistic  astrophysical jets







#Installation
 I strongly encourage to use anaconda and python3
 - https://www.anaconda.com/download/

## download the JeTseT code
- https://gitlab.com/andrea.tramacere/jetset/-/archive/stable/jetset-stable.tar.gz
- untar the  archive  `jetset-stable.tar.gz`
- cd to  the dir `jetset-stable`

## install requirements
    
   - If you use Anaconda: 
     * `while read requirement; do conda install --yes $requirement; done < requirements.txt`
   
   - If you use pip:
     * `pip install -r requirements.txt `


## install jetset
- `python setup.py clean`
- `python setup.py install`

**run all the examples outsied of the installation dir**

requirements installation without pip or conda
---------------------------------
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