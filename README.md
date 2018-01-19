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

python setup.py install