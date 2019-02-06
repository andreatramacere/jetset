<<<<<<< HEAD
![img](./logo/logo_large.png)


A multiwavelenght modeler and fitting framework for the SED of relativistic  astrophysical jets




# Acknowledgements

If you use this code in any kind of scientific publication please cite the following papers:

* `Tramacere A. et al. 2011` http://adsabs.harvard.edu/abs/2011ApJ...739...66T
* `Tramacere A. et al. 2009` http://adsabs.harvard.edu/abs/2009A%26A...501..879T
* `Massaro E. et. al 2006`   http://adsabs.harvard.edu/abs/2006A%26A...448..861M

# Documentation
visit: https://jetset.readthedocs.io/en/latest/

#  Requirements
The following python packages are required:

          python 2.7 or >=3.6 (python 3 is suggested, python 2 should work still fine)
          scipy
          numpy
          astropy
          iminuit (https://github.com/iminuit/iminuit)
         


A C compiler is also necessary, plus the SWIG wrapper generator.

All the dependencies can be installed following the Anaconda method 
(suggested) **OR** the pip method, as described in the following sections.



# Download the JetSeT code
   - Get the source code from: 
     - https://gitlab.com/andrea.tramacere/jetset/-/archive/stable/jetset-stable.tar.gz
   - Uncompress the  archive  `jetset-stable.tar.gz`
   - cd to  the dir `jetset-stable`

# Installation using Anaconda (suggested method)
 - I strongly encourage to use anaconda and python3
   - https://www.anaconda.com/download/

 - #### Download the JetSeT code
   - Get the source code from: 
     - https://gitlab.com/andrea.tramacere/jetset/-/archive/stable/jetset-stable.tar.gz
   - Uncompress the  archive  `jetset-stable.tar.gz`
   - cd to  the dir `jetset-stable`


 - #### Install requirements 
  
    * Linux/MAC : run on the command line
        - `while read requirement; do conda install --yes $requirement; done < requirements.txt`
     
    * Windows   : run on the command line
        - ` FOR /F "delims=~" %f in (requirements.txt) DO conda install --yes "%f"`

 - #### Install JetSeT: 
   run on the command line: 
     * `python setup.py clean`
     * `python setup.py install`

**run all the examples outside of the installation dir**

---------
# Installation using PIP 

 - #### Install requirements 
    - run on the command line: `pip install -r requirements.txt `
    - install SWIG following one of the following methods:
        - SWIG (http://www.swig.org/)
        - on linux Ubuntu:
            - `sudo apt-get install python-dev`
            - `sudo apt-get install swig`
         - on linux Debian:
            - `sudo aptitude install python-dev`
            - `sudo aptitude install swig`
         - on linux Fedora:
            - `sudo yum install python-dev`
            - `sudo yum install swig`
         - on mac:
            - `ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)" < /dev/null 2> /dev/null`
            - `brew install swig`
     

- #### Install JetSeT: 
   run on the command line: 
     * `python setup.py clean`
     * `python setup.py install`


**run all the examples outside of the installation dir**




--------------------------------------------------------------

You can skip the following if your installation was successful


--------------------------------------------------------------

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
 
 
 
=======
# jetset

The code is hosted here: https://gitlab.com/andreatramacere/jetset

The documentation here: https://jetset.readthedocs.io/en/latest/
>>>>>>> 4e79289a1a5884f09a5bc726f765badbf140d02a
