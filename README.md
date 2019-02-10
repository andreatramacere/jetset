[![Build Status](https://travis-ci.org/andreatramacere/jetset.svg?branch=py23)](https://travis-ci.org/andreatramacere/jetset)

![img](./logo/logo_large.png)


JetSeT is  an open source  C/Python   framework  to reproduce radiative and accelerative processes acting in relativistic jets,  
allowing to fit the numerical models to observed data. The main features of this framework are: 

 * handling observed data: re-binning, definition of data sets, bindings to astropy tables and quantities
   definition of complex numerical radiative scenarios: Synchrotron Self-Compton (SSC), external Compton (EC) and EC 
   against the CMB 
 
 * Constraining of the model in the pre-fitting stage, based on accurate  and already published phenomenological trends. 
   In particular, starting from phenomenological parameters, such as spectral indices, peak fluxes and frequencies, and 
   spectral  curvatures, that the code evaluates automatically, the pre-fitting algorithm is able to provide a good 
   starting model,following the phenomenological trends that I have implemented. fitting of multiwavelength SEDs using  
   both frequentist approach (iminuit) and bayesian MCMC sampling (emcee)
 
 * Self-consistent temporal evolution of the plasma under the effect of radiative and accelerative processes, both first  
   order and second order (stochastic acceleration) processes.



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
     - https://gitlab.com/andreatramacere/jetset/-/archive/stable/jetset-stable.tar.gz

     OR

     - https://github.com/andreatramacere/jetset/archive/stable.tar.gz
     
 
### Installation using Anaconda (suggested method)
 - I suggest to use anaconda and python3
 - https://www.anaconda.com/download/
 
 - Uncompress the  archive  `jetset-stable.tar.gz`
 - cd to  the dir `jetset-stable`


 - ##### Install requirements 
  
    * Linux/MAC : run on the command line
        - `while read requirement; do conda install --yes $requirement; done < requirements.txt`
     
    * Windows   : run on the command line
        - ` FOR /F "delims=~" %f in (requirements.txt) DO conda install --yes "%f"`

 - ##### Install JetSeT: 
   run on the command line: 
     * `python setup.py clean`
     * `python setup.py install`

**run all the examples outside of the installation dir**

---------
### Installation using PIP 
 - cd to  the dir `jetset-stable`
 - ##### Install requirements 
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
     

- ##### Install JetSeT: 
   run on the command line: 
     * `python setup.py clean`
     * `python setup.py install`


**run all the examples outside of the installation dir**




# jetset code repoistory

The code is hosted here: 
 - https://gitlab.com/andreatramacere/jetset
 
 OR
 
 -  https://github.com/andreatramacere/jetset
 




# Build documentation

 requires: 
    
 - sphinx
 - pylint
 - sphinx-pyreverse: "https://github.com/alendit/sphinx-pyreverse"
 - nbsphinx: "conda install -c conda-forge nbsphinx"
 - sphinx_rtd_theme: conda install -c anaconda sphinx_rtd_theme 
 - sphinx-bootstrap-theme: 'https://github.com/ryan-roemer/sphinx-bootstrap-theme'
 - sphinx automod: 'https://github.com/astropy/sphinx-automodapi'    
 
 
 

