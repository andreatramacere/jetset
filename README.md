[![pip test](https://github.com/andreatramacere/jetset/actions/workflows/test-pip-workflow.yml/badge.svg?branch=v1.3.x)](https://github.com/andreatramacere/jetset/actions/workflows/test-pip-workflow.yml)
[![conda test](https://github.com/andreatramacere/jetset/actions/workflows/test-conda-workflow.yml/badge.svg?branch=v1.3.x)](https://github.com/andreatramacere/jetset/actions/workflows/test-conda-workflow.yml)
![img](./logo/logo_git.png)


JetSeT is an open source C/Python framework to reproduce radiative and accelerative processes acting in relativistic jets, and galactic objects (beamed and unbeamed), 
allowing to fit the numerical models to observed data. The main features of this framework are:

 * handling observed data: re-binning, definition of data sets, bindings to astropy tables and quantities
   definition of complex numerical radiative scenarios: Synchrotron Self-Compton (SSC), external Compton (EC) and EC 
   against the CMB 
 
 * Constraining of the model in the pre-fitting stage, based on accurate  and already published phenomenological trends. 
   In particular, starting from phenomenological parameters, such as spectral indices, peak fluxes and frequencies, and 
   spectral  curvatures, that the code evaluates automatically, the pre-fitting algorithm is able to provide a good 
   starting model,following the phenomenological trends that I have implemented. fitting of multiwavelength SEDs using  
   both frequentist approach (iminuit) and bayesian MCMC sampling (emcee)
 
 * Self-consistent temporal evolution of the plasma under the effect of radiative, accelerative processes, and adiabatic expansion. Both first order and second order (stochastic acceleration) processes are implemented.



## Acknowledgements

If you use this code in any kind of scientific publication please cite the following papers:

* `Tramacere A. 2020`  https://ui.adsabs.harvard.edu/abs/2020ascl.soft09001T/abstract
* `Tramacere A. et al. 2011` http://adsabs.harvard.edu/abs/2011ApJ...739...66T
* `Tramacere A. et al. 2009` http://adsabs.harvard.edu/abs/2009A%26A...501..879T



# Documentation
visit: https://jetset.readthedocs.io/en/latest/

run the notebook on binder: 
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/andreatramacere/jetset/master)
# Installation 
Read the documentation for further details (e.g. installing form source etc...) [here](https://jetset.readthedocs.io/en/latest/install.html)
## Install  JetSeT from Anaconda 
 
 - create a virtual environment (not necessary, but suggested): 
 
    `conda create --name jetset python=3.10 ipython jupyter`
    
    `conda activate jetset`
     
- install the code:
  
  `conda install -c andreatramacere jetset`
  

## Install  JetSeT from pip 

- create a virtual environment (not necessary, but suggested): 
 
  `pip install virtualenv`
 
  `virtualenv -p python3.10 jetset`

  `source jetset/bin/activate`

  `pip install ipython jupyter`

- install the code:

  `pip install jetset>=1.3`


## Licence

JetSeT is released under a 3-clause BSD  license,  for deatils see
[License](https://github.com/andreatramacere/jetset/blob/master/LICENSE.txt) file 


