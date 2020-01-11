.. JetSeT documentation master file




.. image:: _static/logo_large.png
   :width: 400px



:Author: `Andrea Tramacere <andrea.tramacere@gmail.com>`_




`JetSeT` is an open source  C/Python   framework  to reproduce radiative and accelerative processes acting in relativistic jets,
allowing to fit the numerical models to observed data. The main features of this framework are:

* handling observed data: re-binning, definition of data sets, bindings to astropy tables and quantities definition of complex numerical radiative scenarios: Synchrotron Self-Compton (SSC), external Compton (EC) and EC against the CMB

* Constraining of the model in the pre-fitting stage, based on accurate  and already published phenomenological trends. In particular, starting from phenomenological parameters, such as spectral indices, peak fluxes and frequencies, and
  spectral  curvatures, that the code evaluates automatically, the pre-fitting algorithm is able to provide a good
  starting model,following the phenomenological trends that I have implemented. fitting of multiwavelength SEDs using
  both frequentist approach (iminuit) and bayesian MCMC sampling (emcee)

* Self-consistent temporal evolution of the plasma under the effect of radiative and accelerative processes, both first
  order and second order (stochastic acceleration) processes.

.. Important::
    Acknowledgements: if you use this code in any kind of scientific publication please cite the following papers:

    * `Tramacere A. et al. 2011 <http://adsabs.harvard.edu/abs/2011ApJ...739...66T>`_
    * `Tramacere A. et al. 2009 <http://adsabs.harvard.edu/abs/2009A%26A...501..879T>`_
    * `Massaro E. et. al 2006 <http://adsabs.harvard.edu/abs/2006A%26A...448..861M>`_


.. _user-docs:



Documentation
-------------

.. toctree::
   :maxdepth: 1

   installation <install.rst>
   what's new in jetset 1.1.0  <new_release.rst>
   user guide <user_guide/user_guide.rst>
   code documentation (API) <api/modules.rst>


..   jet model <user_guide/jet_model/Jet_example.rst>
..   model fitting  <user_guide/model_fit/fit_example.rst>
..

..   introduction <intro/intro.rst>
..   tutorial <tutorial/tutorial.rst>








Index
------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

