.. JetSeT documentation master file




.. image:: _static/logo_large.png
   :width: 400px



:Author: `Andrea Tramacere <andrea.tramacere@gmail.com>`_




JetSeT is an open source  C/Python   framework  to reproduce radiative and accelerative processes acting in relativistic jets,
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

    * `Tramacere A. 2020 <https://ui.adsabs.harvard.edu/abs/2020ascl.soft09001T/abstract>`_
    * `Tramacere A. et al. 2011 <http://adsabs.harvard.edu/abs/2011ApJ...739...66T>`_
    * `Tramacere A. et al. 2009 <http://adsabs.harvard.edu/abs/2009A%26A...501..879T>`_
    * `Massaro E. et. al 2006 <http://adsabs.harvard.edu/abs/2006A%26A...448..861M>`_



.. _user-docs:


Documentation
-------------


.. toctree::
   :maxdepth: 1

   installation <install.rst>
   what's new in jetset 1.1.2  <new_release.rst>
   user guide <user_guide/user_guide.rst>
   code documentation (API) <api/modules.rst>


Testers
-------
I would like to thanks the most active testers:

Hubing Xiao, Cosimo Nigro, Vaidehi S. Paliya, Sara Buson, Sonal Patel


License
-------
JetSeT is released under a 3-clause BSD  license - for details see the
`License <https://github.com/andreatramacere/jetset/blob/master/LICENSE.txt>`_ file










