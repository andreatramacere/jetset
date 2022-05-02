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

* Self-consistent temporal evolution of the plasma under the effect of radiative, accelerative processes, and adiabatic expansion. Both first order and second order (stochastic acceleration) processes are implemented.

.. Important::
    Acknowledgements: if you use this code in any kind of scientific publication please cite the following papers:

    * `Tramacere A. 2020 <https://ui.adsabs.harvard.edu/abs/2020ascl.soft09001T/abstract>`_
    * `Tramacere A. et al. 2011 <http://adsabs.harvard.edu/abs/2011ApJ...739...66T>`_
    * `Tramacere A. et al. 2009 <http://adsabs.harvard.edu/abs/2009A%26A...501..879T>`_



.. _user-docs:

.. toctree::
   :maxdepth: 1
   :caption: Documentation:

   installation <install.rst>
   what's new in jetset 1.2.0  <new_release.rst>
   user guide <user_guide/user_guide.rst>
   code documentation (API) <api/modules.rst>
   bibliography <references.rst>



.. nbgallery::
    :caption: New/updated in v1.2.0:
    :name: rst-gallery

    user_guide/documentation_notebooks/depending_pars/depending_pars.rst
    user_guide/documentation_notebooks/composite_model/Composite_model.rst
    user_guide/documentation_notebooks/custom_emitters_distr/custom_emitters.rst
    user_guide/documentation_notebooks/jet_model_phys_SSC/Jet_example_phys_SSC.rst
    user_guide/documentation_notebooks/hadronic_pp_jet/hadornic.rst
    user_guide/documentation_notebooks/temporal_evolution/Temp_Ev_one_zone_only_cooling.rst
    user_guide/documentation_notebooks/temporal_evolution/Temp_Ev_two_zones_acc_and_cooling.rst
    user_guide/documentation_notebooks/temporal_evolution/Temp_Ev_two_zones_acc_and_cooling_adb_exp.rst
    user_guide/documentation_notebooks/phen_constr/SSC_th_bkg.rst
    user_guide/documentation_notebooks/sherpa_plugin/sherpa-plugin-sherpa-interface.rst
    user_guide/documentation_notebooks/sherpa_plugin/sherpa-plugin-jetset-interface.rst
    user_guide/documentation_notebooks/gammapy_plugin/gammapy_plugin.rst

Testers:

I would like to thanks the most active testers:

Hubing Xiao, Cosimo Nigro, Vaidehi S. Paliya, Sara Buson, Sonal Patel


License:

JetSeT is released under a 3-clause BSD  license - for details see the
`License <https://github.com/andreatramacere/jetset/blob/master/LICENSE.txt>`_ file










