.. _jetset_index:




.. image:: _static/logo_large_no_border.png
   :width: 400px



:Author: `Andrea Tramacere <andrea.tramacere@gmail.com>`_

.. _least_squares:  https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html
.. _iminuit: https://scikit-hep.org/iminuit/ 
.. _emcee: https://emcee.readthedocs.io/en/stable/   
.. _astropy: https://www.astropy.org/
.. _gammapy: https://gammapy.org/
.. _scipy: https://scipy.org/
.. _numpy: https://numpy.org/
.. _matplotlib: https://matplotlib.org/stable/

JetSeT is an open source  C/Python   framework  to reproduce radiative and accelerative processes acting in relativistic jets, and galactic objects (beamed and unbeamed),
allowing to fit the numerical models to observed data. The main features of this framework are:

* handling observed data: re-binning, definition of data sets, bindings to astropy tables and quantities definition of complex numerical radiative scenarios: Synchrotron Self-Compton (SSC), external Compton (EC) and EC against the CMB

* Constraining of the model in the pre-fitting stage, based on accurate  and already published phenomenological trends. In particular, starting from phenomenological parameters, such as spectral indices, peak fluxes and frequencies, and
  spectral  curvatures, that the code evaluates automatically, the pre-fitting algorithm is able to provide a good
  starting model,following the phenomenological trends that I have implemented. fitting of multiwavelength SEDs using
  both frequentist approach (`iminuit`_, scipy `least_squares`_) and bayesian MCMC sampling (`emcee`_)

* Self-consistent temporal evolution of the plasma under the effect of radiative, accelerative processes, and adiabatic expansion. Both first order and second order (stochastic acceleration) processes are implemented.

.. Important::
    Acknowledgements: if you use this code in any kind of scientific publication please cite the following papers:

    * `Tramacere A. 2020 <https://ui.adsabs.harvard.edu/abs/2020ascl.soft09001T/abstract>`_
    * `Tramacere A. et al. 2011 <http://adsabs.harvard.edu/abs/2011ApJ...739...66T>`_
    * `Tramacere A. et al. 2009 <http://adsabs.harvard.edu/abs/2009A%26A...501..879T>`_

    Please, consider also citing `astropy`_, `gammapy`_, `scipy`_, `iminuit`_, `emcee`_, and `matplotlib`_, if you use functionalities involving the corresponding package 

.. _user-docs:

.. toctree::
   :maxdepth: 1
   :caption: Documentation:

   installation <install.rst>
   what's new in JetSeT v1.3.0? <change_log_specs/v1.3.0>
   changelog  <changelog.rst>
   user guide <user_guide/user_guide.rst>
   code documentation (API) <api/modules.rst>
   bibliography <references.rst>



.. nbgallery::
    :caption: New/updated in v 1.2.0-1.3.0:
    :name: rst-gallery

    documentation_notebooks/notebooks/depending_pars/depending_pars.ipynb
    user_guide/documentation_notebooks/notebooks/composite_model/Composite_model.rst
    user_guide/documentation_notebooks/notebooks/custom_emitters_distr/custom_emitters.rst
    user_guide/documentation_notebooks/notebooks/jet_model_phys_SSC/Jet_example_phys_SSC.rst
    user_guide/documentation_notebooks/notebooks/hadronic_pp_jet/hadronic.rst
    user_guide/documentation_notebooks/notebooks/temporal_evolution/Temp_Ev_one_zone_only_cooling.rst
    user_guide/documentation_notebooks/notebooks/temporal_evolution/Temp_Ev_two_zones_acc_and_cooling.rst
    user_guide/documentation_notebooks/notebooks/temporal_evolution/Temp_Ev_two_zones_acc_and_cooling_adb_exp.rst
    user_guide/documentation_notebooks/notebooks/phen_constr/SSC_th_bkg.rst
    user_guide/documentation_notebooks/notebooks/sherpa_plugin/sherpa-plugin-sherpa-interface.rst
    user_guide/documentation_notebooks/notebooks/sherpa_plugin/sherpa-plugin-jetset-interface.rst
    user_guide/documentation_notebooks/notebooks/gammapy_plugin/gammapy_plugin.rst
    user_guide/documentation_notebooks/notebooks/galactic/galactic.rst

Testers:

I would like to thank the most active testers:

Hubing Xiao, Cosimo Nigro, Vaidehi S. Paliya, Sara Buson, Sonal Patel, Jayant Abhir, Axel Arbet-Engels, Jessica Luna Cervantes


License:

JetSeT is released under a 3-clause BSD  license - for details see the
`License <https://github.com/andreatramacere/jetset/blob/master/LICENSE.txt>`_ file










