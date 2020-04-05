What's new in version 1.1.2
===========================

In the following are listed the main new features common to 1.1.1 and 1.1.2

Bugfix on MCMC and moved to mcmc v3
-----------------------------------
This release fixes a bug for the MCMC method, and we have moved to emcee v3

EBL absorption
-----------------
EBL models are implemented using a 2D interpolation where the x and y axes represent the redshift and the frequency, and the z axes represents the value of :math:`e^{-\tau}`.
(See this section of the user guide  :ref:`ebl_model`)

Included models are

* Franceschini 2008 :cite:`Franceschini2008`
* Finke 2010 :cite:`Finke2010`
* Dominguez 2011 :cite:`Dominguez2011`



Composite models
-----------------
The `FitModel` is now able to handle additive and multiplicative models, giving to the user the possibility to define the composite model expression using a single instruction e.g.
(See this section of the user guide  :ref:`composite_models`)

.. code-block:: python

   composite_model.composite_expr='(jet_flaring + steady_jet) * Franceschini_2008'






Custom emitters distributions
-----------------------------
The user can define custom distributions of emitting particles. See this section of the user guide  :ref:`custom_emitters_guide`



Improved fit interface
-----------------------------
The new intefrace allows to repeat the fit for a desired numbers of time, in order to get a better convergence (See this section of the user guide  :ref:`model_fitting_1`)

.. code-block:: python

  model_minimizer_lsb=ModelMinimizer('lsb')
  best_fit_lsb=model_minimizer_lsb.fit(fit_model_lsb,sed_data,1E11,1E29,fitname='SSC-best-fit-minuit',repeat=3)










References
----------



.. bibliography:: references.bib
