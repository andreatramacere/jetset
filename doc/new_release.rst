What's new in version 1.1.2
===========================

This release fixes a bug for the MCMC methiod and adds sevaral new features


In the following are listed the main new features common to 1.1.1 and 1.1.2



#. EBL absorption

   EBL models are implemented using a 2D interpolation where the x and y axes represent the redshift and the frequency, and the z axes represents the value of :math:`e^{-\tau}`
   Included models are

   * Franceschini 2008 :cite:`Franceschini2008`
   * Finke 2010 :cite:`Finke2010`
   * Dominguez 2011 :cite:`Dominguez2011`



#. Composite models

   The `FitModel` is now able to handle additive and multiplicative models, giving to the user to define the composite model expression using a single instruction e.g.

   .. code-block:: python

       composite_model.composite_expr='(jet_flaring + steady_jet) * Franceschini_2008'




#. Custom emitters distributions.

   The user can difine custom distribution of emitting particles



#. Improved fit interface

   The new intefrace allows to repeat the fit for a desired numbers of time, in ordert to get a better convergence

   .. code-block:: python

      model_minimizer_lsb=ModelMinimizer('lsb')
      best_fit_lsb=model_minimizer_lsb.fit(fit_model_lsb,sed_data,1E11,1E29,fitname='SSC-best-fit-minuit',repeat=3)


.. important::
    starting from version 1.1.2 the `FitModel` class has been improved to handle parameters from different models, allowing
    to link parameters among different models, or to set the same parameter for different models see the :ref:`composite_models`



.. bibliography:: references.bib
