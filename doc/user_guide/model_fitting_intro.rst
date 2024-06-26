.. _model_fitting_intro:


Model Fitting Introduction
==========================
In this section we provide some guidelines and caveats for mode fitting and minimizers in JetSeT.  In general to perform a model fitting to data,

you can use the :class:`.ModelMinimizer` class from the :mod:`.minimizer` module, or you can use specific plugins as for the case of sherpa and gammapy,  documented here:

- :ref:`sherpa_minimizer_plugin`

- :ref:`sherpa_plugin` 

- :ref:`gammapy_plugin`


.. _least_squares:  https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html
.. _iminuit: https://scikit-hep.org/iminuit/ 
.. _emcee: https://emcee.readthedocs.io/en/stable/   
In the following we will describe the frequentist model fitting using the :class:`.ModelMinimizer`  using the ``minuit``  minimizer (wrapping the `iminuit`_ package ) and the ``lsb`` minimizer (wrapping the scipy `least_squares`_
package),  and the Bayesian approach using the :class:`.McmcSampler` interface to ecee `emcee`_ package. 



Data handling
-------------
Before starting with the description of the model fitting implementation in JeSeT, please read  Handling errors and systematics section in the :ref:`data_format`.

Frequentist model fitting
-------------------------
.. _frequentist_model_fitting:


FitModel creation
^^^^^^^^^^^^^^^^^
The first step consists in creating a fit :class:`.FitModel` object. Assuming you have already a :class:`.Jet`, in this example is the object ``prefit_jet`` (either constrained from data, see :ref:`phenom_constr` and the examples in :ref:`model_fitting_examples`, or any instance of the :class:`.Jet` class)

.. code:: ipython3
    
    from jetset.model_manager import  FitModel
    fit_model_lsb=FitModel( jet=prefit_jet, name='SSC-best-fit-lsb',template=None) 



.. note::
   Starting from JetSeT ``version>=1.1.2`` the implementation of the  composite model  :class:`.FitModel`  for setting  parameters requires to specify the model component.
   and this holds also for the `freeze` method and for setting  `fit_range` intervals, and for the methods relate to parameters setting in general.
   :class:`.FitModel` 



Before proceeding with the minimization, it is important to freeze or free the needed parameters, 
and provide fit boundaries. As default, boundaries will be set to the physical boundaries of the 
parameters, but, depending on your particular source and state, you want to provide tighter boundaries eg:

.. code:: ipython3


    fit_model_minuit.freeze('jet_leptonic','z_cosm')
    fit_model_minuit.freeze('jet_leptonic','R_H')
    fit_model_minuit.freeze('jet_leptonic','R')
    fit_model_minuit.jet_leptonic.parameters.R.fit_range=[5E15,1E17]
    fit_model_minuit.jet_leptonic.parameters.gmin.fit_range=[10,1000]
    fit_model_minuit.jet_leptonic.parameters.gmax.fit_range=[5E5,1E7]
    fit_model_minuit.jet_leptonic.parameters.gamma0_log_parab.fit_range=[1E3,1E5]

    fit_model_minuit.jet_leptonic.parameters.beam_obj.fit_range=[5,50]

.. note:: 
    Setting fit_range can speed up and improve the fit convergence and is advised. Anyhow, the ranges should be used judged by the user each time according to the physics of the particular source. Please, have look at 
    the examples in :ref:`model_fitting_examples` to have a broader view of this topic.

A good strategy is to run first a ``lsb`` fit and then, using the same fit_model, run a fit with ``minuit``. Using minuit we notice that sometimes the fit will converge, but the quality will not be enough (valid==false) to run minos. 
Anyhow, as shown in the MCMC sampling, it still possible to estimate asymmetric errors by means of MCMC sampling.


ModelMinimizer creation
^^^^^^^^^^^^^^^^^^^^^^^
The second step consists in the creation of a  :class:`.ModelMinimizer` object



- for ``lsb`` (wrapping: scipy.optimize.least_squares)
 
  .. code:: ipython3
    
    from jetset.minimizer import ModelMinimizer
    model_minimizer=ModelMinimizer('lsb')

- for ``minuit``:
  
  .. code:: ipython3

        from jetset.minimizer import ModelMinimizer
        model_minimizer=ModelMinimizer('minuit')


  .. note::
    For ``minuit``, starting from JetSeT ``version==1.3.0``, ``simplex`` will be run before ``migrad``, this allows a better sampling of the parameter space. You can skip ``simplex`` as follows:

    .. code:: ipython3

       model_minimizer.minimizer.add_simplex=False

both for the case of ``lsb`` and ``minuit`` specific kw for the minimizer can be accessed by the dictionary:



.. code:: ipython3


    model_minimizer.minimizer.conf_dict


eg: for the case of ``lsb`` that dictionary would be:

.. code:: ipython3


    {'xtol': None,
     'ftol': 1e-08,
     'gtol': 1e-08,
     'jac': '3-point',
     'loss': 'linear',
     'tr_solver': None,
     'f_scale': 1}

.. warning::
    avoid to change those parameters, unless you are completely aware of what you are doing, and you had a full reading the of the ``lsb`` / ``minuit`` documentation.





Running the minimization process
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The third step will consist in actually calling the minimization process

.. code:: ipython3

    best_fit_res=model_minimizer.fit(fit_model=fit_model,
                                     sed_data=sed_data,
                                     nu_fit_start=1E11,
                                     nu_fit_stop=1E29,
                                     max_ev=None,
                                     use_UL=True,
                                     fitname='SSC-best-fit',
                                     repeat=1)

notice that:

-   ``sed_data`` is your :class:`.ObsData` object (see :ref:`data_format` section from more info)


-  ``nu_fit_start`` and ``nu_fit_stop`` correspond to the range of the data, in Hz, used for the model fitting

-  ``use_UL`` if set to ``True`` will take into account also the upper limits, and will skip them in ``use_UL=False``

-   ``repeat``, introduced in version 1.1.2, allows to repeat the fit process, and will provide a better fit convergence. For example, setting ``repeat=3`` the fit process will be repeated 3 times.



See the examples in :ref:`model_fitting_examples`  to inspect the output for each minimizer. The plot for the covariance matrix can be obtained

.. code:: ipython3

    p=model_minimizer.plot_corr_matrix()



Bayesian model fitting with emcee
---------------------------------
.. _bayesian_model_fitting:

Building the mcmc object
^^^^^^^^^^^^^^^^^^^^^^^^
The  :class:`.McmcSampler` interface to emcee `emcee`_ package, in the following we show how to perform a sampling of the
model parameter space, starting from a frequentist best fit result.

it could either a:

-  previously saved instance of a model minimizer 
   
   .. code:: ipython3
   
        from jetset.mcmc import McmcSampler
        model_minimizer = ModelMinimizer.load_model('model_minimizer_minuit.pkl')
        mcmc=McmcSampler(model_minimizer)


- or  any model minimizer instance created in your notebook/script. We create the ``mcmc`` object
   
  .. code:: ipython3
   
        from jetset.mcmc import McmcSampler
        mcmc=McmcSampler(model_minimizer)



Setting the labels
^^^^^^^^^^^^^^^^^^
Now we need to set the parameters to perform the sample: 

- we can use the same free parameters used in the frequentist minimizer object

  .. code:: ipython3

        mcmc.set_labels()


- or define a different sample (or a subsample), defining a dictionary where each key correspond to  the name of the model component, and the value of the key to the list of parameters name for that component.
  For example, assuming that the  :class:`.FitModel` member of the :class:`.ModelMinimizer` object has the ``jet_leptonic`` component  

  .. code:: ipython3
    
    labels=['N','B','beam_obj','s','gamma0_log_parab']
    model_name='jet_leptonic'
    use_labels_dict={model_name:labels}

    mcmc.set_labels(use_labels_dict=use_labels_dict)

  in this case the sampler will use only the ``['N','B','beam_obj','s','gamma0_log_parab']`` par from `model_minimizer.fit_model.jet_leptonic` component


Setting the priors
^^^^^^^^^^^^^^^^^^
Now we need to set the priors. For the current release we are using flat priors centered on the best fit values. 

We provide three different approaches 

- Relative bounds: 
 
  .. code:: ipython3

    mcmc.set_bounds(bound=5.0,bound_rel=True)


  setting ``bound=5.0`` and ``bound_rel=True`` means that: 
  
  - the prior interval will be defined as  ``[best_fit_val - delta_m , best_fit_val + delta_p]``

  - with ``delta_p=delta_m=best_fit_val*bound``

  It is possible to define asymmetric boundaries e.g. ``bound=[2.0,5.0]`` meaning that: 

  - ``delta_p = min(par.best_fit_val*bound[1], par.fit_range_max)``
 
  - ``delta_m = max(par.best_fit_val*bound[0], par.fit_range_min)``

- Absolute bounds:
  
  .. code:: ipython3
    
     mcmc.set_bounds(bound=5.0,bound_rel=False)


  setting ``bound=5.0`` and ``bound_rel=False`` means that: 
  
  - the prior interval will be defined as  ``[best_fit_val - delta_m , best_fit_val + delta_p]``

  - with ``delta_p = delta_m = best_fit_err*bound``

  It is possible to define asymmetric boundaries e.g. ``bound=[2.0,5.0]`` meaning that:
 
  - ``delta_p = par.best_fit_err*bound[1]``
  
  -  ``delta_m = par.best_fit_err*bound[0]``



Running the sampler
^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    mcmc.run_sampler(nwalkers=20, burnin=50,steps=500,progress='notebook')

The number of walkers ``nwalkers``, if not specified, is set to :math:`4 \times` (``number of sampled parameters``), if you pass a lower number, a ``RuntimeError`` will be raised.
Anyhow,  ``nwalkers``, ``burnin`` and ``steps``, should be chosen depending on the particular analysis, and is strongly advised to read the  `emcee`_  documentation.


Please, read the MCMC section of  examples in :ref:`model_fitting_examples` for a showcase of the usage and features of the   :class:`.McmcSampler`.