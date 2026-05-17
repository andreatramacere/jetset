
.. _model_fitting_intro:

.. _least_squares:  https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html
.. _iminuit: https://scikit-hep.org/iminuit/ 
.. _emcee: https://emcee.readthedocs.io/en/stable/   
.. _ultranest: https://johannesbuchner.github.io/UltraNest/index.html
.. _astropy: https://www.astropy.org/
.. _gammapy: https://gammapy.org/
.. _scipy: https://scipy.org/
.. _numpy: https://numpy.org/
.. _matplotlib: https://matplotlib.org/stable/

Model Fitting Introduction
==========================
In this section we provide some guidelines and caveats for mode fitting and minimizers in JetSeT.  In general to perform a model fitting to data, you can use the :class:`.ModelMinimizer` class from the :mod:`.minimizer` module, or you can use specific plugins as for the case of sherpa and gammapy,  documented here:

- :ref:`sherpa_minimizer_plugin`

- :ref:`sherpa_plugin` 

- :ref:`gammapy_plugin`

   

In the following we will describe the 

- frequentist model fitting using the :class:`.ModelMinimizer`  using the ``minuit``  minimizer (wrapping the `iminuit`_ package ) and the ``lsb`` minimizer (wrapping the scipy `least_squares`_ package),

- the Bayesian approach using :class:`.McmcSampler` (backend: `emcee`_) and ``UltraNestSampler`` (backend: `ultranest`_).



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



Bayesian model fitting with emcee and ultranest
------------------------------------------------
.. _bayesian_model_fitting:

JetSeT provides two Bayesian samplers with a common workflow:

- :class:`.McmcSampler` (backend: `emcee <https://emcee.readthedocs.io/en/stable/>`_)
- ``UltraNestSampler`` from :mod:`jetset.mcmc_ultranest` (backend: `UltraNest <https://johannesbuchner.github.io/UltraNest/index.html>`_)

Both samplers start from a :class:`.ModelMinimizer`, use the same parameter/bounds setup, and share the same post-processing helpers
(``sampler_parameters``, ``plot_chain``, ``corner_plot``, ``plot_model``, ``save``/``load``).

Building the sampler object
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Starting from a frequentist best-fit result:
""""""""""""""""""""""""""""""""""""""""""""
.. code:: ipython3

    from jetset.minimizer import ModelMinimizer
    model_minimizer = ModelMinimizer.load_model('model_minimizer_minuit.pkl')


- emcee backend:

  .. code:: ipython3

      from jetset.mcmc import McmcSampler
      mcmc = McmcSampler(model_minimizer)

- ultranest backend:

  .. code:: ipython3

      from jetset.mcmc_ultranest import UltraNestSampler
      mcmc = UltraNestSampler(model_minimizer)


As default, the mcmc sampler will inherit ``nu_fit_start``, ``nu_fit_stop``, and ``use_UL`` from the ``model_minimizer`` object. To change these values, you can use the following code, setting the values according to your choice.

.. code:: ipython3
    
    model_minimizer.prepare_fit(fit_model,
                                sed_data,
                                nu_fit_start=1E7,
                                nu_fit_stop=1E29,
                                use_UL=True)


starting from a plain model
"""""""""""""""""""""""""""
It is also possible to run Bayesian sampling starting from a plain model (without a previous frequentist minimization) by creating a ``ModelMinimizer('mcmc')`` and calling ``prepare_fit``.
This plain-model workflow is suggested mainly for ``UltraNestSampler``.

.. code:: ipython3

    from jetset.minimizer import ModelMinimizer
    model_minimizer = ModelMinimizer('mcmc')
    model_minimizer.prepare_fit(fit_model,
                                sed_data,
                                nu_fit_start=1E7,
                                nu_fit_stop=1E29,
                                use_UL=True)

    from jetset.mcmc_ultranest import UltraNestSampler
    mcmc = UltraNestSampler(model_minimizer)


New interface after the MCMC refactoring
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. important::
   Starting from v1.4.0, :class:`.McmcSampler` does not use ``labels`` and ``use_labels_dict`` anymore.
   All the ``free`` parameters are sampled. To include/exclude parameters, ``freeze/free`` them.

You can inspect sampler parameters with:

.. code:: ipython3

    mcmc.sampler_parameters

To exclude/include parameters from the sampler:

.. code:: ipython3

    mcmc.model.radio_spectrum.parameters.nu_ssa.freeze()
    mcmc.model.radio_spectrum.parameters.nu_ssa.free()


Setting the priors/bounds
^^^^^^^^^^^^^^^^^^^^^^^^^
Priors are uniform within the bounds set for each sampled parameter.

- Global bounds for all currently free parameters:

  .. code:: ipython3

      mcmc.set_bounds(bound=5.0, bound_rel=True)

- Custom bounds for specific parameters:

  .. code:: ipython3

      mcmc.set_bounds(comp_name='jet_leptonic', par_name='N', par_bounds=[1E-5,10])
      mcmc.set_bounds(comp_name='radio_spectrum', par_name='alpha_radio', par_bounds=[0,1])

By default ``preserve_fit_range=True``; this keeps MCMC bounds consistent with the fit ranges set in the frequentist stage.

.. warning::
   For ``UltraNestSampler`` all sampled parameters must have finite and valid ``[min,max]`` bounds; otherwise a ``RuntimeError`` is raised.


Running the sampler
^^^^^^^^^^^^^^^^^^^
- emcee:

  .. code:: ipython3

      mcmc.run_sampler(nwalkers=20, burnin=50, steps=500, progress='notebook')

- ultranest:

  .. code:: ipython3

      mcmc.run_sampler(min_num_live_points=400,
                       dlogz=0.5,
                       frac_remain=0.01,
                       ultranest_output_dir='ultranest_runs/mrk421_plain')

  ``ultranest_output_dir`` controls where UltraNest writes its run files.

  - ultranest with  Open MPI helper:
  
    .. important::
 
      To run ultranest via openmpi, you have to install
        - mpi4py (pip/mamba install mpi4py)
        - OpenMPI (https://docs.open-mpi.org/en/main/installing-open-mpi/quickstart.html#binary-packages)

    .. code:: ipython3

      from jetset.mcmc_ultranest import run_open_mpi
      mcmc = run_open_mpi(mcmc,
                          n_proc=8,
                          ultranest_output_dir='ultranest_runs/mrk421_openmpi')
      
    
    ``ultranest_output_dir`` controls where UltraNest writes its run files, and will be placed under the ``run_mpi`` directory.



Corner plot for composite models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For composite models, ``corner_plot`` can be used in two modes:

- ``per_component=True`` (default): one corner plot per model component.

  .. code:: ipython3

      f = mcmc.corner_plot(per_component=True)

- ``per_component=False``: a single global corner plot including parameters from all sampled components.

  .. code:: ipython3

      f = mcmc.corner_plot(per_component=False)


Caveats and differences: emcee vs ultranest
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- Common caveat:

  In all cases, the quality of parameter bounds is critical for convergence and for physically meaningful posteriors.

- emcee specific caveats:

  - It is strongly recommended to start from a converged frequentist best-fit.
  - ``nwalkers`` default is :math:`4 \times` (number of sampled parameters).
  - ``nwalkers`` must be at least :math:`2 \times` (number of sampled parameters), otherwise a ``RuntimeError`` is raised.
  - ``burnin`` and ``steps`` must be tuned case by case.
  - After the run, you can retune burn-in using autocorrelation time:

    .. code:: ipython3

        mcmc.tune_burnin()

    This recomputes the posterior samples after updating ``burnin`` from ``sampler.get_autocorr_time``.
  - ``tune_burnin`` is meaningful only after ``mcmc.run_sampler(...)`` has been executed.
  - Autocorrelation-time estimates can be unstable for short or non-converged chains; if needed, increase ``steps`` and rerun before trusting the tuned burn-in.
  - Always inspect chains (e.g. with ``mcmc.plot_chain()``) before and after ``tune_burnin()``.
  - For advanced usage and diagnostics, see the emcee autocorrelation documentation:
    `https://emcee.readthedocs.io/en/stable/tutorials/autocorr/ <https://emcee.readthedocs.io/en/stable/tutorials/autocorr/>`_.

- ultranest specific caveats:

  - Starting from a minimized model is not mandatory; a plain-model workflow is supported and can be effective.
  - If you start from a plain model, use physically motivated finite bounds and robust sampling settings (e.g. not quick-test values for ``min_num_live_points``, ``nsteps``, ``min_ess``).
  - There is no burn-in phase (internally ``burnin=0``); ``tune_burnin`` is not part of the ultranest workflow.
  - ``acceptance_fraction`` is not used as a diagnostic (set to ``nan``).
  - UltraNest provides Bayesian evidence estimates:

    .. code:: ipython3

        print(mcmc.logz, mcmc.logzerr)

  - Notebook configurations such as very low ``min_num_live_points``/``nsteps`` are quick-test setups only; for robust analyses use stricter/default settings.
 
    .. important::
    To run ultranest via openmpi, you have to install
        - mpi4py (pip/mamba install mpi4py)
        - OpenMPI (https://docs.open-mpi.org/en/main/installing-open-mpi/quickstart.html#binary-packages)

Please, read the Bayesian sections in :ref:`model_fitting_examples` for end-to-end workflows with both backends.
