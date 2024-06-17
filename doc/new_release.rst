
Changelog
=========


Version 1.3.0
------------- 

New features
^^^^^^^^^^^^

- Introduction of C threads: Starting from version 1.3.0 the :class:`.JetBase` class and all the derived classes, perform the C computation using threads. 
  This increase the computational speed Each time you create a new Jet object, you will get a log noticing how many C threads have been created.
  The number of threads is automatically determined according to the number of cores and threads of your CPU.
  You can revert to single thread, or set a custom number of threads, using the :meth:`.JetBase.set_num_c_threads` and passing the number of threads:

    .. code-block:: python

        from jetset.jet_model import Jet
        my_jet=Jet(electron_distribution='lppl')
        
        my_jet.set_num_c_threads(1)

    You can increase or decrease the number of C threads using the same method.
    It is not advised to exceed the actual number of threads offered by your CPU, 
    since will not increase the speed performance.

- Improved MCMC functionalities
  
  - labels setting and bounds have been placed in specific methods
   
    .. code-block:: python

        from jetset.mcmc import McmcSampler
        model_minimizer = ModelMinimizer.load_model('model_minimizer_minuit.pkl')
        mcmc=McmcSampler(model_minimizer)
        
        labels=['N','B','beam_obj','s','gamma0_log_parab']
        model_name='jet_leptonic'
        use_labels_dict={model_name:labels}

        mcmc.set_labels(use_labels_dict=use_labels_dict)

        mcmc.set_bounds(bound=5.0,bound_rel=True)

  - added option to plot model and residuals, for the mcmc best fit (obtained setting the pars to their  quantiles=0.5 posterior), instead of the frequentist best-fit model, passing ``plot_mcmc_best_fit_model=True``
    
    .. code-block:: python

        p=mcmc.plot_model(sed_data=sed_data,fit_range=[1E11, 2E27],size=100,quantiles=[0.05,0.95], plot_mcmc_best_fit_model=True)

  - customization of lables for plotting

    .. code-block:: python

        mcmc.set_plot_label('N',r'$N$')
        mcmc.set_plot_label('B',r'$B$')
        mcmc.set_plot_label('beam_obj',r'$\delta$')
        mcmc.set_plot_label('s',r'$s$')
        mcmc.set_plot_label('gamma0_log_parab',r'$\gamma_0$')

- Updated EBL models
  
  - Franceschini et al. (2008) [Franceschini2008]_
  
  - Finke et al. (2010) [Finke2010]_ 
  
  - Dominguez et al. (2011) [Dominguez2011]_

  - Dominguez & Saldana-Lopez (2023) [Dominguez2023]_, [Saldana-Lopez2021]_

- Improved EC computation for large angles for anisotropic external fields

- Added convenience methods for conical jet and EC fields:
  
  - method :meth:`.JetBase.make_conical_jet` class will set parameters dependencies to have  conical jet constraining the blob radius

  - method :meth:`.JetBase.set_EC_dependencies` class  will set parameters dependencies to have scaling relations between BLR and DT radius and disk luminosity
  
- Improved dependent parameters: handling of astropy units has been improved  the functional dependency of the parameters

- Improved serialization: saved models will not break if astropy or numba break their interface in future releases

- Heat map for correlation matrix: added convenience method to plot correlation matrix


bug fixing 
^^^^^^^^^^
- fixed typo in `EC_components_list` kw in the method :meth:`.ObsConstrain.constrain_SSC_EC_model`

Version 1.2.2
------------- 



New features
^^^^^^^^^^^^
- Introduction of a class for galactic objects (beamed and unbeamed) .This class is useful for PWN, SNR, or jetted galactic objects. It is still preliminary, and will be improved soon See this section of the user guide  :ref:`galactic_guide`


- quantile confidence range for MCMC plots: the MCMC plot method now accepts quantities for model confidence band eg:

 .. code-block:: python

    mcmc.plot_model(sed_data=sed_data,fit_range=[1E11, 2E27],size=100,quantiles=[0.05,0.95])



bug fixing 
^^^^^^^^^^
- hadronic pp: the grid for the equilibrium evolution was not starting from gmin=1, but from the lowest energy of the secondary
- temporal evolution: in the light curves with crossing time, the `(1+z)` factor term was missing in the times column


Version 1.2.1
-------------

New features
^^^^^^^^^^^^

- Introduction of depending pars: model parameters can be linked via functional dependence.  See this section of the user guide  :ref:`depending_parameters`


- Temporal Evolution: the python interface to perform self-consistent temporal evolution of leptonic emitters under acceleration and cooling has been added. See this section of the user guide  :ref:`temp_ev`

- Hadronic pp emission: the python interface to perform self-consistent temporal evolution of leptonic emitters under acceleration and cooling has been added. See this section of the user guide  :ref:`hadronic_pp_jet_guide`


- Theoretical background for SSC model: a detailed explanation of the theoretical background for SSC/EC model has been added. See this section of the user guide  :ref:`ssc_th_bkg`


- Emitters distributions: the emitters distribution class has be improved. See this section of the user guide  :ref:`custom_emitters_guide`

- EBL absorption: linking of parameters has been updated See this section of the user guide  :ref:`ebl_model`


- Plugins:
 - JetSeT plugins to Sherpa (:ref:`sherpa_plugin`, :ref:`sherpa_minimizer_plugin`) 
 - Gammapy plugin  (:ref:`gammapy_plugin`)




