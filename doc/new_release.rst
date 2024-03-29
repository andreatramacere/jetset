What's new in version 1.2.2
===========================

bug fixing 
----------
- hadronic pp: the grid for the equilibrium evolution was not starting from gmin=1, but from the lowest energy of the secondary
- temporal evolution: in the light curves with crossing time, the `(1+z)` factor term was missing in the times column


Introduction of a class for galactic objects (beamed and unbeamed) 
------------------------------------------------------------------
This class is useful for PWN, SNR, or jetted galactic objects
It is still preliminary, and will be improved soon
See this section of the user guide  :ref:`galactic_guide`


quantile confidence range for MCMC plots 
-----------------------------------------
The MCMC plot method now accepts quantities for model confidence band eg:

.. code-block:: python

    mcmc.plot_model(sed_data=sed_data,fit_range=[1E11, 2E27],size=100,quantiles=[0.05,0.95])


What's new in version 1.2.1
---------------------------

Introduction of depending pars
-----------------------------------
Model parameters can be linked via functional dependence.
See this section of the user guide  :ref:`depending_parameters`


Temporal Evolution
-----------------------------------
The python interface to perform self-consistent temporal evolution of leptonic emitters under
acceleration and cooling has been added.
See this section of the user guide  :ref:`temp_ev`

Hadronic pp emission
-----------------------------------
The python interface to perform self-consistent temporal evolution of leptonic emitters under
acceleration and cooling has been added.
See this section of the user guide  :ref:`hadronic_pp_jet_guide`


Theoretical background for SSC model
------------------------------------
A detailed explanation of the theoretical background for SSC/EC model has been added.
See this section of the user guide  :ref:`ssc_th_bkg`


Emitters distributions
-----------------------------
The emitters distribution has be improved. See this section of the user guide  :ref:`custom_emitters_guide`

EBL absorption
-----------------
linking of parameters has been updated See this section of the user guide  :ref:`ebl_model`


Plugins 
-----------------
 - JetSeT plugins to Sherpa (:ref:`sherpa_plugin`, :ref:`sherpa_minimizer_plugin`) 
 - Gammapy plugin  (:ref:`gammapy_plugin`)




