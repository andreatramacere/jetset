
.. _jet_physical_guide:



physical setup
==============


In this section we describe how  to build a model of jet able to reproduce SSC/EC emission processes, using the :class:`.Jet` class from the :mod:`.jet_model` module. to This class thorough a flexible and intuitive interface allows to access the C numerical code that provides an accurate and fast computation of the synchrotron and inverse Compton processes.  

basic setup
-----------

A jet can be built using the  the :class:`.Jet` class, istanciating a jet object, in the following way:

.. code:: ipython3

    %matplotlib inline
    from jetset.jet_model import Jet
    my_jet=Jet(name='test',electron_distribution='lppl',)

This instruction will create:
    * a ``Jet`` object with ``name`` **test**,
    * using as electron distribution the **lppl** model, that is a log-parabola with a low-energy power-law branch.
    * using as working directory **test_jet_prod**

For a list of possible distribution you can run the command 

.. code:: ipython3

    Jet.available_electron_distributions()


.. parsed-literal::

    lp: log-parabola: $K*\gamma^p$
    pl: powerlaw
    lppl: log-parabola with low-energy powerlaw branch
    lpep: log-parabola defined by peak energy
    plc: powerlaw with cut-off
    bkn: broken powerlaw
    spitkov: 
    lppl_pile_up: 
    bkn_pile_up: 




The parameters of the model are accessible throug the instruction

.. code:: ipython3

    my_jet.show_pars()


.. parsed-literal::

    -------------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                 | Units            | value         | phys. boundaries              | log
    -------------------------------------------------------------------------------------------------------------------
     B                | magnetic_field       | G                | +1.000000e-01 | [+0.000000e+00,No           ] | False 
     N                | electron_density     | cm^-3            | +1.000000e+02 | [+0.000000e+00,No           ] | False 
     R                | region_size          | cm               | +1.569897e+01 | [+0.000000e+00,+3.000000e+01] | True 
     beam_obj         | beaming              |                  | +1.000000e+01 | [+1.000000e+00,No           ] | False 
     gamma0_log_parab | turn-over-energy     | Lorentz-factor   | +1.000000e+04 | [+1.000000e+00,+1.000000e+08] | False 
     gmax             | high-energy-cut-off  | Lorentz-factor   | +1.000000e+06 | [+1.000000e+00,+1.000000e+15] | False 
     gmin             | low-energy-cut-off   | Lorentz-factor   | +2.000000e+00 | [+1.000000e+00,+1.000000e+05] | False 
     r                | spectral_curvature   |                  | +4.000000e-01 | [-1.500000e+01,+1.500000e+01] | False 
     s                | LE_spectral_slope    |                  | +2.000000e+00 | [-1.000000e+01,+1.000000e+01] | False 
     z_cosm           | redshift             |                  | +1.000000e-01 | [+0.000000e+00,No           ] | False 
    -------------------------------------------------------------------------------------------------------------------



Each parameter has default values. All the parameters listed are handled by :class:`.ModelParameterArray`, and each parameter is an instance of the the :class:`.JetParameter`. class


To get a full description of the model you can use the instruction

.. code:: ipython3

    my_jet.show_model()


.. parsed-literal::

    
    -------------------------------------------------------------------------------------------------------------------
    jet model description
    -------------------------------------------------------------------------------------------------------------------
    name: test  
    
    electron distribution:
     type: lppl  
     electron energy grid size:  1001
     gmin grid : 2.000000e+00
     gmax grid : 1.000000e+06
    
    radiative fields:
     seed photons grid size:  100
     IC emission grid size:  50
     source emissivity lower bound :  1.000000e-120
     spectral components:
       name:Sum, state: on
       name:Sync, state: self-abs
       name:SSC, state: on
    
    SED info:
     nu grid size :200
     nu mix (Hz): 1.000000e+06
     nu max (Hz): 1.000000e+30
    
    flux plot lower bound   :  1.000000e-30
    
    -------------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                 | Units            | value         | phys. boundaries              | log
    -------------------------------------------------------------------------------------------------------------------
     B                | magnetic_field       | G                | +1.000000e-01 | [+0.000000e+00,No           ] | False 
     N                | electron_density     | cm^-3            | +1.000000e+02 | [+0.000000e+00,No           ] | False 
     R                | region_size          | cm               | +1.569897e+01 | [+0.000000e+00,+3.000000e+01] | True 
     beam_obj         | beaming              |                  | +1.000000e+01 | [+1.000000e+00,No           ] | False 
     gamma0_log_parab | turn-over-energy     | Lorentz-factor   | +1.000000e+04 | [+1.000000e+00,+1.000000e+08] | False 
     gmax             | high-energy-cut-off  | Lorentz-factor   | +1.000000e+06 | [+1.000000e+00,+1.000000e+15] | False 
     gmin             | low-energy-cut-off   | Lorentz-factor   | +2.000000e+00 | [+1.000000e+00,+1.000000e+05] | False 
     r                | spectral_curvature   |                  | +4.000000e-01 | [-1.500000e+01,+1.500000e+01] | False 
     s                | LE_spectral_slope    |                  | +2.000000e+00 | [-1.000000e+01,+1.000000e+01] | False 
     z_cosm           | redshift             |                  | +1.000000e-01 | [+0.000000e+00,No           ] | False 
    -------------------------------------------------------------------------------------------------------------------
    -------------------------------------------------------------------------------------------------------------------


as you can notice, you can now access further information regarding the model, such as numerical configuration of the grid. These parameters will be discussed 
in the :ref:`jet_numerical_guide' section

setting the parameters
----------------------

assume you want to change some of the parameters in your model, you can use two methods: 

1) using the :class:`.Jet.set_par()` method 

.. code:: ipython3

    my_jet.set_par('B',val=0.2)
    my_jet.set_par('gamma0_log_parab',val=5E3)
    my_jet.set_par('gmin',val=1E2)
    my_jet.set_par('gmax',val=1E8)
    my_jet.set_par('R',val=14.5)
    my_jet.set_par('N',val=1E3)

2) accessing directly the parameter 

.. code:: ipython3

    my_jet.parameters.B.val=0.2
    my_jet.parameters.r.val=0.4

investigating the electron distribution
---------------------------------------

.. code:: ipython3

    my_jet.show_electron_distribution()


.. parsed-literal::

    -------------------------------------------------------------------------------------------------------------------
    electron distribution:
     type: lppl  
     electron energy grid size:  1001
     gmin grid : 2.000000e+00
     gmax grid : 1.000000e+06
    
    -------------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                 | Units            | value         | phys. boundaries              | log
    -------------------------------------------------------------------------------------------------------------------
     N                | electron_density     | cm^-3            | +1.000000e+03 | [+0.000000e+00,No           ] | False 
     gamma0_log_parab | turn-over-energy     | Lorentz-factor   | +5.000000e+03 | [+1.000000e+00,+1.000000e+08] | False 
     gmax             | high-energy-cut-off  | Lorentz-factor   | +1.000000e+08 | [+1.000000e+00,+1.000000e+15] | False 
     gmin             | low-energy-cut-off   | Lorentz-factor   | +1.000000e+02 | [+1.000000e+00,+1.000000e+05] | False 
     r                | spectral_curvature   |                  | +4.000000e-01 | [-1.500000e+01,+1.500000e+01] | False 
     s                | LE_spectral_slope    |                  | +2.000000e+00 | [-1.000000e+01,+1.000000e+01] | False 
    -------------------------------------------------------------------------------------------------------------------


.. code:: ipython3

    p=my_jet.electron_distribution.plot()



.. image:: Jet_example_phys_files/Jet_example_phys_21_0.png


.. code:: ipython3

    p=my_jet.electron_distribution.plot3p()



.. image:: Jet_example_phys_files/Jet_example_phys_22_0.png


.. code:: ipython3

    import numpy as np
    p=None
    for r in np.linspace(0.3,1,10):
        my_jet.parameters.r.val=r
        if p is None:
            p=my_jet.electron_distribution.plot3p()
        else:
            p=my_jet.electron_distribution.plot3p(p)



.. image:: Jet_example_phys_files/Jet_example_phys_23_0.png


using log values for electron distribution parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    my_jet=Jet(name='test',electron_distribution='lppl',electron_distribution_log_values=True)
    my_jet.show_model()


.. parsed-literal::

    
    -------------------------------------------------------------------------------------------------------------------
    jet model description
    -------------------------------------------------------------------------------------------------------------------
    name: test  
    
    electron distribution:
     type: lppl  
     electron energy grid size:  1001
     gmin grid : 2.000000e+00
     gmax grid : 1.000000e+06
    
    radiative fields:
     seed photons grid size:  100
     IC emission grid size:  50
     source emissivity lower bound :  1.000000e-120
     spectral components:
       name:Sum, state: on
       name:Sync, state: self-abs
       name:SSC, state: on
    
    SED info:
     nu grid size :200
     nu mix (Hz): 1.000000e+06
     nu max (Hz): 1.000000e+30
    
    flux plot lower bound   :  1.000000e-30
    
    -------------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                 | Units            | value         | phys. boundaries              | log
    -------------------------------------------------------------------------------------------------------------------
     B                | magnetic_field       | G                | +1.000000e-01 | [+0.000000e+00,No           ] | False 
     N                | electron_density     | cm^-3            | +1.000000e+02 | [+0.000000e+00,No           ] | False 
     R                | region_size          | cm               | +1.569897e+01 | [+0.000000e+00,+3.000000e+01] | True 
     beam_obj         | beaming              |                  | +1.000000e+01 | [+1.000000e+00,No           ] | False 
     gamma0_log_parab | turn-over-energy     | Lorentz-factor   | +4.000000e+00 | [+0.000000e+00,+8.000000e+00] | True 
     gmax             | high-energy-cut-off  | Lorentz-factor   | +6.000000e+00 | [+0.000000e+00,+1.500000e+01] | True 
     gmin             | low-energy-cut-off   | Lorentz-factor   | +3.010300e-01 | [+0.000000e+00,+5.000000e+00] | True 
     r                | spectral_curvature   |                  | +4.000000e-01 | [-1.500000e+01,+1.500000e+01] | False 
     s                | LE_spectral_slope    |                  | +2.000000e+00 | [-1.000000e+01,+1.000000e+01] | False 
     z_cosm           | redshift             |                  | +1.000000e-01 | [+0.000000e+00,No           ] | False 
    -------------------------------------------------------------------------------------------------------------------
    -------------------------------------------------------------------------------------------------------------------



evaluate and plot the model
---------------------------

At this point we can evaluate the emission for this jet model using the
instruction

.. code:: ipython3

    my_jet.eval()

.. code:: ipython3

    my_jet.show_pars()


.. parsed-literal::

    -------------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                 | Units            | value         | phys. boundaries              | log
    -------------------------------------------------------------------------------------------------------------------
     B                | magnetic_field       | G                | +1.000000e-01 | [+0.000000e+00,No           ] | False 
     N                | electron_density     | cm^-3            | +1.000000e+02 | [+0.000000e+00,No           ] | False 
     R                | region_size          | cm               | +1.569897e+01 | [+0.000000e+00,+3.000000e+01] | True 
     beam_obj         | beaming              |                  | +1.000000e+01 | [+1.000000e+00,No           ] | False 
     gamma0_log_parab | turn-over-energy     | Lorentz-factor   | +4.000000e+00 | [+0.000000e+00,+8.000000e+00] | True 
     gmax             | high-energy-cut-off  | Lorentz-factor   | +6.000000e+00 | [+0.000000e+00,+1.500000e+01] | True 
     gmin             | low-energy-cut-off   | Lorentz-factor   | +3.010300e-01 | [+0.000000e+00,+5.000000e+00] | True 
     r                | spectral_curvature   |                  | +4.000000e-01 | [-1.500000e+01,+1.500000e+01] | False 
     s                | LE_spectral_slope    |                  | +2.000000e+00 | [-1.000000e+01,+1.000000e+01] | False 
     z_cosm           | redshift             |                  | +1.000000e-01 | [+0.000000e+00,No           ] | False 
    -------------------------------------------------------------------------------------------------------------------


and plot the corresponding SED:

.. code:: ipython3

    from jetset.plot_sedfit import PlotSED
    my_plot=PlotSED()
    my_jet.plot_model(plot_obj=my_plot)
    my_plot.rescale(y_max=-13,y_min=-17.5,x_min=8)



.. image:: Jet_example_phys_files/Jet_example_phys_32_0.png


alternatively, you can call the ``plot_model`` method without passing a
``Plot`` object

.. code:: ipython3

    my_plot=my_jet.plot_model()
    my_plot.rescale(y_max=-13,y_min=-17.5,x_min=8)



.. image:: Jet_example_phys_files/Jet_example_phys_34_0.png


the ``my_plot`` object returned will be built on the fly by the
``plot_model`` method

if you wanto to have interacitve plot:

1) in a jupyter notebook use:

.. code-block:: no

    %matplotlib notebook


2) in an ipython terminal

.. code-block:: python
    
    from matplotlib import pylab as plt
    plt.ion()

comparing models on the same plot
---------------------------------

to compare the same model after changing a parameter

.. code:: ipython3

    my_jet=Jet(name='test',electron_distribution='lppl',)
    my_jet.set_par('B',val=0.2)
    my_jet.set_par('gamma0_log_parab',val=5E3)
    my_jet.set_par('gmin',val=1E2)
    my_jet.set_par('gmax',val=1E8)
    my_jet.set_par('R',val=14.5)
    my_jet.set_par('N',val=1E3)
    
    my_jet.parameters.gamma0_log_parab.val=1E4
    my_jet.eval()
    my_plot=my_jet.plot_model(label='gamma0_log_parab=1E4',comp='Sum')
    my_jet.set_par('gamma0_log_parab',val=1.0E5)
    my_jet.eval()
    my_plot=my_jet.plot_model(my_plot,label='gamma0_log_parab=1E5',comp='Sum')
    my_plot.rescale(y_max=-13,y_min=-17.5,x_min=8)



.. image:: Jet_example_phys_files/Jet_example_phys_39_0.png


saving a plot
-------------

to save the plot

.. code:: ipython3

    my_plot.save('jet1.png')

saving and lodaing a model
--------------------------

.. code:: ipython3

    my_jet.save_model('test_model.dat')

.. code:: ipython3

    my_jet_new=Jet.load_model('test_model.dat')


.. parsed-literal::

    -------------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                 | Units            | value         | phys. boundaries              | log
    -------------------------------------------------------------------------------------------------------------------
     B                | magnetic_field       | G                | +1.000000e-01 | [+0.000000e+00,No           ] | False 
     N                | electron_density     | cm^-3            | +1.000000e+02 | [+0.000000e+00,No           ] | False 
     R                | region_size          | cm               | +1.569897e+01 | [+0.000000e+00,+3.000000e+01] | True 
     beam_obj         | beaming              |                  | +1.000000e+01 | [+1.000000e+00,No           ] | False 
     gamma0_log_parab | turn-over-energy     | Lorentz-factor   | +1.000000e+04 | [+1.000000e+00,+1.000000e+08] | False 
     gmax             | high-energy-cut-off  | Lorentz-factor   | +1.000000e+06 | [+1.000000e+00,+1.000000e+15] | False 
     gmin             | low-energy-cut-off   | Lorentz-factor   | +2.000000e+00 | [+1.000000e+00,+1.000000e+05] | False 
     r                | spectral_curvature   |                  | +4.000000e-01 | [-1.500000e+01,+1.500000e+01] | False 
     s                | LE_spectral_slope    |                  | +2.000000e+00 | [-1.000000e+01,+1.000000e+01] | False 
     z_cosm           | redshift             |                  | +1.000000e-01 | [+0.000000e+00,No           ] | False 
    -------------------------------------------------------------------------------------------------------------------


switching on/off the particle distribution normalization
--------------------------------------------------------

As default the electron distributions are normalized, i.e. are multiplied by a constant ``N_0``, in such a way that :

:math:`\int_{\gamma_{min}}^{\gamma_{max}} n(\gamma) d\gamma =1`, 

it means the the value `N`, refers to the actual density of emitters.
If you want to chance this behavior, you can start looking at the sate of ``Norm_distr`` flag with the following command

.. code:: ipython3

    my_jet.Norm_distr




.. parsed-literal::

    1



and then you can switch off the normalization withe command

.. code:: ipython3

    my_jet.switch_Norm_distr_OFF()

or set back the normalization on with

.. code:: ipython3

    my_jet.switch_Norm_distr_ON()

setting the particle density from observed Fluxes or Luminosities
-----------------------------------------------------------------

It is possible to set the density of emitting particle starting from some observed luminosity or flux (see the method     :meth:`.Jet.set_N_from_nuFnu`,th:`.Jet.set_N_from_nuLnu`)

.. code:: ipython3

    my_jet=Jet(name='test',electron_distribution='lppl')

this is the initial value of N

.. code:: ipython3

    my_jet.parameters.N.val




.. parsed-literal::

    100.0



we now want to set the value of ``N`` in order that the observed synchrotron flux at a given frequency matches a desired value. 
For example, assume that we wish to set ``N`` in order that  the synchrotron flux at math:`10^{15}` Hz is exactly matching the desired value of :math:`10^{-=14}` ergs cm-2 s-1. We can accomplish this by using the :class:`.Jet.get_par_by_name()` as follows: 

.. code:: ipython3

    
    my_jet.set_N_from_nuFnu(nuFnu_obs=1E-14,nu_obs=1E15)

This is the updated value of ``N``, obtained in order to match the given
flux at the given frequency

.. code:: ipython3

    my_jet.get_par_by_name('N').val




.. parsed-literal::

    249.04461454958587



.. code:: ipython3

    my_jet.parameters.show_pars()


.. parsed-literal::

    -------------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                 | Units            | value         | phys. boundaries              | log
    -------------------------------------------------------------------------------------------------------------------
     N                | electron_density     | cm^-3            | +2.490446e+02 | [+0.000000e+00,No           ] | False 
     gmin             | low-energy-cut-off   | Lorentz-factor   | +2.000000e+00 | [+1.000000e+00,+1.000000e+05] | False 
     gmax             | high-energy-cut-off  | Lorentz-factor   | +1.000000e+06 | [+1.000000e+00,+1.000000e+15] | False 
     s                | LE_spectral_slope    |                  | +2.000000e+00 | [-1.000000e+01,+1.000000e+01] | False 
     r                | spectral_curvature   |                  | +4.000000e-01 | [-1.500000e+01,+1.500000e+01] | False 
     gamma0_log_parab | turn-over-energy     | Lorentz-factor   | +1.000000e+04 | [+1.000000e+00,+1.000000e+08] | False 
     R                | region_size          | cm               | +1.569897e+01 | [+0.000000e+00,+3.000000e+01] | True 
     B                | magnetic_field       | G                | +1.000000e-01 | [+0.000000e+00,No           ] | False 
     beam_obj         | beaming              |                  | +1.000000e+01 | [+1.000000e+00,No           ] | False 
     z_cosm           | redshift             |                  | +1.000000e-01 | [+0.000000e+00,No           ] | False 
    -------------------------------------------------------------------------------------------------------------------


.. code:: ipython3

    my_jet.eval()
    my_plot=my_jet.plot_model(label='set N from F=1E-14')
    my_plot.rescale(y_max=-13,y_min=-17.5,x_min=8)



.. image:: Jet_example_phys_files/Jet_example_phys_63_0.png


as you can see, the synchrotron flux at :math:`10^{15}` Hz is exactly matching the desired value of :math:`10^{-14}` ergs cm-2 s-1.
Alternatively, the value of N  can be obtained using the rest-frame luminosity and  frequency, using the :class:`.Jet.set_N_from_nuLnu()

.. code:: ipython3

    my_jet.set_N_from_nuLnu(L_0=1E43,nu_0=1E15)

where ``L_0`` is the rest-frame luminosity in erg/s at the rest-frame frequency ``nu_0`` in Hz.



## setting the beaming factor

It is possible to set the beaming factor according to the relativistic BulkFactor and viewing angle, this can be done by setting the ``beaming_expr`` kw in the Jet constructor, possible choices are

* `delta` to provide directly the beaming factor (default)
* `bulk_theta` to provide the BulkFactor and the jet  viewing angle 

.. code:: ipython3

    my_jet=Jet(name='test',electron_distribution='lppl',beaming_expr='bulk_theta')

.. code:: ipython3

    my_jet.parameters.show_pars()


.. parsed-literal::

    -------------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                 | Units            | value         | phys. boundaries              | log
    -------------------------------------------------------------------------------------------------------------------
     N                | electron_density     | cm^-3            | +1.000000e+02 | [+0.000000e+00,No           ] | False 
     gmin             | low-energy-cut-off   | Lorentz-factor   | +2.000000e+00 | [+1.000000e+00,+1.000000e+05] | False 
     gmax             | high-energy-cut-off  | Lorentz-factor   | +1.000000e+06 | [+1.000000e+00,+1.000000e+15] | False 
     s                | LE_spectral_slope    |                  | +2.000000e+00 | [-1.000000e+01,+1.000000e+01] | False 
     r                | spectral_curvature   |                  | +4.000000e-01 | [-1.500000e+01,+1.500000e+01] | False 
     gamma0_log_parab | turn-over-energy     | Lorentz-factor   | +1.000000e+04 | [+1.000000e+00,+1.000000e+08] | False 
     R                | region_size          | cm               | +1.569897e+01 | [+0.000000e+00,+3.000000e+01] | True 
     B                | magnetic_field       | G                | +1.000000e-01 | [+0.000000e+00,No           ] | False 
     theta            | jet-viewing-angle    | deg              | +1.000000e-01 | [+0.000000e+00,No           ] | False 
     BulkFactor       | jet-bulk-factor      | Lorentz-factor   | +1.000000e+01 | [+1.000000e+00,No           ] | False 
     z_cosm           | redshift             |                  | +1.000000e-01 | [+0.000000e+00,No           ] | False 
    -------------------------------------------------------------------------------------------------------------------


the actual value of the beaming factor can be obtained using the :meth:`.Jet.get_beaming`

.. code:: ipython3

    my_jet.get_beaming()




.. parsed-literal::

    19.943844732554165



We can change the value of ``theta`` and get the updated value of the beaming factor

.. code:: ipython3

    my_jet.set_par('theta',val=10.)

.. code:: ipython3

    my_jet.get_beaming()




.. parsed-literal::

    4.968041140891955



of course setting ``beaming_expr=delta`` we get the same beaming
expression as in the default case

.. code:: ipython3

    my_jet=Jet(name='test',electron_distribution='lppl',beaming_expr='delta')

.. code:: ipython3

    my_jet.parameters.show_pars()


.. parsed-literal::

    -------------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                 | Units            | value         | phys. boundaries              | log
    -------------------------------------------------------------------------------------------------------------------
     N                | electron_density     | cm^-3            | +1.000000e+02 | [+0.000000e+00,No           ] | False 
     gmin             | low-energy-cut-off   | Lorentz-factor   | +2.000000e+00 | [+1.000000e+00,+1.000000e+05] | False 
     gmax             | high-energy-cut-off  | Lorentz-factor   | +1.000000e+06 | [+1.000000e+00,+1.000000e+15] | False 
     s                | LE_spectral_slope    |                  | +2.000000e+00 | [-1.000000e+01,+1.000000e+01] | False 
     r                | spectral_curvature   |                  | +4.000000e-01 | [-1.500000e+01,+1.500000e+01] | False 
     gamma0_log_parab | turn-over-energy     | Lorentz-factor   | +1.000000e+04 | [+1.000000e+00,+1.000000e+08] | False 
     R                | region_size          | cm               | +1.569897e+01 | [+0.000000e+00,+3.000000e+01] | True 
     B                | magnetic_field       | G                | +1.000000e-01 | [+0.000000e+00,No           ] | False 
     beam_obj         | beaming              |                  | +1.000000e+01 | [+1.000000e+00,No           ] | False 
     z_cosm           | redshift             |                  | +1.000000e-01 | [+0.000000e+00,No           ] | False 
    -------------------------------------------------------------------------------------------------------------------


accessing individual spectral components
----------------------------------------

It is possible to access specific spectral components of our model

.. code:: ipython3

    my_jet=Jet(name='test',electron_distribution='lppl',beaming_expr='bulk_theta')
    my_jet.eval()

We can obtain this information anytime using the :meth:`.Jet.list_spectral_components` method

.. code:: ipython3

    
    my_jet.list_spectral_components()


.. parsed-literal::

    Sum
    Sync
    SSC


the on-screen message is telling us which components have been
evaluated.

and we cann access a specific component using the :meth:`.Jet.get_spectral_component_by_name` method

.. code:: ipython3

    Sync=my_jet.get_spectral_component_by_name('Sync')

OR

.. code:: ipython3

    Sync=my_jet.spectral_components.Sync

and from the ``SED`` object we can extract both the nu and nuFnu array

.. code:: ipython3

    nu_sync=Sync.SED.nu
    nuFnu_sync=Sync.SED.nuFnu

.. code:: ipython3

    print (nuFnu_sync[::10])


.. parsed-literal::

    [1.00000000e-120 1.00000000e-120 1.18346083e-022 1.87412089e-018
     4.45026043e-016 1.78624983e-015 7.07667943e-015 2.69215529e-014
     7.95326288e-014 1.35642311e-013 1.22398936e-013 1.54978292e-014
     4.52069023e-028 1.00000000e-120 1.00000000e-120 1.00000000e-120
     1.00000000e-120 1.00000000e-120 1.00000000e-120 1.00000000e-120]


External Compton
----------------

Broad Line Region
~~~~~~~~~~~~~~~~~

.. code:: ipython3

    my_jet=Jet(name='BLR example',electron_distribution='bkn')
    my_jet.add_EC_component('EC_BLR')
    my_jet.show_model()


.. parsed-literal::

    
    -------------------------------------------------------------------------------------------------------------------
    jet model description
    -------------------------------------------------------------------------------------------------------------------
    name: BLR example  
    
    electron distribution:
     type: bkn  
     electron energy grid size:  1001
     gmin grid : 2.000000e+00
     gmax grid : 1.000000e+06
    
    radiative fields:
     seed photons grid size:  100
     IC emission grid size:  50
     source emissivity lower bound :  1.000000e-120
     spectral components:
       name:Sum, state: on
       name:Sync, state: self-abs
       name:SSC, state: on
       name:EC_BLR, state: on
       name:Disk, state: on
    
    SED info:
     nu grid size :200
     nu mix (Hz): 1.000000e+06
     nu max (Hz): 1.000000e+30
    
    flux plot lower bound   :  1.000000e-30
    
    -------------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                 | Units            | value         | phys. boundaries              | log
    -------------------------------------------------------------------------------------------------------------------
     B                | magnetic_field       | G                | +1.000000e-01 | [+0.000000e+00,No           ] | False 
     L_Disk           | Disk                 | erg/s            | +1.000000e+45 | [+0.000000e+00,No           ] | False 
     N                | electron_density     | cm^-3            | +1.000000e+02 | [+0.000000e+00,No           ] | False 
     R                | region_size          | cm               | +1.569897e+01 | [+0.000000e+00,+3.000000e+01] | True 
     R_BLR_in         | BLR                  | cm               | +1.000000e+18 | [+0.000000e+00,No           ] | False 
     R_BLR_out        | BLR                  | cm               | +2.000000e+18 | [+0.000000e+00,No           ] | False 
     R_H              | Disk                 | cm               | +1.000000e+17 | [+0.000000e+00,No           ] | False 
     R_ext_Sw         | Disk                 | Sw. radii        | +5.000000e+02 | [+0.000000e+00,No           ] | False 
     R_inner_Sw       | Disk                 | Sw. radii        | +3.000000e+00 | [+0.000000e+00,No           ] | False 
     T_Disk           | Disk                 | K                | +1.000000e+05 | [+0.000000e+00,No           ] | False 
     accr_eff         | Disk                 |                  | +1.000000e-01 | [+0.000000e+00,No           ] | False 
     beam_obj         | beaming              |                  | +1.000000e+01 | [+1.000000e+00,No           ] | False 
     gamma_break      | turn-over-energy     | Lorentz-factor   | +1.000000e+04 | [+1.000000e+00,+1.000000e+08] | False 
     gmax             | high-energy-cut-off  | Lorentz-factor   | +1.000000e+06 | [+1.000000e+00,+1.000000e+15] | False 
     gmin             | low-energy-cut-off   | Lorentz-factor   | +2.000000e+00 | [+1.000000e+00,+1.000000e+05] | False 
     p                | LE_spectral_slope    |                  | +2.000000e+00 | [-1.000000e+01,+1.000000e+01] | False 
     p_1              | HE_spectral_slope    |                  | +3.000000e+00 | [-1.000000e+01,+1.000000e+01] | False 
     tau_BLR          | BLR                  |                  | +1.000000e-01 | [+0.000000e+00,+1.000000e+00] | False 
     z_cosm           | redshift             |                  | +1.000000e-01 | [+0.000000e+00,No           ] | False 
    -------------------------------------------------------------------------------------------------------------------
    -------------------------------------------------------------------------------------------------------------------


.. code:: ipython3

    my_jet.set_par('L_Disk',val=1E46)
    my_jet.set_par('gmax',val=1E5)
    my_jet.set_par('gmin',val=2.)
    
    my_jet.set_par('p',val=1.5)
    my_jet.set_par('p_1',val=3.5)
    my_jet.set_par('R',val=15.5)
    my_jet.set_par('B',val=1.0)
    my_jet.set_par('z_cosm',val=0.6)
    my_jet.set_par('beam_obj',val=25)
    my_jet.set_par('gamma_break',val=5E2)
    my_jet.set_N_from_nuLnu(nu_0=1E14,L_0=1E45)

.. code:: ipython3

    my_jet.eval()
    p=my_jet.plot_model()
    p.rescale(y_min=-15,y_max=-10,x_min=8,x_max=27)



.. image:: Jet_example_phys_files/Jet_example_phys_94_0.png


Dusty Torus
~~~~~~~~~~~

.. code:: ipython3

    my_jet.add_EC_component('DT')
    my_jet.show_model()


.. parsed-literal::

    
    -------------------------------------------------------------------------------------------------------------------
    jet model description
    -------------------------------------------------------------------------------------------------------------------
    name: BLR example  
    
    electron distribution:
     type: bkn  
     electron energy grid size:  1001
     gmin grid : 2.000000e+00
     gmax grid : 1.000000e+05
    
    radiative fields:
     seed photons grid size:  100
     IC emission grid size:  50
     source emissivity lower bound :  1.000000e-120
     spectral components:
       name:Sum, state: on
       name:Sync, state: self-abs
       name:SSC, state: on
       name:EC_BLR, state: on
       name:Disk, state: on
       name:DT, state: on
    
    SED info:
     nu grid size :200
     nu mix (Hz): 1.000000e+06
     nu max (Hz): 1.000000e+30
    
    flux plot lower bound   :  1.000000e-30
    
    -------------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                 | Units            | value         | phys. boundaries              | log
    -------------------------------------------------------------------------------------------------------------------
     B                | magnetic_field       | G                | +1.000000e+00 | [+0.000000e+00,No           ] | False 
     L_Disk           | Disk                 | erg/s            | +1.000000e+46 | [+0.000000e+00,No           ] | False 
     N                | electron_density     | cm^-3            | +6.428643e+03 | [+0.000000e+00,No           ] | False 
     R                | region_size          | cm               | +1.550000e+01 | [+0.000000e+00,+3.000000e+01] | True 
     R_BLR_in         | BLR                  | cm               | +1.000000e+18 | [+0.000000e+00,No           ] | False 
     R_BLR_out        | BLR                  | cm               | +2.000000e+18 | [+0.000000e+00,No           ] | False 
     R_DT             | DT                   | cm               | +5.000000e+18 | [+0.000000e+00,No           ] | False 
     R_H              | Disk                 | cm               | +1.000000e+17 | [+0.000000e+00,No           ] | False 
     R_ext_Sw         | Disk                 | Sw. radii        | +5.000000e+02 | [+0.000000e+00,No           ] | False 
     R_inner_Sw       | Disk                 | Sw. radii        | +3.000000e+00 | [+0.000000e+00,No           ] | False 
     T_DT             | DT                   | K                | +1.000000e+02 | [+0.000000e+00,No           ] | False 
     T_Disk           | Disk                 | K                | +1.000000e+05 | [+0.000000e+00,No           ] | False 
     accr_eff         | Disk                 |                  | +1.000000e-01 | [+0.000000e+00,No           ] | False 
     beam_obj         | beaming              |                  | +2.500000e+01 | [+1.000000e+00,No           ] | False 
     gamma_break      | turn-over-energy     | Lorentz-factor   | +5.000000e+02 | [+1.000000e+00,+1.000000e+08] | False 
     gmax             | high-energy-cut-off  | Lorentz-factor   | +1.000000e+05 | [+1.000000e+00,+1.000000e+15] | False 
     gmin             | low-energy-cut-off   | Lorentz-factor   | +2.000000e+00 | [+1.000000e+00,+1.000000e+05] | False 
     p                | LE_spectral_slope    |                  | +1.500000e+00 | [-1.000000e+01,+1.000000e+01] | False 
     p_1              | HE_spectral_slope    |                  | +3.500000e+00 | [-1.000000e+01,+1.000000e+01] | False 
     tau_BLR          | BLR                  |                  | +1.000000e-01 | [+0.000000e+00,+1.000000e+00] | False 
     tau_DT           | DT                   |                  | +1.000000e-01 | [+0.000000e+00,+1.000000e+00] | False 
     z_cosm           | redshift             |                  | +6.000000e-01 | [+0.000000e+00,No           ] | False 
    -------------------------------------------------------------------------------------------------------------------
    -------------------------------------------------------------------------------------------------------------------


.. code:: ipython3

    my_jet.eval()
    p=my_jet.plot_model()
    p.rescale(y_min=-15,y_max=-10,x_min=8,x_max=27)



.. image:: Jet_example_phys_files/Jet_example_phys_97_0.png


.. code:: ipython3

    my_jet.add_EC_component('EC_DT')
    my_jet.eval()
    p=my_jet.plot_model()
    p.rescale(y_min=-15,y_max=-10,x_min=8,x_max=27)



.. image:: Jet_example_phys_files/Jet_example_phys_98_0.png


