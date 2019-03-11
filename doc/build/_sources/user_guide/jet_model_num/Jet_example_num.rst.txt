
.. _jet_numerical_guide:

Numerical setup
===============

Changing the grid size for the electron distribution
----------------------------------------------------

.. code:: ipython3

    from jetset.jet_model import Jet
    my_jet=Jet(name='test',electron_distribution='lppl',)
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


It is possible to change the size of the grid for the electron
distribution. It is worth noting that at lower values of the grid size
the speed will increase, **but it is not recommended to go below 100**.

The actual value of the grid size is returned by the :meth:`.Jet.gamma_grid_size`

.. code:: ipython3

    print (my_jet.gamma_grid_size)

and this value can be changed using the method :meth:`.Jet.set_gamma_grid_size`. In the following we show the result for a grid of size=10, as anticipated the final integration will be not satisfactory

.. code:: ipython3

    my_jet.set_gamma_grid_size(10)
    my_jet.eval()
    sed_plot=my_jet.plot_model()
    sed_plot.rescale(x_min=8,y_min=-20,y_max=-12)

.. code:: ipython3

    my_jet.set_gamma_grid_size(100)
    my_jet.eval()
    sed_plot=my_jet.plot_model()
    sed_plot.rescale(x_min=8,y_min=-20,y_max=-12)

.. code:: ipython3

    my_jet.set_gamma_grid_size(1000)
    my_jet.eval()
    sed_plot=my_jet.plot_model()
    sed_plot.rescale(x_min=8,y_min=-20,y_max=-12)

.. code:: ipython3

    my_jet.set_gamma_grid_size(10000)
    my_jet.eval()
    sed_plot=my_jet.plot_model()
    sed_plot.rescale(x_min=8,y_min=-20,y_max=-12)

Changing the grid size for the seed photons
-------------------------------------------

.. code:: ipython3

    my_jet=Jet(name='test',electron_distribution='lppl',)
    my_jet.show_model()

we can get the current value of the seed photons grid size using attribute :meth:`.Jet.nu_seed_size`

**in the current version there is lit of the size to 1000**

.. code:: ipython3

    print (my_jet.nu_seed_size)

and this value can be changed using the method :meth:`.Jet.set_seed_nu_size`. In the following we show the result for a grid of nu_size=10

.. code:: ipython3

    my_jet.nu_seed_size=10
    my_jet.eval()
    sed_plot=my_jet.plot_model()
    sed_plot.rescale(x_min=8,y_min=-20,y_max=-12)

Changing the grid size for the IC process spectra
-------------------------------------------------

**in the current version there is lit of the size to 1000**

.. code:: ipython3

    my_jet=Jet(name='test',electron_distribution='lppl',)
    my_jet.show_model()

.. code:: ipython3

    print(my_jet.IC_nu_size)

.. code:: ipython3

    my_jet.IC_nu_size=10
    my_jet.eval()
    sed_plot=my_jet.plot_model()
    sed_plot.rescale(x_min=8,y_min=-20,y_max=-12)

.. code:: ipython3

    my_jet.IC_nu_size=100
    my_jet.eval()
    sed_plot=my_jet.plot_model()
    sed_plot.rescale(x_min=8,y_min=-20,y_max=-12)


