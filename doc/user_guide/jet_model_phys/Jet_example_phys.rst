
jet\_model user guide
=====================


In this section we describe how to use the module  :mod:`.jet_model`
to buil a model of jet able to reproduce SSC/EC emission processes.
The :mod:`.jet_model`  allows to build a jet  model  providing an interface 
to call the C numerical code that provides an accurate and fast computation of the synchrotron and inverse Compoton processes.  The python wrappper is  built using SWIG framework. 

basic setup
-----------

A jet can be built using the  the :class:`.Jet` class, istanciating a jet object.

.. code:: ipython3

    %matplotlib inline
    from jetset import jet_model
    my_jet=jet_model.Jet(name='test',electron_distribution='lppl',)

This instruction will create:
    * a ``Jet`` object with ``name`` **test**,
    * using as electron distribution the **lppl** model, that is a log-parabola witha low-energy power-law branch.
    * using as wokring directory **test_jet_prod**

For a list of possible distribution refer to :func:`.jet_model.build_electron_distribution_dic`
    
The parameters of the model are accessible throug the instruction

.. code:: ipython3

    my_jet.parameters.show_pars()


.. parsed-literal::

    --------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                     | Units            | value         | phys. boundaries
    --------------------------------------------------------------------------------------------------------------
     N                | electron_density         | cm^-3            | +1.000000e+02 | [+0.000000e+00,No           ]  
     gmin             | low-energy-cut-off       | Lorentz-factor   | +2.000000e+00 | [+1.000000e+00,No           ]  
     gmax             | high-energy-cut-off      | Lorentz-factor   | +1.000000e+08 | [+1.000000e+00,No           ]  
     s                | LE_spectral_slope        |                  | +2.000000e+00 | [-1.000000e+01,+1.000000e+01]  
     r                | spectral_curvature       |                  | +4.000000e-01 | [-1.500000e+01,+1.500000e+01]  
     gamma0_log_parab | turn-over-energy         | Lorentz-factor   | +1.000000e+04 | [+1.000000e+00,No           ]  
     R                | region_size              | cm               | +3.000000e+15 | [+0.000000e+00,No           ]  
     B                | magnetic_field           | G                | +1.000000e-01 | [+0.000000e+00,No           ]  
     beam_obj         | beaming                  |                  | +1.000000e+01 | [+1.000000e+00,No           ]  
     z_cosm           | redshift                 |                  | +1.000000e-01 | [+0.000000e+00,No           ]  
    --------------------------------------------------------------------------------------------------------------


Each parameter has default values. All the parameters listed are handled by :class:`.ModelParameterArray`, and each parameter is an instance of the the :class:`.JetParameter`.

setting the parameters
----------------------

assume you want to change some of the parameters in your model:

.. code:: ipython3

    my_jet.set_par('B',val=0.2)
    my_jet.set_par('gamma0_log_parab',val=5E3)
    my_jet.set_par('gmin',val=1E2)
    my_jet.set_par('R',val=5E14)
    my_jet.set_par('N',val=1E3)

evaluate and plot the model
---------------------------

At this point we can evaluate the emission for this jet model using the
instruction

.. code:: ipython3

    my_jet.eval()

and plot the corresponding SED:

.. code:: ipython3

    from jetset.plot_sedfit import PlotSED
    my_plot=PlotSED()
    my_jet.plot_model(plot_obj=my_plot)
    my_plot.rescale(y_min=-17.5,x_min=8)



.. image:: Jet_example_phys_files/Jet_example_phys_15_0.png


alternatively, you can call the ``plot_model`` method without passing a
``Plot`` object

.. code:: ipython3

    my_plot=my_jet.plot_model()
    my_plot.rescale(y_min=-17.5,x_min=8)



.. image:: Jet_example_phys_files/Jet_example_phys_17_0.png


the ``my_plot`` objet returned will be built on the fly by the
``plot_model`` method

comparing models on the same plot
---------------------------------

to compare the same model after changing a parameter

.. code:: ipython3

    my_plot=PlotSED()
    
    my_jet.plot_model(my_plot,label='gamma0_log_parab=1E4',comp='Sum')
    my_jet.set_par('gamma0_log_parab',val=1.0E5)
    my_jet.eval()
    my_jet.plot_model(my_plot,label='gamma0_log_parab=1E5',comp='Sum')
    my_plot.rescale(y_min=-17,x_min=8)



.. image:: Jet_example_phys_files/Jet_example_phys_21_0.png


saving a plot
-------------

to save the plot

.. code:: ipython3

    my_plot.save('jet1.png')

saving and loadaing model
-------------------------

.. code:: ipython3

    my_jet.save_model('test_model.dat')

.. code:: ipython3

    my_jet_new=jet_model.Jet.load_model('test_model.dat')


.. parsed-literal::

    Sum
    Sync
    SSC
    -----------------------------------------------------------------------------------------
    model parameters for jet model:
    electron distribution type = lppl  
    --------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                     | Units            | value         | phys. boundaries
    --------------------------------------------------------------------------------------------------------------
     R                | region_size              | cm               | +3.000000e+15 | [+0.000000e+00,No           ]  
     B                | magnetic_field           | G                | +1.000000e-01 | [+0.000000e+00,No           ]  
     beam_obj         | beaming                  |                  | +1.000000e+01 | [+1.000000e+00,No           ]  
     z_cosm           | redshift                 |                  | +1.000000e-01 | [+0.000000e+00,No           ]  
     N                | electron_density         | cm^-3            | +1.000000e+02 | [+0.000000e+00,No           ]  
     gmin             | low-energy-cut-off       | Lorentz-factor   | +2.000000e+00 | [+1.000000e+00,No           ]  
     gmax             | high-energy-cut-off      | Lorentz-factor   | +1.000000e+08 | [+1.000000e+00,No           ]  
     s                | LE_spectral_slope        |                  | +2.000000e+00 | [-1.000000e+01,+1.000000e+01]  
     r                | spectral_curvature       |                  | +4.000000e-01 | [-1.500000e+01,+1.500000e+01]  
     gamma0_log_parab | turn-over-energy         | Lorentz-factor   | +1.000000e+04 | [+1.000000e+00,No           ]  
    --------------------------------------------------------------------------------------------------------------
    -----------------------------------------------------------------------------------------


set N from observed Flux
------------------------

It is possible to set the density of emitting particle starting from some observed Luminosity (see the method     :meth:`.Jet.set_N_from_nuFnu`)

.. code:: ipython3

    my_jet=jet_model.Jet(name='test',electron_distribution='lppl')

this is the initial value of N

.. code:: ipython3

    my_jet.get_par_by_name('N').val




.. parsed-literal::

    100.0



.. code:: ipython3

    
    my_jet.set_N_from_nuFnu(nuFnu_obs=1E-14,nu_obs=1E12)

This is the updated value of N, obtained in order to match the given
flux at the given frequencies

.. code:: ipython3

    my_jet.get_par_by_name('N').val




.. parsed-literal::

    30781.15088279897



.. code:: ipython3

    my_jet.parameters.show_pars()


.. parsed-literal::

    --------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                     | Units            | value         | phys. boundaries
    --------------------------------------------------------------------------------------------------------------
     N                | electron_density         | cm^-3            | +3.078115e+04 | [+0.000000e+00,No           ]  
     gmin             | low-energy-cut-off       | Lorentz-factor   | +2.000000e+00 | [+1.000000e+00,No           ]  
     gmax             | high-energy-cut-off      | Lorentz-factor   | +1.000000e+08 | [+1.000000e+00,No           ]  
     s                | LE_spectral_slope        |                  | +2.000000e+00 | [-1.000000e+01,+1.000000e+01]  
     r                | spectral_curvature       |                  | +4.000000e-01 | [-1.500000e+01,+1.500000e+01]  
     gamma0_log_parab | turn-over-energy         | Lorentz-factor   | +1.000000e+04 | [+1.000000e+00,No           ]  
     R                | region_size              | cm               | +3.000000e+15 | [+0.000000e+00,No           ]  
     B                | magnetic_field           | G                | +1.000000e-01 | [+0.000000e+00,No           ]  
     beam_obj         | beaming                  |                  | +1.000000e+01 | [+1.000000e+00,No           ]  
     z_cosm           | redshift                 |                  | +1.000000e-01 | [+0.000000e+00,No           ]  
    --------------------------------------------------------------------------------------------------------------


.. code:: ipython3

    my_jet.eval()
    my_plot=my_jet.plot_model(label='set N from F=1E-14')
    my_plot.rescale(y_min=-17,x_min=8)



.. image:: Jet_example_phys_files/Jet_example_phys_37_0.png


setting beaming factor
----------------------

It is possible to set the bemaing factor according to the realativistic
BulkFactor and viewing angle, this can be done by setting the
``beaming_expr`` kw in the Jet constructor, possbile choiches are

-  ``delta`` to provide directly the beaming factor (default)
-  ``bulk_theta`` to provide the BulkFactor and the jet viewing angle

.. code:: ipython3

    my_jet=jet_model.Jet(name='test',electron_distribution='lppl',beaming_expr='bulk_theta')

.. code:: ipython3

    my_jet.parameters.show_pars()


.. parsed-literal::

    --------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                     | Units            | value         | phys. boundaries
    --------------------------------------------------------------------------------------------------------------
     N                | electron_density         | cm^-3            | +1.000000e+02 | [+0.000000e+00,No           ]  
     gmin             | low-energy-cut-off       | Lorentz-factor   | +2.000000e+00 | [+1.000000e+00,No           ]  
     gmax             | high-energy-cut-off      | Lorentz-factor   | +1.000000e+08 | [+1.000000e+00,No           ]  
     s                | LE_spectral_slope        |                  | +2.000000e+00 | [-1.000000e+01,+1.000000e+01]  
     r                | spectral_curvature       |                  | +4.000000e-01 | [-1.500000e+01,+1.500000e+01]  
     gamma0_log_parab | turn-over-energy         | Lorentz-factor   | +1.000000e+04 | [+1.000000e+00,No           ]  
     R                | region_size              | cm               | +3.000000e+15 | [+0.000000e+00,No           ]  
     B                | magnetic_field           | G                | +1.000000e-01 | [+0.000000e+00,No           ]  
     theta            | jet-viewing-angle        | deg              | +1.000000e-01 | [+0.000000e+00,No           ]  
     BulkFactor       | jet-bulk-factor          | Lorentz-factor   | +1.000000e+01 | [+1.000000e+00,No           ]  
     z_cosm           | redshift                 |                  | +1.000000e-01 | [+0.000000e+00,No           ]  
    --------------------------------------------------------------------------------------------------------------


the actual value of the beaming factor che be obatained using the :meth:`.Jet.get_beaming`

.. code:: ipython3

    print ('beaming factor=',my_jet.get_beaming())


.. parsed-literal::

    beaming factor= 19.943844732554165


We can change the value of ``theta`` and get the updated value of the beaming factor

.. code:: ipython3

    my_jet.set_par('theta',val=10.)

.. code:: ipython3

    print('beaming factor=',my_jet.get_beaming())


.. parsed-literal::

    beaming factor= 4.968041140891955


of course setting `beaming_expr=delta` we get the same beaming expression as in the default case

.. code:: ipython3

    my_jet=jet_model.Jet(name='test',electron_distribution='lppl',beaming_expr='delta')

.. code:: ipython3

    my_jet.parameters.show_pars()


.. parsed-literal::

    --------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                     | Units            | value         | phys. boundaries
    --------------------------------------------------------------------------------------------------------------
     N                | electron_density         | cm^-3            | +1.000000e+02 | [+0.000000e+00,No           ]  
     gmin             | low-energy-cut-off       | Lorentz-factor   | +2.000000e+00 | [+1.000000e+00,No           ]  
     gmax             | high-energy-cut-off      | Lorentz-factor   | +1.000000e+08 | [+1.000000e+00,No           ]  
     s                | LE_spectral_slope        |                  | +2.000000e+00 | [-1.000000e+01,+1.000000e+01]  
     r                | spectral_curvature       |                  | +4.000000e-01 | [-1.500000e+01,+1.500000e+01]  
     gamma0_log_parab | turn-over-energy         | Lorentz-factor   | +1.000000e+04 | [+1.000000e+00,No           ]  
     R                | region_size              | cm               | +3.000000e+15 | [+0.000000e+00,No           ]  
     B                | magnetic_field           | G                | +1.000000e-01 | [+0.000000e+00,No           ]  
     beam_obj         | beaming                  |                  | +1.000000e+01 | [+1.000000e+00,No           ]  
     z_cosm           | redshift                 |                  | +1.000000e-01 | [+0.000000e+00,No           ]  
    --------------------------------------------------------------------------------------------------------------


accessing individual components
-------------------------------

It is possible to access specific spectral components of oura model

.. code:: ipython3

    my_jet=jet_model.Jet(name='test',electron_distribution='lppl',beaming_expr='bulk_theta')
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

and from the ``SED`` object we can extract both the nu and nuFnu array

.. code:: ipython3

    nu_sync=Sync.SED.nu
    nuFnu_sync=Sync.SED.nuFnu

.. code:: ipython3

    print (nuFnu_sync)


.. parsed-literal::

    [1.00000000e-30 1.00000000e-30 1.00000000e-30 1.00000000e-30
     1.00000000e-30 1.00000000e-30 1.00000000e-30 1.00000000e-30
     1.00000000e-30 1.00000000e-30 1.00000000e-30 1.99195355e-26
     4.32084831e-26 9.42080231e-26 2.07796467e-25 4.64226006e-25
     1.07421835e-24 2.57262866e-24 6.33304102e-24 1.62172378e-23
     4.21830732e-23 1.10565753e-22 2.90749001e-22 7.64703810e-22
     2.01137205e-21 5.29044499e-21 1.39152631e-20 3.66005742e-20
     9.62674977e-20 2.53161731e-19 6.65135791e-19 1.74532879e-18
     4.47162406e-18 1.04798662e-17 2.22510780e-17 3.42183006e-17
     4.73795507e-17 5.98146121e-17 7.08393890e-17 8.27924884e-17
     9.58571035e-17 1.10361592e-16 1.26894305e-16 1.45793212e-16
     1.67429058e-16 1.92253041e-16 2.20741484e-16 2.53430172e-16
     2.90956596e-16 3.34034909e-16 3.83482107e-16 4.40231844e-16
     5.05364252e-16 5.80111171e-16 6.65860835e-16 7.64238867e-16
     8.77094393e-16 1.00644889e-15 1.15473841e-15 1.32470643e-15
     1.51922268e-15 1.74184487e-15 1.99657203e-15 2.28728218e-15
     2.61887108e-15 2.99699329e-15 3.42637539e-15 3.91284208e-15
     4.46398022e-15 5.08436971e-15 5.77822165e-15 6.55504712e-15
     7.41688282e-15 8.35893256e-15 9.39317225e-15 1.05162021e-14
     1.17012783e-14 1.29677483e-14 1.43077405e-14 1.56618661e-14
     1.70669190e-14 1.85140020e-14 1.99149955e-14 2.13199797e-14
     2.27192694e-14 2.40181924e-14 2.52556860e-14 2.64340055e-14
     2.74611617e-14 2.83590869e-14 2.91498397e-14 2.97544176e-14
     3.01736120e-14 3.04551553e-14 3.05414488e-14 3.04103209e-14
     3.01368848e-14 2.96892404e-14 2.90230677e-14 2.82372882e-14
     2.73248535e-14 2.62226222e-14 2.50449532e-14 2.38042485e-14
     2.24239213e-14 2.10225626e-14 1.96144574e-14 1.81447517e-14
     1.66964251e-14 1.52898970e-14 1.38900619e-14 1.25442884e-14
     1.12742894e-14 1.00574503e-14 8.91399490e-15 7.86233347e-15
     6.88692857e-15 5.99001626e-15 5.18464148e-15 4.45908613e-15
     3.80577564e-15 3.23238294e-15 2.72949108e-15 2.28586629e-15
     1.90500228e-15 1.57929305e-15 1.29769920e-15 1.06107838e-15
     8.63339561e-16 6.96082614e-16 5.58275014e-16 4.45484724e-16
     3.52308323e-16 2.76800930e-16 2.16201897e-16 1.67118893e-16
     1.27789365e-16 9.67429184e-17 7.19114364e-17 5.19813320e-17
     3.65724953e-17 2.45791387e-17 1.51568234e-17 8.62852876e-18
     4.37420009e-18 1.73809120e-18 5.58993124e-19 1.38471677e-19
     1.84044040e-20 1.46229352e-21 6.75656391e-23 6.84161407e-25
     2.05558538e-27 1.70711836e-30 1.00000000e-30 1.00000000e-30
     1.00000000e-30 1.00000000e-30 1.00000000e-30 1.00000000e-30
     1.00000000e-30 1.00000000e-30 1.00000000e-30 1.00000000e-30
     1.00000000e-30 1.00000000e-30 1.00000000e-30 1.00000000e-30
     1.00000000e-30 1.00000000e-30 1.00000000e-30 1.00000000e-30
     1.00000000e-30 1.00000000e-30 1.00000000e-30 1.00000000e-30
     1.00000000e-30 1.00000000e-30 1.00000000e-30 1.00000000e-30
     1.00000000e-30 1.00000000e-30 1.00000000e-30 1.00000000e-30
     1.00000000e-30 1.00000000e-30 1.00000000e-30 1.00000000e-30
     1.00000000e-30 1.00000000e-30 1.00000000e-30 1.00000000e-30
     1.00000000e-30 1.00000000e-30 1.00000000e-30 1.00000000e-30
     1.00000000e-30 1.00000000e-30 1.00000000e-30 1.00000000e-30]

