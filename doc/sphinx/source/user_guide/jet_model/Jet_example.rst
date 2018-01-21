
jet\_model module user guide
============================

.. code:: ipython2

    import jetset
    from jetset import jet_model
    myJet=jet_model.Jet(name='test',electron_distribution='lppl')



In this section we describe how to use the module  :mod:`.jet_model`
to buil a model of jet able to reproduce SSC/EC emission processes. 


The :mod:`.jet_model`  allows to build a jet  model  providing an interface 
to call the BlazarSED code. The BlazarSED code is a numerical 
accurate C code, to evaluate SSC/EC emission processes in a relativistic jet. 
The python wrappper is  built using SWIG. 

A jet can be built using the  the :class:`.Jet` class, istanciating a jet object.


import jetset from jetset import jet\_model
myJet=jet\_model.Jet(name='test',electron\_distribution='lppl')

This instruction will create:
    * a ``Jet`` object with ``name`` **test**,
    * using as electron distribution the **lppl** model, that is a log-parabola witha low-energy power-law branch.
    * using as wokring directory **test_dir**
The parameters of the model arre accessible throug the instruction

.. code:: ipython2

    myJet.parameters.show_pars()


.. parsed-literal::

    --------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                     | Units            | value         | phys. boundaries
    --------------------------------------------------------------------------------------------------------------
     gmax             | high-energy-cut-off      | Lorentz-factor   | +1.000000e+08 | [+1.000000e+00,No           ]  
     gmin             | low-energy-cut-off       | Lorentz-factor   | +2.000000e+00 | [+1.000000e+00,No           ]  
     N                | electron_density         | cm^-3            | +1.000000e+02 | [+0.000000e+00,No           ]  
     s                | LE_spectral_slope        |                  | +2.000000e+00 | [-1.000000e+01,+1.000000e+01]  
     r                | spectral_curvature       |                  | +4.000000e-01 | [-1.000000e+01,+1.000000e+01]  
     gamma0_log_parab | turn-over-energy         | Lorentz-factor   | +1.000000e+04 | [+1.000000e+00,No           ]  
     z_cosm           | redshift                 |                  | +1.000000e-01 | [+0.000000e+00,No           ]  
     B                | magnetic_field           | G                | +1.000000e-01 | [+0.000000e+00,No           ]  
     R                | region_size              | cm               | +3.000000e+15 | [+0.000000e+00,No           ]  
     beam_obj         | beaming                  |                  | +1.000000e+01 | [+1.000000e+00,No           ]  
    --------------------------------------------------------------------------------------------------------------


Each parameter has default values. All the parameters listed are handled by
:class:`.ModelParameterArray`, and each parameter is an instance of the the
:class:`.JetParameter`.

At this point one can evaluate the SSC/EC emission for this jet model
using the instruction

.. code:: ipython2

    myJet.eval()


.. parsed-literal::

    ('fill name', 'Sum')
    ('fill name', 'Sync')
    ('fill name', 'SSC')


and plot the corresponding SED:

.. code:: ipython2

    from jetset.plot_sedfit import Plot
    
    myPlot=Plot()
    
    myPlot.add_model_plot(myJet,autoscale=True)
    myPlot.save('jet.png')


.. parsed-literal::

    running PyLab in interactive mode



.. image:: Jet_example_files/Jet_example_11_1.png


To change one of the parameter in the model:

.. code:: ipython2

    myPlot=Plot()
    
    myPlot.add_model_plot(myJet,autoscale=True)
    myJet.set_par('gamma0_log_parab',val=1.0E5)
    
    myJet.eval()
    
    myPlot.add_model_plot(myJet,label='gamma0_log_parab=1E5',autoscale=True)
    
    myPlot.save('jet1.png')


.. parsed-literal::

    running PyLab in interactive mode
    ('fill name', 'Sum')
    ('fill name', 'Sync')
    ('fill name', 'SSC')



.. image:: Jet_example_files/Jet_example_13_1.png


To plot all the components

.. code:: ipython2

    myPlot=Plot()
    for c in myJet.spectral_components: myPlot.add_model_plot(c.SED,autoscale=True)


.. parsed-literal::

    running PyLab in interactive mode



.. image:: Jet_example_files/Jet_example_15_1.png

