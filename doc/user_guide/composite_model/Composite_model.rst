.. _composite_models:


Composite Models
================

.. code:: ipython3

    from jetset.jet_model import Jet
    from jetset.plot_sedfit import PlotSED
    from jetset.model_manager import FitModel


Composite models allow to combine together different models, such as Jet, and templates, including additive or multiplicative models, and give to the user the possibility to define the functional form of the model composition using a very simple and intuitive form such as:

`'jet1+jet2'*Franceschini_2008`

that sums two jet models SEDs, and apply to both of them the `Franceschini_2008` EBL absorption.

Building composite models it is very easy. Composite models are handled by   the :class:`.FitModel` class, as shown by the following examples. 

Combine a Jet model with the EBL model (Multiplicative case)
------------------------------------------------------------

We start by combining a Jet model with the EBL absorption model, i.e. a
multiplicative model. First, we define our Jet model

.. code:: ipython3

    from jetset.jet_model import Jet
    my_jet=Jet(electron_distribution='lppl',name='jet_flaring')

Second, we define the EBL model, and we use in this case the `Franceschini_2008` model ( read the section :ref:`ebl_model`  for more info regarding the EBL models)

.. code:: ipython3

    from jetset.template_2Dmodel import EBLAbsorptionTemplate
    ebl_franceschini=EBLAbsorptionTemplate.from_name('Franceschini_2008')

Now we add the components models to the the :class:`.FitModel` class, using the :class:`.FitModel.add_component()` method 

.. code:: ipython3

    composite_model=FitModel(nu_size=500,name='EBL corrected')
    composite_model.add_component(my_jet)
    composite_model.add_component(ebl_franceschini)



.. parsed-literal::

    /Users/orion/anaconda3/envs/jetset/lib/python3.7/site-packages/jetset-1.1.2-py3.7-macosx-10.9-x86_64.egg/jetset/model_manager.py:160: UserWarning: no cosmology defined, using default FlatLambdaCDM(name="Planck13", H0=67.8 km / (Mpc s), Om0=0.307, Tcmb0=2.725 K, Neff=3.05, m_nu=[0.   0.   0.06] eV, Ob0=0.0483)
      warnings.warn('no cosmology defined, using default %s'%self.cosmo)


the waring message is just telling that you are not passing any specific
cosmology model to the ``FitModel`` class, so it is using a default one

.. code:: ipython3

    composite_model.show_pars()


.. parsed-literal::

        model name          name             par type           units          val      phys. bound. min phys. bound. max  log  frozen
    ----------------- ---------------- ------------------- --------------- ------------ ---------------- ---------------- ----- ------
          jet_flaring             gmin  low-energy-cut-off lorentz-factor* 2.000000e+00     1.000000e+00     1.000000e+09 False  False
          jet_flaring             gmax high-energy-cut-off lorentz-factor* 1.000000e+06     1.000000e+00     1.000000e+15 False  False
          jet_flaring                N    emitters_density         1 / cm3 1.000000e+02     0.000000e+00               -- False  False
          jet_flaring                s   LE_spectral_slope                 2.000000e+00    -1.000000e+01     1.000000e+01 False  False
          jet_flaring                r  spectral_curvature                 4.000000e-01    -1.500000e+01     1.500000e+01 False  False
          jet_flaring gamma0_log_parab    turn-over-energy lorentz-factor* 1.000000e+04     1.000000e+00     1.000000e+09 False  False
          jet_flaring                R         region_size              cm 5.000000e+15     1.000000e+03     1.000000e+30 False  False
          jet_flaring              R_H     region_position              cm 1.000000e+17     0.000000e+00               -- False   True
          jet_flaring                B      magnetic_field               G 1.000000e-01     0.000000e+00               -- False  False
          jet_flaring         beam_obj             beaming Lorentz-factor* 1.000000e+01     1.000000e-04               -- False  False
          jet_flaring           z_cosm            redshift                 1.000000e-01     0.000000e+00               -- False  False
    Franceschini_2008           z_cosm            redshift                 1.000000e+00     0.000000e+00               -- False   True


Since, both the Jet model the EBL share the same parameter, i.e. the
redshift, we link the two parameters

.. code:: ipython3

    composite_model.link_par(par_name='z_cosm',model_name_list=['jet_flaring'],root_model_name='Franceschini_2008')

.. code:: ipython3

    composite_model.show_pars()


.. parsed-literal::

        model name                name                  par type           units          val      phys. bound. min phys. bound. max  log  frozen
    ----------------- --------------------------- ------------------- --------------- ------------ ---------------- ---------------- ----- ------
          jet_flaring                        gmin  low-energy-cut-off lorentz-factor* 2.000000e+00     1.000000e+00     1.000000e+09 False  False
          jet_flaring                        gmax high-energy-cut-off lorentz-factor* 1.000000e+06     1.000000e+00     1.000000e+15 False  False
          jet_flaring                           N    emitters_density         1 / cm3 1.000000e+02     0.000000e+00               -- False  False
          jet_flaring                           s   LE_spectral_slope                 2.000000e+00    -1.000000e+01     1.000000e+01 False  False
          jet_flaring                           r  spectral_curvature                 4.000000e-01    -1.500000e+01     1.500000e+01 False  False
          jet_flaring            gamma0_log_parab    turn-over-energy lorentz-factor* 1.000000e+04     1.000000e+00     1.000000e+09 False  False
          jet_flaring                           R         region_size              cm 5.000000e+15     1.000000e+03     1.000000e+30 False  False
          jet_flaring                         R_H     region_position              cm 1.000000e+17     0.000000e+00               -- False   True
          jet_flaring                           B      magnetic_field               G 1.000000e-01     0.000000e+00               -- False  False
          jet_flaring                    beam_obj             beaming Lorentz-factor* 1.000000e+01     1.000000e-04               -- False  False
          jet_flaring z_cosm(L,Franceschini_2008)            redshift                           --               --               -- False   True
    Franceschini_2008                   z_cosm(R)            redshift                 1.000000e+00     0.000000e+00               -- False   True


As you can see, now the paramter ``z_cosm`` in ``Franceschini_2008`` is
the root parameter (flagged by the R in parenthesis), and the one
belonging to the ``jet_flaring`` component is the linked one (flagged by
the L in parenthesis).

These methods are alternative and equivalent ways to set a parameter in
a composite model

.. code:: ipython3

    composite_model.jet_flaring.parameters.z_cosm.val=0.1
    composite_model.set_par('jet_flaring','z_cosm',0.1)
    composite_model.set_par(my_jet,'z_cosm',0.1)

And now, we can define the functional form of the model composition,
just by writing the mathematical expression as a string, using the model
names reported in the model description table, and that’s it!

.. code:: ipython3

    composite_model.show_model_components()


.. parsed-literal::

    
    -------------------------------------------------------------------------------------------------------------------
    Composite model description
    -------------------------------------------------------------------------------------------------------------------
    name: EBL corrected  
    type: composite_model  
    components models:
     -model name: jet_flaring model type: jet
     -model name: Franceschini_2008 model type: table2D
    
    -------------------------------------------------------------------------------------------------------------------


.. code:: ipython3

    composite_model.composite_expr='jet_flaring*Franceschini_2008'

.. code:: ipython3

    composite_model.jet_flaring.IC_nu_size=150
    composite_model.eval()
    p=composite_model.plot_model()
    p.rescale(y_max=-12)



.. image:: Composite_model_files/Composite_model_22_0.png


Sum two jets (steady and flaring) and apply the EBL to both (Multiplicative and additive)
-----------------------------------------------------------------------------------------

Assume that now we want to sum to jet models (a steady and flaring
component) and apply to both of them the EBL absorption.

.. code:: ipython3

    composite_model=FitModel(nu_size=500,name='EBL corrected flaring+steady')
    composite_model.add_component(my_jet)
    composite_model.add_component(ebl_franceschini)


.. parsed-literal::

    /Users/orion/anaconda3/envs/jetset/lib/python3.7/site-packages/jetset-1.1.2-py3.7-macosx-10.9-x86_64.egg/jetset/model_manager.py:160: UserWarning: no cosmology defined, using default FlatLambdaCDM(name="Planck13", H0=67.8 km / (Mpc s), Om0=0.307, Tcmb0=2.725 K, Neff=3.05, m_nu=[0.   0.   0.06] eV, Ob0=0.0483)
      warnings.warn('no cosmology defined, using default %s'%self.cosmo)


.. code:: ipython3

    steady_jet=Jet(electron_distribution='plc',name='steady_jet')
    composite_model.add_component(steady_jet)
    composite_model.show_model_components()


.. parsed-literal::

    
    -------------------------------------------------------------------------------------------------------------------
    Composite model description
    -------------------------------------------------------------------------------------------------------------------
    name: EBL corrected flaring+steady  
    type: composite_model  
    components models:
     -model name: jet_flaring model type: jet
     -model name: Franceschini_2008 model type: table2D
     -model name: steady_jet model type: jet
    
    -------------------------------------------------------------------------------------------------------------------


.. code:: ipython3

    composite_model.link_par(par_name='z_cosm',model_name_list=['steady_jet'],root_model_name='Franceschini_2008') 

.. code:: ipython3

    composite_model.show_pars()


.. parsed-literal::

        model name                name                  par type           units          val      phys. bound. min phys. bound. max  log  frozen
    ----------------- --------------------------- ------------------- --------------- ------------ ---------------- ---------------- ----- ------
          jet_flaring                        gmin  low-energy-cut-off lorentz-factor* 2.000000e+00     1.000000e+00     1.000000e+09 False  False
          jet_flaring                        gmax high-energy-cut-off lorentz-factor* 1.000000e+06     1.000000e+00     1.000000e+15 False  False
          jet_flaring                           N    emitters_density         1 / cm3 1.000000e+02     0.000000e+00               -- False  False
          jet_flaring                           s   LE_spectral_slope                 2.000000e+00    -1.000000e+01     1.000000e+01 False  False
          jet_flaring                           r  spectral_curvature                 4.000000e-01    -1.500000e+01     1.500000e+01 False  False
          jet_flaring            gamma0_log_parab    turn-over-energy lorentz-factor* 1.000000e+04     1.000000e+00     1.000000e+09 False  False
          jet_flaring                           R         region_size              cm 5.000000e+15     1.000000e+03     1.000000e+30 False  False
          jet_flaring                         R_H     region_position              cm 1.000000e+17     0.000000e+00               -- False   True
          jet_flaring                           B      magnetic_field               G 1.000000e-01     0.000000e+00               -- False  False
          jet_flaring                    beam_obj             beaming Lorentz-factor* 1.000000e+01     1.000000e-04               -- False  False
          jet_flaring z_cosm(L,Franceschini_2008)            redshift                           --               --               -- False   True
    Franceschini_2008                   z_cosm(R)            redshift                 1.000000e-01     0.000000e+00               -- False   True
           steady_jet                        gmin  low-energy-cut-off lorentz-factor* 2.000000e+00     1.000000e+00     1.000000e+09 False  False
           steady_jet                        gmax high-energy-cut-off lorentz-factor* 1.000000e+06     1.000000e+00     1.000000e+15 False  False
           steady_jet                           N    emitters_density         1 / cm3 1.000000e+02     0.000000e+00               -- False  False
           steady_jet                           p   LE_spectral_slope                 2.000000e+00    -1.000000e+01     1.000000e+01 False  False
           steady_jet                   gamma_cut    turn-over-energy lorentz-factor* 1.000000e+04     1.000000e+00     1.000000e+09 False  False
           steady_jet                           R         region_size              cm 5.000000e+15     1.000000e+03     1.000000e+30 False  False
           steady_jet                         R_H     region_position              cm 1.000000e+17     0.000000e+00               -- False   True
           steady_jet                           B      magnetic_field               G 1.000000e-01     0.000000e+00               -- False  False
           steady_jet                    beam_obj             beaming Lorentz-factor* 1.000000e+01     1.000000e-04               -- False  False
           steady_jet z_cosm(L,Franceschini_2008)            redshift                           --               --               -- False   True


.. code:: ipython3

    composite_model.steady_jet.IC_nu_size=150


.. code:: ipython3

    composite_model.composite_expr=composite_model.composite_expr='(jet_flaring + steady_jet) * Franceschini_2008'

.. code:: ipython3

    composite_model.eval()
    p=composite_model.plot_model()
    p.rescale(y_max=-12)



.. image:: Composite_model_files/Composite_model_31_0.png

