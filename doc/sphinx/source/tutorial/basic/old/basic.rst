basic concepts
=====================


work place
-------------------------
The first step is to import the package, and to set a work place:

.. literalinclude:: ../../../../test/test_tutorial_basic.py
	:lines: 1-6

the latter command will set the output directory to `./` and will add a flag `TEST` to the
output products

data format
-------------------------

The SED data can be stored in ASCII file in a  quite flexible way, but some requirements are needed:

 - you must provide at least two columns for frequencies and fluxes
 - frequencies are in *Hz*
 - fluxes are in  *cgs*, .....

The header of the file can contain some meta-data that are sourced when the data are loaded.
The meta-data available are :

 - z : redhsift
 - resframe: restframe of the data `src` or `obs` 
 - data_scale: scale of the data `lin-lin` or `log-log`
 - dataType: structure of the comumns with the SED data


the meaning of these meta-data is explained in detail in :class:`BlazarSEDFit.data_loader.ObsData` class 
documentation. The meta-data can be included in the header with  line like:

.. code::

	# md meta-data-name meta-data-value

A typical structure of SED data file, including meta-data declaration  is the following:

.. literalinclude:: ../../../../BlazarSEDFit/test_data/SEDs_data/SED_MW_Mrk421.dat
   :lines: 1-20
   :emphasize-lines: 1-5
 





loading data
-------------------------

The most effective way to import the SED data is to create an object 
instance of :class:`BlazarSEDFit.data_loader.ObsData` class 
(see the documentation for the :doc:`data_loader <../../modules_doc/data_loader>` module)
The package provides some test SEDs, accessible as follows:

.. literalinclude:: ../../../../examples/tutorial/test_tutorial_basic.py
	:lines: 1-8
	
producing the following output:
	
.. code-block:: python

	['/Users/orion/.local/lib/python2.6/site-packages/BlazarSEDFit/test_data/SEDs_data/SED_MW_Mrk421.dat',
	 '/Users/orion/.local/lib/python2.6/site-packages/BlazarSEDFit/test_data/SEDs_data/SED_MW_Mrk501.dat']
	


to load the SED of  Mrk 421, the first one in the list:

.. literalinclude:: ../../../../examples/tutorial/test_tutorial_basic.py
	:lines: 9-12
	

.. code:: 

	=============================================================================================

	*** getting meta-data from file header
	=============================================================================================
	
	=============================================================================================
	
	*** loading data ***
	---> loading data for file=/Users/orion/Library/Python/2.7/lib/python/site-packages/BlazarSEDFit/test_data/SEDs_data/SED_MW_Mrk501.dat
	---> found these col ID and names: [2, 1, 0] ['dy', 'y', 'x']
	---> z=3.360000e-02
	---> restframe=obs
	---> obj_name=Mrk-501
	---> data len=120
	---> filtering for dupl entries
	---> remove duplicate entries
	---> data len after filtering dupl entries=120
	---> filtering for UL
	=============================================================================================
	
As you can see the all the meta-data have been properly sourced from the SED file header. You also get information on 
the lenght of the data, before and after elimination of duplicated entries, and upper limits
These meta-data are parameters needed by the 
:class:`BlazarSEDFit.data_loader.ObsData` constructor. 
	
	
	
plotting data
-------------------------

We can now plot our SED using the :class:`BlazarSEDFit.plot_sedfit.Plot` class 
(see the documentation for the :doc:`plot_sedfit <../../modules_doc/plot_sedfit>` module)

.. literalinclude:: ../../../../examples/tutorial/test_tutorial_basic.py
	:lines: 13-16

That will produce

.. image::  ../../../../examples/tutorial/first-Trial/SED_data.png
   :align:   center    
   :width: 85 %
  

grouping data
-----------------------------------
As you can see, due to the overlapping of different instruments and to different time snapshots, 
some points have multiple values. Although this is not a problem for the fit process, you might 
want to rebin your data.
This can be obtained with the following command:

.. literalinclude:: ../../../../examples/tutorial/test_tutorial_basic.py
	:lines: 19-20



handling errors and systematics
-----------------------------------
Another important issues when dealing with fitting of data, is the proper handling of errors.
Typically one might need to add systematics for different reasons:
 
 - data are not really simultaneous, and you want to add systematics to take this into account
 - data (typically IR up to UV), might have very small errors compared to those at higher energies.
   This might bias the minimizer to accomodate the parameters  in order to fit 'better' the low
   frequencies branch.
  
For these reasons the package offer the possibility to add systematics 
  

.. literalinclude:: ../../../../examples/tutorial/test_tutorial_basic.py
	:lines: 22-27

with this command we add 10% systematics for data between :math:`10^{6}<\nu<10^{19}` Hz

.. image::  ../../../../examples/tutorial/first-Trial/SED_data_rebinned.png
   :align:   center    
   :width: 85 %


sed shaping
=================


spectral indices
----------------
.. literalinclude:: ../../../../examples/tutorial/test_tutorial_basic.py
	:lines: 31-41


.. image::  ../../../../examples/tutorial/first-Trial/SED_indices_rebinned.png
   :align:   center    
   :width: 85 %

loglog polyfit
-------------------------
.. literalinclude:: ../../../../examples/tutorial/test_tutorial_basic.py
	:lines: 45-59
	

.. image::  ../../../../examples/tutorial/first-Trial/SED_shaped_rebinned.png
   :align:   center    
   :width: 85 %

sed shaping report
------------------
.. literalinclude:: ../../../../examples/tutorial/test_tutorial_basic.py
	:lines: 60-61

.. code::
	
	=============================================================================================

	*** SEDShape values ***
	---> spectral inidces values
	---> name = radio            range=[6.000 ,10.000] log(Hz)  photon.val=-1.295207e+00, err=1.359751e-01
	---> name = radio_mm         range=[10.000,11.000] log(Hz)  photon.val=-1.455018e+00, err=5.280900e-02
	---> name = mm_IR            range=[10.300,13.700] log(Hz)  photon.val=-1.296277e+00, err=3.749587e-02
	---> name = IR_Opt           range=[12.300,14.700] log(Hz)  photon.val=-2.087455e+00, err=5.433977e-01
	---> name = Opt_UV           range=[14.000,16.000] log(Hz)  photon.val=-2.665891e+00, err=1.419430e-01
	---> name = BBB              range=[14.800,16.200] log(Hz)  photon.val=-2.282189e+00, err=5.738888e-01
	---> name = UV_X             range=[15.000,17.500] log(Hz)  photon.val=-1.873128e+00, err=7.268880e-03
	---> name = X                range=[16.000,19.000] log(Hz)  photon.val=-2.111490e+00, err=3.364662e-02
	---> name = Fermi            range=[22.380,25.380] log(Hz)  photon.val=-1.847682e+00, err=1.485163e-02
	
	---> peak values
	---> sync       nu_p=+1.703050e+01 (err=+9.434834e-02)  nuFnu_p=-1.030007e+01 (err=+1.891358e-02) curv.=-6.398249e-02 (err=+7.839578e-03)
	
	---> IC         nu_p=+2.534625e+01 (err=+9.936733e-02)  nuFnu_p=-1.058254e+01 (err=+2.259612e-02) curv.=-1.167252e-01 (err=+1.187810e-02)
	=============================================================================================
		
		
constraining the model
=======================


SSC constraining 
----------------

.. literalinclude:: ../../../../examples/tutorial/test_tutorial_basic.py
	:lines: 65-73

.. image::  ../../../../examples/tutorial/first-Trial/obs_constr_lppl.png
   :align:   center    
   :width: 85 %
   
the report for the constrained model will reported in table:

.. code::
	
	show pars
	-----------------------------------------------------------------------------------------
	model parameters for jet model:
	
	electron distribution type = lppl
	--------------------------------------------------------------------------------------------------------------
	model parameters:
	 Name             | Type                     | Units            | value         | phys. boundaries
	--------------------------------------------------------------------------------------------------------------
	 z_cosm           | redshift                 |                  | +3.360000e-02 | [+0.000000e+00,No           ]
	 B                | magnetic_field           | G                | +1.000000e-01 | [+0.000000e+00,No           ]
	 R                | region_size              | cm               | +4.245923e+15 | [+0.000000e+00,No           ]
	 beam_obj         | beaming                  |                  | +2.500000e+01 | [+1.000000e+00,No           ]
	 gmax             | high-energy-cut-off      | Lorentz-factor   | +1.642067e+06 | [+1.000000e+00,No           ]
	 gmin             | low-energy-cut-off       | Lorentz-factor   | +7.474642e+01 | [+1.000000e+00,No           ]
	 N                | electron_density         | cm^-3            | +2.166529e+02 | [+0.000000e+00,No           ]
	 s                | LE_spectral_slope        |                  | +2.244597e+00 | [-1.000000e+01,+1.000000e+01]
	 r                | spectral_curvature       |                  | +3.199125e-01 | [-1.000000e+01,+1.000000e+01]
	 gamma0_log_parab | turn-over-energy         | Lorentz-factor   | +7.222931e+03 | [+1.000000e+00,No           ]
	--------------------------------------------------------------------------------------------------------------
	-----------------------------------------------------------------------------------------

to get any time the report for the model:

.. code ::

	jet_model.show_pars()

Least-squares curve fitting
=============================

.. literalinclude:: ../../../../examples/tutorial/test_tutorial_basic.py
	:lines: 77-91

.. image::  ../../../../examples/tutorial/first-Trial/SSC_best_fit_lppl.png
   :align:   center    
   :width: 85 %

the fit information will be reported:

.. code ::
	
	=============================================================================================

	*** start fit process ***
	initial pars:
	--------------------------------------------------------------------------------------------------------------
	model parameters:
	 Name             | Type                     | Units            | value         | phys. boundaries
	--------------------------------------------------------------------------------------------------------------
	 z_cosm           | redshift                 |                  | +3.360000e-02 | [+0.000000e+00,No           ]
	 B                | magnetic_field           | G                | +1.000000e-01 | [+0.000000e+00,No           ]
	 R                | region_size              | cm               | +4.245923e+15 | [+0.000000e+00,No           ]
	 beam_obj         | beaming                  |                  | +2.500000e+01 | [+1.000000e+00,No           ]
	 gmax             | high-energy-cut-off      | Lorentz-factor   | +1.642067e+06 | [+1.000000e+00,No           ]
	 gmin             | low-energy-cut-off       | Lorentz-factor   | +7.474642e+01 | [+1.000000e+00,No           ]
	 N                | electron_density         | cm^-3            | +2.166529e+02 | [+0.000000e+00,No           ]
	 s                | LE_spectral_slope        |                  | +2.244597e+00 | [-1.000000e+01,+1.000000e+01]
	 r                | spectral_curvature       |                  | +3.199125e-01 | [-1.000000e+01,+1.000000e+01]
	 gamma0_log_parab | turn-over-energy         | Lorentz-factor   | +7.222931e+03 | [+1.000000e+00,No           ]
	 L_host           | nuFnu-scale              | erg cm^-2 s^-1   | -1.003477e+01 | [-2.000000e+01,+2.000000e+01]
	--------------------------------------------------------------------------------------------------------------
	**************************************************************************************************
	Fit report
	
	Model:
	--------------------------------------------------------------------------------------------------------------
	model parameters:
	 Name             | Type                     | Units            | value         | phys. boundaries
	--------------------------------------------------------------------------------------------------------------
	 z_cosm           | redshift                 |                  | +3.360000e-02 | [+0.000000e+00,No           ]
	 B                | magnetic_field           | G                | +4.189317e-02 | [+0.000000e+00,No           ]
	 R                | region_size              | cm               | +7.407974e+15 | [+0.000000e+00,No           ]
	 beam_obj         | beaming                  |                  | +2.500000e+01 | [+1.000000e+00,No           ]
	 gmax             | high-energy-cut-off      | Lorentz-factor   | +1.642069e+06 | [+1.000000e+00,No           ]
	 gmin             | low-energy-cut-off       | Lorentz-factor   | +6.073412e+01 | [+1.000000e+00,No           ]
	 N                | electron_density         | cm^-3            | +1.695550e+02 | [+0.000000e+00,No           ]
	 s                | LE_spectral_slope        |                  | +2.252950e+00 | [-1.000000e+01,+1.000000e+01]
	 r                | spectral_curvature       |                  | +2.729864e-01 | [-1.000000e+01,+1.000000e+01]
	 gamma0_log_parab | turn-over-energy         | Lorentz-factor   | +1.291370e+04 | [+1.000000e+00,No           ]
	 L_host           | nuFnu-scale              | erg cm^-2 s^-1   | -1.003477e+01 | [-2.000000e+01,+2.000000e+01]
	--------------------------------------------------------------------------------------------------------------
	
	converged=2
	calls=69
	mesg=The relative error between two consecutive iterates is at most 0.000000
	dof=23
	chisq=54.145014, chisq/red=2.354131 null hypothesis sig=0.000256
	
	best fit pars
	---------------------------------------------------------------------------------------------------
	best-fit parameters:
	  Name            | best-fit value| best-fit err  | start value   | fit boundaries
	---------------------------------------------------------------------------------------------------
	 z_cosm           | Frozen        | Frozen        | +3.360000e-02 | [+0.000000e+00,No           ]
	 B                | +4.189317e-02 | +1.418345e-02 | +1.000000e-01 | [+0.000000e+00,No           ]
	 R                | +7.407974e+15 | +2.811495e+15 | +4.245923e+15 | [+0.000000e+00,No           ]
	 beam_obj         | Frozen        | Frozen        | +2.500000e+01 | [+1.000000e+00,No           ]
	 gmax             | +1.642069e+06 | +8.184449e+05 | +1.642067e+06 | [+1.000000e+00,No           ]
	 gmin             | +6.073412e+01 | +3.026311e+01 | +7.474642e+01 | [+1.000000e+00,No           ]
	 N                | +1.695550e+02 | +1.430634e+02 | +2.166529e+02 | [+0.000000e+00,No           ]
	 s                | +2.252950e+00 | +8.596979e-02 | +2.244597e+00 | [-1.000000e+01,+1.000000e+01]
	 r                | +2.729864e-01 | +1.059576e-01 | +3.199125e-01 | [-1.000000e+01,+1.000000e+01]
	 gamma0_log_parab | +1.291370e+04 | +1.398577e+04 | +7.222931e+03 | [+1.000000e+00,No           ]
	 L_host           | Frozen        | Frozen        | -1.003477e+01 | [-1.225748e+01,-8.257480e+00]
	---------------------------------------------------------------------------------------------------
	1e+11 1e+28
	=============================================================================================

to get any time the  fir report:

.. code ::

	best_fit.show_report()
	
