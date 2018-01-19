Interactive session
=========================

.. content::

	
initial setup
-------------------------
.. code::

	SED_file='M87_SED.txt'
	workplace=SEDFit.set_workplace(out_dir='./',flag='TEST')

loading data:
-------------------------
The most effective way to import the SED data is to make an object 
instance of :class:`BlazarSEDFit.data_loader.ObsData` class. 
In the interactive session this object can be created with a call like this: 


.. code::

	mySEDdata=SEDFit.ObsData(dataFile=SED_file)  

in this case the header of the file stores all the information
to correctly initialize the ObsData instance:

.. code ::

	# z  0.0308
	# restframe  obs
	# dataScale  lin-lin
	# dataType x,y,dy
	# obj_name     J1104+3812,Mrk 421
	# Frequency [Hz]  EnergyFlux [erg/cm2/s]
	#  Xval           Yval       YvalError
	2.299540e+09 1.340900e-14 3.910000e-16
	2.639697e+09 1.793088e-14 3.231099e-26
	4.799040e+09 2.313600e-14 2.400000e-16
	4.805039e+09 1.773414e-14 1.773414e-15
	......
	
the first header line sets the source redshift. The second the rest frame 
of the data, the third is  

thi	
	
	mySEDdata.add_systematic(0.01)
	mySEDdata.bin_data(bin_width=0.2)
	    
	#SEDdata.add_systematic(0.15)
	mySEDdata.set_error(0.1)
	
	#FermiPlanck
	#SEDdata=SEDFit.ObsData(dataType='x,y,dx,dy,DataSet',dataFile=SED_file,dataSetFilter=[-1,0,1,2])
	
	
	myPlot=SEDFit.Plot(mySEDdata,interactive=True,workplace=workplace,engine='PyLab')
	myPlot.add_data_plot(mySEDdata)
	print mySEDdata.obj_name
	#raw_input("press enter...")
	
	SEDShape=SEDFit.SEDShape(mySEDdata)
	
	SEDShape.eval_indices()
	    
	SEDShape.sync_fit()
	
	SEDShape.IC_fit()
	
	myPlot.add_model_plot(SEDShape.sync)
	myPlot.add_model_plot(SEDShape.IC)
	#raw_input("press enter...")
	
	SEDShape.show_values()
	
	SED_obspar=SEDFit.SED_obspar(1.5,3.5,0.01,'bkn',86400*30,1.5E8,SEDShape=SEDShape,workplace=workplace) 
	
	SEDModel=SED_obspar.constrain_SSC_model()
	myPlot.add_model_plot(SEDModel.SED)
	
	#raw_input("press enter...")
	#fit_par_array.show_pars()
	SEDModel.fit_pars.set('z_cosm','freeze')
	#fit_par_array.get_par('p').freeze()
	#fit_par_array.get_par('B').freeze()
	#
	#   
	
