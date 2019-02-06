#!/usr/bin/env python

import BlazarSEDFit as SEDFit

SEDFit.set_workplace(out_dir='first-Trial',flag='Mrk501')

SEDFit.test_SEDs

mySED=SEDFit.test_SEDs[1]

mySEDdata=SEDFit.ObsData(data_file=mySED)

myPlot=SEDFit.Plot(mySEDdata,interactive=True)

myPlot.add_data_plot(mySEDdata,autoscale=True)

myPlot.save('SED_data.png')

mySEDdata.group_data(bin_width=0.2)

myPlot.save('SED_data_rebinned.png')

mySEDdata.add_systematics(0.1,[10.**6,10.**19])

#mySEDdata.set_error(0.1)

myPlot.add_data_plot(mySEDdata,label='grouped+syst')

myPlot.save('SED_data_rebinned.png')

SEDShape=SEDFit.SEDShape(mySEDdata)

SEDShape.eval_indices()

myPlot=SEDFit.Plot(mySEDdata,interactive=True)

myPlot.add_data_plot(mySEDdata,autoscale=True)

for model in SEDShape.index_models:

    myPlot.add_model_plot(model,label=model.name,line_style='--')

for index in SEDShape.index_models:
    print "=====>",index.name, SEDShape.indices.get_by_name(index.name).val.sed,SEDShape.indices.get_by_name(index.name).val.spectral
    
myPlot.save('SED_indices_rebinned.png')

myPlot=SEDFit.Plot(mySEDdata,interactive=True)

myPlot.add_data_plot(mySEDdata,autoscale=True)

SEDShape.sync_fit(check_host_gal_template=True)

myPlot.add_model_plot(SEDShape.sync_fit,label='sync, poly-fit')

myPlot.add_model_plot(SEDShape.host_gal,label='host-gal')

myPlot.add_model_plot(SEDShape.sync_fit_model,label='sync+host, poly-fit')

SEDShape.IC_fit()

myPlot.add_model_plot(SEDShape.IC,'IC, poly-fit')

SEDShape.show_values()

myPlot.save('SED_shaped_rebinned.png')

fit_Plot=SEDFit.Plot(mySEDdata,interactive=True)
    
fit_Plot.add_data_plot(mySEDdata,autoscale=True)

SED_obspar=SEDFit.ObsConstrain(beaming=25,B_range=[0.01,0.1],distr_e='lppl',t_var_sec=3*86400,nu_cut_IR=3.0E7,SEDShape=SEDShape) 

jet_model=SED_obspar.constrain_SSC_model()

fit_Plot.add_model_plot(jet_model.SED,label='obs-constr-lppl')

myPlot.save('obs_constr_lppl.png')
    
SEDModel=SEDFit.FitModel( jet=jet_model, name='SSC-best-fit',  template=SEDShape.host_gal)

SEDModel.set('z_cosm','frozen')

SEDModel.set('beam_obj','frozen')

SEDModel.set('nuFnu_p_host','frozen')

#SEDModel.set('L_host',fit_range=[-10.5,-9.5])

SEDModel.show_pars()
    
best_fit=SEDFit.fit_SED(SEDModel,mySEDdata,10.0**8 ,10**28.0,fitname='SSC-best-fit-lppl')

fit_Plot.add_model_plot(SEDModel,label='SSC-best-fit')

fit_Plot.add_residual_plot(SEDModel,autoscale=True)

myPlot.save('SSC_best_fit_lppl.png')




def test_best_fit(model):


    fit_Plot=SEDFit.Plot(mySEDdata,interactive=True)
    
    fit_Plot.add_data_plot(mySEDdata,autoscale=True)
    
    SED_obspar=SEDFit.ObsConstrain([15,25],[0.01,0.1],model,86400,5.0E10,SEDShape=SEDShape) 
    
    jet_model=SED_obspar.constrain_SSC_model()
    
    fit_Plot.add_model_plot(jet_model.SED,label='obs-constr-%s'%model)
    
    myPlot.save('obs_constr_%s.png'%model)
        
    SEDModel=SEDFit.FitModel( jet=jet_model, name='SSC-best-fit',  template=SEDShape.host)
    
    SEDModel.set('z_cosm','frozen')
    
    SEDModel.set('beam_obj','frozen')
    
    SEDModel.set('n',fit_range=[-10.5,-9.5])
   
    SEDModel.show_pars()
        
    best_fit=SEDFit.fit_SED(SEDModel,mySEDdata,10.0**11,10**27.5,fitname='SSC-best-fit-%s'%model)
    
    fit_Plot.add_model_plot(SEDModel,label='SSC-best-fit')
    
    fit_Plot.add_residual_plot(SEDModel,autoscale=True)

    myPlot.save('SSC_best_fit_%s.png'%model)

model_list=['bkn','plc']

#for model in model_list:

    #test_best_fit(model)
