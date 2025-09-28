import  pytest
import astropy.units as u
import  numpy as np
from jetset.test_data_helper import  test_SEDs
from jetset.data_loader import ObsData,Data
from jetset.test_data_helper import  test_SEDs
from jetset.sed_shaper import  SEDShape
from .radio_plugin import RadioSpectrum
from jetset.gammapy_plugin import GammapyJetsetModelFactory
from gammapy.datasets import FluxPointsDataset,Datasets
from gammapy.estimators import FluxPoints
from jetset.obs_constrain import ObsConstrain
from jetset.model_manager import  FitModel
from gammapy.modeling.models import SkyModel
from gammapy.modeling import Fit

def test_composite_model_gp():
    data=Data.from_file(test_SEDs[1])
    sed_data=ObsData(data_table=data)
    sed_data.group_data(bin_width=0.1)
    sed_data.add_systematics(0.1,[10.**6,10.**29])
    my_shape=SEDShape(sed_data)
    my_shape.eval_indices(minimizer='lsb',silent=True)


    mm,best_fit=my_shape.sync_fit(check_host_gal_template=False,
                    Ep_start=None,
                    minimizer='lsb',
                    silent=True,
                    fit_range=[10.,21.])


    my_shape.IC_fit(fit_range=[23.,29.],minimizer='minuit',silent=True)




    sed_obspar=ObsConstrain(beaming=25,
                            B_range=[0.001,0.1],
                            distr_e='lppl',
                            t_var_sec=3*86400,
                            nu_cut_IR=1E12,
                            SEDShape=my_shape)


    jet=sed_obspar.constrain_SSC_model(electron_distribution_log_values=False,silent=True)


    jet.parameters.z_cosm.freeze()
    jet.parameters.R_H.freeze()
    jet.parameters.R.freeze()
    jet.parameters.gmin.freeze()

    jet.parameters.gmax.fit_range=[1E5,1E7]
    jet.parameters.s.fit_range=[1,3]
    jet.parameters.r.fit_range=[0,5]
    jet.parameters.B.fit_range=[1E-4,1]
    jet.parameters.N.fit_range=[1E-3,10]
    jet.parameters.gamma0_log_parab.fit_range=[1E3,1E5]
    jet.parameters.beam_obj.fit_range=[5,50]

    jet.add_user_par(name='B0',units='G',val=1E3,val_min=0,val_max=None)
    jet.add_user_par(name='R0', units='cm', val=5E13, val_min=0, val_max=None)
    jet.add_user_par(name='m_B', val=1, val_min=1, val_max=2)
    jet.parameters.R0.frozen=True
    jet.parameters.B0.frozen=True

    def par_func(R0,B0,R_H,m_B):
        return B0*np.power((R0/R_H),m_B)

    jet.make_dependent_par(par='B', depends_on=['B0', 'R0', 'R_H','m_B'], par_expr=par_func)

    jet.add_user_par(name='theta_open',val=3,units='deg',val_min=1,val_max=5)
    jet.make_dependent_par(par='R', depends_on=['R_H','theta_open'],par_expr='np.tan(np.radians(theta_open))*R_H')
    jet.parameters.R_H.free()
    jet.parameters.R_H.val=5E17

    composite_model=FitModel(jet=jet, name='gammapy',template=None)
    composite_model.add_component(RadioSpectrum())
    composite_model.radio_spectrum.parameters.nu_cut.val=2E10
    composite_model.radio_spectrum.parameters.nu_ssa.val=1E8

    composite_model.radio_spectrum.parameters.alpha_radio.val=-0.3
    composite_model.radio_spectrum.parameters.nuFnu_p.val=5E-14
    composite_model.radio_spectrum.parameters.alpha_radio.frozen=True
    composite_model.radio_spectrum.parameters.nu_ssa.frozen=True
    composite_model.jet_leptonic.parameters.m_B.frozen=True
    composite_model.composite_expr='jet_leptonic+radio_spectrum'


    gammapy_jet_model=GammapyJetsetModelFactory(composite_model)
    gammapy_jet_model.parameters.to_table()



    composite_model.show_model()
    fp=FluxPoints.from_table(sed_data.gammapy_table,sed_type='e2dnde', format='gadf-sed')
    
    sky_model = SkyModel(name="SSC model Mrk 421", spectral_model=gammapy_jet_model)

    gammapy_jet_model.evaluate()
    
    datasets = Datasets()
    E_min_fit = (1e9 * u.Hz).to("eV", equivalencies=u.spectral())
    fp=FluxPoints.from_table(sed_data.gammapy_table,sed_type='e2dnde', format='gadf-sed')
    dataset_mrk421 = FluxPointsDataset(data=fp,models=sky_model)

    #this workaround was needed with version 1.2
    dataset_mrk421.mask_fit= dataset_mrk421.data.energy_ref >= E_min_fit
    dataset_mrk421.mask_fit=dataset_mrk421.mask_fit.reshape(dataset_mrk421.mask_safe.shape)

    datasets = Datasets(dataset_mrk421)
    datasets.models=sky_model

    

    #conf_dict=dict(tol=1E-8)

    fitter = Fit(backend='scipy')

    results = fitter.run(datasets=datasets)






 