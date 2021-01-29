import pytest

#@pytest.fixture
#def plot():
#   input = False
#   return input


def hadronic_model(plot=True):
    from jetset.jet_model import Jet
    j = Jet(proton_distribution='plc')
    j.parameters.gmin.val = 2
    j.parameters.gmax.val = 1E8
    j.parameters.NH_pp.val = 1E10
    j.parameters.N.val = 1E1
    j.parameters.B.val = 80

    j.parameters.p.val = 2.5
    j.eval()
    j.show_model()
    if plot is True:
        j.plot_model()



def custom_emitters(plot=True):
    from jetset.jet_model import Jet

    from jetset.jet_emitters import EmittersDistribution
    import numpy as np
    def distr_func_bkn(gamma_break, gamma, s1, s2):
        return np.power(gamma, -s1) * (1. + (gamma / gamma_break)) ** (-(s2 - s1))

    n_e = EmittersDistribution('custom_bkn',spectral_type='bkn')
    n_e.add_par('gamma_break', par_type='turn-over-energy', val=1E3, vmin=1., vmax=None, unit='lorentz-factor')
    n_e.add_par('s1', par_type='LE_spectral_slope', val=2.5, vmin=-10., vmax=10, unit='')
    n_e.add_par('s2', par_type='LE_spectral_slope', val=3.2, vmin=-10., vmax=10, unit='')
    n_e.set_distr_func(distr_func_bkn)
    n_e.parameters.show_pars()
    n_e.parameters.s1.val = 2.0
    n_e.parameters.s2.val = 3.5
    if plot is True:
        n_e.plot()

    my_jet= Jet(electron_distribution=n_e)
    my_jet.Norm_distr = True
    my_jet.parameters.N.val = 5E4
    my_jet.eval()
    diff= np.fabs(np.trapz(n_e.n_gamma_e, n_e.gamma_e )- my_jet.parameters.N.val)
    print('diff',diff)
    assert ( diff<1E-3)

def data(plot=True):
    from jetset.data_loader import ObsData, Data
    from jetset.test_data_helper import test_SEDs

    data = Data.from_file(test_SEDs[0])
    sed_data = ObsData(data_table=data)
    sed_data.show_data_sets()
    sed_data.filter_data_set('-1', exclude=True)
    sed_data.show_data_sets()
    sed_data.reset_data()
    sed_data.show_data_sets()
    sed_data.filter_data_set('-1', exclude=False)
    sed_data.show_data_sets()
    sed_data.reset_data()

    data=Data.from_file(test_SEDs[2])
    sed_data=ObsData(data_table=data)
    sed_data.group_data(bin_width=0.2)

    sed_data.add_systematics(0.1,[10.**6,10.**29])
    if plot is True:
        p=sed_data.plot_sed()

    return sed_data

def spectral_indices(sed_data,plot=True):
    from jetset.sed_shaper import SEDShape

    my_shape = SEDShape(sed_data)
    my_shape.eval_indices(silent=True)
    if plot is True:
        p = my_shape.plot_indices()
        p.rescale(y_min=-15, y_max=-6)

    return my_shape

def sed_shaper(my_shape, plot=True):

    mm, best_fit = my_shape.sync_fit(check_host_gal_template=True,
                                     Ep_start=None,
                                     minimizer='lsb',
                                     silent=True,
                                     fit_range=[10, 21])

    best_fit.show_report()
    best_fit.save_report('synch_shape_fit_rep.txt')

    mm, best_fit= my_shape.IC_fit(fit_range=[23, 29], minimizer='minuit',silent=True)

    if plot is True:
        p = my_shape.plot_shape_fit()
        p.rescale(y_min=-15)

    best_fit.show_report()
    best_fit.save_report('IC_shape_fit_rep.txt')
    my_shape.save_values('sed_shape_values.ecsv')
    return my_shape

def model_constr(my_shape, plot=True):
    from jetset.obs_constrain import ObsConstrain

    sed_obspar = ObsConstrain(beaming=25,
                              B_range=[0.001, 0.1],
                              distr_e='lppl',
                              t_var_sec=3 * 86400,
                              nu_cut_IR=1E11,
                              SEDShape=my_shape)

    prefit_jet = sed_obspar.constrain_SSC_model(electron_distribution_log_values=False,silent=True)
    prefit_jet.save_model('prefit_jet_gal_templ.pkl')

    return prefit_jet

def model_fit_lsb(sed_data,my_shape, plot=True):
    from jetset.minimizer import fit_SED,ModelMinimizer
    from jetset.model_manager import FitModel
    from jetset.jet_model import Jet

    jet_lsb = Jet.load_model('prefit_jet_gal_templ.pkl')
    jet_lsb.set_gamma_grid_size(200)

    fit_model_lsb = FitModel(jet=jet_lsb, name='SSC-best-fit-lsb', template=my_shape.host_gal)
    fit_model_lsb.freeze('jet_leptonic','z_cosm')
    fit_model_lsb.freeze('jet_leptonic','R_H')
    fit_model_lsb.jet_leptonic.parameters.beam_obj.fit_range = [5, 50]
    fit_model_lsb.jet_leptonic.parameters.R.fit_range = [10 ** 15.5, 10 ** 17.5]
    fit_model_lsb.jet_leptonic.parameters.gmax.fit_range = [1E4, 1E8]
    fit_model_lsb.host_galaxy.parameters.nuFnu_p_host.frozen = False
    fit_model_lsb.host_galaxy.parameters.nu_scale.frozen = True

    model_minimizer_lsb, best_fit_lsb = fit_SED(fit_model_lsb, sed_data, 10.0 ** 11, 10 ** 29.0,
                                                fitname='SSC-best-fit-lsb', minimizer='lsb',silent=True)

    best_fit_lsb.save_report('best-fit-minuit-report.txt')
    best_fit_lsb.show_report()
    fit_model_lsb.save_model('fit_model_lsb.pkl')
    fit_model_lsb_new = FitModel.load_model('fit_model_lsb.pkl')

    model_minimizer_lsb.save_model('model_minimizer_lsb.pkl')
    model_minimizer_lsb_new=ModelMinimizer.load_model('model_minimizer_lsb.pkl')


    return jet_lsb, model_minimizer_lsb_new, fit_model_lsb_new


def model_fit_minuit(sed_data,my_shape, plot=True):
    from jetset.minimizer import fit_SED
    from jetset.model_manager import FitModel
    from jetset.jet_model import Jet

    jet_minuit = Jet.load_model('prefit_jet_gal_templ.pkl')
    jet_minuit.set_gamma_grid_size(200)

    fit_model_minuit = FitModel(jet=jet_minuit, name='SSC-best-fit-minuit', template=my_shape.host_gal)
    fit_model_minuit.freeze('jet_leptonic','z_cosm')
    fit_model_minuit.freeze('jet_leptonic','R_H')
    fit_model_minuit.jet_leptonic.parameters.beam_obj.fit_range = [5, 50]
    fit_model_minuit.jet_leptonic.parameters.R.fit_range = [10 ** 15.5, 10 ** 17.5]
    fit_model_minuit.host_galaxy.parameters.nuFnu_p_host.frozen = False
    fit_model_minuit.host_galaxy.parameters.nu_scale.frozen = True

    model_minimizer_minuit, best_fit_minuit = fit_SED(fit_model_minuit, sed_data, 10.0 ** 11, 10 ** 29.0,
                                                      fitname='SSC-best-fit-minuit', minimizer='minuit',silent=True)
    best_fit_minuit.show_report()
    best_fit_minuit.save_report('best-fit-minuit-report.txt')
    fit_model_minuit.save_model('fit_model_minuit.pkl')

    return jet_minuit,model_minimizer_minuit

def test_emcee(jet_bf_model,model_minimizer):
    from jetset.mcmc import McmcSampler

    mcmc = McmcSampler(model_minimizer)
    mcmc.run_sampler(nwalkers=150, burnin=10, steps=50,threads=2)



def test_build_bessel():
    from jetset.jet_model import Jet
    Jet().eval()

@pytest.mark.users
def test_jet(plot=True):
    print('--------> test_jet',plot)
    from jetset.jet_model import Jet
    j=Jet()
    j.eval()
    j.energetic_report()

    if plot is True:
        j.plot_model()
        j.emitters_distribution.plot()
        j.emitters_distribution.plot2p()
        j.emitters_distribution.plot3p()
        j.emitters_distribution.plot3p(energy_unit='eV')
        j.emitters_distribution.plot3p(energy_unit='erg')
    j.save_model('test_jet.pkl')
    j_new=Jet.load_model('test_jet.pkl')


def test_full(plot=False):
    from jetset.plot_sedfit import plt
    plt.ioff()
    test_jet(plot)
    hadroinc_model(plot)
    sed_data=data()
    print('done')
    my_shape=spectral_indices(sed_data)
    print('done')
    my_shape=sed_shaper(my_shape)
    print('done')
    prefit_jet=model_constr(my_shape)
    print('done')
    jet_lsb, model_minimizer_lsb, fit_model_lsb =model_fit_lsb(sed_data,my_shape)
    jet_minuit,model_minimizer_minuit=model_fit_minuit(sed_data,my_shape)

def test_short(plot=False):
    #from jetset.plot_sedfit import  plt
    #plt.ioff()
    test_jet(plot)
    sed_data = data(plot)
    print('done')
    my_shape = spectral_indices(sed_data,plot)
    print('done')
    my_shape = sed_shaper(my_shape,plot)
    print('done')
    prefit_jet = model_constr(my_shape,plot)
    print('done')
    jet_lsb, model_minimizer_lsb,fit_model_lsb = model_fit_lsb(sed_data, my_shape,plot)
    print('TEST PASSED: OK')


@pytest.mark.users
def test_users():
    test_short(plot=False)