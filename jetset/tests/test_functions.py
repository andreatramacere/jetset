import sys




def data():
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
    p=sed_data.plot_sed()

    return sed_data

def spectral_indices(sed_data):
    from jetset.sed_shaper import SEDShape





    my_shape = SEDShape(sed_data)
    my_shape.eval_indices(silent=True)
    p = my_shape.plot_indices()
    p.rescale(y_min=-15, y_max=-6)

    return my_shape

def sed_shaper(my_shape):

    mm, best_fit = my_shape.sync_fit(check_host_gal_template=True,
                                     Ep_start=None,
                                     minimizer='lsb',
                                     silent=True,
                                     fit_range=[10, 21])

    best_fit.show_report()
    best_fit.save_report('synch_shape_fit_rep.txt')

    mm, best_fit= my_shape.IC_fit(fit_range=[23, 29], minimizer='minuit')
    p = my_shape.plot_shape_fit()
    p.rescale(y_min=-15)
    best_fit.show_report()
    best_fit.save_report('IC_shape_fit_rep.txt')
    my_shape.save_values('sed_shape_values.ecsv')
    return my_shape

def model_constr(my_shape):
    from jetset.obs_constrain import ObsConstrain

    sed_obspar = ObsConstrain(beaming=25,
                              B_range=[0.001, 0.1],
                              distr_e='lppl',
                              t_var_sec=3 * 86400,
                              nu_cut_IR=1E11,
                              SEDShape=my_shape)

    prefit_jet = sed_obspar.constrain_SSC_model(electron_distribution_log_values=False)
    prefit_jet.save_model('prefit_jet_gal_templ.dat')

    return prefit_jet

def model_fit_lsb(sed_data,my_shape):
    from jetset.minimizer import fit_SED
    from jetset.model_manager import FitModel
    from jetset.jet_model import Jet

    jet_lsb = Jet.load_model('prefit_jet_gal_templ.dat')
    jet_lsb.set_gamma_grid_size(200)

    fit_model_lsb = FitModel(jet=jet_lsb, name='SSC-best-fit-lsb', template=my_shape.host_gal)
    fit_model_lsb.freeze('z_cosm')
    fit_model_lsb.freeze('R_H')
    fit_model_lsb.parameters.beam_obj.fit_range = [5, 50]
    fit_model_lsb.parameters.R.fit_range = [10 ** 15.5, 10 ** 17.5]
    fit_model_lsb.parameters.gmax.fit_range = [1E4, 1E8]
    fit_model_lsb.parameters.nuFnu_p_host.frozen = False
    fit_model_lsb.parameters.nu_scale.frozen = True

    model_minimizer_lsb, best_fit_lsb = fit_SED(fit_model_lsb, sed_data, 10.0 ** 11, 10 ** 29.0,
                                                fitname='SSC-best-fit-lsb', minimizer='lsb')

    best_fit_lsb.save_report('best-fit-minuit-report.txt')
    fit_model_lsb.save_model('fit_model_lsb.dat')

def model_fit_minuit(sed_data,my_shape):
    from jetset.minimizer import fit_SED
    from jetset.model_manager import FitModel
    from jetset.jet_model import Jet

    jet_minuit = Jet.load_model('prefit_jet_gal_templ.dat')
    jet_minuit.set_gamma_grid_size(200)

    fit_model_minuit = FitModel(jet=jet_minuit, name='SSC-best-fit-minuit', template=my_shape.host_gal)
    fit_model_minuit.freeze('z_cosm')
    fit_model_minuit.freeze('R_H')
    fit_model_minuit.parameters.beam_obj.fit_range = [5, 50]
    fit_model_minuit.parameters.R.fit_range = [10 ** 15.5, 10 ** 17.5]
    fit_model_minuit.parameters.nuFnu_p_host.frozen = False
    fit_model_minuit.parameters.nu_scale.frozen = True

    model_minimizer_minuit, best_fit_minuit = fit_SED(fit_model_minuit, sed_data, 10.0 ** 11, 10 ** 29.0,
                                                      fitname='SSC-best-fit-minuit', minimizer='minuit')

    best_fit_minuit.save_report('best-fit-minuit-report.txt')
    fit_model_minuit.save_model('fit_model_minuit.dat')


def test_build_bessel():
    from jetset.jet_model import Jet
    Jet().eval()


def test_jet():
    from jetset.jet_model import Jet
    j=Jet()
    j.eval()
    j.energetic_report()
    j.plot_model()
    j.emitters_distribution.plot()
    j.emitters_distribution.plot2p()
    j.emitters_distribution.plot3p()
    j.emitters_distribution.plot3p(energy_unit='eV')
    j.emitters_distribution.plot3p(energy_unit='erg')

def test_full():
    from jetset.plot_sedfit import plt
    plt.ioff()
    test_jet()
    sed_data=data()
    print('done')
    my_shape=spectral_indices(sed_data)
    print('done')
    my_shape=sed_shaper(my_shape)
    print('done')
    prefit_jet=model_constr(my_shape)
    print('done')
    model_fit_lsb(sed_data,my_shape)
    model_fit_minuit(sed_data,my_shape)

def test_short():
    from jetset.plot_sedfit import  plt
    plt.ioff()
    test_jet()
    sed_data = data()
    print('done')
    my_shape = spectral_indices(sed_data)
    print('done')
    my_shape = sed_shaper(my_shape)
    print('done')
    prefit_jet = model_constr(my_shape)
    print('done')
    model_fit_lsb(sed_data,my_shape)
