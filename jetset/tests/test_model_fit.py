
def test_model_fit_minuit(sed_data=None, prefit_jet=None, plot=True):
    from jetset.minimizer import fit_SED,ModelMinimizer
    from jetset.model_manager import FitModel
    from jetset.jet_model import Jet

    if sed_data is None:
        from .test_data import   test_data_loader
        sed_data=test_data_loader(plot=plot)

    if prefit_jet is None:
        from .test_phenom_constr import test_model_constr
        prefit_jet, my_shape = test_model_constr()

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
    best_fit_minuit.save_report('best-fit-minuit-report.pkl')
    best_fit_minuit.bestfit_table.write('best-fit-minuit-report.ecsv')
    model_minimizer_minuit.save_model('model_minimizer_minuit.pkl')

    model_minimizer_minuit = ModelMinimizer.load_model('model_minimizer_lsb.pkl')
    fit_model_minuit.save_model('fit_model_minuit.pkl')

    return jet_minuit,model_minimizer_minuit


def test_model_fit_lsb(sed_data=None, prefit_jet=None, plot=True):
    from jetset.minimizer import fit_SED,ModelMinimizer
    from jetset.model_manager import FitModel
    from jetset.jet_model import Jet

    if sed_data is None:
        from .test_data import test_data_loader
        sed_data = test_data_loader(plot=plot)

    if prefit_jet is None:
        from .test_phenom_constr import test_model_constr
        prefit_jet, my_shape= test_model_constr()

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

    best_fit_lsb.save_report('best-fit-minuit-report.pkl')
    best_fit_lsb.show_report()
    fit_model_lsb.save_model('fit_model_lsb.pkl')
    fit_model_lsb_new = FitModel.load_model('fit_model_lsb.pkl')


    model_minimizer_lsb.save_model('model_minimizer_lsb.pkl')
    model_minimizer_lsb_new=ModelMinimizer.load_model('model_minimizer_lsb.pkl')

    return jet_lsb, model_minimizer_lsb_new, fit_model_lsb_new
