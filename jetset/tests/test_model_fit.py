import  pytest

def test_prepare_fit(sed_data=None, prefit_jet=None, plot=True,sed_number=1):
    from jetset.jet_model import Jet

    if sed_data is None:
        from .test_data import test_data_loader
        sed_data = test_data_loader(plot=plot, sed_number=sed_number)

    if prefit_jet is None:
        from .test_phenom_constr import test_model_constr
        prefit_jet, my_shape = test_model_constr(sed_data=sed_data)

    if hasattr(my_shape, 'host_gal'):
        template = my_shape.host_gal
    else:
        template = None

    jet = Jet.load_model('prefit_jet.pkl')
    jet.set_gamma_grid_size(200)

    return template,jet,sed_data



def test_model_fit(sed_data=None, prefit_jet=None, plot=True,sed_number=1,minimizer='lsb'):
    from jetset.minimizer import fit_SED,ModelMinimizer
    from jetset.model_manager import FitModel

    template, jet,sed_data=test_prepare_fit(sed_data=sed_data,prefit_jet=prefit_jet,plot=plot,sed_number=sed_number)

    fit_model = FitModel(jet=jet, name='SSC-best-fit-minuit', template=template)
    fit_model.freeze('jet_leptonic','z_cosm')
    fit_model.freeze('jet_leptonic','R_H')
    fit_model.jet_leptonic.parameters.beam_obj.fit_range = [5, 50]
    fit_model.jet_leptonic.parameters.R.fit_range = [10 ** 15.5, 10 ** 17.5]
    if minimizer == 'minuit':
        fit_model.jet_leptonic.parameters.gmin.fit_range = [2, 200]
        fit_model.jet_leptonic.parameters.gmax.fit_range = [1E5, 1E7]
        fit_model.jet_leptonic.parameters.B.fit_range = [1E-3,2]

    if template is not None:
        fit_model.host_galaxy.parameters.nuFnu_p_host.frozen = False
        fit_model.host_galaxy.parameters.nu_scale.frozen = True

    model_minimizer = ModelMinimizer(minimizer)
    best_fit = model_minimizer.fit(fit_model, sed_data, 10 ** 11., 10 ** 29.0,
                                                 fitname='SSC-best-fit-minuit', repeat=1)
    best_fit.show_report()
    best_fit.save_report('best-fit-%s-report.pkl'%minimizer)
    best_fit.bestfit_table.write('best-fit-%s-report.ecsv'%minimizer)
    model_minimizer.save_model('model_minimizer_%s.pkl'%minimizer)

    model_minimizer = ModelMinimizer.load_model('model_minimizer_%s.pkl'%minimizer)
    fit_model.save_model('fit_model_%s.pkl'%minimizer)

    return fit_model, model_minimizer,sed_data

def test_model_fit_dep_pars(sed_data=None, prefit_jet=None, plot=True,sed_number=1,minimizer='lsb'):
    from jetset.minimizer import fit_SED,ModelMinimizer
    from jetset.model_manager import FitModel

    template, jet,sed_data=test_prepare_fit(sed_data=sed_data,prefit_jet=prefit_jet,plot=plot,sed_number=sed_number)

    fit_model = FitModel(jet=jet, name='SSC-best-fit-minuit', template=template)
    fit_model.freeze('jet_leptonic','z_cosm')
    fit_model.jet_leptonic.parameters.beam_obj.fit_range = [5, 50]
    fit_model.jet_leptonic.parameters.R_H.val = 5E17
    fit_model.jet_leptonic.parameters.R_H.frozen = False
    fit_model.jet_leptonic.parameters.R_H.fit_range = [1E15, 1E19]
    fit_model.jet_leptonic.parameters.R.fit_range = [10 ** 15.5, 10 ** 17.5]
    fit_model.jet_leptonic.add_user_par(name='B0', units='G', val=1E3, val_min=0, val_max=None)
    fit_model.jet_leptonic.add_user_par(name='R0', units='cm', val=5E13, val_min=0, val_max=None)
    fit_model.jet_leptonic.parameters.R0.frozen = True
    fit_model.jet_leptonic.parameters.B0.frozen = True

    fit_model.jet_leptonic.make_dependent_par(par='B', depends_on=['B0', 'R0', 'R_H'], par_expr='B0*(R0/R_H)')

    if minimizer == 'minuit':
        fit_model.jet_leptonic.parameters.gmin.fit_range = [2, 200]
        fit_model.jet_leptonic.parameters.gmax.fit_range = [1E5, 1E7]
        fit_model.jet_leptonic.parameters.B.fit_range = [1E-3,2]

    if template is not None:
        fit_model.host_galaxy.parameters.nuFnu_p_host.frozen = False
        fit_model.host_galaxy.parameters.nu_scale.frozen = True

    model_minimizer = ModelMinimizer(minimizer)
    best_fit = model_minimizer.fit(fit_model, sed_data, 10 ** 11., 10 ** 29.0,
                                                 fitname='SSC-best-fit-minuit', repeat=1)
    best_fit.show_report()
    best_fit.save_report('best-fit-%s-report.pkl'%minimizer)
    best_fit.bestfit_table.write('best-fit-%s-report.ecsv'%minimizer)
    model_minimizer.save_model('model_minimizer_%s.pkl'%minimizer)

    model_minimizer = ModelMinimizer.load_model('model_minimizer_%s.pkl'%minimizer)
    fit_model.save_model('fit_model_%s.pkl'%minimizer)

    return fit_model, model_minimizer,sed_data

