import  pytest

def test_spectral_indices(sed_data=None,plot=True,sed_number=None):
    from jetset.sed_shaper import SEDShape
    if sed_data is None:
        from .test_data import test_data_loader
        sed_data=test_data_loader(sed_number=sed_number)

    my_shape = SEDShape(sed_data)
    my_shape.eval_indices(silent=True)
    if plot is True:
        p = my_shape.plot_indices()
        p.rescale(y_min=-15, y_max=-6)

    return my_shape

def test_sed_shaper(my_shape=None, plot=True, sed_number=None,sed_data=None,check_host_gal_template=False):
    if my_shape is None:
        my_shape=test_spectral_indices(plot=plot,sed_number=sed_number,sed_data=sed_data)

    mm, best_fit = my_shape.sync_fit(check_host_gal_template=check_host_gal_template,
                                     Ep_start=None,
                                     minimizer='lsb',
                                     silent=True,
                                     fit_range=[10, 21])

    best_fit.show_report()
    best_fit.save_report('synch_shape_fit_rep.pkl')
    best_fit.bestfit_table.write('IC_shape_fit_rep.ecsv')
    mm, best_fit= my_shape.IC_fit(fit_range=[23, 29], minimizer='minuit',silent=True)

    if plot is True:
        p = my_shape.plot_shape_fit()
        p.rescale(y_min=-15)

    best_fit.show_report()
    best_fit.save_report('IC_shape_fit_rep.pkl')
    best_fit.bestfit_table.write('IC_shape_fit_rep.ecsv')
    my_shape.save_values('sed_shape_values.ecsv')
    return my_shape

def test_model_constr(my_shape=None, plot=True,sed_number=None,sed_data=None,check_host_gal_template=False):
    from jetset.obs_constrain import ObsConstrain

    if my_shape is None:
        my_shape=test_sed_shaper(plot=plot,sed_number=sed_number,sed_data=sed_data,check_host_gal_template=check_host_gal_template)

    sed_obspar = ObsConstrain(beaming=25,
                              B_range=[0.001, 0.1],
                              distr_e='lppl',
                              t_var_sec=3 * 86400,
                              nu_cut_IR=1E11,
                              SEDShape=my_shape)

    prefit_jet = sed_obspar.constrain_SSC_model(electron_distribution_log_values=False,silent=True)
    prefit_jet.save_model('prefit_jet.pkl')

    return prefit_jet, my_shape

