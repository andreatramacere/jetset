import  pytest
from .base_class import TestBase
from .test_data import load_data_tool
class TestPhenomenologyConstr(TestBase):

    def integration_suite(self,sed_number=1,plot=False):
        prepare_asset(sed_number=sed_number,plot=plot)
        
        sed_data=load_data_tool(sed_number=sed_number,plot=False)
        my_shape=spectral_indices_tool(sed_data=sed_data,plot=plot,sed_number=sed_number)
        my_shape=sed_shaper_tool(my_shape=my_shape,sed_data=sed_data,sed_number=sed_number)
        prefit_jet, my_shape,template=model_constr_tool(my_shape=my_shape,sed_data=sed_data,plot=plot,sed_number=sed_number)
        
        return dict(jet=prefit_jet, shape=my_shape,template=template,sed_data=sed_data)

    def test_spectral_indices(self,sed_data=None,plot=True,sed_number=1):
        _=spectral_indices_tool(sed_data=sed_data,plot=plot,sed_number=sed_number)

    def test_sed_shaper(self,my_shape=None, plot=True, sed_number=1,sed_data=None,check_host_gal_template=False):
        _=sed_shaper_tool(my_shape=my_shape, plot=plot, sed_number=sed_number,sed_data=sed_data,check_host_gal_template=check_host_gal_template)

    def test_model_constr(self,my_shape=None, plot=True,sed_number=1,sed_data=None,check_host_gal_template=False):
        _=model_constr_tool(my_shape=my_shape, plot=plot,sed_number=sed_number,sed_data=sed_data,check_host_gal_template=check_host_gal_template)


def prepare_asset(sed_number=1,plot=False):
    sed_data=load_data_tool(sed_number=sed_number,plot=False)
    my_shape=spectral_indices_tool(sed_data=sed_data,plot=plot,sed_number=sed_number)
    my_shape=sed_shaper_tool(my_shape=my_shape,sed_data=sed_data,sed_number=sed_number)
    prefit_jet, my_shape,template=model_constr_tool(my_shape=my_shape,sed_data=sed_data,plot=plot,sed_number=sed_number)
    
    return dict(jet=prefit_jet, shape=my_shape,template=template,sed_data=sed_data)


def spectral_indices_tool(sed_data=None,plot=True,sed_number=1):
    from jetset.sed_shaper import SEDShape
    from .test_data import TestData
    if sed_data is None:
        sed_data=load_data_tool(sed_number=sed_number,plot=False)

    my_shape = SEDShape(sed_data)
    my_shape.eval_indices(silent=True)
    if plot is True:
        p = my_shape.plot_indices()
        p.setlim(y_min=1E-15, y_max=1E-6)

    return my_shape



def sed_shaper_tool(my_shape=None, plot=True, sed_number=1,sed_data=None,check_host_gal_template=False):
    if my_shape is None:
        my_shape=spectral_indices_tool(plot=plot,sed_number=sed_number,sed_data=sed_data)
    
    if sed_number == 2:
        check_host_gal_template=True
    
    print('==> check_host_gal_template',check_host_gal_template,sed_number)
    mm, best_fit = my_shape.sync_fit(check_host_gal_template=check_host_gal_template,
                                    Ep_start=None,
                                    minimizer='lsb',
                                    silent=True,
                                    fit_range=[10, 21])

    best_fit.show_report()
    best_fit.save_report('synch_shape_fit_rep.pkl')
    best_fit.bestfit_table.write('IC_shape_fit_rep.ecsv',overwrite=True)
    mm, best_fit= my_shape.IC_fit(fit_range=[23, 29], minimizer='minuit',silent=True)

    if plot is True:
        p = my_shape.plot_shape_fit()
        p.setlim(y_min=1E-15)

    best_fit.show_report()
    best_fit.save_report('IC_shape_fit_rep.pkl')
    best_fit.bestfit_table.write('IC_shape_fit_rep.ecsv',overwrite=True)
    my_shape.save_values('sed_shape_values.ecsv')
    return my_shape


def model_constr_tool(my_shape=None, plot=True,sed_number=1,sed_data=None,check_host_gal_template=False):
    from jetset.obs_constrain import ObsConstrain

    if sed_number == 2:
        check_host_gal_template=True
    print('==> check_host_gal_template', check_host_gal_template, sed_number)
    if my_shape is None:
        my_shape=sed_shaper_tool(plot=plot,sed_number=sed_number,sed_data=sed_data,check_host_gal_template=check_host_gal_template)

    sed_obspar = ObsConstrain(beaming=25,
                            B_range=[0.001, 0.1],
                            distr_e='lppl',
                            t_var_sec=3 * 86400,
                            nu_cut_IR=1E11,
                            SEDShape=my_shape)

    prefit_jet = sed_obspar.constrain_SSC_model(electron_distribution_log_values=False,silent=True)
    prefit_jet.save_model('prefit_jet.pkl')
    template=None
    if sed_number == 2:
        template = my_shape.host_gal
    return prefit_jet, my_shape, template