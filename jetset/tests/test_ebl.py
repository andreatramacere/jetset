import pytest

def test_ebl(plot=True):
    import numpy as np
    import matplotlib.pyplot as plt

    from jetset.template_2Dmodel import EBLAbsorptionTemplate
    ebl_dominguez = EBLAbsorptionTemplate.from_name('Dominguez_2010')
    ebl_finke = EBLAbsorptionTemplate.from_name('Finke_2010')
    ebl_franceschini = EBLAbsorptionTemplate.from_name('Franceschini_2008')

    z = 0.1
    nu = np.logspace(23, 30, 100)
    ebl_dominguez.parameters.z_cosm.val = z
    ebl_dominguez.eval(nu=nu)
    ebl_finke.parameters.z_cosm.val = z
    ebl_finke.eval(nu=nu)
    ebl_franceschini.parameters.z_cosm.val = z
    ebl_franceschini.eval(nu=nu)
    if plot is True:
        p = ebl_dominguez.plot_model()
        ebl_finke.plot_model(p)
        ebl_franceschini.plot_model(p)
        p.rescale(y_min=-10, x_max=29)

    nu = 1E26
    z_range = np.linspace(0.001, 1, 100)
    y_fr = np.zeros(z_range.size)
    y_fi = np.zeros(z_range.size)
    y_do = np.zeros(z_range.size)
    for ID, z in enumerate(z_range):
        ebl_franceschini.parameters.z_cosm.val = z
        ebl_finke.parameters.z_cosm.val = z
        ebl_dominguez.parameters.z_cosm.val = z
        y_fr[ID] = ebl_franceschini.eval(nu=nu, get_model=True)
        y_fi[ID] = ebl_finke.eval(nu=nu, get_model=True)
        y_do[ID] = ebl_dominguez.eval(nu=nu, get_model=True)

    if plot is True:
        plt.plot(z_range, y_fr, label='%s' % ebl_franceschini.name)
        plt.plot(z_range, y_fi, label='%s' % ebl_finke.name)
        plt.plot(z_range, y_do, label='%s' % ebl_dominguez.name)

        plt.xlabel('z')
        plt.ylabel(r'$exp^{-\tau}$')
        plt.legend()
        plt.semilogy()
        plt.title(r'$\nu=%1.1E Hz$' % nu)

    z_range = np.linspace(0.001, 1, 100)
    y_fr = np.zeros(z_range.size)
    y_fi = np.zeros(z_range.size)
    y_do = np.zeros(z_range.size)
    nu = 1E27
    for ID, z in enumerate(z_range):
        ebl_franceschini.parameters.z_cosm.val = z
        ebl_finke.parameters.z_cosm.val = z
        ebl_dominguez.parameters.z_cosm.val = z
        y_fr[ID] = ebl_franceschini.eval(nu=nu, get_model=True)
        y_fi[ID] = ebl_finke.eval(nu=nu, get_model=True)
        y_do[ID] = ebl_dominguez.eval(nu=nu, get_model=True)

    if plot is True:
        plt.plot(z_range, y_fr, label='%s' % ebl_franceschini.name)
        plt.plot(z_range, y_fi, label='%s' % ebl_finke.name)
        plt.plot(z_range, y_do, label='%s' % ebl_dominguez.name)

        plt.xlabel('z')
        plt.ylabel(r'$exp^{-\tau}$')
        plt.legend()
        plt.semilogy()
        plt.title(r'$\nu=%1.1E Hz$' % nu)


def test_ebl_jet(plot=True,):
    from jetset.jet_model import Jet
    from jetset.template_2Dmodel import EBLAbsorptionTemplate
    from jetset.model_manager import FitModel

    my_jet = Jet(electron_distribution='lppl', name='jet_flaring')
    my_jet.parameters.z_cosm.val = 0.01

    ebl_franceschini = EBLAbsorptionTemplate.from_name('Franceschini_2008')

    composite_model = FitModel(nu_size=500, name='EBL corrected')
    composite_model.add_component(my_jet)
    composite_model.add_component(ebl_franceschini)

    composite_model.show_pars()

    composite_model.link_par(par_name='z_cosm', from_model='Franceschini_2008', to_model='jet_flaring')
    v=0.03001
    my_jet.parameters.z_cosm.val = v
    assert (composite_model.Franceschini_2008.parameters.z_cosm.val==v)
    assert (composite_model.Franceschini_2008.parameters.z_cosm.linked==True)

    composite_model.composite_expr = '%s*%s'%(my_jet.name,ebl_franceschini.name)
    composite_model.eval()

    if plot is True:
        composite_model.plot_model()

    composite_model.save('ebl_jet.pkl')
    my_jet=composite_model.load('ebl_jet.pkl')
    my_jet.eval()

def test_ebl_jet_fit(plot=True):
    from .test_phenom_constr import test_model_constr

    prefit_jet, my_shape = test_model_constr()
    prefit_jet.show_model()

    from jetset.template_2Dmodel import EBLAbsorptionTemplate
    ebl_franceschini = EBLAbsorptionTemplate.from_name('Franceschini_2008')
    ebl_franceschini.show_model()

    from jetset.model_manager import FitModel
    composite_model = FitModel(nu_size=500, name='EBL corrected', template=my_shape.host_gal)
    composite_model.add_component(prefit_jet)
    composite_model.add_component(ebl_franceschini)
    composite_model.link_par(par_name='z_cosm', from_model='Franceschini_2008', to_model=prefit_jet.name)

    composite_model.composite_expr = '(%s+host_galaxy)*Franceschini_2008'%prefit_jet.name

    assert (composite_model.Franceschini_2008.parameters.z_cosm.val == prefit_jet.parameters.z_cosm.val)
    assert (composite_model.Franceschini_2008.parameters.z_cosm.linked is True)

    composite_model.show_model()
    composite_model.eval()

    if plot is True:
        composite_model.plot_model()

    from jetset.minimizer import ModelMinimizer


    composite_model.freeze(prefit_jet, 'R_H')
    composite_model.freeze(prefit_jet, 'z_cosm')
    composite_model.jet_leptonic.parameters.beam_obj.fit_range = [5, 50]
    composite_model.jet_leptonic.parameters.R.fit_range = [10 ** 15.5, 10 ** 17.5]
    composite_model.jet_leptonic.parameters.gmax.fit_range = [1E4, 1E8]

    composite_model.host_galaxy.parameters.nuFnu_p_host.frozen = False
    composite_model.host_galaxy.parameters.nu_scale.frozen = True
    composite_model.jet_leptonic.nu_size = 200
    composite_model.jet_leptonic.IC_nu_size = 100
    model_minimizer_lsb = ModelMinimizer('lsb')
    best_fit_lsb = model_minimizer_lsb.fit(composite_model, my_shape.sed_data, 1E11, 1E29, fitname='SSC-best-fit-lsb', repeat=1)
    best_fit_lsb.save_report('best-fit-minuit-report.pkl')
    best_fit_lsb.save_model('fit_model_minuit.pkl')