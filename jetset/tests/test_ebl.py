import pytest
from .base_class import TestBase

class TestEBL(TestBase):

    def integration_suite(self,plot=False):
        self.test_ebl(plot=plot)
        self.test_ebl_jet(plot=plot)

    def test_ebl(self,plot=True):
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
            p.setlim(y_min=1E-10, x_max=1E29)

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


    def test_ebl_jet(self,plot=True):
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
        assert (composite_model.Franceschini_2008.parameters.z_cosm.val == composite_model.jet_flaring.parameters.z_cosm.val)
        composite_model.composite_expr = '%s*%s'%(my_jet.name,ebl_franceschini.name)
        composite_model.eval()

        if plot is True:
            composite_model.plot_model()

        composite_model.save_model('ebl_jet.pkl')
        new_composite_model=FitModel.load_model('ebl_jet.pkl')
        v=2.0
        new_composite_model.jet_flaring.parameters.z_cosm.val=v
        assert (new_composite_model.Franceschini_2008.parameters.z_cosm.val == v)
        assert (new_composite_model.Franceschini_2008.parameters.z_cosm.linked == True)

    def test_ebl_jet_fit(self,plot=True,sed_number=2,minimizer='lsb'):
        from .test_model_fit import TestModelFit
        

        template, jet,sed_data = TestModelFit().prepare_model(sed_number=sed_number)

        from jetset.template_2Dmodel import EBLAbsorptionTemplate
        ebl_franceschini = EBLAbsorptionTemplate.from_name('Franceschini_2008')
        ebl_franceschini.show_model()

        from jetset.model_manager import FitModel
        composite_model = FitModel(nu_size=500, name='EBL corrected', template=template)
        composite_model.add_component(jet)
        composite_model.add_component(ebl_franceschini)
        composite_model.link_par(par_name='z_cosm', from_model='Franceschini_2008', to_model=jet.name)

        if template is not None:
            composite_model.composite_expr = '(%s+host_galaxy)*Franceschini_2008'%jet.name
        else:
            composite_model.composite_expr = '(%s)*Franceschini_2008' % jet.name

        assert (composite_model.Franceschini_2008.parameters.z_cosm.val == composite_model.jet_leptonic.parameters.z_cosm.val)
        assert (composite_model.Franceschini_2008.parameters.z_cosm.linked is True)

        composite_model.show_model()
        composite_model.eval()

        if plot is True:
            composite_model.plot_model()

        from jetset.minimizer import ModelMinimizer

        composite_model.freeze(jet, 'R_H')
        composite_model.freeze(jet, 'z_cosm')

        if minimizer == 'minuit':
            composite_model.jet_leptonic.parameters.gmin.fit_range = [2, 200]
            composite_model.jet_leptonic.parameters.gmax.fit_range = [1E5, 1E7]

        if template is not None:
            composite_model.host_galaxy.parameters.nuFnu_p_host.frozen = False
            composite_model.host_galaxy.parameters.nu_scale.frozen = True

        composite_model.jet_leptonic.parameters.beam_obj.fit_range = [5, 50]
        composite_model.jet_leptonic.parameters.R.fit_range = [10 ** 15.5, 10 ** 17.5]
        composite_model.jet_leptonic.parameters.gmax.fit_range = [1E4, 1E8]


        composite_model.jet_leptonic.nu_size = 200
        composite_model.jet_leptonic.IC_nu_size = 100

        model_minimizer = ModelMinimizer(minimizer)
        best_fit = model_minimizer.fit(composite_model, sed_data, 10 ** 11., 10 ** 29.0,
                                    fitname='SSC-best-fit-minuit', repeat=1)
        best_fit.show_report()
        best_fit.save_report('best-fit-%s-report.pkl' % minimizer)
        best_fit.bestfit_table.write('best-fit-%s-report.ecsv' % minimizer,overwrite=True)
        model_minimizer.save_model('model_minimizer_%s.pkl' % minimizer)

        model_minimizer = ModelMinimizer.load_model('model_minimizer_%s.pkl' % minimizer)
        composite_model.save_model('fit_model_%s.pkl' % minimizer)