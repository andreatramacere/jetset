import  pytest
from .base_class import TestBase


class TestModelFit(TestBase):

    def test(self,plot=False):
        self.integration_suite(plot=plot,sed_number=1)

    def integration_suite(self,plot=False,sed_number=1,phenom_dict=None,use_dep_pars=False,use_ebl=False,skip_minuit=False):
        if sed_number is not None:
            template,jet,sed_data=self.prepare_model(plot=plot,sed_number=sed_number)
        elif phenom_dict is not None:
            template = phenom_dict['template']
            jet = phenom_dict['jet']
            sed_data = phenom_dict['sed_data']
        else:
            raise RuntimeError('please provide either sed_number or  phenom_dict')

        fit_model=self.prepare_fit_model(jet=jet,use_dep_pars=use_dep_pars,use_ebl=use_ebl,template=template)
        fit_model, model_minimizer,sed_data=self.model_fit(fit_model=fit_model,
                                                            sed_data=sed_data, 
                                                            minimizer='lsb')


        if skip_minuit is False:
            fit_model, model_minimizer,sed_data=self.model_fit(fit_model=fit_model,
                                                                sed_data=sed_data, 
                                                                minimizer='minuit')

        return dict(fit_model=fit_model, model_minimizer=model_minimizer,sed_data=sed_data)

    def prepare_model(self, plot=True,sed_number=1):
        from .test_phenom_constr import TestPhenomenologyConstr
        t=TestPhenomenologyConstr()
        phenom_dict=t.integration_suite(sed_number=sed_number,plot=plot)
        prefit_jet=phenom_dict['jet']
        sed_data=phenom_dict['sed_data']
        template=phenom_dict['template']
        prefit_jet.set_gamma_grid_size(200)
        return template,prefit_jet,sed_data

    def prepare_fit_model(self,jet,use_dep_pars=False,use_ebl=False,template=None,minimizer='lsb'):
        from jetset.model_manager import FitModel

    
        fit_model = FitModel(jet=jet, name='SSC-best-fit-minuit', template=template)

        if use_ebl is True:
            from jetset.template_2Dmodel import EBLAbsorptionTemplate
            ebl_franceschini = EBLAbsorptionTemplate.from_name('Franceschini_2008')
            ebl_franceschini.show_model()
            fit_model.add_component(ebl_franceschini)
            fit_model.link_par(par_name='z_cosm', from_model='Franceschini_2008', to_model=jet.name)

        fit_model.freeze('jet_leptonic','z_cosm')
        fit_model.jet_leptonic.parameters.beam_obj.fit_range = [5, 50]
        fit_model.jet_leptonic.parameters.R_H.val = 5E17
        fit_model.jet_leptonic.parameters.R_H.frozen = False
        fit_model.jet_leptonic.parameters.R_H.fit_range = [1E15, 1E19]
        fit_model.jet_leptonic.parameters.R.fit_range = [10 ** 15.5, 10 ** 17.5]
        
        if use_dep_pars is True:
            if fit_model.parameters.get_par_by_name('jet_leptonic','B0') is None:
                fit_model.jet_leptonic.add_user_par(name='B0', units='G', val=1E3, val_min=0, val_max=None)
            
            if fit_model.parameters.get_par_by_name('jet_leptonic','R0') is None:
                fit_model.jet_leptonic.add_user_par(name='R0', units='cm', val=5E13, val_min=0, val_max=None)
            
            fit_model.jet_leptonic.parameters.R0.frozen = True
            fit_model.jet_leptonic.parameters.B0.frozen = True

            if fit_model.jet_leptonic.parameters.B._is_dependent  is False:
                fit_model.jet_leptonic.make_dependent_par(par='B', depends_on=['B0', 'R0', 'R_H'], par_expr='B0*(R0/R_H)')


        if use_ebl is True:
            if template is not None:
                fit_model.composite_expr = '(%s+host_galaxy)*Franceschini_2008'%jet.name
            else:
                fit_model.composite_expr = '(%s)*Franceschini_2008' % jet.name
        else:
            if template is not None:
                fit_model.composite_expr = '(%s+host_galaxy)'%jet.name
          

       

        if template is not None:
            fit_model.host_galaxy.parameters.nuFnu_p_host.frozen = False
            fit_model.host_galaxy.parameters.nu_scale.frozen = True

        
        return fit_model
        

    def model_fit(self,fit_model,sed_data, plot=False,minimizer='lsb'):
        
        from jetset.minimizer import ModelMinimizer
        
        model_minimizer = ModelMinimizer(minimizer)
        
        if minimizer == 'minuit':
            fit_model.jet_leptonic.parameters.gmin.fit_range = [2, 300]
            fit_model.jet_leptonic.parameters.gmax.fit_range = [1E5, 1E7]

        best_fit = model_minimizer.fit(fit_model, sed_data, 10 ** 11., 10 ** 29.0,
                                                    fitname='SSC-best-fit-%s'%minimizer, repeat=1)
        best_fit.show_report()
        best_fit.save_report('best-fit-%s-report.pkl'%minimizer)
        best_fit.bestfit_table.write('best-fit-%s-report.ecsv'%minimizer,overwrite=True)
        model_minimizer.save_model('model_minimizer_%s.pkl'%minimizer)

        model_minimizer = ModelMinimizer.load_model('model_minimizer_%s.pkl'%minimizer)
        fit_model.save_model('fit_model_%s.pkl'%minimizer)

        return fit_model, model_minimizer,sed_data

   