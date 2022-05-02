import pytest
from astropy.table import Table
from .base_class import TestBase


class TestCompositeModel(TestBase):

    def integration_suite(self,plot=False):
        self._all(plot=plot)

    def test_composite_model_pars(self,plot=False):
        from jetset.model_manager import  FitModel

        from jetset.template_2Dmodel import EBLAbsorptionTemplate
        ebl_franceschini=EBLAbsorptionTemplate.from_name('Franceschini_2008')

        from jetset.jet_model import Jet
        jet=Jet()

        jet.set_gamma_grid_size(100)
        fit_model=FitModel(jet=jet, name='test')

        fit_model.add_component(ebl_franceschini)
        fit_model.link_par(par_name='z_cosm', from_model='Franceschini_2008', to_model='jet_leptonic')
        fit_model.composite_expr='(jet_leptonic)*Franceschini_2008'
        fit_model.freeze('jet_leptonic','gmin')
        fit_model.parameters._build_par_table()
        fit_model.save_model('fit_model.pkl')
    
        _fit_model=FitModel.load_model('fit_model.pkl')
        _fit_model.jet_leptonic.parameters.gmin.frozen=False
        _fit_model.parameters._build_par_table()
        _t=_fit_model.parameters._par_table
        _t.sort(_t.colnames)
        

        _fit_model_clone=fit_model.clone()
        _fit_model_clone.jet_leptonic.parameters.gmin.frozen=False
        _fit_model_clone.parameters._build_par_table()
        _t_clone=_fit_model_clone.parameters._par_table
        _t_clone.sort(_t.colnames)
    

        assert(all([_t[i]['frozen']==_t_clone[i]['frozen'] for i in range(len(_t))]))
    
        assert(all([_t[i]['val']==_t_clone[i]['val'] for i in range(len(_t)) if _t.mask[i] is False]))
        
        fit_model.jet_leptonic.parameters.gmin.frozen=False
        fit_model.parameters._build_par_table()
        t=fit_model.parameters._par_table
        t.sort(_t.colnames)
        
        
        assert(all([_t[i]['frozen']==t[i]['frozen'] for i in range(len(_t))]))
        assert(all([_t[i]['val']==t[i]['val'] for i in range(len(_t)) if _t.mask[i] is False]))
        
        for _m in [fit_model,_fit_model,_fit_model_clone]:
            p_m=_m.parameters.get_par_by_name('jet_leptonic','z_cosm')
            p_l=_m.parameters.get_par_by_name('Franceschini_2008','z_cosm')
            assert( p_m._is_dependent == False )
            assert( p_l._is_dependent == True )
            p_l_m=p_l._master_pars[0]
            assert(p_l_m.name == p_m.name)