
__author__ = "Andrea Tramacere"


import warnings

import numpy as np

from  . import minimizer



from .model_parameters import  CompositeModelParameterArray

from .spectral_shapes import  SED
   
from .base_model import  Model

from .plot_sedfit import  PlotSED

from .utils import  clean_var_name

from .jet_model import Jet

from .cosmo_tools import  Cosmo

import  dill as pickle

__all__=['FitModel']

class CompositeModelContainer(object):

    def __init__(self):
        self._components_list=[]
        self._components_value=[]
        self._components_value_dict = {}
        self.parameters=CompositeModelParameterArray()

    def add_component(self, model_comp,fit_model):
        try:
            assert (model_comp.name not in [m_comp.name for m_comp in self._components_list])
        except Exception as e:
            raise RuntimeError('model name:', model_comp.name, 'already assigned',e)

        try:
            assert (model_comp not in self._components_list)
        except Exception as e:
            raise RuntimeError('model:', model_comp, 'already added',e)

        self._components_list.append(model_comp)
        self._components_value.append(None)

        self._components_value_dict[model_comp.name] = self._components_value[-1]

        self.parameters.add_model_parameters(model_comp)
        setattr(self, clean_var_name(model_comp.name), model_comp)
        setattr(fit_model, clean_var_name(model_comp.name), model_comp)

    @property
    def components_list(self):
        return self._components_list

    def get_model_by_name(self, model_name,get_idx=False):
        model = None
        idx=None
        for ID,m in enumerate(self._components_list):
            if m.name == model_name:
                model = m
                idx=ID
        if model is None:
            warnings.warn('Model',model_name,'not present')
        if get_idx is False:
            return model
        else:
            return model,idx

    def del_component(self,model_name,fit_model):
        m,ID=self.get_model_by_name(model_name,get_idx=True)
        if m is not None:
            _p = self._components_list.pop(ID)
            _p = self._components_value.pop(ID)
            _p = self.parameters.del_model_parameters(m)
            del self._components_value_dict[m]

        delattr(self,model_name)
        delattr(fit_model, model_name)

    def show_pars(self):
        self.parameters.show_pars()

    def show_model(self):
        for c in self._components_list:
            #print()
            c.show_model()
            #print()





class FitModel(Model):
    """
    """
    
    def __init__(self,
                 elec_distr=None,
                 jet=None,
                 name='no-name',
                 out_dir=None,
                 flag=None,
                 template=None,
                 loglog_poly=None,
                 analytical=None,
                 nu_size=200,
                 cosmo=None,
                 composite_expr=None,
                 **keywords):

        super(FitModel,self).__init__( model_type='composite_model', **keywords)
       
        if  jet is not None and elec_distr is not None:
            #!! warning or error?
            raise RuntimeError("you can't provide both elec_distr and jet, only one")

        self.sed_data=None
        self.nu_min_fit=1E6
        self.nu_max_fit=1E30
        self.name=name
        self.SED=SED(name=self.name)
        self.nu_min=1E6
        self.nu_max=1E30
        self.nu_size=nu_size
        self.flux_plot_lim=1E-30
        self.components=CompositeModelContainer()
        self.parameters=self.components.parameters
        self.composite_expr = composite_expr

        if elec_distr is not None:
            jet=Jet(cosmo=cosmo,name=flag, electron_distribution=elec_distr, jet_workplace=None)
            self.add_component(jet)
        
        if jet is not None:
            self.add_component(jet)

        if cosmo is None:
          self.cosmo=Cosmo()
          warnings.warn('no cosmology defined, using default %s'%self.cosmo)
        else:
            self.cosmo = cosmo

        if template is not None:
            self.add_component(template)

        if loglog_poly is not None:
            self.add_component(loglog_poly)

        if analytical is not None:
            self.add_component(analytical)

    def plot_model(self,plot_obj=None,clean=False,sed_data=None,frame='obs',skip_components=False,label=None,skip_sub_components=False, density=False):
        plot_obj=self._set_up_plot(plot_obj,sed_data,frame,density)

        if clean is True:
            plot_obj.clean_model_lines()

        if skip_components is False:
            line_style = '--'

            for mc in self.components._components_list:
                comp_label = mc.name
                if hasattr(mc,'SED'):
                    try:
                        plot_obj.add_model_plot(mc.SED, line_style=line_style,label=comp_label,flim=self.flux_plot_lim, density=density, frame=frame)
                    except Exception as e:
                        try:
                            mc.eval()
                            plot_obj.add_model_plot(mc.SED, line_style=line_style, label=comp_label,
                                                    flim=self.flux_plot_lim, density=density, frame=frame)
                        except Exception as e:
                            raise RuntimeError('for model', mc.name, "problem with plotting SED", e)


                if skip_sub_components is False:
                    if hasattr(mc,'spectral_components_list'):
                        for c in mc.spectral_components_list:
                            comp_label = c.name
                            if comp_label!='Sum':
                                if hasattr(c, 'SED'):
                                    try:
                                        plot_obj.add_model_plot(c.SED, line_style=line_style, label='  -%s'%comp_label, flim=self.flux_plot_lim, density=density, frame=frame)
                                    except Exception as e:
                                        try:
                                            #print('==> reval',mc.name)
                                            mc.eval()
                                            _c=mc.spectral_components.get_spectral_component_by_name(c.name)
                                            plot_obj.add_model_plot(_c.SED, line_style=line_style, label='  -%s'%comp_label, flim=self.flux_plot_lim, density=density, frame=frame)
                                        except Exception as e:
                                            raise RuntimeError('for model', mc.name, "spectral component",c.name, "problem with plotting SED", e)
        line_style = '-'
        if label is None:
            label=self.name

        plot_obj.add_model_plot(self.SED, line_style=line_style, label=label, flim=self.flux_plot_lim,fit_range=[self.nu_min_fit,self.nu_max_fit], density=density, frame=frame  )
        plot_obj.add_model_residual_plot(data=sed_data, model=self, fit_range=[self.nu_min_fit, self.nu_max_fit])

        #if frame == 'src' and sed_data is not None:
        #    sed_data.z = z_sed_data

        return plot_obj



    def set_nu_grid(self,nu_min=None,nu_max=None,nu_size=None):
        if nu_size is not None:
            self.nu_size=nu_size
        
        if nu_min is not None:
            self.nu_min=nu_min
        
        if nu_max is not None:
            self.nu_max=nu_max
        
        for model_comp in self.components._components_list:
        
            if nu_size is not None:
                model_comp.nu_size=nu_size
            
            if nu_min is not None:
                model_comp.nu_min=nu_min
            
            if nu_max is not None:
                model_comp.nu_max=nu_max

    def set(self,model,par_name, *args, **kw):
        self.parameters.set(model, par_name, *args, **kw)

    def set_par(self,model,par_name,val):
        self.parameters.set(model, par_name, val=val)

    def get(self,model,par_name,*args):
        self.parameters.get(model,par_name,*args)

    def get_par_by_name(self,model,par_name):
        return self.parameters.get_par_by_name(model, par_name)

    def freeze(self,model,par_name):
        self.parameters.freeze(model,par_name)

    def free(self,model,par_name):
        self.parameters.free(model,par_name)

    def free_all(self,):
        self.parameters.free_all()

    def freeze_all(self,):
        self.parameters.freeze_all()

    def add_component(self,m):
        self.components.add_component(m,self)

    def del_component(self,m):
        self.components.del_component(m,self)

    @property
    def composite_expr(self):
        return self._composite_expr

    def link_par(self,par_name,from_model,to_model):
        if isinstance(from_model, list) is False:
            from_model = [from_model]
        self.parameters.link_par(par_name,from_model, to_model)

    @composite_expr.setter
    def composite_expr(self,expr_string):
        if expr_string is None:
            self._composite_expr = expr_string
        else:
            try:
                for key, val in self.components._components_value_dict.items():
                    exec(key + '= 1.0')
                eval(expr_string)
            except Exception as e:
                raise RuntimeError('function string not valid',e)

        self._composite_expr=expr_string

    def _eval_composite_func(self,loglog):
        #transform each key into a local var
        for key, val in self.components._components_value_dict.items():
            if loglog is True:
                exec(key + '=np.power(10.,**val)')
            else:
                exec(key + '=val')

        return eval(self.composite_expr)

    def _eval_model(self, lin_nu, log_nu, loglog,fill_SED):
        lin_model = np.zeros(lin_nu.size)
        log_model = None
        for model_comp in self.components._components_list:
            model_comp.cosmo=self.cosmo
            #print('--> eval composite component', model_comp.name)
            if loglog is False:
                self.components._components_value_dict[model_comp.name]=model_comp.eval(nu=lin_nu, fill_SED=fill_SED, get_model=True, loglog=loglog)
            else:
                self.components._components_value_dict[model_comp.name] = model_comp.eval(nu=log_nu, fill_SED=fill_SED,get_model=True, loglog=loglog)
            if self.composite_expr is None:
                if loglog is False:
                    lin_model += self.components._components_value_dict[model_comp.name]
                else:
                    lin_model += np.power(10.,self.components._components_value_dict[model_comp.name])

        if self.composite_expr is not None:
            lin_model = self._eval_composite_func(loglog)

        #lin_model[lin_model< self.flux_plot_lim]=self.flux_plot_lim
        if loglog is True:
            log_model=np.log10(lin_model)

        return lin_model, log_model

    def eval(self,nu=None,fill_SED=True,get_model=False,loglog=False,label=None,phys_output=False):

        out_model= None
        #print('--> model mananger eval 1')
        lin_nu, log_nu = self._prepare_nu_model(nu, loglog)
        #print('--> model mananger eval 2',lin_nu[0],log_nu[0] )
        lin_model,log_model  = self._eval_model(lin_nu, log_nu, loglog, fill_SED)

        if fill_SED is True:
            self._fill(lin_nu, lin_model)

        if get_model is True:
            if loglog is True:
                out_model = log_model
            else:
                out_model = lin_model

        return out_model

    @classmethod
    def load_model(cls, file_name):
         c = pickle.load(open(file_name, "rb"))
         return cls._build_model(c)
    
    @staticmethod
    def _build_model(c):
        try:
            ml=c.components.components_list[::]
            for m in ml:
                try:
                    c.del_component(m.name)
                except:
                    pass
            for m in ml:
                c.add_component(m)

            for p in c.parameters.par_array:
                if p._linked is True:
                    p._linked = False
                    p._is_dependent = False
                    #print(p.name,p._root_par,[p.model],p._linked_root_model,p.immutable)
                    c.parameters.link_par(p._root_par.name,[p.model.name],p._linked_root_model.name)

            for m in c.components.components_list:
                if isinstance(m,Jet):
                    m._fix_par_dep_on_load()
            if isinstance(c, Model):
                c.eval()
                return c
            else:
                raise RuntimeError('The model you loaded is not valid please check the file name')

        except Exception as e:
            raise RuntimeError(e)

    def clone(self):
        return self._build_model(pickle.loads(pickle.dumps(self, protocol=pickle.HIGHEST_PROTOCOL)))

    def show_model_components(self):
        print("")
        print('-'*80)

        print("Composite model description")
        print('-'*80)
        print("name: %s  " % (self.name))
        print("type: %s  " % (self.model_type))
        print("components models:")
        for m in self.components._components_list:
            print(' -model name:', m.name, 'model type:', m.model_type)
        print('')
        print('-'*80)

    def show_model(self):
        print("")
        print('-'*80)

        print("Composite model description")
        print('-'*80)
        print("name: %s  " % (self.name))
        print("type: %s  " % (self.model_type))
        print("components models:")
        for m in self.components._components_list:
            print(' -model name:',m.name,'model type:', m.model_type)
        print('')
        print('-'*80)

        print("individual component description")

        self.components.show_model()
        print('-'*80)

