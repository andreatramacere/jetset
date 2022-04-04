
__author__ = "Andrea Tramacere"


import numpy as np
import json
import dill as pickle
import warnings
import inspect

from .model_parameters import ModelParameterArray, ModelParameter
from .spectral_shapes import SED
from .data_loader import  ObsData
from .utils import  get_info
from .plot_sedfit import  PlotSED


from .cosmo_tools import  Cosmo

__all__=['Model','MultiplicativeModel']




class Model(object):
    
    
    def __init__(self,name='no-name',nu_size=200,model_type='base_model',scale='lin-lin',cosmo=None,nu_min=None,nu_max=None):
        
        self.model_type=model_type
        
        self.name=name

        self.SED = SED(name=self.name)
    
        self.parameters = ModelParameterArray()    
        
        self._scale=scale
        
        self.nu_size=nu_size
        
        self.nu_min=nu_min
        
        self.nu_max=nu_max

        self.flux_plot_lim = 1E-30

        if cosmo is None:
          self.cosmo=Cosmo()

        self._set_version(v=None)
        
    @property
    def version(self):
        return self._version

    def _set_version(self, v=None):
        if v is None:
            self._version = get_info()['version']
        else:
            self._version = v

    def _prepare_nu_model(self,nu,loglog):

        if nu is None:
            x1 = np.log10(self.nu_min)
            x2 = np.log10(self.nu_max)
            lin_nu = np.logspace(x1, x2, self.nu_size)
            log_nu = np.log10(lin_nu)
        else:

            if np.shape(nu) == ():
                nu = np.array([nu])

            if loglog is True:
                lin_nu = np.power(10., nu)
                log_nu = nu
            else:
                log_nu = np.log10(nu)
                lin_nu = nu

        return lin_nu,log_nu

    def _eval_model(self,lin_nu,log_nu,loglog):
        log_model=None

        if loglog is False:
            lin_model = self.lin_func(lin_nu)
        else:
            if hasattr(self, 'log_func'):
                log_model = self.log_func(log_nu)
                lin_model = np.power(10., log_model)
            else:
                lin_model = self.lin_func(lin_nu)
                lin_model[lin_model<0.]=self.flux_plot_lim
                log_model = np.log10(lin_model)

        return lin_model,log_model

    def _fill(self, lin_nu, lin_model):

        if hasattr(self,'SED'):
            self.SED.fill(nu=lin_nu, nuFnu=lin_model)
            z=self.get_par_by_type('redshift')

            if z is not None:
                z=z.val

            if z is None and hasattr(self, 'get_redshift'):
                z=self.get_redshift()
                z = z

            #print('--> fill z ',self.name,self.name,z)
            if z is not None:
                if hasattr(self,'get_DL_cm'):
                    dl = self.get_DL_cm('redshift')
                else:
                    dl = self.cosmo.get_DL_cm(z)
                self.SED.fill_nuLnu( z =z, dl = dl)
        else:
            warnings.warn('model',self.name,'of type',type(self),'has no SED member')

        if hasattr(self, 'spectral_components_list'):
            for i in range(len(self.spectral_components_list)):
                self.spectral_components_list[i].fill_SED(lin_nu=lin_nu,skip_zeros=False)

    def eval(self, fill_SED=True, nu=None, get_model=False, loglog=False, label=None, **kwargs):

        out_model = None
        #print('--> base model 1', nu[0])
        lin_nu,log_nu=self._prepare_nu_model(nu,loglog)
        #print('--> base model 2', lin_nu[0], log_nu[0],'loglog',loglog)
        lin_model,log_model  = self._eval_model(lin_nu,log_nu,loglog)

        if fill_SED is True:
            self._fill(lin_nu, lin_model)

        if get_model is True:

            if loglog is True:
                out_model = log_model
            else:
                out_model = lin_model

        return out_model


    def _set_up_plot(self,plot_obj,sed_data,frame,density):
        if plot_obj is None:
            plot_obj=PlotSED( frame = frame, density=density,sed_data=sed_data)

        z_sed_data = None
        if frame == 'src' and sed_data is not None:
            z_sed_data = sed_data.z
            if self.get_par_by_type('redshift') is not None:
                sed_data.z = self.get_par_by_type('redshift').val

        #if sed_data is not None:
        #    plot_obj.add_data_plot(sed_data)


        if frame == 'src' and z_sed_data is not None:
            sed_data.z = z_sed_data

        return plot_obj


    def plot_model(self,plot_obj=None,clean=False,sed_data=None,frame='obs',skip_components=False,label=None,line_style='-', density=False):

        plot_obj=self._set_up_plot(plot_obj,sed_data,frame,density)

        if clean is True:
            plot_obj.clean_model_lines()

        if label is None:
            label = self.name

        if hasattr(self,'SED'):
            plot_obj.add_model_plot(self.SED, line_style=line_style,label =label,flim=self.flux_plot_lim,density=density, frame=frame)



        if skip_components is False:
            if hasattr(self,'spectral_components_list'):
                for c in self.spectral_components_list:
                    #print('--> c name', c.name)
                    comp_label = c.name
                    line_style = '--'
                    if comp_label!='Sum':
                        if hasattr(c, 'SED'):
                            plot_obj.add_model_plot(c.SED, line_style=line_style, label='  -%s'%comp_label, flim=self.flux_plot_lim, density=density, frame=frame)

        line_style = '-'


        #plot_obj.add_model_residual_plot(data=sed_data, model=self,fit_range=np.log10([self.nu_min_fit,self.nu_max_fit]) )

        return plot_obj

    def set_nu_grid(self,nu_min=None,nu_max=None,nu_size=None):
        if nu_size is not None:
            self.nu_size=nu_size
        
        if nu_min is not None:
            self.nu_min=nu_min
        
        if nu_max is not None:
            self.nu_max=nu_max


    def lin_func(self,lin_nu):
        return np.ones(lin_nu.size) * self.flux_plot_lim
    
    def log_func(self,log_nu):
        return np.log10(self.lin_func(np.power(10,log_nu)))



    def get_residuals(self, data, log_log=False,filter_UL=True):
        if isinstance(data,ObsData):
            data=data.data

        model = self.eval(nu=data['nu_data'], fill_SED=False, get_model=True, loglog=False)

        if filter_UL ==True:
            msk=data['UL']==False
        else:
            msk=np.ones(data.size,dt=True)

        residuals = (data['nuFnu_data'] - model) /  data['dnuFnu_data']

        nu_residuals=data['nu_data']


        if log_log == False:
            return nu_residuals[msk], residuals[msk]
        else:
            return  np.log10(nu_residuals[msk]),  residuals[msk]

    def save_model(self, file_name):

        pickle.dump(self, open(file_name, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)


    @classmethod
    def load_model(cls, file_name):
        try:
            c = pickle.load(open(file_name, "rb"))
            c._fix_par_dep_on_load()
            if isinstance(c, Model):
                c.eval()
                return c
            else:
                raise RuntimeError('The model you loaded is not valid please check the file name')

        except Exception as e:
            raise RuntimeError(e)

    def _fix_par_dep_on_load(self,):
        for p in self.parameters.par_array:
            #print("==> 1", p.name, p._master_par_list, p._depending_par_expr,p.model.name)
            if p._is_dependent is True and p._linked is False:
                self._is_dependent = False
                #print("==> 2", p.name, p._master_par_list, p._depending_par_expr)
                self.make_dependent_par(p.name, p._master_par_list, p._depending_par_expr)

    #def _set_pars_dep(self):
    #    for p in self.parameters.par_array:
    #        if

    def clone(self):
        return  pickle.loads(pickle.dumps(self))

    def show_model(self):
        print("")
        print('-'*80)

        print("model description")
        print('-' * 80)
        print("name: %s  " % (self.name))
        print("type: %s  " % (self.model_type))

        print('')

        print('-'*80)
        self.parameters.show_pars()

        print('-'*80)

    def show_pars(self, sort_key='par type'):
        return self.parameters.show_pars(sort_key=sort_key)

    def show_best_fit_pars(self):
        self.parameters.show_best_fit_pars()

    def set_par(self,par_name,val):
        """
        shortcut to :class:`ModelParametersArray.set` method
        set a parameter value

        :param par_name: (srt), name of the parameter
        :param val: parameter value

        """

        self.parameters.set(par_name, val=val)



    def get_par_by_type(self,par_type):
        """

        get parameter by type

        """
        for param in self.parameters.par_array:
            if param.par_type==par_type:
                return param

        return None

    def get_par_by_name(self,par_name):
        """

        get parameter by type

        """
        for param in self.parameters.par_array:
            if param.name==par_name:
                return param

        return None

    def dep_func_get_default_args(self, par_func):
        signature = inspect.signature(par_func)
        d = []
        for k, v in signature.parameters.items():
            #print('==> par',k,v)
            p = self.get_par_by_name(k)
            if p is not None:
                d.append(k)
            else:
                raise RuntimeError('argument', k, 'is not valid, should be a model parameter name')
        return d

    def _test_par_expr(self,master_par_list,par_expr):
        if type(par_expr) == str:
            for p_name in master_par_list:
                exec(p_name + '= 1')
            try:
                eval(par_expr)
                pass
            except:
                raise RuntimeError('the parameter expression is not valid')

    def make_dependent_par(self, par, depends_on, par_expr):
        master_par_list = depends_on

        dep_par=self.parameters.get_par_by_name(par)

        if dep_par.name in master_par_list:
            raise RuntimeError("depending parameter:", dep_par.name, "can't be in master par list",master_par_list)

        self._test_par_expr(master_par_list,par_expr)

        dep_par.freeze()
        dep_par._is_dependent = True
        dep_par.par_expr = par_expr
        dep_par._func = dep_par._eval_par_func
        dep_par._master_par_list=master_par_list
        for p in master_par_list:
            m = self.parameters.get_par_by_name(p)
            dep_par._add_master_par(m)
            m._add_depending_par(dep_par)

        for p in master_par_list:

            m = self.parameters.get_par_by_name(p)
            if m._is_dependent is False:
                m.val=m.val
        if isinstance(par_expr,str):
            _par_expr=par_expr
        else:
            _par_expr=inspect.getsource(par_expr)
        print('==> par', dep_par.name, 'is now depending on', master_par_list, f'according to expr:{dep_par.name} =\n{_par_expr}'.format(dep_par.name,_par_expr))

    def add_user_par(self,name,val,units='',val_min=None,val_max=None):
        self.parameters.add_par(ModelParameter(name=name,units=units,val=val,val_min=val_min,val_max=val_max,par_type='user_defined'))


class MultiplicativeModel(Model):

    def __init__(self, name='no-name', nu_size=100, model_type='multiplicative_model', scale='lin-lin'):
        super(MultiplicativeModel, self).__init__(name=name, nu_size=nu_size, model_type=model_type,scale=scale)
        delattr(self,'SED')


