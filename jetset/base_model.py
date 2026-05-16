"""Base model abstractions and shared utilities used across JetSeT models."""


__author__ = "Andrea Tramacere"


import numpy as np
import copy
import dill as pickle
import warnings
import inspect
import numbers
from astropy.table import Table
import warnings

from .model_parameters import ModelParameterArray, ModelParameter
from .spectral_shapes import SED
from .data_loader import  ObsData
from .utils import  get_info
from .plot_sedfit import  PlotSED
from .utils import check_frame,unexpected_behaviour

from .cosmo_tools import  Cosmo

__all__=['Model','MultiplicativeModel']




class Model(object):
    
    
    """Base class for analytical and numerical SED models."""
    def __init__(self,name='no_name',
                 nu_size=200,
                 model_type='base_model',
                 scale='lin-lin',
                 cosmo=None,
                 nu_min=None,
                 nu_max=None):
        
        """Initialize a model container.

        Parameters
        ----------
        name : str, optional
            Model name.
        nu_size : int, optional
            Number of frequencies used when building an internal evaluation grid.
        model_type : str, optional
            Descriptive model type label.
        scale : str, optional
            Preferred plotting/evaluation scale.
        cosmo : Cosmo, optional
            Cosmology helper. If omitted, a default :class:`Cosmo` instance is used.
        nu_min : float, optional
            Minimum frequency of the model grid in Hz.
        nu_max : float, optional
            Maximum frequency of the model grid in Hz.
        """
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
        else:
            self.cosmo=cosmo

        self._set_version(v=None)
        
    @property
    def version(self):
        """Return the package version used to create this model.

        Returns
        -------
        str
            Version string.
        """
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

        """Evaluate model fluxes.

        Parameters
        ----------
        fill_SED : bool, optional
            If ``True``, update the model SED object after evaluation.
        nu : array-like or float, optional
            Frequency grid in Hz (or log10(Hz) if ``loglog`` is ``True``).
            If omitted, the internal ``nu_min``/``nu_max``/``nu_size`` grid is used.
        get_model : bool, optional
            If ``True``, return evaluated values.
        loglog : bool, optional
            If ``True``, treat input/output frequency and model values in log10 space.
        label : str, optional
            Reserved for subclasses/plotting integrations.
        **kwargs
            Extra keyword arguments accepted for subclass compatibility.

        Returns
        -------
        ndarray or None
            Evaluated model array when ``get_model`` is ``True``, otherwise ``None``.
        """
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


        if frame == 'src' and z_sed_data is not None:
            sed_data.z = z_sed_data

        return plot_obj


    def plot_model(self,plot_obj=None,clean=False,sed_data=None,frame='obs',skip_components=False,label=None,line_style='-', density=False):

        """Plot model SED and optional spectral components.

        Parameters
        ----------
        plot_obj : PlotSED, optional
            Existing plotting object. If omitted, a new one is created.
        clean : bool, optional
            If ``True``, clear existing model lines from ``plot_obj``.
        sed_data : ObsData, optional
            Optional observed data used by the plotting helper.
        frame : {'obs', 'src'}, optional
            Output frame for plotting.
        skip_components : bool, optional
            If ``True``, do not plot individual spectral components.
        label : str, optional
            Label for the model curve.
        line_style : str, optional
            Matplotlib line style.
        density : bool, optional
            If ``True``, use density representation in plotting helper.

        Returns
        -------
        PlotSED
            Plot object with model curves added.
        """
        plot_obj=self._set_up_plot(plot_obj,sed_data,frame,density)

        if clean is True:
            plot_obj.clean_model_lines()

        if label is None:
            label = self.name

        if hasattr(self,'SED'):
            plot_obj.add_model_plot(self.SED, line_style=line_style,label =label,flim=self.flux_plot_lim, frame=frame)



        if skip_components is False:
            if hasattr(self,'spectral_components_list'):
                for c in self.spectral_components_list:
                    #print('--> c name', c.name)
                    comp_label = c.name
                    line_style = '--'
                    if comp_label!='Sum':
                        if hasattr(c, 'SED'):
                            plot_obj.add_model_plot(c.SED, line_style=line_style, label='  -%s'%comp_label, flim=self.flux_plot_lim, frame=frame)

        line_style = '-'


        return plot_obj

    def set_nu_grid(self,nu_min=None,nu_max=None,nu_size=None):
        """Set model frequency-grid settings.

        Parameters
        ----------
        nu_min : float, optional
            Minimum frequency in Hz.
        nu_max : float, optional
            Maximum frequency in Hz.
        nu_size : int, optional
            Number of samples in the grid.
        """
        if nu_size is not None:
            self.nu_size=nu_size
        
        if nu_min is not None:
            self.nu_min=nu_min
        
        if nu_max is not None:
            self.nu_max=nu_max


    def lin_func(self,lin_nu):
        """Return model values in linear space.

        Parameters
        ----------
        lin_nu : ndarray
            Frequency array in Hz.

        Returns
        -------
        ndarray
            Model values in linear ``nuFnu`` units.
        """
        return np.ones(lin_nu.shape) * self.flux_plot_lim
    
    def log_func(self,log_nu):
        """Return model values in log10 space.

        Parameters
        ----------
        log_nu : ndarray
            Frequency array in log10(Hz).

        Returns
        -------
        ndarray
            Model values in log10 space.
        """
        return np.log10(self.lin_func(np.power(10,log_nu)))



    def get_residuals(self, data, log_log=False,filter_UL=True):
        """Compute residuals between observed data and model prediction.

        Parameters
        ----------
        data : ObsData or table-like
            Input data table containing ``nu_data``, ``nuFnu_data``,
            ``dnuFnu_data`` and ``UL`` columns.
        log_log : bool, optional
            If ``True``, return frequency axis in log10(Hz).
        filter_UL : bool, optional
            If ``True``, exclude upper-limit points.

        Returns
        -------
        tuple of ndarray
            ``(nu_axis, residuals)``.
        """
        if isinstance(data,ObsData):
            data=data.data

        model = self.eval(nu=data['nu_data'], fill_SED=False, get_model=True, loglog=False)

        if filter_UL ==True:
            msk=data['UL']==False
        else:
            msk=np.ones(len(data),dtype=bool)

        residuals = (data['nuFnu_data'] - model) /  data['dnuFnu_data']

        nu_residuals=data['nu_data']


        if log_log == False:
            return nu_residuals[msk], residuals[msk]
        else:
            return  np.log10(nu_residuals[msk]),  residuals[msk]

    def save_model(self, file_name):
        """Serialize the model to disk.

        Parameters
        ----------
        file_name : str
            Output pickle file path.
        """
        pickle.dump(self, open(file_name, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

    def save_model(self, file_name):
        """Serialize the model to disk.

        Parameters
        ----------
        file_name : str
            Output pickle file path.
        """
        pickle.dump(self, open(file_name, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)


    @classmethod
    def load_model(cls, file_name_or_obj,from_string=False):
        """Load a serialized model.

        Parameters
        ----------
        file_name_or_obj : str or bytes or file-like
            Serialized model source accepted by the pickle loader.
        from_string : bool, optional
            If ``True``, interpret ``file_name_or_obj`` as in-memory content.

        Returns
        -------
        Model
            Reconstructed and evaluated model instance.

        Raises
        ------
        RuntimeError
            If deserialization fails or object type is not valid.
        """
        try:
            c=cls._load_pickle(file_name_or_obj,from_string=from_string)
            c._fix_par_dep_on_load(verbose=True)
            if isinstance(c, Model):
                c.eval()
                return c
            else:
                raise RuntimeError('The model you loaded is not valid please check the file name')

        except Exception as e:
            raise RuntimeError(e)

    def _fix_par_dep_on_load(self,verbose=True):
        for p in self.parameters.par_array:
            _master_par_names_list=[]
            if hasattr(p,'_depending_pars'):
                p._depending_pars=[]
            if hasattr(p,'_master_pars'):
                for _p in p._master_pars:
                    if hasattr(_p,'name'):
                        if _p.name not in _master_par_names_list:
                            _master_par_names_list.append(_p.name)
                    else:
                        warnings.warn('found master par without name')
                p._master_pars=[]
                p._master_par_list=[]
                p._master_par_names_list=_master_par_names_list

        for p in self.parameters.par_array:
            if p._is_dependent is True and p._linked is False:
                _depending_par_expr=copy.deepcopy(p._depending_par_expr)
                self.make_dependent_par(p.name, depends_on=p._master_par_names_list, par_expr=_depending_par_expr,set_par_expr_source_code=True,verbose=verbose)
                if hasattr(p,'_master_par_names_list'):
                    delattr(p,'_master_par_names_list')

            
    @staticmethod
    def _load_pickle(file_name_or_obj,from_string=False):
        if from_string:
            c = pickle.loads(file_name_or_obj)
        else:
            c = pickle.load(open(file_name_or_obj, "rb"))
        return c
    
    def clone(self):
        """Return a deep clone of the model via in-memory serialization.

        Returns
        -------
        Model
            Cloned model instance.
        """
        return self.load_model(pickle.dumps(self, protocol=pickle.HIGHEST_PROTOCOL),from_string=True)

    def show_model(self):
        """Print model summary and parameters."""
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
        """Display model parameters.

        Parameters
        ----------
        sort_key : str, optional
            Column used to sort displayed parameters.

        Returns
        -------
        object
            Output returned by ``ModelParameterArray.show_pars``.
        """
        return self.parameters.show_pars(sort_key=sort_key)

    def show_best_fit_pars(self):
        """Display best-fit parameter values."""
        self.parameters.show_best_fit_pars()

    def set_par(self,par_name,val):
        """Set a parameter value by name.

        Parameters
        ----------
        par_name : str
            Parameter name.
        val : float or int
            New parameter value.
        """

        self.parameters.set(par_name, val=val)



    def get_par_by_type(self,par_type):
        """Return first parameter matching a parameter type.

        Parameters
        ----------
        par_type : str
            Parameter type label.

        Returns
        -------
        ModelParameter or None
            Matching parameter, if found.
        """
        for param in self.parameters.par_array:
            if param.par_type==par_type:
                return param

        return None

    def get_par_by_name(self,par_name):
        """Return parameter by name.

        Parameters
        ----------
        par_name : str
            Parameter name.

        Returns
        -------
        ModelParameter or None
            Matching parameter, if found.
        """
        for param in self.parameters.par_array:
            if param.name==par_name:
                return param

        return None

    def dep_func_get_default_args(self, par_func):
        """Validate dependency-function arguments against model parameters.

        Parameters
        ----------
        par_func : callable
            Function used as dependency expression.

        Returns
        -------
        list of str
            Ordered list of parameter names accepted by ``par_func``.

        Raises
        ------
        RuntimeError
            If ``par_func`` uses argument names not present in the model.
        """
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
            _par_values = {p_name: 1 for p_name in master_par_list}
            try:
                eval(par_expr, {"np": np, "__builtins__": __builtins__}, _par_values)
                pass
            except:
                raise RuntimeError('the parameter expression is not valid')

    def make_dependent_par(self, par, depends_on, par_expr,verbose=True,set_par_expr_source_code=True,master_pars=None):
        #print("\n  ===> make par: ",par, "depending on : ",depends_on, " START")
        """Make dependent par.
        
        Parameters
        ----------
        par : object
            Parameter object or parameter name.
        depends_on : object
            Names of master parameters used by the dependency.
        par_expr : object
            Expression defining a dependent parameter.
        verbose : bool, optional
            If ``True``, print additional information.
        set_par_expr_source_code : bool, optional
            If ``True``, store source code for the dependency expression.
        master_pars : object, optional
            Master-parameter objects used by the dependency.
        """
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
            try:
                m = self.parameters.get_par_by_name(p)
                dep_par._add_master_par(m,verbose=verbose)
                m._add_depending_par(dep_par)
            except Exception as e:
                message='problem with parameter name: %s'%p
                message+='\nexception:%s'%str(e)
                raise RuntimeError(message)
              
        for p in master_par_list:

            m = self.parameters.get_par_by_name(p)
            #trigger parameter update 
            if m._is_dependent is False:
                try:
                    m.val=m.val
                except:
                    pass
        if set_par_expr_source_code is True:
            dep_par._set_par_expr_source_code()
            if verbose is True:
                dep_par.par_expression_source_code
        
    
    def add_user_par(self,name,val,units='',val_min=None,val_max=None):
        """Add a user-defined parameter to the model.

        Parameters
        ----------
        name : str
            Parameter name.
        val : float
            Initial value.
        units : str, optional
            Parameter units label.
        val_min : float, optional
            Lower physical bound.
        val_max : float, optional
            Upper physical bound.
        """
        self.parameters.add_par(ModelParameter(name=name,units=units,val=val,val_min=val_min,val_max=val_max,par_type='user_defined'))



    def set_fit_range(self,down_tol=0.1,up_tol=100):
        """Set fit ranges for numeric parameters.

        Parameters
        ----------
        down_tol : float, optional
            Multiplicative factor for lower fit bound.
        up_tol : float, optional
            Relative factor for upper fit bound.
        """
        for p in self.parameters.par_array:
            if isinstance(p.val, numbers.Number):
                if p.val_min is not None:
                    p.fit_range_min=max(p.val_min,p.val_lin*down_tol)
                else:
                    p.fit_range_min = p.val_lin*down_tol
                if p.val_max is not None:
                    p.fit_range_max=min(p.val_max,p.val_lin *(1+up_tol))
                else:
                    p.fit_range_max = p.val_lin*(1+up_tol)


    def build_table(self, restframe='obs'):

        """Build SED table for the current model state.

        Parameters
        ----------
        restframe : {'obs', 'src'}, optional
            Frame used to build frequency/flux columns.
        """
        _names = ['nu']
        _cols=[]
        if hasattr(self,'SED'):
            check_frame(restframe)
            if restframe=='obs':
                 _cols.append(self.SED.nu)
            elif restframe=='src':
                _cols.append(self.SED.nu_src)
            else:
                unexpected_behaviour()

            
            if restframe == 'obs':
                _names.append('nuFnu')
                _cols.append(self.SED.nuFnu)
            else:
                _names.append('nuLnu_src')
                _cols.append(self.SED.nuLnu_src)
        
            _meta=dict(model_name=self.name)       
            _meta['restframe']= restframe
            self._SED_table = Table(_cols, names=_names,meta=_meta)
        else:
            self._SED_table = None

   
    def sed_table(self, restframe='obs'):
        """Return SED table, evaluating the model if needed.

        Parameters
        ----------
        restframe : {'obs', 'src'}, optional
            Frame used to build frequency/flux columns.

        Returns
        -------
        astropy.table.Table or None
            SED table for the current model state.
        """
        try:
            self.build_table(restframe=restframe)
        except:
            self.eval()
            self.build_table(restframe=restframe)
                 
        return  self._SED_table
   

class MultiplicativeModel(Model):

    """Model subclass for multiplicative spectral components."""
    def __init__(self, name='no_name', nu_size=100, model_type='multiplicative_model', scale='lin-lin'):
        """Initialize a multiplicative model.

        Parameters
        ----------
        name : str, optional
            Model name.
        nu_size : int, optional
            Number of evaluation frequencies.
        model_type : str, optional
            Model type label.
        scale : str, optional
            Preferred plotting/evaluation scale.
        """
        super(MultiplicativeModel, self).__init__(name=name, nu_size=nu_size, model_type=model_type,scale=scale)
        delattr(self,'SED')
