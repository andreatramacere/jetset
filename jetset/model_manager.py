"""Composite model containers and orchestration helpers for fitting."""


__author__ = "Andrea Tramacere"


import warnings

import numpy as np

#from  . import minimizer



from .model_parameters import  CompositeModelParameterArray

from .spectral_shapes import  SED
   
from .base_model import  Model

#from .plot_sedfit import  PlotSED

from .utils import  clean_var_name

from .jet_model import Jet

from .cosmo_tools import  Cosmo

#import  dill as pickle

__all__=['FitModel']

class CompositeModelContainer(object):

    """Container that manages component models inside a composite fit model.

    Notes
    -----
    Tracks component registration order, merged parameter arrays, and cached
    per-component evaluated values used during composite model evaluation.
    """
    def __init__(self):
        """Initialize container state for model components.

        Notes
        -----
        Tracks component objects, evaluated values, and merged parameters.
        """
        self._components_list=[]
        self._components_value=[]
        self._components_value_dict = {}
        self.parameters=CompositeModelParameterArray()

    def add_component(self, model_comp,fit_model):
        """Add component.
        
        Parameters
        ----------
        model_comp : object
            Component model instance to add/remove/query.
        fit_model : object
            Model instance used for fitting.
        """
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
        """Components list.
        
        Returns
        -------
        object
            Requested value.
        """
        return self._components_list

    def get_model_by_name(self, model_name,get_idx=False):
        """Return model by name.
        
        Parameters
        ----------
        model_name : object
            Name of the model/component.
        get_idx : bool, optional
            If ``True``, also return the component index.
        
        Returns
        -------
        object
            Requested value.
        """
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
        """Del component.
        
        Parameters
        ----------
        model_name : object
            Name of the model/component.
        fit_model : object
            Model instance used for fitting.
        """
        m,ID=self.get_model_by_name(model_name,get_idx=True)
        if m is not None:
            _ = self._components_list.pop(ID)
            _ = self._components_value.pop(ID)
            _ = self.parameters.del_model_parameters(m)
            del self._components_value_dict[m.name]

        delattr(self,model_name)
        delattr(fit_model, model_name)

    def show_pars(self):
        """Display parameters for all registered components.

        Notes
        -----
        Delegates formatting to the internal parameter-array helper.
        """
        self.parameters.show_pars()

    def show_model(self):
        """Display each registered component summary.

        Notes
        -----
        Calls ``show_model`` on each component in insertion order.
        """
        for c in self._components_list:
            #print()
            c.show_model()
            #print()





class FitModel(Model):
    """Composite spectral model used for fitting observational datasets.

    Notes
    -----
    Combines one or more component models (for example, jets and templates),
    exposes a unified parameter interface, and evaluates either summed
    components or user-defined composite expressions.
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

        
        """Create a new `FitModel` instance.
        
        Parameters
        ----------
        elec_distr : object, optional
            Electron-distribution object used to initialize a jet.
        jet : object, optional
            Jet model instance.
        name : str, optional
            Name identifier.
        out_dir : object, optional
            Output directory path.
        flag : object, optional
            Optional flag/name suffix used in object naming.
        template : object, optional
            Template-model component instance.
        loglog_poly : object, optional
            Log-log polynomial component instance.
        analytical : object, optional
            Analytical component instance.
        nu_size : int, optional
            Number of points for frequency grids.
        cosmo : object, optional
            Cosmology helper used for frame/luminosity conversions.
        composite_expr : object, optional
            Composite expression used to combine component outputs.
        **keywords : dict
            Additional keyword arguments.
        """
        super(FitModel,self).__init__(model_type='composite_model',
                                      nu_size=nu_size,
                                      name=name,
                                      cosmo=cosmo,
                                       **keywords)
       
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
        self.spectral_components_table_list=[]
        
        if elec_distr is not None:
            jet=Jet(cosmo=cosmo,name=flag, electron_distribution=elec_distr, jet_workplace=None)
            self.add_component(jet)
        
        if jet is not None:
            self.add_component(jet)

        if cosmo is None:
            if jet is not None:
                  self.cosmo=jet.cosmo
                  m='no cosmology defined, using the one from jet %s'%str(jet.cosmo)
            else:
                 self.cosmo = Cosmo()
                 m='no cosmology defined, using %s'%str(self.cosmo )
            warnings.warn(m)
        else:
            self.cosmo = cosmo

        if template is not None:
            self.add_component(template)

        if loglog_poly is not None:
            self.add_component(loglog_poly)

        if analytical is not None:
            self.add_component(analytical)

    def plot_model(self,plot_obj=None,clean=False,sed_data=None,frame='obs',skip_components=False,label=None,skip_sub_components=False, only_components=False, density=False):
        """Plot the composite model, optional components, and residuals.

        Parameters
        ----------
        plot_obj : PlotSED, optional
            Existing plot object to update. If ``None``, a new one is created.
        clean : bool, optional
            If ``True``, clear previously drawn model lines before plotting.
        sed_data : ObsData, optional
            Observational SED data. Used to initialize the plot helper and to
            compute residuals. If ``None``, residual points are not added.
        frame : {'obs', 'src'}, optional
            Frame used to display model curves.
        skip_components : bool, optional
            If ``True``, do not plot top-level component models.
        label : str, optional
            Label for the total composite model curve. Defaults to ``self.name``.
        skip_sub_components : bool, optional
            If ``True``, do not plot sub-components for each component model.
        only_components : bool, optional
            If ``True``, skip plotting the total composite curve (``self.SED``)
            and plot only component/sub-component curves (subject to skip flags).
        density : bool, optional
            Passed to the plotting helper when creating a new plot object.

        Returns
        -------
        PlotSED
            Plot object with model curves (and residuals when ``sed_data`` is
            provided).
        """
        plot_obj=self._set_up_plot(plot_obj,sed_data,frame,density)

        if clean is True:
            plot_obj.clean_model_lines()

        if skip_components is False:
            line_style = '--'

            for mc in self.components._components_list:
              
                comp_label = mc.name
                if hasattr(mc,'SED'):
                    try:
                        plot_obj.add_model_plot(mc.SED, line_style=line_style,label=comp_label,flim=self.flux_plot_lim, frame=frame)
                    except Exception as e:
                        try:
                            mc.eval()
                            plot_obj.add_model_plot(mc.SED, line_style=line_style, label=comp_label,
                                                    flim=self.flux_plot_lim, frame=frame)
                        except Exception as e:
                            raise RuntimeError('for model', mc.name, "problem with plotting SED", e)


                if skip_sub_components is False:
                    if hasattr(mc,'spectral_components_list'):
                        for c in mc.spectral_components_list:
                            comp_label = c.name
                            if comp_label!='Sum':
                                if hasattr(c, 'SED'):
                                    try:
                                        plot_obj.add_model_plot(c.SED, line_style=line_style, label='  -%s'%comp_label, flim=self.flux_plot_lim, frame=frame)
                                    except Exception as e:
                                        try:
                                            #print('==> reval',mc.name)
                                            mc.eval()
                                            _c=mc.spectral_components.get_spectral_component_by_name(c.name)
                                            plot_obj.add_model_plot(_c.SED, line_style=line_style, label='  -%s'%comp_label, flim=self.flux_plot_lim, frame=frame)
                                        except Exception as e:
                                            raise RuntimeError('for model', mc.name, "spectral component",c.name, "problem with plotting SED", e)
        line_style = '-'
        if label is None:
            label=self.name

        if not only_components:
            plot_obj.add_model_plot(self.SED, line_style=line_style, label=label, flim=self.flux_plot_lim,fit_range=[self.nu_min_fit,self.nu_max_fit], frame=frame  )
        plot_obj.add_model_residual_plot(data=sed_data, model=self, fit_range=[self.nu_min_fit, self.nu_max_fit])

        #if frame == 'src' and sed_data is not None:
        #    sed_data.z = z_sed_data

        return plot_obj



    def set_nu_grid(self,nu_min=None,nu_max=None,nu_size=None):
        """Set nu grid.
        
        Parameters
        ----------
        nu_min : object, optional
            Minimum frequency in Hz.
        nu_max : object, optional
            Maximum frequency in Hz.
        nu_size : object, optional
            Number of points for frequency grids.
        """
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
        """Set.
        
        Parameters
        ----------
        model : object
            Model instance.
        par_name : object
            Parameter name.
        *args : tuple
            Additional positional arguments.
        **kw : dict
            Additional keyword-value mapping.
        """
        self.parameters.set(model, par_name, *args, **kw)

    def set_par(self,model,par_name,val):
        """Set par.
        
        Parameters
        ----------
        model : object
            Model instance.
        par_name : object
            Parameter name.
        val : object
            Value to assign.
        """
        self.parameters.set(model, par_name, val=val)

    def get(self,model,par_name,*args):
        """Get.
        
        Parameters
        ----------
        model : object
            Model instance.
        par_name : object
            Parameter name.
        *args : tuple
            Additional positional arguments.
        """
        self.parameters.get(model,par_name,*args)

    def get_par_by_name(self,model,par_name):
        """Return par by name.
        
        Parameters
        ----------
        model : object
            Model instance.
        par_name : object
            Parameter name.
        
        Returns
        -------
        object
            Requested value.
        """
        return self.parameters.get_par_by_name(model, par_name)

    def freeze(self,model,par_name):
        """Freeze.
        
        Parameters
        ----------
        model : object
            Model instance.
        par_name : object
            Parameter name.
        """
        self.parameters.freeze(model,par_name)

    def free(self,model,par_name):
        """Free.
        
        Parameters
        ----------
        model : object
            Model instance.
        par_name : object
            Parameter name.
        """
        self.parameters.free(model,par_name)

    def free_all(self,):
        """Set all model parameters to free state.

        Notes
        -----
        Delegates to ``CompositeModelParameterArray.free_all``.
        """
        self.parameters.free_all()

    def freeze_all(self,):
        """Freeze all model parameters.

        Notes
        -----
        Delegates to ``CompositeModelParameterArray.freeze_all``.
        """
        self.parameters.freeze_all()

    def add_component(self,m):
        """Add component.
        
        Parameters
        ----------
        m : object
            Model/component object.
        """
        self.components.add_component(m,self)

    def del_component(self,m):
        """Del component.
        
        Parameters
        ----------
        m : object
            Model/component object.
        """
        self.components.del_component(m,self)

    @property
    def composite_expr(self):
        """Composite expr.
        
        Returns
        -------
        object
            Requested value.
        """
        return self._composite_expr

    def link_par(self,par_name,from_model,to_model):
        """Link par.
        
        Parameters
        ----------
        par_name : object
            Parameter name.
        from_model : object
            Source model/component for parameter linking.
        to_model : object
            Destination model/component for parameter linking.
        """
        if isinstance(from_model, list) is False:
            from_model = [from_model]
        self.parameters.link_par(par_name,from_model, to_model)

    @composite_expr.setter
    def composite_expr(self,expr_string):
        """Composite expr.
        
        Parameters
        ----------
        expr_string : object
            String expression for composite-model evaluation.
        """
        if expr_string is None:
            self._composite_expr = expr_string
        else:
            try:
                _components_namespace = {key: 1.0 for key in self.components._components_value_dict}
                eval(expr_string, {"np": np, "__builtins__": __builtins__}, _components_namespace)
            except Exception as e:
                raise RuntimeError('function string not valid',e)

        self._composite_expr=expr_string

    def _eval_composite_func(self,loglog):
        #transform each key into a local var
        _components_namespace = {}
        for key, val in self.components._components_value_dict.items():
            if loglog is True:
                _components_namespace[key] = np.power(10., val)
            else:
                _components_namespace[key] = val

        return eval(self.composite_expr, {"np": np, "__builtins__": __builtins__}, _components_namespace)

    def _eval_model(self, lin_nu, log_nu, loglog,fill_SED):
        lin_model = np.zeros(lin_nu.shape)
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

        """Evaluate model output.
        
        Parameters
        ----------
        nu : object, optional
            Frequency values in Hz.
        fill_SED : bool, optional
            If ``True``, store evaluated values into SED containers.
        get_model : bool, optional
            If ``True``, return model values.
        loglog : bool, optional
            If ``True``, operate in log10 space.
        label : object, optional
            Label used in output or plots.
        phys_output : bool, optional
            If ``True``, return physical-output quantities when available.
        
        Returns
        -------
        object
            Computed value.
        """
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
    def load_model(cls, file_name_or_obj, from_string=False):
         """Load object state from disk.
         
         Parameters
         ----------
         file_name_or_obj : object
             Serialized model path or already-open object.
         from_string : bool, optional
             If ``True``, deserialize from an in-memory string payload.
         
         Returns
         -------
         object
             Loaded object.
         """
         c = cls._load_pickle(file_name_or_obj,from_string=from_string)
         return cls._build_model(c)
    
    @staticmethod
    def _build_model(c):
        try:
            ml=c.components.components_list[::]
            for m in ml:
                try:
                    #print("===> m.name",m.name)
                    c.del_component(m.name)
                except Exception as e:
                    raise RuntimeError('for model',m.name,e)
            for m in ml:
                c.add_component(m)

            for p in c.parameters.par_array:
                if p._linked is True:
                    p._linked = False
                    p._is_dependent = False
                    #print(p.name,p._root_par,[p.model],p._linked_root_model,p.immutable)
                    c.parameters.link_par(p._root_par.name,[p.model.name],p._linked_root_model.name)

            #for m in c.components.components_list:
            #    if isinstance(m,Jet):
            #        m._fix_par_dep_on_load()
            if isinstance(c, Model):
                c.eval()
                return c
            else:
                raise RuntimeError('The model you loaded is not valid please check the file name')

        except Exception as e:
            raise RuntimeError(e)


    def set_fit_range(self,down_tol=0.1,up_tol=100):
        """Set fit range.
        
        Parameters
        ----------
        down_tol : float, optional
            Lower tolerance for parameter scans/intervals.
        up_tol : int, optional
            Upper tolerance for parameter scans/intervals.
        """
        for m in self.components.components_list:
            m.set_fit_range(down_tol=down_tol,up_tol=up_tol)
                

    #def clone(self):
    #    return self.load_model(pickle.loads(pickle.dumps(self, protocol=pickle.HIGHEST_PROTOCOL)))

    def show_model_components(self):
        """Print a concise overview of component models.

        Notes
        -----
        This view does not print individual parameter tables.
        """
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
        """Print full composite-model summary and component details.

        Notes
        -----
        Includes the per-component ``show_model`` output.
        """
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

    def sed_tables_dict(self, restframe='obs'):
        """Sed tables dict.
        
        Parameters
        ----------
        restframe : str, optional
            Target rest frame for output.
        
        Returns
        -------
        object
            Computed value.
        """
        self.eval()
        self._sed_tables_dict={}
        
        self._sed_tables_dict[self.name]=self.sed_table(restframe=restframe)
        for comp in self.components.components_list:
            if hasattr(comp,'sed_table'):
                if comp.sed_table(restframe=restframe) is not None:
                    self._sed_tables_dict[comp.name]=comp.sed_table(restframe=restframe)
        
        return self._sed_tables_dict
