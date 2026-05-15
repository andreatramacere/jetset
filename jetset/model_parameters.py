__author__ = "Andrea Tramacere"

import  numpy as np
import warnings
from astropy.table import Table,vstack,MaskedColumn
from astropy import  units as u
from .utils import clean_var_name
from functools import wraps
import ast
import inspect


__all__=['ModelParameter','ModelParameterArray','Value']

class SettingDependentParError(Exception):
    """custom error for dependent error"""
    def __init__(self, message):
        """Create a new `SettingDependentParError` instance.
        
        Parameters
        ----------
        message : object
            Parameter controlling message.
        """
        super().__init__(message)
        

def is_notebook():
    """Is notebook.
    
    Returns
    -------
    object
        Computed result.
    """
    try:
        from IPython import get_ipython
        if "IPKernelApp" not in get_ipython().config:  # pragma: no cover
            raise ImportError("console")
            return False
    except:
        return False
    else:  # pragma: no cover
        return True

def _show_table(t):
    if is_notebook():
        try:
            from IPython.display import display
            display(t.show_in_notebook(show_row_index=False, display_length=100, backend='classic'))
            #display(t.show_in_notebook(show_row_index=False, display_length=100, auto_fit_columns=True  ))
        except Exception as e:
            try:
                from IPython.display import display
                display(t)
            except Exception as e:
                t.pprint_all()
    else:
        t.pprint_all()
class Value(object):
    """Numeric parameter value with unit and log/linear views.

    Notes
    -----
    Centralizes conversions between stored value, linear value, and log value,
    while tracking dimensionality and Astropy-compatible units.
    """

    def __init__(self,val,units,islog=False):
        """Create a new `Value` instance.
        
        Parameters
        ----------
        val : object
            Value to assign.
        units : object
            Parameter controlling units.
        islog : bool, optional
            Parameter controlling islog.
        """
        self.val=val
        self.islog=islog
        self.units=units
        self._adimensional=False

    def __repr__(self):
        return '%s'%self._val

    @property
    def islog(self):
        """Islog.
        
        Returns
        -------
        object
            Requested value.
        """
        return  self._islog

    @islog.setter
    def islog(self,val):
        """Islog.
        
        Parameters
        ----------
        val : object
            Value to assign.
        """
        self._islog=val

    @property
    def val(self):
        """Val.
        
        Returns
        -------
        object
            Requested value.
        """
        return self._val

    @val.setter
    def val(self, val):
        """Val.
        
        Parameters
        ----------
        val : object
            Value to assign.
        """
        self._val = val

    @property
    def lin(self):
        """Lin.
        
        Returns
        -------
        object
            Requested value.
        """
        if self._val is None:
            return None
        if self._islog is True:
            return 10**self._val
        else:
            return self._val

    @property
    def log(self):
        """Log.
        
        Returns
        -------
        object
            Requested value.
        """
        if self._val is None:
            return None
        if self._islog is True:
            return self._val
        else:
            return np.log10(self._val)

    @property
    def units(self):
        """Units.
        
        Returns
        -------
        object
            Requested value.
        """
        return self._units

    @units.setter
    def units(self, p_unit,verbose=False):
        """Units.
        
        Parameters
        ----------
        p_unit : object
            Parameter controlling p unit.
        verbose : bool, optional
            Parameter controlling verbose.
        """
        try:
            self._units = u.Unit(p_unit)
            if self._units == '':
                self._adimensional = True
            else:
                self._adimensional = False
            #print(units,type(self._units))
        except Exception as e:
            if verbose is True:
                print('units not valid for astropy units ',p_unit)

            if hasattr(self,'_units'):
                if '*' not in str(self._units):
                    self._units = p_unit + '*'
            else:
                if '*' not in p_unit:
                    self._units = p_unit + '*'
                else:
                    self._units = p_unit

            self._adimensional = True
            
            

class ModelParameter(object):
    """
    This class is the base class for models parameters. The following keywords 
    arguments can be passed to the constructor
    
    Parameters
    ----------      
    name          :  (str)    parameter name (default='unk')
    val           :  (float)  parameter value 
    par_type      :  (str)    parameter type (default='unk')
    units         :  (str)    units of the parameter, (default='No')
    val_min       :  (float)  minimum physical value 
    val_max       :  (float)  maximum physical value
    val_start     :  (float)  starting value
    val_last_call :  (float)  last call value
    fit_range_min :  (float)  minimum boundary value for the fit
    fit_range_max :  (float)  maximum boundary value for the fit  
    frozen        :  (bool)   boolean flag for frozen parameter (default=False)
    log           :  (bool)   boolean flag for log-scale value (default=False)
                    
    Attributes
    ----------
    name          :  (str)         parameter name (default='unk')                      
    val           :  (float)       parameter value                                    
    par_type      :  (str)         parameter type (default='unk')                     
    units         :  (str)         units of the parameter, (default='No')             
    val_min       :  (float)       minimum physical value                             
    val_max       :  (float)       maximum physical value                             
    val_start     :  (float)       starting value                                     
    val_last_call :  (float)       last call value                                    
    fit_range_min :  (float)       minimum boundary value for the fit                 
    fit_range_max :  (float)       maximum boundary value for the fit                 
    frozen        :  (bool)        boolean flag for frozen parameter (default=False)  
    log           :  (bool)        boolean flag for log-scale value (default=False)    
 
    """
    
    
    
    def __init__(self, **keywords):


        """Create a new `ModelParameter` instance.
        
        Parameters
        ----------
        **keywords : dict
            Parameter controlling keywords.
        """
        self.allowed_keywords={'name':None}
        self.allowed_keywords['val']=None
        self.allowed_keywords['par_type']=None
        self.allowed_keywords['allowed_par_types']=None
        self.allowed_keywords['units']='No'
        self.allowed_keywords['val_min']=None
        self.allowed_keywords['val_max']=None
        self.allowed_keywords['val_start']=None
        self.allowed_keywords['val_last_call']=None
        self.allowed_keywords['fit_range_min']=None
        self.allowed_keywords['fit_range_max']=None
        self.allowed_keywords['fit_range']=None
        self.allowed_keywords['best_fit_val']=None
        self.allowed_keywords['best_fit_err']=None
        self.allowed_keywords['frozen']=False
        self.allowed_keywords['log']=False
        self.allowed_keywords['allowed_values'] = None
        self.allowed_keywords['_is_dependent'] = False
        self.allowed_keywords['_depending_pars'] = []
        self.allowed_keywords['_func'] = None
        self.allowed_keywords['_master_pars'] = []
        self.allowed_keywords['_root_par'] = None
        self.allowed_keywords['_linked'] = False
        self.allowed_keywords['_linked_root_model'] = None
        self.allowed_keywords['_master_par_list'] = []
        self.allowed_keywords['_depending_par_expr'] = None
        self.allowed_keywords['_par_expr_text'] = None
        self._linked = False
        self._linked_root_model = None
        self._root_par = None
        _v = None
        _l = False
        _units= None
        self.par_type=None

        if '_is_dependent' in keywords.keys():
            self._is_dependent =  keywords['_is_dependent']
        else:
            self._is_dependent = self.allowed_keywords['_is_dependent']

        if '_depending_pars' in keywords.keys():
            self._depending_pars = keywords['_depending_pars']
        else:
            self._depending_pars = self.allowed_keywords['_depending_pars']

        if '_func' in keywords.keys():
            self._func = keywords['_func']
        else:
            self._func = self.allowed_keywords['_func']

        if '_master_pars' in keywords.keys():
            self._master_pars = keywords['_master_pars']
        else:
            self._master_pars = self.allowed_keywords['_master_pars']

        if 'val' in keywords.keys():
            _v = keywords['val']

        if 'log' in keywords.keys():
            _l = keywords['log']
        
        if 'units' in keywords.keys():
           _units= keywords['units']
        
        if '_par_expr_text' in keywords.keys():
            if keywords['_par_expr_text'] is None:
                pass
            else:
                self._par_expr_text = keywords['_par_expr_text']

        self._val = Value(val=_v, islog=_l,units=_units)

        self.hidden=False
        
        for kw in self.allowed_keywords.keys():
            if kw == 'val':
                pass
            elif kw == 'log':
                pass
            elif kw == 'units':
                pass
            elif kw == 'par_type':
                pass
            else:
                setattr(self,kw,self.allowed_keywords[kw])


        self.set(**keywords)




    @property
    def linked(self):
        """Linked.
        
        Returns
        -------
        object
            Requested value.
        """
        return self._linked

    @property
    def islog(self):
        """Islog.
        
        Returns
        -------
        object
            Requested value.
        """
        return self._val.islog

    @property
    def val(self ):
        """Val.
        
        Returns
        -------
        object
            Requested value.
        """
        return self._val.val

    @property
    def val_lin(self):
        """Val lin.
        
        Returns
        -------
        object
            Requested value.
        """
        return self._val.lin

    @property
    def val_log(self):
        """Val log.
        
        Returns
        -------
        object
            Requested value.
        """
        return self._val.log

    @val.setter
    def val(self,val):
        """Val.
        
        Parameters
        ----------
        val : object
            Value to assign.
        """
        self.set(val=val)

    @property
    def units(self):
        """Units.
        
        Returns
        -------
        object
            Requested value.
        """
        return self._val.units

    @units.setter
    def units(self,val):
        """Units.
        
        Parameters
        ----------
        val : object
            Value to assign.
        """
        self._val.units=val

    @property
    def adimensional(self):
        """Adimensional.
        
        Returns
        -------
        object
            Requested value.
        """
        if hasattr(self._val,'_adimensional'):
            return self._val._adimensional
        else:
            return False

    def to(self, units):
        """To.
        
        Parameters
        ----------
        units : object
            Parameter controlling units.
        
        Returns
        -------
        object
            Computed result.
        """
        try:
            return (self.val*self.units).to(units)
        except Exception as e:
            print(e)
            message='the units of this parameter:'+ str(self.units) +', are not convertible'
            warnings.warn(message)
    @property
    def fit_range(self):

        """Fit range.
        
        Returns
        -------
        object
            Requested value.
        """
        return [self.fit_range_min,self.fit_range_max]

    @fit_range.setter
    def fit_range(self, fit_range=None):
        """Fit range.
        
        Parameters
        ----------
        fit_range : object, optional
            Range for fit.
        """
        if fit_range is None:
            fit_range=[None,None]
        if isinstance(fit_range,tuple):
            pass
        elif isinstance(fit_range,list):
            pass
        else:
            raise RuntimeError('fit_range bust me list or tuple with length=2')
        if len(fit_range)!=2:
            raise RuntimeError('fit_range bust me list or tuple with length=2')

        self.fit_range_min=fit_range[0]
        self.fit_range_max=fit_range[1]
        if fit_range[0] is not None and self.val<fit_range[0]:
            raise RuntimeError('par',self.name, 'value',self.val,'< fit range min',fit_range[0])
        if fit_range[1] is not None and self.val>fit_range[1]:
            raise RuntimeError('par',self.name, 'value',self.val,'> fit range max',fit_range[1])

    def reset_dependencies(self):
        """Reset dependencies."""
        for mp in self._master_pars:
            if hasattr(mp, '_depending_pars'):
                if self in mp._depending_pars:
                    mp._depending_pars.remove(self)
        self._linked = False
        self._linked_root_model = None
        self._func=None
        self._is_dependent = False
        self._master_par_list=[]
        #NOTE: probably I have commented out self._depending_pars=[]
        #NOTE  since I don't want to remove the depending pars
        #NOTE: probably is to be done only on model loading?
        self._depending_pars=[]
        self._master_pars=[]
        self.par_expr=None
       
        

    def _add_depending_par(self,par):
        if par not in self._depending_pars:
            self._depending_pars.append(par)
      
    def _add_master_par(self,par,verbose=False):
        if par not in self._master_pars :
            if verbose:
                print("adding par:",par.name,"to ",self.name)
            self._master_pars.append(par)

    @property
    def par_expression_source_code(self):
        """Par expression source code.
        
        Returns
        -------
        object
            Requested value.
        """
        if hasattr(self,'_par_expr_text'):
            pass
        else:
            self._set_par_expr_source_code()
        print('==> par', self.name, 'is depending on', [_p.name for _p in self._master_pars],  f'according to expr:   {self.name} =\n{self._par_expr_text}'.format(self.name,self._par_expr_text))

    def _set_par_expr_source_code(self):
        if isinstance(self.par_expr,str):
            _par_expr_text=self.par_expr
        else:
            try:
                _par_expr_text=inspect.getsource(self.par_expr)
            except Exception as e:
                print('the source code of the function was not accessible due to the following error: ',e)
                if hasattr(self,'_par_expr_text'):
                    print('recovering saved code')
                    _par_expr_text=self._par_expr_text
                else:
                    _par_expr_text=None
                    
        self.set(_par_expr_text=_par_expr_text,skip_dep_par_warning=True)

    @property
    def par_expr(self):
        """Par expr.
        
        Returns
        -------
        object
            Requested value.
        """
        return self._depending_par_expr

    @par_expr.setter
    def par_expr(self, expr_string):
        """Par expr.
        
        Parameters
        ----------
        expr_string : object
            Parameter controlling expr string.
        """
        self._depending_par_expr = expr_string

    def _eval_par_func(self):
        #transform par name and value into a local var
        #print('working on ',self.name)
        #TODO:THIS HOLDS ONLY FOR NUMPY <1.22, should be removed
        warnings.filterwarnings('ignore', message='invalid value encountered in reciprocal*')
       
        if type(self._depending_par_expr) == str:
            _par_values = {}
            for _user_par_ in self._master_pars:
                if _user_par_.adimensional:
                    _par_values[_user_par_.name] = _user_par_.val_lin
                else:
                    _par_values[_user_par_.name] = _user_par_.val_lin*u.Unit(str(_user_par_.units))
            res = eval(self._depending_par_expr, {"np": np, "u": u, "__builtins__": __builtins__}, _par_values)
        elif callable(self._depending_par_expr) is True:
            _par_values={}
            for ID, _user_par_ in enumerate(self._master_pars):
                if _user_par_.adimensional:
                    _par_values[_user_par_.name] = _user_par_.val_lin
                else:
                    _par_values[_user_par_.name] = _user_par_.val_lin*u.Unit(str(_user_par_.units))
            res=self._depending_par_expr(**_par_values)
     
        _unit = None
       
        if hasattr(res,'unit'):
            _unit=res.unit
 
        if hasattr(res,'value'):
            res=res.value
            

        if self.islog is True:
            res=np.log10(res)
       

        return res


    def set(self, *args, skip_dep_par_warning=False, **keywords):
        """
        sets a parameter value checking for physical boundaries
        
        Parameters: keywords of the constructor
         
        """
        keys = keywords.keys()
        if self.immutable is False or skip_dep_par_warning is True:
            pass
        else:
            #warnings.warn('\n\n *** you are trying to set a dependent parameter:%s *** \n'%self.name)
            raise SettingDependentParError('\n\n *** you are trying to set a dependent parameter:%s *** \n'%self.name)
            #return


        for kw in keys:

            if kw in  self.allowed_keywords.keys() :
                if kw == 'val':

                    if self.allowed_values is not None:

                        if keywords[kw] not in self.allowed_values:
                            raise RuntimeError('parameter  %s' %(self.name), 'the value', keywords[kw] , 'is not in the allowed list',self.allowed_values)

                    self._val.val = keywords[kw]
                    #excluding both [] and None
                    if self._depending_pars:
                        for p in self._depending_pars:
                            _dep_func = p._func
                            if callable(_dep_func) is False:
                                raise RuntimeError('depending parameter', p.name, 'has no callable _func')
                            p.set(val=_dep_func(),skip_dep_par_warning=True)
                elif kw == 'log':
                    self._val.islog = keywords[kw]

                elif kw== 'units':
                    self._val.units = keywords[kw]

                elif kw == 'par_type':
                    if self.allowed_par_types is not None:

                        if  keywords[kw] not in self.allowed_par_types:
                            msg = "parameter %s the  type %s is not allowed" % (self.name,keywords[kw]) + "\n please choose among %s" % self.allowed_par_types
                            raise ValueError("%s" % msg)
                    setattr(self, kw, keywords[kw])
                else:
                    setattr(self,kw,keywords[kw])
            else:
                
                print ("wrong keyword=%s, not in%s "%(kw, self.allowed_keywords.keys()))
                
                raise ValueError

        if self.fit_range is not None:
            self.fit_range_max=self.fit_range[1]
            self.fit_range_min=self.fit_range[0]
 
        self.val_last_call=self.val    
        
        if self.val_min is not None:
            if self.val<self.val_min:
                raise RuntimeError("par=%s  = %e out of boundary=%e"%(self.name,self.val,self.val_min))
         
        if self.val_max is not None:
            if self.val>self.val_max:
                raise RuntimeError("par=%s  = %e out of boundary=%e"%(self.name,self.val,self.val_max))

        if self.val_min is not None and self.fit_range_min is None:
            self.fit_range_min=self.val_min
        
        if self.val_max is not None and self.fit_range_max is None:
            self.fit_range_max=self.val_max
    
    
    
    def get(self,*args ):
        
       
        """Get.
        
        Parameters
        ----------
        *args : tuple
            Additional positional arguments.
        
        Returns
        -------
        object
            Computed result.
        """
        for arg in args:
            
            if arg in  self.allowed_keywords.keys():
                
                return getattr(self,arg)
                    
            else:
                
                print ("wrong arg=%s, not in%s "%(arg, self.allowed_keywords.keys()))
                
                raise ValueError

    @property
    def immutable(self):
        """Immutable.
        
        Returns
        -------
        object
            Requested value.
        """
        return self._is_dependent or self._linked

    @property
    def frozen(self):
        """Frozen.
        
        Returns
        -------
        object
            Requested value.
        """
        return self._frozen

    @frozen.setter
    def frozen(self,v,skip_dep_par_warning=False):
        """Frozen.
        
        Parameters
        ----------
        v : object
            Parameter controlling v.
        skip_dep_par_warning : bool, optional
            If ``True``, skip dep par warning.
        """
        if self.immutable is True and skip_dep_par_warning is True:
            raise RuntimeError('frozen state of linked/dependent parameter:',self.name , 'can not be changed, please update your script')

        if v not in [True,False]:
            raise RuntimeError('par',self.name,'only True or False are allowed')
        else:
            self._frozen = v



    def freeze(self):
        """
        freezes a parameter
        """
        self.frozen=True

    
        
    def free(self):
        """
        make a parameter free
        """ 
        self.frozen=False



    def get_fit_initial_value(self):
        """
        Gives the initial fit value of the parameter

        Returns
        -------
        val_start : value    
            the parameter initial value
        """
        return self.val_start
    
    def set_fit_initial_value(self,val):
        """
        Sets the initial fit value of the parameter

       
        """
        self.val_start=val
    
    
    def show(self):
        """
        Prints the description of a parameter
        """
        print (self.get_description())
    
    def show_best_fit(self):
        """
        Prints the best-fit description of a parameter
        """
        print (self.get_bestfit_description())
    
    def get_description(self,nofields=False):
        """
        gives the value of each member of the  :class:`ModelParameter` objects, except the best-fit values
 
        Returns
        -------
        descr : (str)
            a string describing all the parameter values, except the best-fit values
        """

        if self.val_min is None:
            val_min='No'
        else:
            val_min=str("%+e"%self.val_min)
            
        if self.val_max is None:
            val_max='No'
        else:
            val_max=str("%+e"%self.val_max)
        
        if self.val is None:
            val='No'
        else:
            val=str("%+e"%self.val)
        
        if nofields==False:
            descr= "name = %-16s  type = %-20s  units = %-16s  val = %s  phys-bounds = [%-13s,%-13s] islog = %s  frozen= %s "%(self.name, self.par_type ,
                                                                                   self.units, val, val_min, val_max,self.islog,self.frozen)
        else:
            descr= " %-16s | %-20s | %-16s | %s | [%-13s,%-13s] | %s | %s"%(self.name, self.par_type ,
                                                                                   self.units, val, val_min, val_max,self.islog,self.frozen)
            
        return descr


    def get_bestfit_description(self,nofields=False):
        """
        gives the value of each member of the  :class:`ModelParameter` objects, suited for  the best-fit values
        
        Returns
        -------
        descr : (str)
            a string describing all the parameter values, suited for  the best-fit values
        """
        
        if self.val_start is None:
            val_start='No'
        else:
            val_start=str("%+e"%self.val_start)
        
        if self.fit_range_min is None:
            fit_range_min='No'
        else:
            fit_range_min=str("%+e"%self.fit_range_min)
            
        if self.fit_range_max is None:
            fit_range_max='No'
        else:
            fit_range_max=str("%+e"%self.fit_range_max)
        
        if self.best_fit_val is None:
            best_fit_val='No'
        else:
            best_fit_val=str("%+e"%self.best_fit_val)

        if hasattr(self,'err_p') and  hasattr(self,'err_m'):
            best_fit_err_m=str("%+e"%self.err_m)
            best_fit_err_p=str("%+e"%self.err_p)
        else:
            best_fit_err_m='#'
            if self.best_fit_err is None:
                best_fit_err_p='No'
            else:
                best_fit_err_p=str("%+e"%self.best_fit_err)
        
        if self.frozen==True:
            best_fit_val='Frozen'
            best_fit_err_p='Frozen'
            best_fit_err_m='Frozen'
        if nofields==False:
            descr= "name = %-16s best-fit val=%-13s  best-fit err p=%-13s best-fit err m=%-13s start-val=%-13s   fit-bounds=[%-13s,%-13s]"%(self.name,best_fit_val,best_fit_err_p,best_fit_err_m,val_start,fit_range_min,fit_range_max)
        else:
            descr= " %-16s | %-13s | %-13s | %-13s | %-13s | [%-13s,%-13s]"%(self.name,best_fit_val,best_fit_err_p,best_fit_err_m,val_start,fit_range_min,fit_range_max)


        return descr

    def identity_func(self):
        """Identity func.
        
        Returns
        -------
        object
            Computed result.
        """
        return self._root_par.val


    def make_log(self):
        if self.islog:
            pass
        else:
            
            if self.val_max is not None:
               self.val_max=self._handle_zero_in_log_pars(self.val_max)
            if self.val_min is not None:
                self.val_min=self._handle_zero_in_log_pars(self.val_min)
           
            
            if self.val_start is not None:
                self.val_start=self._handle_zero_in_log_pars(self.val_start)
            if self.val_last_call is not None:
                self.val_last_call=self._handle_zero_in_log_pars(self.val_last_call)
            if self.fit_range_min is not None:
                self.fit_range_min=self._handle_zero_in_log_pars(self.fit_range_min)
            if self.fit_range_max is not None:
                self.fit_range_max=self._handle_zero_in_log_pars(self.fit_range_max)

            lin_val=self.val_lin
            #NOTE: get lin val before setting par log
            self._val.islog=True
            self.set(val=self._handle_zero_in_log_pars(lin_val))
    def _handle_zero_in_log_pars(self,v):
        if v<=0:
            return -200
        else:
            return np.log10(v)


def create_a_function( **kwargs):

    """Create a function.
    
    Parameters
    ----------
    **kwargs : dict
        Additional keyword arguments.
    
    Returns
    -------
    object
        Computed result.
    """
    def function_template(**kwargs):
        """Function template.
        
        Parameters
        ----------
        **kwargs : dict
            Additional keyword arguments.
        """
        return

    return function_template


class CompositeModelParameterArray(object):
    """Unified parameter interface across multiple model components.

    Notes
    -----
    Collects per-component parameter arrays, supports cross-component linking,
    and builds merged parameter/best-fit tables for composite models.
    """

    def __repr__(self):
        return str(self.show_pars())

    def __init__(self,):

        """Create a new `CompositeModelParameterArray` instance."""
        self.all_frozen = False
        self._parameters = []
        self._model_comp = []

    def reset_dependencies(self):
        """Reset dependencies."""
        for p in self.par_array:
            p.reset_dependencies()




    def link_par(self,par_name,model_name_list,root_model_name):
        """Link par.
        
        Parameters
        ----------
        par_name : object
            Parameter controlling par name.
        model_name_list : object
            List of model name.
        root_model_name : object
            Parameter controlling root model name.
        """
        m_root=self.get_model_by_name(root_model_name)
        p_root=m_root.get_par_by_name(par_name)

        for m_name in model_name_list:
            dep_model=self.get_model_by_name(m_name)

            if dep_model is not None:
                dep_par = dep_model.get_par_by_name(par_name)
                if dep_par is not None:
                    if p_root == dep_par:
                        raise RuntimeError(" root and linked parameter can't be the same")
                    if dep_par.immutable is True:
                        raise RuntimeError(" this parameter is already linked or dependent ")
                    if m_root==dep_par.model:
                        raise RuntimeError(" linked and root model must be different")
                    
                    dep_par._root_par = p_root
                    dep_par._func=dep_par.identity_func
                    dep_par._linked = True
                    dep_par._is_dependent = True
                    dep_par._linked_root_model=m_root
                    dep_par._add_master_par(p_root)
                    p_root._add_depending_par(dep_par)
                    dep_par.freeze()
                     
                    p_root.set(val=p_root.val,skip_dep_par_warning=True)

            else:
                 self._handle_missing_component_error(m_name)



    def _handle_missing_component_error(self, model_name):

        if hasattr(model_name, 'name'):
            name = model_name.name
        elif type(model_name) == str:
            name = model_name
        else:
            name = 'passed'
        s = 'Model component'+ name +' not present'
        raise RuntimeError(s)

    def add_model_parameters(self, model):
        """Add model parameters.
        
        Parameters
        ----------
        model : object
            Model instance.
        """
        try:
            assert (model.name not in [m.name for m in self._model_comp])
        except:
            raise RuntimeError('model name:', model.name, 'already assigned')


        model.parameters.model=model
        for p in model.parameters.par_array:
            p.model=model

        self._parameters.append(model.parameters)
        self._model_comp.append(model)

    def del_model_parameters(self, model_name):
        """Del model parameters.
        
        Parameters
        ----------
        model_name : object
            Parameter controlling model name.
        """
        m,ID=self.get_model_by_name(model_name,get_idx=True)
        if m is not None:
            #_p = self._comp_par_array.pop(ID)
            _p = self._model_comp.pop(ID)
            _p = self._parameters.pop(ID)
        else:
             self._handle_missing_component_error(model_name)


    def get_model_by_name(self, model_name,get_idx=False):
        """Return model by name.
        
        Parameters
        ----------
        model_name : object
            Parameter controlling model name.
        get_idx : bool, optional
            Index/identifier for get idx.
        
        Returns
        -------
        object
            Requested value.
        """
        try:
            if isinstance(model_name, str):
                pass
            else:
                model_name = model_name.name

            selected_model = None
            idx=None
            for ID,m in enumerate(self._model_comp):
                if m.name == model_name:
                    selected_model = m
                    idx=ID


            if selected_model is None:
                self._handle_missing_component_error(model_name)
            if get_idx is False:
                return selected_model
            else:
                return selected_model,idx
        except Exception as e:
           message = str(e)
           message += '\n'
           message += 'Starting from veriosn 1.2.0, FitModel is a CompositeModel, hence to set parameters ' \
                      'you have to pass as first paramter the model name or model object of the corresponing parameter e.g. \n' \
                      '''   
                            fit_model.set_par('model-name',value) 
                            OR 
                            fit_model.set_par(jet,value) 
                      '''
    def get_par_by_name(self, model_name,par_name):
        """Return par by name.
        
        Parameters
        ----------
        model_name : object
            Parameter controlling model name.
        par_name : object
            Parameter controlling par name.
        
        Returns
        -------
        object
            Requested value.
        """
        p=None
        m=self.get_model_by_name(model_name)

        if m is None:
            self._handle_missing_component_error(model_name)
        else:
            p=m.parameters.get_par_by_name(par_name)

        return p


    def _build_par_table(self):
        _l=[]
        for ID, mc in enumerate(self._model_comp):
            self._parameters[ID]._build_par_table()
            if len(self._parameters[ID]._par_table) >= 1:
     
                _l.append(self._parameters[ID]._par_table)
        self._par_table=vstack(_l)

    def _build_best_fit_par_table(self):
        _l = []
        for ID, mc in enumerate(self._model_comp):
            self._parameters[ID]._build_best_fit_par_table()
            if len(self._parameters[ID]._best_fit_par_table) >= 1:
                
                _l.append(self._parameters[ID]._best_fit_par_table)
        self._best_fit_par_table = vstack(_l)

    @property
    def par_array(self):
        """Par array.
        
        Returns
        -------
        object
            Requested value.
        """
        pa=[]
        for p in  self._parameters:
            pa.extend(p.par_array)

        return pa

    @property
    def par_table(self):
        """Par table.
        
        Returns
        -------
        object
            Requested value.
        """
        self._build_par_table()
        return self._par_table

    @property
    def best_fit_par_table(self):
        """Best fit par table.
        
        Returns
        -------
        object
            Requested value.
        """
        self._build_best_fit_par_table()
        return self._best_fit_par_table

    def show_pars(self, getstring=False, sort_key=None):
        """Display pars.
        
        Parameters
        ----------
        getstring : bool, optional
            Parameter controlling getstring.
        sort_key : object, optional
            Parameter controlling sort key.
        
        Returns
        -------
        object
            Computed result.
        """
        self._build_par_table()
        if sort_key is not None:
            self.par_table.sort(sort_key)

        if getstring == True:
            return self.par_table.pformat_all()
        else:
            #self.par_table.pprint_all()
            _show_table(self.par_table)

    def show_best_fit_pars(self, getstring=False):
        """Display best fit pars.
        
        Parameters
        ----------
        getstring : bool, optional
            Parameter controlling getstring.
        
        Returns
        -------
        object
            Computed result.
        """
        self._build_best_fit_par_table()
        if getstring == True:
            return self._best_fit_par_table.pformat_all()
        else:
            #return self._best_fit_par_table
            #.pprint_all()
            _show_table(self._best_fit_par_table)

    #@compositr_parameter_setter
    def freeze(self, model_name,par_name):
        """Freeze.
        
        Parameters
        ----------
        model_name : object
            Parameter controlling model name.
        par_name : object
            Parameter controlling par name.
        """
        self.set(model_name,par_name, 'frozen')

    #@compositr_parameter_setter
    def free(self,model_name, par_name):
        """Free.
        
        Parameters
        ----------
        model_name : object
            Parameter controlling model name.
        par_name : object
            Parameter controlling par name.
        """
        self.set(model_name,par_name, 'free')

    #@compositr_parameter_setter
    def set(self, model_name, par_name, *args, **kw):
        """Set.
        
        Parameters
        ----------
        model_name : object
            Parameter controlling model name.
        par_name : object
            Parameter controlling par name.
        *args : tuple
            Additional positional arguments.
        **kw : dict
            Parameter controlling kw.
        """
        m=self.get_model_by_name(model_name)
        if m is not None:
            m.parameters.set(par_name, *args, **kw)
        else:
            self._handle_missing_component_error(model_name)

    #@compositr_parameter_setter
    def set_par(self,model_name, par_name, val):
        """Set par.
        
        Parameters
        ----------
        model_name : object
            Parameter controlling model name.
        par_name : object
            Parameter controlling par name.
        val : object
            Value to assign.
        """
        m=self.get_model_by_name(model_name)
        if m is not None:
            m.parameters.set(par_name, val=val)
        else:
            self._handle_missing_component_error(model_name)

    #@compositr_parameter_setter
    def get(self, model_name,par_name, field_name, *args, **kw):
        #print('-->', par_name, field_name)
        """Get.
        
        Parameters
        ----------
        model_name : object
            Parameter controlling model name.
        par_name : object
            Parameter controlling par name.
        field_name : object
            Parameter controlling field name.
        *args : tuple
            Additional positional arguments.
        **kw : dict
            Parameter controlling kw.
        
        Returns
        -------
        object
            Computed result.
        """
        m = self.get_model_by_name(model_name)
        if m is not None:
            pass
        else:
            self._handle_missing_component_error(model_name)
        #print('-->',par_name,field_name)
        return m.parameters.get(par_name, field_name, *args, **kw )

    #@compositr_parameter_setter
    def get_val(self, model_name,par_name):
        """Return val.
        
        Parameters
        ----------
        model_name : object
            Parameter controlling model name.
        par_name : object
            Parameter controlling par name.
        
        Returns
        -------
        object
            Requested value.
        """
        m = self.get_model_by_name(model_name)
        if m is not None:
            pass
        else:
            self._handle_missing_component_error(model_name)

        return m.parameters.get(par_name, 'val')

    def freeze_all(self):
        """Freeze all."""
        self.all_frozen = True
        for p_arr in self.par_array:
           p_arr.freeze()

    def free_all(self):
        """Free all."""
        self.all_frozen = False
        for p_arr in self.par_array:
            p_arr.free()


class ModelParameterArray(object):
    """
    This class provide and interface to handle an array of :class:`ModelParameter` objects.
    
    
    Attributes
    ----------
    par_array :  list 
        list of :class:`ModelParameter` objects
    """   

    def __repr__(self):
        return str(self.show_pars())



    def __init__(self,model=None):
            
        """
        Constructor 
        """
        
        self.par_array=[]
        self.all_frozen=False

        self.properties={}
        self.model=model
        self._numeric_fields = ['val', 'phys. bound. min', 'phys. bound. max','bestfit val','err +','err -','start val','fit range min','fit range max']


    def reset_dependencies(self):
        """Reset dependencies."""
        for p in self.par_array:
            p.reset_dependencies()

    def add_par(self,par):
        """
        adds a new :class:`ModelParameter` object  to the `par_array`
        """
        try:
            assert (isinstance(par,ModelParameter))
        except:
            raise RuntimeError('parameter is not an instance of',type(ModelParameter))

        try:
            assert (par.name not in [p.name for p in self.par_array])
        except:
            raise RuntimeError('parameter name:',par.name,'already assigned')

        par.model = self.model
        self.par_array.append(par)

        setattr(self,clean_var_name(par.name), par)
        self.properties[par.name]=par
        par.model=self.model
    @property
    def names(self):
        """Names.
        
        Returns
        -------
        object
            Requested value.
        """
        return [p.name for p in self.par_array]

    def _build_par_table(self,names_list=None):
        _model_name = []
        _name=[]
        _type=[]
        _unit=[]
        _val=[]
        _bound_min=[]
        _bound_max=[]
        _islog=[]
        _frozen=[]
        _fields=[_model_name,_name,_type,_unit,_val,_bound_min,_bound_max,_islog,_frozen]
        _names=['model name','name','par type','units','val','phys. bound. min','phys. bound. max','log','frozen']


        if self.model is  None:
           _fields.pop(0)
           _names.pop(0)

        for par in self.par_array:

            append=False

            if names_list is not None:
                if par.name in names_list:
                    append=True
            else:
                append = True

            if type(par.val) is str:
                try:
                    ast.literal_eval(par.val)
                except:
                    append= False

            if par.hidden is True:
                append=False

            if append is True:
                if self.model is not None:
                    _model_name.append(self.model.name)
                else:
                    _model_name.append('no_name')

                if par._linked is True and par._linked_root_model != self.model:
                    _p_name = par.name+'(L,%s)'%par._linked_root_model.name

                    _type.append(par.par_type)
                    _unit.append(par.units)
                    _val.append(None)
                    _bound_min.append(None)
                    _bound_max.append(None)
                    _islog.append(par.islog)
                    _frozen.append(par.frozen)
                else:
                    if par._linked is False and par._depending_pars!=[]:
                        _p_name = par.name + '(M)'
                    else:
                        _p_name = par.name

                    _type.append(par.par_type)
                    _unit.append(par.units)
                    _val.append(par.val)
                    _bound_min.append(par.val_min)
                    _bound_max.append(par.val_max)
                    _islog.append(par.islog)
                    _frozen.append(par.frozen)

                if par._is_dependent is True and par._linked is False:
                    for p in par._master_pars:
                        _p_name = '*'+ par.name + '(D,%s)' % p.name
                    if par._linked_root_model is not None:
                        if par._linked_root_model != self.model:
                            _p_name +=  par._linked_root_model.name

                _name.append(_p_name)

        t=Table(_fields,names=_names,masked=False)

        self._fromat_column_entry(t)

        self._par_table= t


    def _fromat_column_entry(self, t):
        for n in self._numeric_fields:
            if n in t.colnames:

                try:
                    if None in t[n].data:

                        t[n] = MaskedColumn(t[n].data, name=n, dtype=np.float64, mask=t[n].data==None)
                    else:
                        t[n] = MaskedColumn(t[n].data, name=n, dtype=np.float64)
                    t[n].format = '%e'
                except:

                    for ID,v in enumerate(t[n].data):
                        try:
                            c=ast.literal_eval(t[n].data[ID])
                            if type(c) == int:
                                t[n].data[ID] = '%d' % c
                            else:
                                t[n].data[ID] = '%e' % c
                        except:
                           pass




    def _build_best_fit_par_table(self, names_list=None):

        _name=[]
        _model_name = []
        _best_fit_val=[]
        _best_fit_err_p=[]
        _best_fit_err_m=[]
        _val_start=[]
        _fit_range_min=[]
        _fit_range_max=[]
        _frozen=[]
        _val=[]
        _fields=[_model_name,_name,_val,_best_fit_val,_best_fit_err_p,_best_fit_err_m,_val_start,_fit_range_min,_fit_range_max,_frozen]
        _names=['model name','name','val','bestfit val','err +','err -','start val','fit range min','fit range max','frozen']

        if self.model is None:
           _fields.pop(0)
           _names.pop(0)

        for par in self.par_array:

            append=False

            if names_list is not None:
                if par.name in names_list:
                    append=True
            else:
                append = True

            if par.hidden is True:
                append=False

            if type(par.val_start) is str:
                try:
                    ast.literal_eval(par.val_start)
                except:
                    append= False

            if append:
                if self.model is not None:
                    _model_name.append(self.model.name)
                else:
                    _model_name.append('no_name')

                if par._linked is True and par._linked_root_model != self.model:
                    _p_name = par.name + '(L,%s)' % par._linked_root_model.name

                    _val_start.append(None)
                    _best_fit_val.append(None)
                    _best_fit_err_p.append(None)
                    _best_fit_err_m.append(None)

                    _fit_range_min.append(par.fit_range_min)
                    _fit_range_max.append(par.fit_range_max)

                else:
                    if par._linked is False and par._depending_pars!=[]:
                        _p_name = par.name + '(M)'
                    else:
                        _p_name = par.name
                    if par.frozen == True:
                        best_fit_val = None
                        best_fit_err_p = None
                        best_fit_err_m = None
                    else:
                        best_fit_val=par.best_fit_val
                        if hasattr(par, 'err_p') and hasattr(par, 'err_m'):
                            best_fit_err_m =  par.err_m
                            best_fit_err_p =  par.err_p
                        else:
                            best_fit_err_m = None
                            if par.best_fit_err is None:
                                best_fit_err_p = None
                            else:
                                best_fit_err_p =  par.best_fit_err

                    _val_start.append(par.val_start)
                    _best_fit_val.append(best_fit_val)
                    _best_fit_err_p.append(best_fit_err_p)
                    _best_fit_err_m.append(best_fit_err_m)

                    _fit_range_min.append(par.fit_range_min)
                    _fit_range_max.append(par.fit_range_max)

                if par._is_dependent is True and par._linked is False:
                    for p in par._master_pars:
                        _p_name = '*' + par.name + '(D,%s)' % p.name
                    if par._linked_root_model is not None:
                        if par._linked_root_model != self.model:
                            _p_name += par._linked_root_model.name

                _name.append(_p_name)
                _frozen.append(par.frozen)
                _val.append(par.val)

        t = Table(_fields, names=_names, masked=False)

        self._fromat_column_entry(t)

        self._best_fit_par_table= t



    def __setattr__(self, name, value):
        if "properties" in self.__dict__ and name in self.properties:
            raise AttributeError('this member is protected, use del_par,add_par, to add/remove, and set() to set values or .val attribute')
        else:
            self.__dict__[name] = value


    def del_par(self,par):
        
        """Del par.
        
        Parameters
        ----------
        par : object
            Parameter controlling par.
        """
        self.par_array.remove(par)
        delattr(self,par.name)
        self.properties.pop(par.name)

    def get_par_by_name(self,name, verbose=False):
        """
        selects a parameter by name
       
        Parameters
        ----------
        name : (str) parameter name
        
        Returns
        -------
        item : the :class:`ModelParameter` element of the `par_array` with the corresponding name
        
        """
        for pi in range(len(self.par_array)):
            if self.par_array[pi].name==name:
                return self.par_array[pi]
        else:
            if verbose:
                print ("no par with name %s found"%name)
                print ("pars in array are:")
                self.show_pars()
            
            return None
                
    def get_par_by_type(self,par_type):
        """

        get parameter by type

        """
        for param in self.par_array:
            if param.par_type==par_type:
                return param

        return None
    

    def get_pars_by_type(self,par_type):
        """

        get parameter by type

        """
        pars=[]
        for param in self.par_array:
            if param.par_type==par_type:
                pars.append(param)

        return pars
    
    @property
    def par_table(self):
        """Par table.
        
        Returns
        -------
        object
            Requested value.
        """
        self._build_par_table()
        return self._par_table

    @property
    def best_fit_par_table(self):
        """Best fit par table.
        
        Returns
        -------
        object
            Requested value.
        """
        self._build_best_fit_par_table()
        return self._best_fit_par_table

    def show_pars(self,getstring=False,names_list=None,sort_key=None):

        """Display pars.
        
        Parameters
        ----------
        getstring : bool, optional
            Parameter controlling getstring.
        names_list : object, optional
            List of names.
        sort_key : object, optional
            Parameter controlling sort key.
        
        Returns
        -------
        object
            Computed result.
        """
        self._build_par_table(names_list=names_list)
        if sort_key is not None:
            self.par_table.sort(sort_key)

        if getstring==True:
            return self.par_table.pformat_all()
        else:
            _show_table(self.par_table)
 
    
    
    def show_best_fit_pars(self,getstring=False):

        """Display best fit pars.
        
        Parameters
        ----------
        getstring : bool, optional
            Parameter controlling getstring.
        
        Returns
        -------
        object
            Computed result.
        """
        self._build_best_fit_par_table()
        if getstring == True:
            return self.best_fit_par_table.pformat_all()
        else:
            return self.best_fit_par_table


    def set(self,par_name,*args, **keywords):
        
        """
        sets to a given value a given parameter
        
        Parameters
        ----------
        
        par_name     : (str) name of the parameter 
        keywords     :  keywords to set the value or the range of the parameter 
            
         
         
           
        Examples
        --------
        if parameters is a :class:`ModelParameterArray` object:
        
        .. code::
     
            parameters.set('R',val=1E16) 
            parameters.set('R',fit_range=[1E16,1E17])
        
        

        """
        
        
        par=self.get_par_by_name(par_name)

        args_list=['free','frozen']

        if  args!=():
            for  arg in args:
            
                if arg not in args_list:
                    print ("argument: %s,  not in allowed args"%arg)
                    print ("allowed args= ",args_list)
                    return
                
                 
                if arg=='frozen':                
                    par.freeze()
                
                
                if  arg=='free':
                    par.free()

                
            
        
        if keywords!={}:
        
            if par is None:
                raise RuntimeWarning('parameter %s is not present in the model'%par_name)
            else:
                par.set( **keywords)
        
        
        
    def get(self,par_name,arg):
        
        """
        gets the argument of a given parameter
        
        Parameters
        ----------
        
        par_name     : (str) name of the parameter 
        arg     :  keyword
            
         
         
           
        Examples
        --------
        if parameters is a :class:`ModelParameterArray` object:
        
        .. code::
     
            parameters.get('R') 
            parameters.get('frozen')
        
        

        """
        
        
        par=self.get_par_by_name(par_name)
        
        return par.get(arg)
    
    def freeze_all(self):
        """Freeze all."""
        self.all_frozen=True
        for pi in range(len(self.par_array)):
            self.par_array[pi].freeze()
                
        
    def free_all(self):
        """Free all."""
        self.all_frozen = False
        for pi in range(len(self.par_array)):
            self.par_array[pi].free()
        

    def _serialize_pars(self):
        _par_keys=['val','val_min','val_max','val_start','val_last_call','fit_range_min','fit_range_max','best_fit_val',
                   'best_fit_err','frozen','allowed_values','_linked','_is_dependent','_func','_master_pars',
                   '_linked_root_model','_depending_pars','_root_par','','_master_par_list','_depending_par_expr','_par_expr_text','units','par_type','log']
        _par_dict = {}
        for par in self.par_array:
            _val_dict={}

            for k in _par_keys:
                if hasattr(par,k):
                    if k=='units':
                        _val_dict[k]=str(getattr(par,k))
                    else:
                        _val_dict[k]=getattr(par,k)

            if par.islog:
                _val_dict['log']=True
            else:
                _val_dict['log']=False
            
            _par_dict[par.name] = _val_dict
           

        return _par_dict

    def _decode_pars(self,_par_dict):
        for p_name in _par_dict.keys():
            p=self.get_par_by_name(p_name)
            if p._is_dependent is False:
                self.set(p_name,**_par_dict[p_name])
