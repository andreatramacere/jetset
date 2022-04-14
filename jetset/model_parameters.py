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


def is_notebook():
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
            display(t.show_in_notebook(show_row_index=False, display_length=100))
        except:
            try:
                from IPython.display import display
                display(t)
            except:
                t.pprint_all()
    else:
        t.pprint_all()

class Value(object):

    def __init__(self,val,units,islog=False):
        self.val=val
        self.islog=islog
        self.units=units

    def __repr__(self):
        return '%s'%self._val

    @property
    def islog(self):
        return  self._islog

    @islog.setter
    def islog(self,val):
        self._islog=val

    @property
    def val(self):
        return self._val

    @val.setter
    def val(self, val):
        self._val = val

    @property
    def lin(self):
        if self._val is None:
            return None
        if self._islog is True:
            return 10**self._val
        else:
            return self._val

    @property
    def log(self):
        if self._val is None:
            return None
        if self._islog is True:
            return self._val
        else:
            return np.log10(self._val)

    @property
    def units(self):
        return self._units

    @units.setter
    def units(self, p_unit,verbose=False):
        try:
            self._units = u.Unit(p_unit)
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

        self._val = Value(val=_v, islog=_l,units=_units)

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
        return self._linked

    @property
    def islog(self):
        return self._val.islog

    @property
    def val(self ):
        return self._val.val

    @property
    def val_lin(self):
        return self._val.lin

    @property
    def val_log(self):
        return self._val.log

    @val.setter
    def val(self,val):
        self.set(val=val)

    @property
    def units(self):
        return self._val.units

    @units.setter
    def units(self,val):
        self._val.units=val

    def to(self, units):
        if isinstance(self.units,u.Unit):
            return self.val*self.units.to(units)
        else:
            message='the units of this parameter:'+self.units+', are not convertible'
            warnings.warn(message)
    @property
    def fit_range(self):

        return [self.fit_range_min,self.fit_range_max]

    @fit_range.setter
    def fit_range(self, fit_range=None):
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

        #self._fit_range = fit_range
        self.fit_range_min=fit_range[0]
        self.fit_range_max=fit_range[1]
        if fit_range[0] is not None and self.val<fit_range[0]:
            raise RuntimeError('par',self.name, 'value',self.val,'< fit range min',fit_range[0])
        if fit_range[1] is not None and self.val>fit_range[1]:
            raise RuntimeError('par',self.name, 'value',self.val,'> fit range max',fit_range[1])

    def reset_dependencies(self):
        self._linked = False
        self._linked_root_model = None
        self._func=None
        self._is_dependent = False
        self._master_pars=[]
        self._depending_pars=[]

    def _add_depending_par(self,par):
        if par not in self._depending_pars:
            self._depending_pars.append(par)

    def _add_master_par(self,par):
        if par not in self._master_pars:
            self._master_pars.append(par)

    #def make_dependent_par(self, master_par, func, root_model=None):

    #    if self == master_par:
    #        raise RuntimeError(" root and linked parameter can't be the same")
    #    self._is_dependent = True
    #    self._func = func
    #    self._master_pars = master_pars
    #    master_par._add_depending_par(self)

    #    self.freeze()
    #    if root_model is not None:
    #        self._linked_root_model = root_model
    #    self.set(val=self._func(self._master_par.val), skip_dep_par_warning=True)


    # def get_default_args(self, par_expr):
    #     #signature = inspect.signature(func)
    #     #d={}
    #     #for k, v in signature.parameters.items():
    #     #    if isinstance(v.default,ModelParameter):
    #     #        d[k]=v.default
    #     #    else:
    #     #        raise RuntimeError('argument',k,'is not valid, should be a model parameter')
    #     #return d
    #     pass

    @property
    def print_par_expr(self):
        if isinstance(self.par_expr,str):
            _par_expr=self.par_expr
        else:
            _par_expr=inspect.getsource(self.par_expr)
        print('==> par', self.name, 'is depending on', [_p.name for _p in self._master_pars],  f'according to expr:   {self.name} =\n{_par_expr}'.format(self.name,_par_expr))
        

    @property
    def par_expr(self):
        return self._depending_par_expr

    @par_expr.setter
    def par_expr(self, expr_string):
        self._depending_par_expr = expr_string

    def _eval_par_func(self):
        #transform par name and value into a local var
        if type(self._depending_par_expr) == str:
            _par_values= [None]*len(self._master_pars)
            for ID, _user_par_ in enumerate(self._master_pars):
                _par_values[ID] = _user_par_.val_lin
                #print('==> _eval_par_func',_user_par_.name,_par_values[ID])
                exec(_user_par_.name + '=_par_values[ID]')
            res = eval(self._depending_par_expr)
        elif callable(self._depending_par_expr) is True:
            _par_values={}
            for ID, _user_par_ in enumerate(self._master_pars):
                _par_values[_user_par_.name] = _user_par_.val_lin
            res=self._depending_par_expr(**_par_values)

        if self.islog is True:
            res=np.log10(res)

        return res
        #return eval(self.par_expr)


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
            raise RuntimeError('\n\n *** you are trying to set a dependent parameter:%s *** \n'%self.name)
            #return


        for kw in keys:

            if kw in  self.allowed_keywords.keys() :
                if kw == 'val':

                    if self.allowed_values is not None:

                        if keywords[kw] not in self.allowed_values:
                            raise RuntimeError('parameter  %s' %(self.name), 'the value', keywords[kw] , 'is not in the allowed list',self.allowed_values)

                    self._val.val = keywords[kw]
                    if self._depending_pars is not []:
                        for p in self._depending_pars:
                            #print("==> setting dep par",p.name, 'to',p._func(),'p=',p,'master',self,'name',self.name)
                            p.set(val=p._func(),skip_dep_par_warning=True)

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
        
    
        
        #if self.frozen==True:
        #    self.val=self.val_start
        
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
        
       
        for arg in args:
            
            if arg in  self.allowed_keywords.keys():
                
                return getattr(self,arg)
                    
            else:
                
                print ("wrong arg=%s, not in%s "%(arg, self.allowed_keywords.keys()))
                
                raise ValueError

    @property
    def immutable(self):
        return self._is_dependent or self._linked

    @property
    def frozen(self):
        return self._frozen

    @frozen.setter
    def frozen(self,v,skip_dep_par_warning=False):
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
        return self._root_par.val


# class LinkedParameter(ModelParameter):
#
#     def __init__(self,p_name,m_list):
#         super(LinkedParameter,self).__init__()
#         self.m_list = m_list
#         self.name = p_name
#
#     #def set_par(self, par_name, val):
#     #    for model in self.m_list:
#
#     #        pm = model.get_par_by_name(self.name)
#     #        pm.set(par_name, val = val)
#
#     #def set(self, *args, **kw):
#     #    self.set_par(*args, **kw)
#
# class ModelLinkedParameter(object):
#     def __init__(self,name,p):
#
#         self.name = name
#         self.parameters = ModelParameterArray()
#         self.paramters.add_par(p)
#         self.paramters.model=self





def compositr_parameter_setter(method):
    @wraps(method)
    def func_wrapper(self, model_name, *args, **kwargs):
        print('--> model_name',args,kwargs)
        try:
            if isinstance(model_name,str):
                pass
            else:
                model_name=model_name.name
            print('--> model_name', model_name, args,kwargs)
            return method(self, *args, **kwargs)
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




           raise RuntimeError(message)

    return func_wrapper

def create_a_function( **kwargs):

    def function_template(**kwargs):
        return

    return function_template


class CompositeModelParameterArray(object):

    def __repr__(self):
        return str(self.show_pars())

    def __init__(self,):

        self.all_frozen = False
        self._parameters = []
        self._model_comp = []

    def reset_dependencies(self):
        for p in self.par_array:
            p.reset_dependencies()




    def link_par(self,par_name,model_name_list,root_model_name):
        m_root=self.get_model_by_name(root_model_name)
        p_root=m_root.get_par_by_name(par_name)

        for m_name in model_name_list:
            dep_model=self.get_model_by_name(m_name)

            if dep_model is not None:
                dep_par = dep_model.get_par_by_name(par_name)
                if dep_par is not None:
                    print('==> par:',dep_par.name, 'from model:',  dep_model.name, 'linked to same parameter in model', m_root.name )
                    if p_root == dep_par:
                        raise RuntimeError(" root and linked parameter can't be the same")
                    if dep_par.immutable is True:
                        raise RuntimeError(" this parameter is already linked or dependent ")
                    if m_root==dep_par.model:
                        raise RuntimeError(" linked and root model must be different")
                    #p_root._linked_models.append(dep_model)
                    #exec(dep_par.name+'= dep_par')
                    #identity_func=lambda p_root=p_root: p_root.val_lin
                    dep_par._root_par = p_root
                    #dep_par.make_dependent_par(dep_par.identity_func)
                    dep_par._func=dep_par.identity_func
                    dep_par._linked = True
                    dep_par._is_dependent = True
                    dep_par._linked_root_model=m_root
                    dep_par._add_master_par(p_root)
                    p_root._add_depending_par(dep_par)
                    dep_par.freeze()
                    p_root.val=p_root.val
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
        m,ID=self.get_model_by_name(model_name,get_idx=True)
        if m is not None:
            #_p = self._comp_par_array.pop(ID)
            _p = self._model_comp.pop(ID)
            _p = self._parameters.pop(ID)
        else:
             self._handle_missing_component_error(model_name)


    def get_model_by_name(self, model_name,get_idx=False):
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
                #t=copy.copy(self._parameters[ID]._par_table)
                #for c in t.columns:
                 #   t[c] = t[c].astype(np.object)
                _l.append(self._parameters[ID]._par_table)
        #print('--> _l',_l)
        self._par_table=vstack(_l)

    def _build_best_fit_par_table(self):
        _l = []
        for ID, mc in enumerate(self._model_comp):
            self._parameters[ID]._build_best_fit_par_table()
            if len(self._parameters[ID]._best_fit_par_table) >= 1:
                #t = copy.copy(self._parameters[ID]._best_fit_par_table)
                #for c in t.columns:
                #    t[c] = t[c].astype(np.object)
                _l.append(self._parameters[ID]._best_fit_par_table)
        self._best_fit_par_table = vstack(_l)

    @property
    def par_array(self):
        pa=[]
        for p in  self._parameters:
            pa.extend(p.par_array)

        return pa

    @property
    def par_table(self):
        self._build_par_table()
        return self._par_table

    @property
    def best_fit_par_table(self):
        self._build_best_fit_par_table()
        return self._best_fit_par_table

    def show_pars(self, getstring=False, sort_key=None):
        self._build_par_table()
        if sort_key is not None:
            self.par_table.sort(sort_key)

        if getstring == True:
            return self.par_table.pformat_all()
        else:
            #self.par_table.pprint_all()
            _show_table(self.par_table)

    def show_best_fit_pars(self, getstring=False):
        self._build_best_fit_par_table()
        if getstring == True:
            return self._best_fit_par_table.pformat_all()
        else:
            #return self._best_fit_par_table
            #.pprint_all()
            _show_table(self._best_fit_par_table)

    #@compositr_parameter_setter
    def freeze(self, model_name,par_name):
        self.set(model_name,par_name, 'frozen')

    #@compositr_parameter_setter
    def free(self,model_name, par_name):
        self.set(model_name,par_name, 'free')

    #@compositr_parameter_setter
    def set(self, model_name, par_name, *args, **kw):
        m=self.get_model_by_name(model_name)
        if m is not None:
            m.parameters.set(par_name, *args, **kw)
        else:
            self._handle_missing_component_error(model_name)

    #@compositr_parameter_setter
    def set_par(self,model_name, par_name, val):
        m=self.get_model_by_name(model_name)
        if m is not None:
            m.parameters.set(par_name, val=val)
        else:
            self._handle_missing_component_error(model_name)

    #@compositr_parameter_setter
    def get(self, model_name,par_name, field_name, *args, **kw):
        #print('-->', par_name, field_name)
        m = self.get_model_by_name(model_name)
        if m is not None:
            pass
        else:
            self._handle_missing_component_error(model_name)
        #print('-->',par_name,field_name)
        return m.parameters.get(par_name, field_name, *args, **kw )

    #@compositr_parameter_setter
    def get_val(self, model_name,par_name):
        m = self.get_model_by_name(model_name)
        if m is not None:
            pass
        else:
            self._handle_missing_component_error(model_name)

        return m.parameters.get(par_name, 'val')

    def freeze_all(self):
        self.all_frozen = True
        for p_arr in self._parameters:
            for pi in range(len(p_arr)):
                self.par_array[pi].freeze()

    def free_all(self):
        self.all_frozen = False
        for p_arr in self._parameters:
            for pi in range(len(p_arr)):
                self.par_array[pi].free()


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

    #def __str__(self):
    #    return str(self.show_pars())

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
        return [p.name for p in self.par_array]

    def _build_par_table(self,names_list=None):
        #, skip_hidden = False):
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

        #_val = np.array(_val, dtype=np.object)
        t=Table(_fields,names=_names,masked=False)

        self._fromat_column_entry(t)

        self._par_table= t


    def _fromat_column_entry(self, t):
        for n in self._numeric_fields:
            if n in t.colnames:

                try:
                    if None in t[n].data:

                        t[n] = MaskedColumn(t[n].data, name=n, dtype=np.float, mask=t[n].data==None)
                    else:
                        t[n] = MaskedColumn(t[n].data, name=n, dtype=np.float)
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
        # skip_hidden=False):

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

        #_val_start = np.array(_val_start, dtype=np.object)
        #_best_fit_val = np.array(_val_start, dtype=np.object)

        # set false to avoid string in hidden
        t = Table(_fields, names=_names, masked=False)

        self._fromat_column_entry(t)

        self._best_fit_par_table= t



    def __setattr__(self, name, value):
        if "properties" in self.__dict__ and name in self.properties:
            raise AttributeError('this member is protected, use del_par,add_par, to add/remove, and set() to set values or .val attribute')
        else:
            self.__dict__[name] = value


    def del_par(self,par):
        
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

    @property
    def par_table(self):
        self._build_par_table()
        return self._par_table

    @property
    def best_fit_par_table(self):
        self._build_best_fit_par_table()
        return self._best_fit_par_table

    def show_pars(self,getstring=False,names_list=None,sort_key=None):

        self._build_par_table(names_list=names_list)
        if sort_key is not None:
            self.par_table.sort(sort_key)

        if getstring==True:
            return self.par_table.pformat_all()
        else:
            _show_table(self.par_table)
            #return self.par_table
            #.pprint_all()
        
    
    
    def show_best_fit_pars(self,getstring=False):

        self._build_best_fit_par_table()
        if getstring == True:
            return self.best_fit_par_table.pformat_all()
        else:
            return self.best_fit_par_table
            #.pprint_all()


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
        
        #print "ModelParamterArray in model parameters",args,keywords
        
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
        
        #print "ModelParamterArray in model parameters",args,keywords
        
        par=self.get_par_by_name(par_name)
        
        return par.get(arg)
    
    def freeze_all(self):
        self.all_frozen=True
        for pi in range(len(self.par_array)):
            self.par_array[pi].freeze()
                
        
    def free_all(self):
        self.all_frozen = False
        for pi in range(len(self.par_array)):
            self.par_array[pi].free()
        

    def _serialize_pars(self):
        _par_keys=['val','val_min','val_max','val_start','val_last_call','fit_range_min','fit_range_max','best_fit_val',
                   'best_fit_err','frozen','allowed_values','_linked','_is_dependent','_func','_master_pars',
                   '_linked_root_model','_depending_pars','_root_par','','_master_par_list','_depending_par_expr','units','par_type']
        _par_dict = {}
        for par in self.par_array:
            _val_dict={}

            for k in _par_keys:
                if hasattr(par,k):
                    _val_dict[k]=getattr(par,k)

            _par_dict[par.name] = _val_dict

        return _par_dict

    def _decode_pars(self,_par_dict):
        for p_name in _par_dict.keys():
            self.set(p_name,**_par_dict[p_name])

