from __future__ import absolute_import, division, print_function

from builtins import (bytes, str, open, super, range,
                      zip, round, input, int, pow, object, map, zip)
import ast
import  numpy as np
import warnings
from astropy.table import Table,vstack,MaskedColumn
from astropy import  units as u
from .utils import clean_var_name
import  copy
__all__=['ModelParameter','ModelParameterArray','Value','LinkedParameter']





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
    def units(self, units):
        try:
            self._units = u.Unit(units)
            #print(units,type(self._units))
        except Exception as e:
            #print('*',units)
            self._units = units + '*'


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
      


        _v = None
        _l = False
        _units= None
        
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

        self._linked = False
        self._linked_models = []
        self._linked_root_model = None
        self._root = False
        self.set(**keywords)

    @property
    def root(self):
        return self._root

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



    def set(self, *args, **keywords):
        """
        sets a parameter value checking for physical boundaries
        
        Parameters: keywords of the constructor
         
        """

        keys = keywords.keys()

        for kw in keys:
            
            if kw in  self.allowed_keywords.keys() :
                if kw == 'val':

                    if self.allowed_values is not None:

                        if keywords[kw] not in self.allowed_values:
                            raise RuntimeError('parameter  %s' %(self.name), 'the value', keywords[kw] , 'is not in the allowed list',self.allowed_values)

                    self._val.val = keywords[kw]

                elif kw == 'log':
                    self._val.islog = keywords[kw]

                elif kw== 'units':

                    self._val.units = keywords[kw]

                elif kw == 'par_type':

                    if self.allowed_par_types is not None:

                        if  keywords['kw'] not in self.allowed_par_types:
                            msg = "parameter%s the  type %s is not allowed" % (self.name,keywords['kw']) + "\n please choose among %s" % self.allowed_par_types
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
            
    
    
    
    def freeze(self):
        """
        freezes a paramter 
        """
        self.frozen=True

    
        
    def free(self):
        """
        make a paramter free
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
            descr= "name = %-16s  type = %-20s  units = %-16s  val = %s  phys-bounds = [%-13s,%-13s] islog = %s  froze= %s "%(self.name, self.par_type ,
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

class LinkedParameter(ModelParameter):

    def __init__(self,p_name,m_list):
        super(LinkedParameter,self).__init__()
        self.m_list = m_list
        self.name = p_name

    #def set_par(self, par_name, val):
    #    for model in self.m_list:

    #        pm = model.get_par_by_name(self.name)
    #        pm.set(par_name, val = val)

    #def set(self, *args, **kw):
    #    self.set_par(*args, **kw)

class ModelLinkedParameter(object):
    def __init__(self,name,p):

        self.name = name
        self.paramters = ModelParameterArray()
        self.paramters.add_par(p)
        self.paramters.model=self

class CompositeModelParameterArray(object):

    def __repr__(self):
        str(self.show_pars())

    def __init__(self,):

        self.all_frozen = False
        self._parameters = []
        self._model_comp = []




    def link_par(self,par_name,model_name_list,root_model_name):
        #m_name=''
        m_root=self.get_model_by_name(root_model_name)
        p_root=m_root.get_par_by_name(par_name)
        p_root._root = True
        p_root._linked = True
        p_root._linked_root_model = m_root
        for m_name in model_name_list:
            m=self.get_model_by_name(m_name)
            p_root._linked_models.append(m)
            if m is not None:
                p = m.get_par_by_name(par_name)
                m.parameters.del_par(p)
                m.parameters.add_par(p_root)
            else:
                warnings.warn('model with',m_name,'not found in composite model')
            #print('del -->',m.name, p ,p.name)
            #self.hide_model_parameter(m.name, p.name)
            #p._linked = True
            #print('m_name -->', m_name)

        #m_name = m_name[:-1]
        #print('m name. -->', m_name)
        #p=LinkedParameter(p_name, m_list)
        #p=copy.deepcopy(p)
        #p._linked = False
        #m=ModelLinkedParameter(m_name,p)
        #self._model_comp.append(m)
        #self._parameters.append(m.paramters)

    def add_model_parameters(self, model):
        try:
            assert (model.name not in [m.name for m in self._model_comp])
        except:
            raise RuntimeError('model name:', model.name, 'already assigned')


        model.parameters.model=model

        self._parameters.append(model.parameters)
        self._model_comp.append(model)

    def del_model_parameters(self, model_name):
        m,ID=self.get_model_by_name(model_name,get_idx=True)
        if m is not None:
            #_p = self._comp_par_array.pop(ID)
            _p = self._model_comp.pop(ID)
            _p = self._parameters.pop(ID)
        else:
            warnings.warn('Model', model_name, 'not present')

    #def del_model_parameter(self, model_name,par_name):
    #    m, ID = self.get_model_by_name(model_name, get_idx=True)

    #    p = m.get_par_by_name(par_name)

    #    if len(self._parameters[ID].par_array)<1:
    #        _ = self._model_comp.pop(ID)



    def get_model_by_name(self, model_name,get_idx=False):
        selected_model = None
        idx=None
        for ID,m in enumerate(self._model_comp):
            if m.name == model_name:
                selected_model = m
                idx=ID


        if selected_model is None:
            warnings.warn('Model',model_name,'not present')
        if get_idx is False:
            return selected_model
        else:
            return selected_model,idx

    def get_par_by_name(self, model_name,par_name):
        p=None
        m=self.get_model_by_name(model_name)

        if m is None:
            warnings.warn('Model', model_name, 'not present')
        else:
            p=m.parameters.get_par_by_name(par_name)

        return p


    def _build_par_table(self):
        _l=[]
        for ID, mc in enumerate(self._model_comp):
            self._parameters[ID]._build_par_table()
            if len(self._parameters[ID]._par_table) >= 1:
                _l.append( self._parameters[ID]._par_table)
        #print('--> _l',_l)
        self._par_table=vstack(_l)

    def _build_best_fit_par_table(self):
        _l = []
        for ID, mc in enumerate(self._model_comp):
            self._parameters[ID]._build_best_fit_par_table()
            if len(self._parameters[ID]._best_fit_par_table) >= 1:
                _l.append( self._parameters[ID]._best_fit_par_table)
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
            self.par_table.pprint_all()

    def show_best_fit_pars(self, getstring=False):
        self._build_best_fit_par_table()
        if getstring == True:
            return self._best_fit_par_table.pformat_all()
        else:
            self._best_fit_par_table.pprint_all()

    def freeze(self, model_name,par_name):
        self.set(model_name,par_name, 'frozen')

    def free(self,model_name, par_name):
        self.set(model_name,par_name, 'free')

    def set(self, model_name, par_name, *args, **kw):
        m=self.get_model_by_name(model_name)
        if m is not None:
            m.parameters.set(par_name, *args, **kw)
        else:
            warnings.warn('Model', model_name, 'not present')

    def set_par(self,model_name, par_name, val):
        m=self.get_model_by_name(model_name)
        if m is not None:
            m.parameters.set(par_name, val=val)
        else:
            warnings.warn('Model', model_name, 'not present')

    def get(self, model_name,par_name, *args):
        m = self.get_model_by_name(model_name)
        if m is not None:
            pass
        else:
            warnings.warn('Model', model_name, 'not present')
        return m.parameters.get(par_name, *args)

    def get_val(self, model_name,par_name):
        m = self.get_model_by_name(model_name)
        if m is not None:
            pass
        else:
            warnings.warn('Model', model_name, 'not present')

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
        str(self.show_pars())

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


    def add_par(self,par):
        """
        adds a new :class:`ModelParameter` object  to the `par_array`
        """
        try:
            assert (isinstance(par,ModelParameter))
        except:
            raise RuntimeError('parameter is not an istance of',type(ModelParameter))

        try:
            assert (par.name not in [p.name for p in self.par_array])
        except:
            raise RuntimeError('parameter name:',par.name,'already assigned')

        self.par_array.append(par)

        setattr(self,clean_var_name(par.name), par)
        self.properties[par.name]=par


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
                    if par._linked is True and par._linked_root_model == self.model:
                        _p_name = par.name + '(R)'
                    else:
                        _p_name = par.name
                    _type.append(par.par_type)
                    _unit.append(par.units)
                    _val.append(par.val)
                    _bound_min.append(par.val_min)
                    _bound_max.append(par.val_max)
                    _islog.append(par.islog)
                    _frozen.append(par.frozen)

                _name.append(_p_name)

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
        _fields=[_model_name,_name,_best_fit_val,_best_fit_err_p,_best_fit_err_m,_val_start,_fit_range_min,_fit_range_max,_frozen]
        _names=['model name','name','bestfit val','err +','err -','start val','fit range min','fit range max','frozen']

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
                    if par._linked is True and par._linked_root_model == self.model:
                        _p_name = par.name + '(R)'
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

                _name.append(_p_name)
                _frozen.append(par.frozen)



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
            self.par_table.pprint_all()
        
    
    
    def show_best_fit_pars(self,getstring=False):

        self._build_best_fit_par_table()
        if getstring == True:
            return self._best_fit_par_table.pformat_all()
        else:
            self._best_fit_par_table.pprint_all()


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
        _par_keys=['val','val_min','val_max','val_start','val_last_call','fit_range_min','fit_range_max','best_fit_val','best_fit_err','frozen','allowed_values']
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