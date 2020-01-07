from __future__ import absolute_import, division, print_function

from builtins import (bytes, str, open, super, range,
                      zip, round, input, int, pow, object, map, zip)

import  numpy as np
from astropy.table import Table
from astropy import  units as u

__all__=['ModelParameter','ModelParameterArray','Value']

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

        """
        Constructor
        """
        self.allowed_keywords={'name':'unk'}
        self.allowed_keywords['val']=None
        self.allowed_keywords['par_type']='unk'
        self.allowed_keywords['units']='No'
        self.allowed_keywords['val_min']=None
        self.allowed_keywords['val_max']=None
        self.allowed_keywords['val_start_fit']=None
        self.allowed_keywords['val_last_call']=None
        self.allowed_keywords['fit_range_min']=None
        self.allowed_keywords['fit_range_max']=None
        self.allowed_keywords['fit_range']=None
        self.allowed_keywords['best_fit_val']=None
        self.allowed_keywords['best_fit_err']=None
        self.allowed_keywords['frozen']=False
        self.allowed_keywords['log']=False
        self.allowed_keywords['allowed_values'] = None
      
        #self._skip_kw=['val']
        #defualt

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
            else:
                setattr(self,kw,self.allowed_keywords[kw])




        #parsing user keywords
        self.set(**keywords)
        #self._skip_kw = []
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
        #print "ModelParamter in model parameters",args,keywords
        
        keys = keywords.keys()

        for kw in keys:
            
            if kw in  self.allowed_keywords.keys() :
                if kw == 'val':
                    #print('->',self.allowed_values )
                    if self.allowed_values is not None:

                        if keywords[kw]  not in self.allowed_values:
                            raise RuntimeError('parameter  %s' %(self.name), 'the value', keywords[kw] , 'is not in the allowed list',self.allowed_values)

                    self._val.val = keywords[kw]

                elif kw == 'log':
                    self._val.islog = keywords[kw]

                elif kw== 'units':

                    self._val.units = keywords[kw]
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

    def __init__(self):
            
        """
        Constructor 
        """
        
        self.par_array=[]
        self.all_frozen=False

        self.properties={}

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

        setattr(self,par.name, par)
        self.properties[par.name]=par


    def _build_par_table(self,names_list=None):
        _name=[]
        _type=[]
        _unit=[]
        _val=[]
        _bound_min=[]
        _bound_max=[]
        _islog=[]
        _frozen=[]
        _fields=[_name,_type,_unit,_val,_bound_min,_bound_max,_islog,_frozen]
        _names=['name','par type','units','val','phys. bound. min','phys. bound. max','log','frozen']
        for par in self.par_array:

            append=False

            if names_list is not None:
                if par.name in names_list:
                    append=True
            else:
                append = True

            if append:
                _name.append(par.name)
                _type.append(par.par_type)
                _unit.append(par.units)
                _val.append(par.val)
                _bound_min.append(par.val_min)
                _bound_max.append(par.val_max)
                _islog.append(par.islog)
                _frozen.append(par.frozen)
        #print(len(_fields),len(_names))

        self._par_table= Table(_fields,names=_names)

    
    def _build_best_fit_par_table(self,names_list=None):
        _name=[]
        _best_fit_val=[]
        _best_fit_err_p=[]
        _best_fit_err_m=[]
        _val_start=[]
        _fit_range_min=[]
        _fit_range_max=[]
        _frozen=[]
        _fields=[_name,_best_fit_val,_best_fit_err_p,_best_fit_err_m,_val_start,_fit_range_min,_fit_range_max,_frozen]
        _names=['name','bestfit val','err +','err -','start val','fit range min','fit range max','frozen']

        for par in self.par_array:

            append=False

            if names_list is not None:
                if par.name in names_list:
                    append=True
            else:
                append = True

            if append:
                _name.append(par.name)

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

                _frozen.append(par.frozen)

        #print(len(_fields),len(_names))
        #for ID,n in enumerate(_names):
        #    print(ID,_fields[ID],n)

        self._best_fit_par_table= Table(_fields,names=_names)
    


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
                
        
    @property
    def par_table(self):
        self._build_par_table()
        return self._par_table

    def show_pars(self,getstring=False,names_list=None,sort_key=None):
        """
        shows the information for all the items in the :attr:`~ModelParameterArray.par_array`
        
        """
        
        #text=[]
        #text.append( "-------------------------------------------------------------------------------------------------------------------------")
        #text.append( "model parameters:")
        #text.append( " Name             | Type                 | Units            | value         | phys. boundaries              | log | frozen")
        #text.append( "-------------------------------------------------------------------------------------------------------------------------")
       
        #for par in self.par_array:
        #    if names_list is not None:
        #        if par.name in names_list:
        #            text.append( par.get_description(nofields=True))
        #    else:
        #        text.append(par.get_description(nofields=True))
        #
        #text.append( "-------------------------------------------------------------------------------------------------------------------------")

        self._build_par_table(names_list)
        if sort_key is not None:
            self.par_table.sort(sort_key)

        if getstring==True:
            return self.par_table.pformat_all()
        else:
            self.par_table.pprint_all()
        
    
    
    def show_best_fit_pars(self,getstring=False):
        """
        shows the best-fit information for all the items in the :py:attr:`par_array`
        """
        
        #text=[]
        #text.append("-------------------------------------------------------------------------------------------------------------------")
        #text.append("best-fit parameters:")
        #text.append("  Name            | best-fit value| best-fit err +| best-fit err -|start value   | fit boundaries")
        #text.append("-------------------------------------------------------------------------------------------------------------------")

       
        #for par in self.par_array:
        #    text.append( par.get_bestfit_description(nofields=True))

        #text.append("-------------------------------------------------------------------------------------------------------------------")
        

        #if getstring==True:
        #    return text
        #else:
        #    for line in text:
        #        print (line)

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
        
        