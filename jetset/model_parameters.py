from __future__ import absolute_import, division, print_function

from builtins import (bytes, str, open, super, range,
                      zip, round, input, int, pow, object, map, zip)

import  numpy as np
"""
Module: model_parameters
===================================================================

This module contains the base classes necessary to handle models paramters/



Classes and Inheritance Structure
----------------------------------------------
.. inheritance-diagram:: BlazarSEDFit.model_parameters
    
Classes relations
----------------------------------------------

.. figure::  classes_model_parameters.png
   :align:   center    


      
.. autosummary::
   ModelParameter
   ModelParameterArray

Module API
-------------------------------------------------------------------

"""


__all__=['ModelParameter','ModelParameterArray','Value']

class Value(object):

    def __init__(self,val,islog=False):
        self._val=val
        self._islog=islog

    def __repr__(self):
        return '%s'%self._val

    @property
    def islog(self):
        return  self._islog

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
        
      
        #self._skip_kw=['val']
        #defualt

        _v = None
        _l = False

        if 'val' in keywords.keys():
            _v = keywords['val']

        if 'log' in keywords.keys():
            _l = keywords['log']

        self._val = Value(val=_v, islog=_l)

        for kw in self.allowed_keywords.keys():
            if kw == 'val':
                self.val = keywords[kw]
            if kw == 'log':
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
        self._val.val=val

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
                    self.val = keywords[kw]
                if kw == 'log':
                    pass
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
            descr= "name = %-16s  type = %-20s  units = %-16s  val = %s  phys-bounds = [%-13s,%-13s] islog = %s   "%(self.name, self.par_type ,
                                                                                   self.units, val, val_min, val_max,self.islog)
        else:
            descr= " %-16s | %-20s | %-16s | %s | [%-13s,%-13s] | %s "%(self.name, self.par_type ,
                                                                                   self.units, val, val_min, val_max,self.islog)
            
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
        
        if self.best_fit_err is None:
            best_fit_err='No'
        else:
            best_fit_err=str("%+e"%self.best_fit_err)
        
        if self.frozen==True:
            best_fit_val='Frozen'
            best_fit_err='Frozen'
       
        if nofields==False:
            descr= "name = %-16s best-fit val=%-13s  best-fit err=%-13s start-val=%-13s   fit-bounds=[%-13s,%-13s]"%(self.name,best_fit_val,best_fit_err,val_start,fit_range_min,fit_range_max)
        else:
            descr= " %-16s | %-13s | %-13s | %-13s | [%-13s,%-13s]"%(self.name,best_fit_val,best_fit_err,val_start,fit_range_min,fit_range_max)

        return descr
        

class ModelParameterArray(object):
    """
    This class provide and interface to handle an array of :class:`ModelParameter` objects.
    
    
    Attributes
    ----------
    par_array :  list 
        list of :class:`ModelParameter` objects
    """   
    
    def __init__(self):
            
        """
        Constructor 
        """
        
        self.par_array=[]
        self.all_frozen=False
        
            
    def add_par(self,par):
        """
        adds a new :class:`ModelParameter` object  to the `par_array`
        """
        self.par_array.append(par)
    
    def del_par(self,par):
        
        self.par_array.remove(par)
        
    
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
                
        
    
    
    def show_pars(self,getstring=False):
        """
        shows the information for all the items in the :attr:`~ModelParameterArray.par_array`
        
        """
        
        text=[]
        text.append( "-------------------------------------------------------------------------------------------------------------------")
        text.append( "model parameters:")
        text.append( " Name             | Type                 | Units            | value         | phys. boundaries              |log ")
        text.append( "-------------------------------------------------------------------------------------------------------------------")
       
        for par in self.par_array:
            text.append( par.get_description(nofields=True))
        
        text.append( "-------------------------------------------------------------------------------------------------------------------")
        
        

        if getstring==True:
            return text
        else:
            for line in text:
                print (line)
        
    
    
    def show_best_fit_pars(self,getstring=False):
        """
        shows the best-fit information for all the items in the :py:attr:`par_array`
        """
        
        text=[]
        text.append("---------------------------------------------------------------------------------------------------")
        text.append("best-fit parameters:")
        text.append("  Name            | best-fit value| best-fit err  | start value   | fit boundaries")
        text.append("---------------------------------------------------------------------------------------------------")

       
        for par in self.par_array:
            text.append( par.get_bestfit_description(nofields=True))

        text.append("---------------------------------------------------------------------------------------------------")
        
    
        if getstring==True:
            return text
        else:
            for line in text:
                print (line)
        
        
    
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
                
        

        
        