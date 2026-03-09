

__author__ = "Andrea Tramacere"

from .model_parameters import ModelParameterArray, ModelParameter
from .base_model import Model

from .spectral_shapes import SED

from numpy import log10,power,sqrt,shape,zeros

import numpy as np

from numpy import polyfit,polyval,polyder


__all__=['find_max_cubic','LogCubic','LogLinear','LogLogModel','LogParabolaEp','LogParabolaPL','PolyParameter']

class PolyParameter(ModelParameter):
    """
    This class is a subclass of the :class:`.ModelParameter` class,
    extending the base class to loglog polynomial parameters.
    """
    def __init__(self,polymodel,**keywords):
        
        """Create a new `PolyParameter` instance.
        
        Parameters
        ----------
        polymodel : object
            Parameter controlling polymodel.
        **keywords : dict
            Parameter controlling keywords.
        """
        self.polymodel=polymodel

        self.allowed_par_types=['curvature','peak freq','peak flux','third-degree','spectral-slope','flux-const','turn-over freq']
        

        super(PolyParameter,self).__init__(  **keywords)

        if 'val' in keywords.keys():
            val=keywords['val']
            self.assign_val(self.name,val)
            
        
        
    def set(self,**keywords):
        """Set.
        
        Parameters
        ----------
        **keywords : dict
            Parameter controlling keywords.
        """
        super(PolyParameter,self).set(**keywords )
        
        """
        sets a parameter value checking for physical boundaries 
        """
        if 'val' in keywords.keys():
            self.assign_val(self.name,keywords['val']) 
        
    
    def assign_val(self,name,val):

        """Assign val.
        
        Parameters
        ----------
        name : object
            Name identifier.
        val : object
            Value to assign.
        """
        setattr(self.polymodel,name,val)



class LogLogModel(Model):
    """Base class for spectral models defined in log-log frequency space.

    Notes
    -----
    Subclasses implement ``log_func``; this base class converts it to linear
    ``nuFnu`` values through :meth:`lin_func`.
    """
    
    def __init__(self,nu_size=100, **keywords):
       
        """Create a new `LogLogModel` instance.
        
        Parameters
        ----------
        nu_size : int, optional
            Frequency/energy control value for nu size.
        **keywords : dict
            Parameter controlling keywords.
        """
        super(LogLogModel,self).__init__(**keywords)
        self.model_type='LogLogModel'
    
    
    def lin_func(self,nu):
        """Lin func.
        
        Parameters
        ----------
        nu : object
            Frequency values in Hz.
        
        Returns
        -------
        object
            Computed result.
        """
        nu_log=np.log10(nu)
        return np.power(10.0,self.log_func(nu_log))
         
       
       

class LogLinear(LogLogModel):
    """
    Class to handle log-linear model
    """
    
    def __init__(self,nu_size=100,**keywords):
        """Create a new `LogLinear` instance.
        
        Parameters
        ----------
        nu_size : int, optional
            Number of points for frequency grids.
        **keywords : dict
            Additional keyword arguments.
        """
     
        super(LogLinear,self).__init__(  **keywords)
        
        self.name='LogLinear'
            
        self.SED = SED(name=self.model_type)
    
        self.parameters = ModelParameterArray()
      
        self.alpha=-1.0
        self.parameters.add_par(PolyParameter(self,name='alpha',par_type='spectral-slope',val=-1.0,val_min=-10.,val_max=10.,units=''))
        
       
        self.K=-10
        self.parameters.add_par(PolyParameter(self,name='K',par_type='flux-const',val=-10.0,val_min=None,val_max=None,units='erg cm^-2 s^-1',log=True))
        
         
        
    def log_func(self,log_nu):
    
            #x_log=log_nu
            #print(x_log)
            #print self.Ep,self.Sp,self.b,self.c,x_log
            
            """Log func.
            
            Parameters
            ----------
            log_nu : object
                Frequency/energy control value for log nu.
            
            Returns
            -------
            object
                Computed result.
            """
            return self.K + self.alpha*(log_nu)
    
    
    
    



class LogParabolaEp(LogLogModel):
    
    """
    Class to handle log-parabolic model
    """
    
    def __init__(self,nu_size=100, **keywords):
        """Create a new `LogParabolaEp` instance.
        
        Parameters
        ----------
        nu_size : int, optional
            Number of points for frequency grids.
        **keywords : dict
            Additional keyword arguments.
        """
        
        super(LogParabolaEp,self).__init__(  **keywords)
        
        self.name='LogParabolaEp'
            
        self.SED = SED(name=self.model_type)
    
        self.parameters = ModelParameterArray()
      
        self.b=-1.0
        self.parameters.add_par(PolyParameter(self,name='b',par_type='curvature',val=-1.0,val_min=-10.,val_max=10.,units=''))
        
        self.Ep=14.0
        self.parameters.add_par(PolyParameter(self,name='Ep',par_type='peak freq',val=14.0,val_min=0.,val_max=30.,units='Hz',log=True))
        
        self.Sp=-10
        self.parameters.add_par(PolyParameter(self,name='Sp',par_type='peak flux',val=-10.0,val_min=-30.,val_max=0.,units='erg cm^-2 s^-1',log=True))
        
    def log_func(self,log_nu):
    
            """Log func.
            
            Parameters
            ----------
            log_nu : object
                Frequency/energy control value for log nu.
            
            Returns
            -------
            object
                Computed result.
            """
            x_log=log_nu-self.Ep

            
            return self.Sp + self.b*(x_log*x_log)
    
    
   




class LogParabolaPL(LogLogModel):
    """
    Class to handle a log-par + pl model
    """
    
    def __init__(self,nu_size=100,**keywords):
        """Create a new `LogParabolaPL` instance.
        
        Parameters
        ----------
        nu_size : int, optional
            Number of points for frequency grids.
        **keywords : dict
            Additional keyword arguments.
        """
        
        super(LogParabolaPL,self).__init__(  **keywords)
        
        self.name='LogParabolaPL'
            
        self.SED = SED(name=self.model_type)
    
        self.parameters = ModelParameterArray()
      
        self.b=-1.0
        self.parameters.add_par(PolyParameter(self,name='b',par_type='curvature',val=-1.0,val_min=-10.,val_max=10.,units=''))

        self.E0=14.0
        self.parameters.add_par(PolyParameter(self,name='E0',par_type='turn-over freq',val=14.0,val_min=0.,val_max=30.,units='Hz',log=True))
        
        self.alpha=0.5
        self.parameters.add_par(PolyParameter(self,name='alpha',par_type='spectral-slope',val=0.5,val_min=-10.,val_max=10.,units=''))

        self.K=-10
        self.parameters.add_par(PolyParameter(self,name='K',par_type='flux-const',val=-10.0,val_min=-30.,val_max=0.,units='erg cm^-2 s^-1',log=True))
        
    
    def log_func(self,log_nu):
        
        
        #TODO fix this to numpy array function
        """Log func.
        
        Parameters
        ----------
        log_nu : object
            Frequency/energy control value for log nu.
        
        Returns
        -------
        object
            Computed result.
        """
        if shape(log_nu)==():
            return self.composite_func(log_nu)

        else:
            y_log=zeros(log_nu.size)
            for i in range(log_nu.size):
                y_log[i]=self.composite_func(log_nu[i])

            return y_log
   
    def composite_func(self,log_nu):     
        
        """Composite func.
        
        Parameters
        ----------
        log_nu : object
            Frequency/energy control value for log nu.
        
        Returns
        -------
        object
            Computed result.
        """
        x_log=log_nu-self.E0
        
        if log_nu >=  self.E0:
            
            return self.K + x_log*(self.alpha + self.b*x_log)

        else:    

            return  self.K + x_log*self.alpha 



    


class LogCubic(LogLogModel):
    
    """
    Class to handle log-cubic model
    """
    
    def __init__(self,nu_size=100, **keywords):
        """Create a new `LogCubic` instance.
        
        Parameters
        ----------
        nu_size : int, optional
            Number of points for frequency grids.
        **keywords : dict
            Additional keyword arguments.
        """
        
        
        super(LogCubic,self).__init__(  **keywords)
        
        
        
        self.name='LogCubic'
            
        self.SED = SED(name=self.model_type)
    
        self.parameters = ModelParameterArray()
      
        self.b=-1.0
        self.parameters.add_par(PolyParameter(self,name='b',par_type='curvature',val=-1.0,val_min=-10.,val_max=10.,units=''))
        
        self.c=0.0
        self.parameters.add_par(PolyParameter(self,name='c',par_type='third-degree',val=-1.0,val_min=-10.,val_max=10.,units=''))

        self.Ep=14.0
        self.parameters.add_par(PolyParameter(self,name='Ep',par_type='peak freq',val=14.0,val_min=0.,val_max=30.,units='Hz',log=True))
        
        self.Sp=-10
        self.parameters.add_par(PolyParameter(self,name='Sp',par_type='peak flux',val=-10.0,val_min=-30.,val_max=0.,units='erg cm^-2 s^-1',log=True))
        
        

    def log_func(self,log_nu):
    
    
            """Log func.
            
            Parameters
            ----------
            log_nu : object
                Frequency/energy control value for log nu.
            
            Returns
            -------
            object
                Computed result.
            """
            x_log=log_nu-self.Ep
            
            #print self.Ep,self.Sp,self.b,self.c,x_log
            
            return self.Sp + self.b*(x_log*x_log)+ self.c*(x_log*x_log*x_log)
            
    
    

def find_max_cubic(x_log,y_log,x_range=None):
    
    """Find max cubic.
    
    Parameters
    ----------
    x_log : object
        Parameter controlling x log.
    y_log : object
        Parameter controlling y log.
    x_range : object, optional
        Range for x.
    
    Returns
    -------
    object
        Computed result.
    """
    if x_range is not None:
        msk = x_log >x_range[0]
        msk*= x_log <x_range[1]
        
        x_log=x_log[msk]
        y_log=y_log[msk]
    
    p=polyfit(x_log,y_log,3)
    
    der=polyder(p,1)
    
    delta=der[1]*der[1] -(4*der[0]*der[2])
    if delta>=0:
        #print p,der
        root1=(-der[1]+sqrt(delta))/(2*der[0])
        root2=(-der[1]-sqrt(delta))/(2*der[0])
        #print root1,root2
        if der[0]>0:
            xp=min(root1,root2)
        else:
            xp=max(root1,root2)
        

        return  xp
    else:
        print("!!! no maxima found for cubic fit")

        return None
