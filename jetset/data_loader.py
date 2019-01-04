"""
Moudule: data_loader
===================================================================

This module contains all the classes necessary to load SED data
from a file.  The most effective way to import the SED data is to create an object 
instance of :class:`ObsData` class.  

Classes and Inheritance Structure
-------------------------------------------------------------------

.. inheritance-diagram:: BlazarSEDFit.data_loader
   


Classes relations
----------------------------------------------

.. figure::  classes_data_loader.png    
   :align:   center     


Summary
---------
.. autosummary::
    ObsData
    

Module API
-------------------------------------------------------------------
"""


from __future__ import absolute_import, division, print_function

from builtins import (bytes, str, open, super, range,
                      zip, round, input, int, pow, object, map, zip)

__author__ = "Andrea Tramacere"


import numpy as np
from .cosmo_tools import Cosmo
#from poly_fit import filter_interval
from astropy.table  import  Table

from .output import section_separator
import os

__all__=['get_data_set_msk','get_freq_range_msk','lin_to_log','log_to_lin','ObsData']




class ObsData(object):
    """
    This class provides a  powerful interface  to load observational data stored in a file.
    The following parameters set the corresponding class members: 
 
    :param data_file: (str), path to the data file
    
    :param  z: (float), redisfhit, if not provided it is looked up in the  file meta-data
        
    :param obj_name: (str), if not provided  it is looked up in the  file meta-data
        
    :param restframe: (str) restframe of the data, possible values are ``'src'`` or  ``'obs'``, 
        if not provided it is looked up in the  file meta-data
        
    :param col_types: (str), string specifying the way data are organized in columns in the input file. 
        This string is passed to the method  :func:`set_data_cols`. If not provided
        it is looked up in the  file meta-data. See the corresponding :func:`set_data_cols` 
        documentation for details.
    
    :param col_nums: (str), string describing the corresponding column data number
        
    :param data_scale: (str) ``'lin-lin'`` or ``'log-log'``,it is looked up in the  file meta-data. 
        This parameter allows to specify      if the data in the file are stored as in log-log or lin-lin scale
        
    :param data_set_filter:  a filter to filter the SED data according to the data_set value,  eg: 
        ``'data_set_filter='mw-1'``, will filter data with ``data_set=='mw-1'``,        
        ``'data_set_filter=['mw-1','mw-2']'`` will filter data with ``data_set=='mw-1'`` or ``data_set=='mw-2'``
    
        
        .. note::
            
            Historical data will be only plotted but will not be used in the fit.
    
   
    **Class Members storing the SED data:**

     data used to fit 
     
     - nu_data
     - nuFnu_data 
     - dnu_data
     - dnuFnu_data   
     - T_data 
     - data_set 
     
     upper limits
     
     - nu_data_UL
     - UL_data
   

    
    
    The private method   :func:`_load_data`  populates the data members starting from the data in the file. Data are properly transformed according
    to the ``restframe`` and ``data_scale`` values.
    
      
    
    **Examples**
    

    The following lines shows and example of how to embed **meta-data** in the header of the SED data file,
    just adding in the header of the file a line starting with ``#`` and followed by the 
    indetifier ``md``, the meta-data name and value :
    
    .. literalinclude:: ../../../BlazarSEDFit/test_data/SEDs_data/SED_MW_Mrk421.dat 
       :lines: 1-20
        
    
    Assuming that the data file path has been stored in ``SED_file=/path/to/file/file.txt``, 
    the data can be imported as follows:    
    
    .. code::
    
        from BlazarSEDFit.data_loader import ObsData
        mySEDdata=ObsData(data_table=SED_file)
    
    that is completely equivalent to:
    
    .. code::
    
        mySEDdata=ObsData(data_table=SED_file,col_types='x,y,dy,data_set',z=0.0308,data_scale='lin-lin')
        
    
    
    """
    
    def __init__(self, data_table=None,dupl_filter=False, data_set_filter=None,UL_filtering=False,UL_value=None,**keywords):
        
        """
            
        """
        
        self.z=None
        self.data=None
        self.obj_name=None
        self.restframe=None
        self.col_types=None
        self.col_nums=None
        self.data_table=data_table
        self.data_scale=None

        if data_set_filter is not None:
            self.data_set_filter=data_set_filter.split(',')
        else:
            self.data_set_filter=['No']
            
        self.cosmo_eval=Cosmo(units='cm')
        self.UL_value=UL_value
        self.UL_filtering=UL_filtering
        self.zero_error_replacment=0.2
        self.facke_error=0.2


        
        
        #------------------------------
        #builds a dictionary to bounds
        #keyword to class members
        #and to set allowed values
        #------------------------------
        allowed_keywords={'z':None}
        allowed_keywords['obj_name']=None
        allowed_keywords['restframe']=['obs','src']
        allowed_keywords['data_scale']=['lin-lin','log-log']
        allowed_keywords['data_table']=None
        allowed_keywords['col_types']=None
        allowed_keywords['col_nums']=None

        #loops over keywords
        #and set values
        #values overwrite values from file metadata
        keys = sorted(keywords.keys())
        #print keys
        for kw in keys:
            #print "kw ",kw,keywords[kw]
            if kw in  allowed_keywords.keys():
              
                if allowed_keywords[kw] is not None:
                    #check that the kyword value is correct
                    if keywords[kw] not in allowed_keywords[kw]:                         
                        print ("keyword=%s, has wrong value=%s, allowed are %s"%(kw,keywords[kw], allowed_keywords[kw]))
                        raise ValueError
                        
                    
                setattr(self,kw,keywords[kw])
                    
            else:
                
                print ("wrong keyword=%s, not in%s "%(kw, allowed_keywords.keys()))
                
                raise ValueError

        #print('col_types a', self.col_types)
        if self.data_table is None:

            print("you must provide a valid path for the astropy Table or  an astropy Table object")

            return
        else:
            self._load_data(self.data_table)

        self.set_md()

        print(section_separator)

        #self.set_data_cols(self.col_types,self.col_nums)

        #print('col_types b', self.col_types)
       
        
        
        #calls method to load data
        self._build_data(dupl_filter=dupl_filter)




    def _build_data(self,dupl_filter=False):
        """
        private method to load and build the SED data     
        
        :param dupl_filter: keyword to perfrom filtering of duplicated entries
        :type dupl_filter: bool 
        :ivar dt: numpy dtype  for the SED data
         
        
        
        **Data Processing**
        
        - checks if data errors are provided 
        - separates historical from simultaneosu (used in the fit) data
        - filters upper limits 
        - removes duplicate entries
        - performs restframe transformation
        - performs `lin-lin`, `log-log` transformations
        

       
       
        """

        sed_dt = [('nu_data', 'f8')]
        sed_dt.append(('dnu_data', 'f8'))
        sed_dt.append(('nuFnu_data', 'f8'))
        sed_dt.append(('dnuFnu_data', 'f8'))

        sed_dt.append(('nu_data_log', 'f8'))
        sed_dt.append(('dnu_data_log', 'f8'))
        sed_dt.append(('nuFnu_data_log', 'f8'))
        sed_dt.append(('dnuFnu_data_log', 'f8'))

        sed_dt.append(('dnuFnu_facke', 'f8'))
        sed_dt.append(('dnuFnu_facke_log', 'f8'))

        sed_dt.append(('UL', 'bool'))
        sed_dt.append(('zero_error', 'bool'))
        sed_dt.append(('T_start', 'f8'))
        sed_dt.append(('T_stop', 'f8'))
        sed_dt.append(('data_set', 'S16'))
        # print ('ciccio',sed_dt)
        # _sed_dt=[]
        # for a in sed_dt:
        #    s=(str(a[0]),str(a[1]) )
        #    _sed_dt.append(s)
        # print('ciccio', _sed_dt)

        self.dt = np.dtype(sed_dt)

        #self.data = None
        self.data = np.zeros(len(self.data_table), dtype=self.dt)

        self._col_dict = {'x': 'nu_data'}
        self._col_dict['y'] = 'nuFnu_data'
        self._col_dict['dx'] = 'dnu_data'
        self._col_dict['dy'] = 'dnuFnu_data'
        self._col_dict['T_start'] = 'T_start'
        self._col_dict['T_stop'] = 'T_stop'
        self._col_dict['data_set'] = 'data_set'
        self._col_dict['UL'] = 'UL'

        self._log_col_dict = {'x': 'nu_data_log'}
        self._log_col_dict['y'] = 'nuFnu_data_log'
        self._log_col_dict['dx'] = 'dnu_data_log'
        self._log_col_dict['dy'] = 'dnuFnu_data_log'
        self._log_col_dict['T_start'] = 'T_start'
        self._log_col_dict['T_stop'] = 'T_stop'
        self._log_col_dict['UL'] = 'UL'
        
        if self.data_scale=='lin-lin':
            self._file_col_dict= self._col_dict
            
        
        elif self.data_scale=='log-log':
             self._file_col_dict= self._log_col_dict
        
        else:
            raise RuntimeError('data_scale not specified')

        self.col_types = []
        for n in self.data_table.colnames:
            self.col_types.append(n)
            #print ('-->n',self._file_col_dict,n,self._file_col_dict[n])
            self.data[self._file_col_dict[n]]=self.data_table[n]

        self.col_nums=len(self.col_types)

        print  (section_separator)
        
        
       

        print ("*** loading data ***")
        print ("---> loading data for file=%s" % self.data_table)
        print ("---> found these col ID=%s and names=%s:"%(self.col_nums,self.col_types))
        print ("---> z=%e"%self.z)
        print ("---> restframe=%s"%self.restframe)
        print ("---> obj_name=%s "%self.obj_name)
        print ("---> data_scale=%s "%self.data_scale)
        

        print('col_types',self.col_types)
       
        

        print("---> data len=%d" % len(self.data))

        #-------------------------------------------------------------------------
        # duplicate entries 
        if dupl_filter==True:
            print ("---> filtering for dupl entries")
            self.remove_dupl_entries(self.data)
       
        #-------------------------------------------------------------------------          

        # sets UL entries
        self.set_UL(self.UL_value)

        # set error if not present:
        self.set_facke_error(self.facke_error)





        self.set_zero_error()

        if self.UL_filtering==True:
            self.filter_UL(self.UL_value)
        
        
        if self.data_set_filter!=['No']:
            self.filter_data_set( self.data_set_filter)
       
        #print "3,",self.data['nu_data_log']
        #self.set_error(self.zero_error_replacment, data_msk=self.data['zero_error'])
        #print "4,",self.data['nu_data_log']
        self._set_data_frame_and_scale()
        #print ("5,",self.data['nu_data_log'])
        
       
        
        if self.data['dnuFnu_data'] is None and self.data_scale=='lin-lin':
        
            self.data['dnuFnu_data']=self.data['Fnu_data']*self.facke_error
            
            print ("Warning: error were not provided ")
            print ("         assigning %f            ")
            print ("         set error with .set_error"%self.facke_error)
            #print self.data['dnuFnu_data']
        
        
        if self.data['dnuFnu_data_log'] is None and self.data_scale=='log-log':
        
            self.data['dnuFnu_data_log']=np.ones(self.data['nu_data_log'].size)*self.facke_error
            
            print ("Warning: error were not provided ")
            print ("         assigning %f           ")
            print ("         set error with .set_error"%self.facke_error)
            #print self.data['dnuFnu_data_log']

    
        
        #print self.data_set_filter
        print("---> final data len",self.data.size)
        print  (section_separator)
    
    
    def _load_data(self,data_table):
        if isinstance(data_table, Table):
            self.data_table = data_table
        else:
            pass

        """
        method to load the data from the file 
        """

        self.data_table=Table.read(self.data_table,format='ascii.ecsv')




        

    def _set_data_frame_and_scale(self):
        
        
        #-------------------------------------------------------------------------          
        # handles data cosmological conversion
        # and axis transformation

        if self.restframe=='obs':
            nu_conv_factor=1
            Lum_conv_factor=1
        
        elif self.restframe=='src':
            DL=self.cosmo_eval.DL(self.z)
            print ("--->DL=%e"%DL)
            
            #!! usa le funzioni in frame_converter
            nu_conv_factor=1.0/(1+self.z)
            Lum_conv_factor=1.0/(np.pi*4.0*DL*DL)

        
        #conv_fac=np.log(10)
        
        
  
        #msk=np.invert(self.data['UL'])
        #*np.invert(self.data['zero_error'])
        if self.data_scale=='lin-lin':
            
            if 'dy' in self.col_types :
                self.data['nuFnu_data_log'],self.data['dnuFnu_data_log']=self.lin_to_log(val=self.data['nuFnu_data'], err=self.data['dnuFnu_data'])
            else:
                self.data['nuFnu_data_log']=self.lin_to_log(val=self.data['nuFnu_data'])
                
            if 'dx' in self.col_types:
                self.data['nu_data_log'],self.data['dnu_data_log']=self.lin_to_log(val=self.data['nu_data'], err=self.data['dnu_data'])
            else:
                self.data['nu_data_log']=self.lin_to_log(val=self.data['nu_data'])
            
            self.data['nu_data_log']+= np.log10(nu_conv_factor)
            self.data['nuFnu_data_log']+= np.log10(Lum_conv_factor)
            
            
            
        if self.data_scale=='log-log':
            #print "AAA"
            #print self.data['nuFnu_data_log'],self.data['dnuFnu_data_log']
            if 'dy' in self.col_types:
                self.data['nuFnu_data'],self.data['dnuFnu_data']=self.log_to_lin(log_val=self.data['nuFnu_data_log'], log_err=self.data['dnuFnu_data_log'])
            else:
                self.data['nuFnu_data']=self.log_to_lin(log_val=self.data['nuFnu_data_log'])
                
            if 'dx' in self.col_types:
                self.data['nu_data'],self.data['dnu_data']=self.log_to_lin(log_val=self.data['nu_data_log'], log_err=self.data['dnu_data_log'])
            else:
                self.data['nu_data']=self.log_to_lin(log_val=self.data['nu_data_log'])
            
            
            #print self.data['nuFnu_data'],self.data['dnuFnu_data']
            self.data['nu_data']*=nu_conv_factor
            self.data['nuFnu_data']*=Lum_conv_factor
        
        

        
        
    def set_md(self):


        md_dic={}
        md_dic['z']=None
        md_dic['file_name']=None
        md_dic['obj_name']=None
        md_dic['restframe']=None
        md_dic['data_scale']=None
        md_dic['col_types']=None
        md_dic['col_nums']=None

        for k in md_dic.keys():
            if k in self.data_table.meta.keys():
                md_dic[k]=self.data_table.meta[k]
                setattr(self, k, md_dic[k])

        if self.z is not None:
           self.z=float(self.z)
    
    
    
    
    def filter_data_set(self,filters,exclude=False):

        filters=filters.split(',')

        if exclude==False:
            print ("---> filtering for fit data_set==",filters)
        else:
            print ("---> filtering for fit data_set!=",filters)

        
        msk=np.ones( self.data['nu_data'].size, dtype=bool)
        for filter in filters:
        
            msk1=np.char.decode(self.data['data_set']) == filter
            msk=msk*msk1
        if exclude==True:
            msk=np.invert(msk)
        
        self.data=self.data[msk]
        
        
        print ("---> data len after filtering=%d"%len(self.data['nu_data']))
        

    def set_facke_error(self,):

        if 'dy' not in self.col_types:
            self.set_error(self.facke_error)

        
    def set_UL(self,val=None):
        
        self.UL_value=val

        if self.UL_value is not None:
            print("---> setting  UL")
            if 'dy' in self.col_types and self.data_scale=='lin-lin':
                print ("---> Settin  UL for val",val)
                error_array=self.data['dnuFnu_data']


            elif 'dy' in self.col_types and self.data_scale=='log-log':

                error_array=self.data['dnuFnu_data_log']

            else:
                error_array=None

            if error_array is not None:
                self.data['UL']=error_array<val


        
        
    def set_zero_error(self,val=0.2,replace_zero=True):
        
        self.zero_error_replacment=val
        print("---> replacing zero error with relative error ", val)
        if 'dy' in self.col_types and self.data_scale=='lin-lin':

            error_array=self.data['dnuFnu_data']
                       
        
        elif 'dy' in self.col_types and self.data_scale=='log-log':

            error_array=self.data['dnuFnu_data_log']
        
        else:
            error_array=None

        if error_array is not None:
            #works only on data that are not UL
            self.data['zero_error']=error_array<=0.0
            #self.data['zero_error']*=~self.data['UL']

            if replace_zero == True:
                #self.set_error(self.zero_error_replacment, data_msk=self.data['zero_error'])
                if self.data_scale=='lin-lin':
                    self.data['dnuFnu_data'][self.data['zero_error']]=self.zero_error_replacment*self.data['nuFnu_data'][self.data['zero_error']]

                if self.data_scale=='log-log':
                    self.data['dnuFnu_data_log'][self.data['zero_error']]=self.zero_error_replacment/np.log(10)

        
        
    def filter_UL(self,val=None):
        """
        remove the  upper limits points from  from data
        
        :param val: minimum value to set the upper limit. **As default, negative errors indicates upper limits, hence val=0.**
        :Retruns msk: a boolean array to mask the upper limits, i.e. all the data points with negative errors.
        
        """
        
#         
#         if 'dy' in self.col_types and self.data_scale=='lin-lin':
#             print "---> filtering for UL"
#             error_array=self.data['dnuFnu_data']
#             
#             msk= error_array<val
#             if len(msk)>0:
#                 self.ULnu_data=self.data['nu_data'][msk]
#                 self.UL_data=self.data['nuFnu_data'][msk]
#         
#         if 'dy' in self.col_types and self.data_scale=='log-log':
#             print "---> filtering for UL"
#             erro_array=self.data['dnuFnu_data_log']
#             msk= error_array<val
#             if len(msk)>0:
#                 self.UL_nu_data_log=self.data['nu_data_log'][msk]
#                 self.UL_data_log=self.data['nuFnu_data_log'][msk]
#             
#         
#         
        if val is None:
            val=self.UL_value
            
        print ("---> filtering  UL for val",val)
        if 'dy' in self.col_types and self.data_scale=='lin-lin':   
            error_array=self.data['dnuFnu_data']
        
           
        
        if 'dy' in self.col_types and self.data_scale=='log-log':
            error_array=self.data['dnuFnu_data_log']
        
        
        msk= error_array<val
        msk=np.invert(msk)
        
        self.data=self.data[msk]
        print ("---> data len after filtering=%d"%len(self.data))
        
     
    def filter_time(self,T_min=None,T_max=None,exclude=False):
        """
        filter the data, keeping all the data with 
        T_min <T< T_max if exclude=False (defualt).
        The opposite if exclude=True
        
        
        :param T_min: lower limit of the range (MJD)
        :type T_min: float
        :param T_max: upper limit of the range (MJD)
        :type T_max: float
        """
        
        
        msk1=np.ones(self.data['nu_data'].size,dtype=bool)
        msk2=np.ones(self.data['nu_data'].size,dtype=bool)
        if T_min is not None:
            msk1=self.data['T_start']>=T_min
        
        if T_max is not None:
            msk2=self.data['T_stop']<=T_max
            
        msk=msk1*msk2
        if exclude==True:
            msk=np.invert(msk)
        
        self.data=self.data[msk]
        print ("---> data len after filtering=%d"%len(self.data) )
        
        
    def filter_freq(self,nu_min=None,nu_max=None,exclude=False):
        """
        filter the data, keeping all the data with 
        nu_min <nu< nu_max if exclude=False (defualt).
        The opposite if exclude=True
        
        
        **both nu_max and nu_min are in linear scale**
        
        :param nu_min: lower limit of the range (linear scale, in Hz)
        :type nu_min: float
        :param nu_max: upper limit of the range (linear scale, in Hz)
        :type nu_max: float
        """
        
        msk1=np.ones(self.data['nu_data'].size,dtype=bool)
        msk2=np.ones(self.data['nu_data'].size,dtype=bool)
        if nu_min is not None:
            msk1=self.data['nu_data']>=nu_min
        
        if nu_max is not None:
            msk2=self.data['nu_data']<=nu_max
            
        msk=msk1*msk2
        
        if exclude==True:
            msk=np.invert(msk)
        self.data=self.data[msk]

        print ("---> data len after filtering=%d"%len(self.data))
    
    def reset_data(self):
        self._build_data()
       
  
  
    def remove_dupl_entries(self,data):
        """
         remove duplicate entries
        
        :param data: (array) 2-dim array storing the the table of the data.
        :Returns msk:  a boolean array to mask the duplicated entries
        
        .. note::
        
            One entry is flagged as duplicated as each comlum in a row of the data table
            is equal to the corresponding elements of another row
        """
        print ("---> remove duplicate entries")
        msk=np.array(np.ones(len(data)),dtype=bool)

        for i in range(1,len(data)):
            for j in range(1,len(data[:i])):
                test= data[i]==data[j]
                if test.all():
                    #print i,test,data[i],data[j]
                    msk[i]=False
               
        data=data[msk]
            
            
            
    def group_data(self,N_bin=None,bin_width=None):
        
        """
        function to perform a spectral group of the data   

        :param N_bin: (int)
        :param bin_width: (float) logarthmic
        
        .. note::
    
            To perform  a rebinning of the data has to be provided either ``N_bin`` or ``bin_width``.    
        """
        
        print (section_separator)
        

        if N_bin is not None and bin_width is  not None:
            print ("you can provide either N_bin or bin_width")
            raise ValueError 
        elif N_bin is None and bin_width is None:
            print ("you must provide either N_bin or bin_width")
            raise ValueError
        
        xmin=self.data['nu_data_log'].min()*0.99 
        xmax=self.data['nu_data_log'].max()*1.01
        
      
        if N_bin is None:
            N_bin=int((xmax-xmin)/bin_width)

        if bin_width is None:
            bin_width=(xmax-xmin)/N_bin


        bin_grid=np.linspace(xmin,xmax,N_bin)
        
        
        print ("***  binning data  ***")
        print ("---> N bins=",N_bin)
        print ("---> bin_widht=",bin_width)
  
        self.data_reb=np.zeros(bin_grid.size,dtype=self.dt)

    
        x_bin=np.zeros(bin_grid.size)
        y_bin=np.zeros(bin_grid.size)
        dx_bin=np.zeros(bin_grid.size)
        dy_bin=np.zeros(bin_grid.size)    

        #gives the id of element falling in each bin
        bin_elements_id=np.digitize(self.data['nu_data_log'],bin_grid)
        
        
        for id in range(bin_grid.size):
            msk1=[bin_elements_id==id][0]
            #Remove UL
            msk2=np.invert(self.data['UL'])
            msk=msk1*msk2
            
            if msk.any()==True>0:
                if len(self.data['dnuFnu_data_log'][msk])>1:
                    w=1.0/self.data['dnuFnu_data_log'][msk]
                    #print "w",self.data['dnuFnu_data'][msk]
                    w=w*w
                    
                    #w=1/sig_i^2
                    y_bin[id],sum_w=np.average(self.data['nuFnu_data_log'][msk],axis=0,weights=w,returned=True)
                    V1=sum_w
                    V2=w*w
                    V2=V2.sum()
                    V3=w*(self.data['nuFnu_data_log'][msk]-y_bin[id])*(self.data['nuFnu_data_log'][msk]-y_bin[id])
                    V3=V3.sum()
                    
                    if V3==0 or V1*V1-V2==0:
                    
                        #print"weighted average not possible for bin=",id
                        #print"V3 V1*V1-V2", V1,V1*V1-V2
                        #print"using err=sqrt(1/sum(w))"
                        dy_bin[id]=np.sqrt(1.0/V1)
                    
                    else:
                        dy_bin[id]=np.sqrt((V1/(V1*V1-V2))*V3)
                    
                    dx_bin[id]=bin_width/2
                    x_bin[id]=bin_grid[id]-dx_bin[id]
                    #print"x", x_bin[id],len(w)
                else:
                    #print "xxx"
                    x_bin[id]=self.data['nu_data_log'][msk]
                    y_bin[id]=self.data['nuFnu_data_log'][msk]
                    dx_bin[id]=bin_width/2
                    dy_bin[id]=self.data['dnuFnu_data_log'][msk]
                    
        self.data_reb['nu_data_log']=x_bin
        self.data_reb['nuFnu_data_log']=y_bin
        self.data_reb['dnu_data_log']=dx_bin
        self.data_reb['dnuFnu_data_log']=dy_bin
        
        #remove empty bins
        msk=[self.data_reb['nu_data_log']!=0]
        
        self.data_reb=self.data_reb[msk]
        
        
         
        self.data_reb['nuFnu_data'],self.data_reb['dnuFnu_data']=self.log_to_lin(log_val=self.data_reb['nuFnu_data_log'], log_err=self.data_reb['dnuFnu_data_log'])
          
        self.data_reb['nu_data'],self.data_reb['dnu_data']=self.log_to_lin(log_val=self.data_reb['nu_data_log'], log_err=self.data_reb['dnu_data_log'])
        
        
        self.data=self.data_reb
        
        self.set_facke_error(self.facke_error)
            
        #print  self.data['nu_data']
        #print "!!!!!!! Time and data_set must be handled somehow"

        print (section_separator)
    
    def add_systematics(self,syst,nu_range=None,data_set=None):
        """
        add systematics to errors
        
        :param syst: (float) systematic value (fractional)
        :param nu_range:  array_like of floats, [nu_min,nu_max], optional, range of frequencies to apply sistematics
        """
        
        if nu_range is not None and data_set  is not  None:
            print ("!!! error, either you provide a range of frequencies or a data_set")
            return
        
        if self.data['dnuFnu_data'] is None:
            self.data['dnuFnu_data']=np.zeros(self.data['nuFnu_data'].size)

        msk=None
        if data_set is None:
            if nu_range is None:
                #for log errors
                msk=None
                return
            else:
                msk=get_freq_range_msk(self.data['nu_data'],nu_range)
        else:
              msk=get_data_set_msk(self.data,data_set)
        
        if msk is not None:
            self.data['dnuFnu_data'][msk]=np.sqrt(self.data['dnuFnu_data'][msk]*self.data['dnuFnu_data'][msk]+ (self.data['nuFnu_data'][msk]*self.data['nuFnu_data'][msk]*syst*syst))
            self.data['nuFnu_data_log'][msk],self.data['dnuFnu_data_log'][msk]=self.lin_to_log(val=self.data['nuFnu_data'][msk], err=self.data['dnuFnu_data'][msk])

        else:
            self.data['dnuFnu_data']=np.sqrt(self.data['dnuFnu_data']*self.data['dnuFnu_data'] + (self.data['nuFnu_data']*self.data['nuFnu_data']*syst*syst))
            self.data['nuFnu_data_log'],self.data['dnuFnu_data_log']=self.lin_to_log(val=self.data['nuFnu_data'], err=self.data['dnuFnu_data'])
    
    
    
    def set_error(self,error_value,nu_range=None,data_set=None,data_msk=None):
        """
         set all the paramters to same error
            
        :param error_value: float, value of the error (fractional error)
        :param nu_range:  array_like of floats, [nu_min,nu_max], optional, range of frequencies to apply the error value
        """
        #print self.data['dnuFnu_data']           
        
        if nu_range is not None and data_set is not None:
            print ("!!! error, either you provide a range of frequencies or a data_set")
            return
        
        msk=None     
        if data_set is None:
            if nu_range is None:
                #for log errors
                msk=None
            else:
                msk=get_freq_range_msk(self.data['nu_data'],nu_range)
            
        else:
            msk=get_data_set_msk(self.data,data_set)
                
        if data_msk is not None:
            
            if msk is not None:
                msk=msk*data_msk
            else:
                msk=data_msk
        
        if msk is not None:
            self.data['dnuFnu_data'][msk]=error_value*self.data['nuFnu_data'][msk]
            self.data['nuFnu_data_log'][msk],self.data['dnuFnu_data_log'][msk]=self.lin_to_log(val=self.data['nuFnu_data'][msk], err=self.data['dnuFnu_data'][msk])
        
        else:
            self.data['dnuFnu_data']=error_value*self.data['nuFnu_data']
            self.data['nuFnu_data_log'],self.data['dnuFnu_data_log']=self.lin_to_log(val=self.data['nuFnu_data'], err=self.data['dnuFnu_data'])

       

    def set_facke_error(self,val):
        """
        Sets the value for the facke error
        """
        self.facke_error=val
        self.data['dnuFnu_facke_log']=np.ones(self.data['nu_data_log'].size)*self.facke_error
        self.data['dnuFnu_facke']=self.data['nuFnu_data']*self.facke_error
        
    
    def get_data_points(self,log_log=False,skip_UL=False):
        """
        Gives data point
        """ 
        if    skip_UL==True:
            msk=self.data['UL']==False
        else:
            msk=np.ones(self.data['nu_data_log'].size,dtype=bool)
            
        if   log_log==True:
            return self.data['nu_data_log'][msk], self.data['nuFnu_data_log'][msk], self.data['dnu_data_log'][msk], self.data['dnuFnu_data_log'][msk]
        else:
            return self.data['nu_data'][msk] , self.data['nuFnu_data'][msk] , self.data['dnu_data'][msk] , self.data['dnuFnu_data'][msk] 
        
    
    def show_data_sets(self):
        shown=[]
        for entry in self.data['data_set']:
            if entry not in shown:
                shown.append(entry)
                print (entry)
                
    
    def get_data_sets(self):
        shown=[]
        for entry in np.unique(self.data['data_set']):
            if entry not in shown:
                shown.append(entry.decode('UTF-8'))
        
        return shown
    
    
    def plot_time_spans(self,save_as=None):
        import pylab as plt
        
        fig=plt.figure(figsize=(12,9))
        ax1 = fig.add_subplot(111)
        ax1.set_xlabel('MJD')
        data_sets=self.get_data_sets()
        y=0
        line_style='-'
        for data_set in data_sets:
            print (data_set)
            
            #try:
                
            T1,T2,dT,n=self.get_time_span(data_set=data_set)
            print(T1,T2,dT,n)
            if T1!=-1:
                ax1.plot([T1,T2],[y,y],label=data_set,lw=3,marker='o')
                x=(T1+T2)/2
                ax1.text(x,y+0.3,data_set+' (%d)'%n)
                y=y+1
                #print("---->", dT,T1,T2,y)
            #except Exception as e:
            #    print ('Exception',e)
            #    print ("no Time span for data_set",data_set)
            
        ax1.set_ylim(-0.5,y+0.5)
        fig.show()
        if save_as is not None:
            fig.savefig(save_as)
       
        
       
        
    
    def get_time_span(self,data_set=None):
        return self.find_time_span(data_set=data_set,get_values=True)
    
    def show_time_span(self,data_set=None):
        self.find_time_span(data_set=data_set,silent=False)
    
    def find_time_span(self,data_set=None,silent=True,get_values=False):
        """
        returns Tstart, Tstop, and Delta T for the full data set (if no dat_set 
        is provided), or for a specific data_set
        """
        
        if data_set is None:
            m1= self.data['T_start']!=0
            T1=self.data['T_start'][m1].min()
            T2=self.data['T_stop'][m1].max()
            DT=T2-T1
            nT=self.data['T_start'][m1].size
            
        elif  data_set  in self.data['data_set'].astype(str):
        
            m1= self.data['T_start']!=0
            m2=self.data['data_set'].astype(str)==data_set
            if self.data['data_set'][m2].size>0:
                try:
                    T1=self.data['T_start'][m1*m2].min()
                    T2=self.data['T_stop'][m1*m2].max()
                    DT=T2-T1
                    nT=self.data['T_start'][m1*m2].size
                except:
                    T1=-1
                    T2=-1
                    DT=-1
                    nT=0
                    print ('something wrong with T_start and T_stop columns, check the values please, for data_set=',data_set)
        else:
            print ("no data found for this selection, data_set= ",data_set)
            if data_set not in self.data['data_set']:
                print ("the data_set %s is not present in the data"%data_set)
                print ("possible data_set: ")
                print (self.show_data_sets())
                    
            T1=-1
            T2=-1
            DT=-1
            nT=0
        
        if silent==False:
            print ("T_start=%f T_stop=%f DT=%f, number of points=%d"%(T1,T2,DT,nT))
        if get_values==True: 
            return T1,T2,DT,nT
            
            
    def lin_to_log(self,val=None,err=None):
        
       return lin_to_log(val=val,err=err)
       
    
    def log_to_lin(self,log_val=None,log_err=None):
        
        return log_to_lin(log_val=log_val,log_err=log_err)
        



def lin_to_log(val=None,err=None):
       
       conv_fac=np.log(10)
       
       ret_val=[]
       
       if val is not None:
           log_val=np.log10(val)
           ret_val.append(log_val)
       if err is not None:
           log_err=err/(val*conv_fac)
           ret_val.append(log_err)

       if len(ret_val)==1:
           ret_val=ret_val[0]
       
       return ret_val
   
   
def log_to_lin(log_val=None,log_err=None):
       
       conv_fac=np.log(10)
       
       ret_val=[]
       
       if log_val is not None:
           val=np.power(10.,log_val)
           
           ret_val.append(val)

       if log_err is not None:
           err=log_err*val*conv_fac
           #print err
           ret_val.append(err)

       if len(ret_val)==1:
           ret_val=ret_val[0]
       
       return ret_val


def get_freq_range_msk(x,x_range):
    msk1=x>=x_range[0]
    msk2=x<=x_range[1]
    return msk1*msk2

def get_data_set_msk(x,data_set):
    return x['data_set']==data_set