__author__ = "Andrea Tramacere"


import numpy as np
import copy


from astropy.table  import  Table,Column
from astropy import  units as u
from astropy.units import cds
import  pickle
from .plot_sedfit import PlotSED, plt

from .output import section_separator
from .utils import *
from .cosmo_tools import Cosmo
from .frame_converter import convert_nu_to_src, convert_nuFnu_to_nuLnu_src

cds.enable()

__all__=['get_data_set_msk','get_freq_range_msk','lin_to_log','log_to_lin','ObsData','Data']


class Data(object):
    """
    Class to store obs data
    """
    def __init__(self,
                 data_table=None,
                 n_rows=None,
                 meta_data=None,
                 import_dictionary=None,
                 cosmo=None):
        """_summary_

        Parameters
        ----------
        data_table : astropy table, optional
            _description_, by default None
        n_rows : int, optional
            the number of rows for an empty table if no data_table is passed, by default None
        meta_data : dictionary, optional
            the dictionary for the metadata, by default None
        import_dictionary : dictionary, optional
            a dictionary mapping the column renaming if needed for the data_table, by default None
        cosmo : jetset.cosmo_tools.Cosmo, optional
            the cosmology object, by default None
        """

        if cosmo is None:

            self.cosmo = Cosmo()
        else:

            self.cosmo = cosmo

        self._necessary_meta = ['data_scale', 'z', 'restframe']
        self._allowed_meta={}
        self._allowed_meta['z']= None
        self._allowed_meta['UL_CL'] = None
        self._allowed_meta['restframe'] = ['obs','src']
        self._allowed_meta['data_scale'] = ['lin-lin','log-log']
        self._allowed_meta['obj_name'] = None

        self._names = ['x', 'dx', 'y', 'dy', 'T_start', 'T_stop', 'UL', 'data_set']
        self._dt = ('f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'bool', 'S16')

        self._units = [u.Hz, u.Hz,(u.erg / (u.cm ** 2 * u.s)) ,(u.erg / (u.cm ** 2 * u.s)), cds.MJD, cds.MJD, None, None]

       
        if isinstance(data_table,str):
            data_table = Table.read(data_table, guess=True)

        if data_table is None:
            self._build_empty_table(n_rows,meta_data=meta_data)

        else:
            if import_dictionary is not None:
                for k in import_dictionary.keys():
                    data_table.rename_column(k, import_dictionary[k])

            self._table = copy.deepcopy(data_table)

        if meta_data is not None:
            self._table.meta = {}
            for k in meta_data.keys():

                self.set_meta_data(k,meta_data[k])
        else:
            for k in self._table.meta.keys():
                self.set_meta_data(k,self._table.meta[k])

        self._check_table()
        self._check_frame_and_scale()
        self._convert_units()

        #self._table.sort('x')


    @property
    def table(self):
        return self._table

    @property
    def metadata(self):
        return self._table.meta

    @classmethod
    def from_file(cls,data_table,format='ascii.ecsv',import_dictionary=None,guess=None):
        return cls(data_table= Table.read(data_table, format=format,guess=guess),import_dictionary=import_dictionary)

    def save_file(self, name, format='ascii.ecsv'):
        self._table.write(name,format=format,overwrite=True)



    def _check_frame_and_scale(self):
        if self.metadata['data_scale'] == 'log-log':
            self._table['x'], self._table['dx'] = log_to_lin(self._table['x'], self._table['dx'])
            self._table['y'], self._table['dy'] = log_to_lin(self._table['y'], self._table['dy'])

            self.metadata['data_scale']='lin-lin'

        if self.metadata['restframe'] == 'src':
            if self._table['y'].unit.is_equivalent('erg/s') and  self._table['dy'].unit.is_equivalent('erg/s'):
                pass
            else:
                raise  RuntimeError('when importing src frame table, units for luminosities have to be in erg/s')



            _c=self.cosmo.get_DL_cm(self.metadata['z'])
            _c=1.0/(4*np.pi*_c*_c)
            self._table['y']  = self._table['y'] * _c
            self._table['dy'] = self._table['dy']  * _c
            self._table['y'].unit = 'erg/(cm2 s)'
            self._table['dy'].unit = 'erg/(cm2 s)'
            self.metadata['restframe'] = 'obs'

    def _check_table(self):

        for ID,_cn in enumerate(self._names):
            if _cn not in self._table.colnames:
                self._table.add_column(index=ID,col=Column(name=_cn,dtype=self._dt[ID],unit=self._units[ID],data=np.zeros(len(self._table))))

        for _n in  self._necessary_meta:
            if _n not in  self._table.meta:
                raise RuntimeError('meta data',_n,'not defined, but it is necessary')


    def _convert_units(self):
        for ID,_cn in enumerate(self._names):
            if hasattr(self._table[_cn],'unit'):
                if self._table[_cn].unit is not None and self._units[ID] is not None:
                    try:
                        self._table[_cn]=self._table[_cn].to(self._units[ID])
                    except Exception as e:
                        try:
                            self._table[_cn] = self._table[_cn].to(self._units[ID], equivalencies=u.spectral())
                        except Exception as e:
                            raise RuntimeError('Unit conversion problem for ',self._table[_cn].unit,'to',self._units[ID],repr(e))


    def set_meta_data(self,m,v):
        if m not in self._allowed_meta:
            raise RuntimeError('meta data ',m,'not in allowed',self._allowed_meta.keys())

        if  self._allowed_meta[m] is not None:
            if v not in self._allowed_meta[m]:
                raise RuntimeError('meta data value ', v, 'not in allowed', self._allowed_meta[m])


        self._table.meta[m]=v



    def set_field(self,field,value,unit=None):

        if unit is None:
            if hasattr(self._table[field], 'unit'):
                unit = self._table[field].unit

        self._table[field] = value
        if unit is not None:
            self._table[field].unit = unit

        

    def _build_empty_table(self,n_rows,meta_data=None):


        self._table = Table(np.zeros((n_rows, len(self._names))), names=self._names, dtype=self._dt)
        for ID,c in enumerate(self._table.columns):
            if self._units[ID] is not None:
                #print(ID,c, units[ID])
                self._table[c]*=self._units[ID]

        self._table.meta['z'] = 0
        self._table.meta['UL_CL'] = 0.95
        self._table.meta['restframe'] = 'obs'
        self._table.meta['data_scale'] = 'lin-lin'
        self._table.meta['obj_name'] = 'new-src'

        if meta_data is not None:
            for k in self._table.meta.keys():
                if k in meta_data.keys():
                    self._table.meta[k]=meta_data[k]



    @classmethod
    def from_asdc(cls, asdc_sed_file, obj_name, z, restframe, data_scale):
        with open(asdc_sed_file, 'r') as f:
            lines = f.readlines()
        # print(len(lines),type(lines),lines)
        for l in lines[:]:
            if l.startswith('#'):
                lines.remove(l)

        UL = np.zeros(len(lines), dtype=np.bool)

        for ID, l in enumerate(lines):
            t = l.strip().split(';')

            if len(t) > 1:
                # print(t[1])

                if 'UPPER LIMIT' in t[1]:
                    UL[ID] = True
                    lines[ID] = t[0]

        d = np.genfromtxt(lines)
        d = np.column_stack((d, UL))

        data_table = Table(d, names=['x', 'dx', 'y', 'dy', 'T_start', 'T_stop', 'UL'])
        data_table['x'] = data_table['x'] * u.Hz
        data_table['dx'] = data_table['dx'] * u.Hz
        data_table['y'] = data_table['y'] * (u.erg / (u.cm ** 2 * u.s))
        data_table['dy'] = data_table['dy'] * (u.erg / (u.cm ** 2 * u.s))
        data_table['T_start'] = data_table['T_start'] * cds.MJD
        data_table['T_stop'] = data_table['T_stop'] * cds.MJD
        data_table['UL']=np.array(data_table['UL'], dtype=np.bool)
        data_table.meta['z'] = z
        data_table.meta['restframe'] = restframe
        data_table.meta['data_scale'] = data_scale
        data_table.meta['obj_name'] = obj_name
        return cls(data_table=data_table)


class ObsData(object):
    """ObsData class

    This class provides a  powerful interface  to load observational data stored in a file.
    The following parameters set the corresponding class members
    """


    def __init__(self,
                 cosmo=None,
                 data_table=None,
                 dupl_filter=False,
                 data_set_filter=None,
                 UL_filtering=False,
                 UL_value=None,
                 UL_CL=0.95,

                 **keywords):
        
        """

        Parameters
        ----------
        cosmo
        data_table
        dupl_filter
        data_set_filter
        UL_filtering
        UL_value
        UL_CL
        keywords
        """
        
        self.z=None
        self.data=None
        self.obj_name=None
        self.restframe=None
        self.col_types=None
        self.col_nums=None
        self.data_scale=None
        self.UL_CL=UL_CL

        if data_set_filter is not None:
            self.data_set_filter=data_set_filter.split(',')
        else:
            self.data_set_filter=['No']

        if cosmo is None:
            if hasattr(data_table,'cosmo'):
                self.cosmo=data_table.cosmo
            else:
                self.cosmo=Cosmo()
        else:
            self.cosmo=cosmo

        self.UL_value=UL_value
        self.UL_filtering=UL_filtering
        self.zero_error_replacment=0.2
        self.fake_error=0.2



        if  hasattr(data_table,'table'):
            _t=data_table.table
        else:
            _t=data_table

        self._input_data_table=None

        if isinstance(_t,Table):
            self._input_data_table = _t
        else:
            raise RuntimeError('table is not an astropy Table')

        self.allowed_keywords = {'z': None}
        self.allowed_keywords['obj_name'] = None
        self.allowed_keywords['restframe'] = ['obs', 'src']
        self.allowed_keywords['data_scale'] = ['lin-lin', 'log-log']
        self.allowed_keywords['file_name'] = None
        self.allowed_keywords['UL_CL'] = None
        self.allowed_keywords['col_types'] = None
        self.allowed_keywords['col_nums'] = None


        _skip = ['data_table', 'n_rows']

        if self._input_data_table is not None:
            self._set_kw(self._input_data_table.meta)

        self._set_kw(keywords,_skip)

        self._build_data(dupl_filter=dupl_filter)
        self.data.sort('nu_data')

    def _set_kw(self,keywords,skip=[]):
        keys = sorted(keywords.keys())
        # print keys

        for kw in keys:
            if kw not in skip:
                if kw in self.allowed_keywords.keys():

                    if self.allowed_keywords[kw] is not None:
                        # check that the kyword value is correct
                        if keywords[kw] not in self.allowed_keywords[kw]:
                            print("keyword=%s, has wrong value=%s, allowed are %s" % (
                            kw, keywords[kw], self.allowed_keywords[kw]))
                            raise ValueError

                    setattr(self, kw, keywords[kw])

                else:

                    print("wrong keyword=%s, not in%s " % (kw, self.allowed_keywords.keys()))

                    raise ValueError

        if self.z is not None:
            self.z = float(self.z)

        if self.UL_CL is not None:
            self.UL_CL = float(self.UL_CL)

    @property
    def metadata(self,skip=['col_types','col_nums']):
        _d={}
        for k in self.allowed_keywords.keys():
            if hasattr(self,k) and k not in skip:
                _d[k]=getattr(self,k)
        return _d

    def _build_empty_table(self,n_rows=None):
        sed_dt = [('nu_data', 'f8')]
        sed_dt.append(('dnu_data', 'f8'))
        sed_dt.append(('nuFnu_data', 'f8'))
        sed_dt.append(('dnuFnu_data', 'f8'))

        sed_dt.append(('nu_data_log', 'f8'))
        sed_dt.append(('dnu_data_log', 'f8'))
        sed_dt.append(('nuFnu_data_log', 'f8'))
        sed_dt.append(('dnuFnu_data_log', 'f8'))

        sed_dt.append(('dnuFnu_fake', 'f8'))
        sed_dt.append(('dnuFnu_fake_log', 'f8'))

        sed_dt.append(('UL', 'bool'))
        sed_dt.append(('zero_error', 'bool'))
        sed_dt.append(('T_start', 'f8'))
        sed_dt.append(('T_stop', 'f8'))
        sed_dt.append(('data_set', 'S16'))

        self.dt = np.dtype(sed_dt)

        if n_rows is None:
            n_rows=len(self._input_data_table)

        return Table(rows=np.zeros((n_rows, len(sed_dt))), dtype=[d[1] for d in sed_dt],
                          names=[d[0] for d in sed_dt])

    @staticmethod
    def load(file_name,):
        _obj=pickle.load(open(file_name, "rb"))
        _obj.data.sort('nu_data')
        return _obj

    def save(self, file_name):
        self.data.sort('nu_data')
        pickle.dump(self, open(file_name, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)


    @property
    def table(self):
        return self.data

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


        self.data=self._build_empty_table()

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
        self._log_col_dict['data_set'] = 'data_set'
        self._log_col_dict['UL'] = 'UL'
        
        if self.data_scale=='lin-lin':
            self._file_col_dict= self._col_dict
            
        
        elif self.data_scale=='log-log':
             self._file_col_dict= self._log_col_dict
        
        else:
            raise RuntimeError('data_scale not specified')

        self.col_types = []
        for n in self._input_data_table.colnames:
            if n  in self._file_col_dict.keys():
                self.col_types.append(n)

                self.data[self._file_col_dict[n]]=self._input_data_table[n]

        self.col_nums=len(self.col_types)



        #-------------------------------------------------------------------------
        # duplicate entries 
        if dupl_filter==True:
            #print ("---> filtering for dupl entries")
            self.remove_dupl_entries(self.data)
       
        #-------------------------------------------------------------------------          

        # sets UL entries
        self.set_UL(self.UL_value)

        # set error if not present:
        self.set_fake_error(self.fake_error)





        self.set_zero_error()

        if self.UL_filtering==True:
            self.filter_UL(self.UL_value)
        
        
        if self.data_set_filter!=['No']:
            self.filter_data_set( self.data_set_filter)
       

        self._set_data_frame_and_scale()

       
        
        if self.data['dnuFnu_data'] is None and self.data_scale== 'lin-lin':
        
            self.data['dnuFnu_data']= self.data['Fnu_data'] * self.fake_error
            
            print ("Warning: error were not provided ")
            print ("         assigning %f            ")
            print ("         set error with .set_error"%self.fake_error)
            #print self.data['dnuFnu_data']
        
        
        if self.data['dnuFnu_data_log'] is None and self.data_scale== 'log-log':
        
            self.data['dnuFnu_data_log']= np.ones(self.data['nu_data_log'].size) * self.fake_error
            
            print ("Warning: error were not provided ")
            print ("         assigning %f           ")
            print ("         set error with .set_error"%self.fake_error)

    








    def _set_data_frame_and_scale(self):
        
        
        #-------------------------------------------------------------------------          
        # handles data cosmological conversion
        # and axis transformation

        check_frame(self.restframe)

        if self.restframe=='obs':
            self.nu_conv_factor=1
            self.Lum_conv_factor=1
        
        elif self.restframe=='src':
            #TODO this must be the same as in jetset
            DL=self.cosmo.get_DL_cm(self.z)
            print ("--->DL=%e"%DL)
            
            #!! usa le funzioni in frame_converter
            self.nu_conv_factor=1.0/(1+self.z)
            self.Lum_conv_factor=1.0/(np.pi*4.0*DL*DL)
        else:
            unexpected_behaviour()

        if self.data_scale=='lin-lin':
            
            if 'dy' in self.col_types :
                self.data['nuFnu_data_log'], self.data['dnuFnu_data_log']=self.lin_to_log(val=self.data['nuFnu_data'], err=self.data['dnuFnu_data'])
            else:
                self.data['nuFnu_data_log']=self.lin_to_log(val=self.data['nuFnu_data'])
                
            if 'dx' in self.col_types:
                self.data['nu_data_log'], self.data['dnu_data_log']=self.lin_to_log(val=self.data['nu_data'], err=self.data['dnu_data'])
            else:
                self.data['nu_data_log']=self.lin_to_log(val=self.data['nu_data'])
            
            self.data['nu_data_log']+= np.log10(self.nu_conv_factor)
            self.data['nuFnu_data_log']+= np.log10(self.Lum_conv_factor)
            
            
            
        if self.data_scale=='log-log':

            if 'dy' in self.col_types:
                self.data['nuFnu_data'], self.data['dnuFnu_data']=self.log_to_lin(log_val=self.data['nuFnu_data_log'], log_err=self.data['dnuFnu_data_log'])
            else:
                self.data['nuFnu_data']=self.log_to_lin(log_val=self.data['nuFnu_data_log'])
                
            if 'dx' in self.col_types:
                self.data['nu_data'], self.data['dnu_data']=self.log_to_lin(log_val=self.data['nu_data_log'], log_err=self.data['dnu_data_log'])
            else:
                self.data['nu_data']=self.log_to_lin(log_val=self.data['nu_data_log'])
            
            
            #print self.data['nuFnu_data'],self.data['dnuFnu_data']
            self.data['nu_data']*=self.nu_conv_factor
            self.data['nuFnu_data']*=self.Lum_conv_factor
        
        

    def _set_md_from_data_table(self,name=None,val=None):




        md_dic = {}
        #md_dic['z'] = None
        #md_dic['file_name'] = None
        #md_dic['obj_name'] = None
        #md_dic['restframe'] = None
        #md_dic['data_scale'] = None
        #md_dic['col_types'] = None
        #md_dic['col_nums'] = None

        if name is None:

            for k in self.allowed_keywords.keys():
                if k in self._input_data_table.meta.keys():
                    #md_dic[k]=self._input_data_table.meta[k]
                    setattr(self, k, self._input_data_table.meta[k])
                    #print(k,md_dic[k])
        else:
            if name in md_dic.keys():
                setattr(self, name, val)
            else:
                raise RuntimeError('meta name',name,'not in allowed',md_dic.keys())

        if self.z is not None:
           self.z=float(self.z)

        if self.UL_CL is not None:
            self.UL_CL = float(self.UL_CL)
    
    
    
    def filter_data_set(self,filters,exclude=False,silent=False):

        filters=filters.split(',')

        if silent is False:
            if exclude==False:
                print ("---> including  only data_set/s",filters)
            else:
                print ("---> excluding  data_set/s",filters)

        
        msk=np.ones(self.data['nu_data'].size, dtype=bool)

        if exclude == True:
            msk = np.ones(self.data['nu_data'].size, dtype=bool)

        if exclude == False:
            msk = np.zeros(self.data['nu_data'].size, dtype=bool)

        for filter in filters:
            #print ('filter',filter)
            msk1= self.data['data_set'] == filter
            if exclude == True:
                msk1 = np.invert(msk1)
                msk = np.logical_and(msk, msk1)
            else:
                msk=np.logical_or(msk,msk1)
            if silent is False:
                print('filter', filter, np.sum(msk))

        
        self.data=self.data[msk]
        if silent is False:
            print("---> data sets left after filtering",self.show_data_sets())
            print ("---> data len after filtering=%d" % len(self.data['nu_data']))
        

    def set_fake_error(self,):

        if 'dy' not in self.col_types:
            self.set_error(self.fake_error)

        
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
                self.data['UL']= error_array < val


        
        
    def set_zero_error(self,val=0.2,replace_zero=True):
        
        self.zero_error_replacment=val
        #print("---> replacing zero error with relative error ", val)
        if 'dy' in self.col_types and self.data_scale=='lin-lin':

            error_array=self.data['dnuFnu_data']
                       
        
        elif 'dy' in self.col_types and self.data_scale=='log-log':

            error_array=self.data['dnuFnu_data_log']
        
        else:
            error_array=None

        if error_array is not None:
            #works only on data that are not UL
            self.data['zero_error']= error_array <= 0.0
            #self.data['zero_error']*=~self.data['UL']

            if replace_zero == True:
                #self.set_error(self.zero_error_replacment, data_msk=self.data['zero_error'])
                if self.data_scale=='lin-lin':
                    self.data['dnuFnu_data'][self.data['zero_error']]= self.zero_error_replacment * self.data['nuFnu_data'][self.data['zero_error']]

                if self.data_scale=='log-log':
                    self.data['dnuFnu_data_log'][self.data['zero_error']]= self.zero_error_replacment / np.log(10)

        
        
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
        print ("---> data len after filtering=%d" % len(self.data))
        
     
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
        
        
        msk1=np.ones(self.data['nu_data'].size, dtype=bool)
        msk2=np.ones(self.data['nu_data'].size, dtype=bool)
        if T_min is not None:
            msk1= self.data['T_start'] >= T_min
        
        if T_max is not None:
            msk2= self.data['T_stop'] <= T_max
            
        msk=msk1*msk2
        if exclude==True:
            msk=np.invert(msk)
        
        self.data=self.data[msk]
        print ("---> data len after filtering=%d" % len(self.data))
        
        
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
        
        msk1=np.ones(self.data['nu_data'].size, dtype=bool)
        msk2=np.ones(self.data['nu_data'].size, dtype=bool)
        if nu_min is not None:
            msk1= self.data['nu_data'] >= nu_min
        
        if nu_max is not None:
            msk2= self.data['nu_data'] <= nu_max
            
        msk=msk1*msk2
        
        if exclude==True:
            msk=np.invert(msk)
        self.data=self.data[msk]

        print ("---> data len after filtering=%d" % len(self.data))
    
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
                    #print i,tests,data[i],data[j]
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

        #if nu_min is None:
        xmin= self.data['nu_data_log'].min() * 0.99
        #else:
        #   xmin= np.log10(nu_min)

        #if nu_max is None:
        xmax= self.data['nu_data_log'].max() * 1.01
        #else:
        #xmin = np.log10(nu_max)


        if N_bin is None:
            N_bin=int((xmax-xmin)/bin_width)

        if bin_width is None:
            bin_width=(xmax-xmin)/N_bin


        bin_grid=np.linspace(xmin,xmax,N_bin)
        
        
        print ("***  binning data  ***")
        print ("---> N bins=",N_bin)
        print ("---> bin_widht=",bin_width)
  
        self.data_reb=self._build_empty_table(n_rows=bin_grid.size)

    
        x_bin=np.zeros(bin_grid.size)
        y_bin=np.zeros(bin_grid.size)
        dx_bin=np.zeros(bin_grid.size)
        dy_bin=np.zeros(bin_grid.size)    

        #gives the id of element falling in each bin
        bin_elements_id=np.digitize(self.data['nu_data_log'], bin_grid)
        
        
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
                    y_bin[id],sum_w=np.average(self.data['nuFnu_data_log'][msk], axis=0, weights=w, returned=True)
                    V1=sum_w
                    V2=w*w
                    V2=V2.sum()
                    V3= w * (self.data['nuFnu_data_log'][msk] - y_bin[id]) * (self.data['nuFnu_data_log'][msk] - y_bin[id])
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
        
        #set original units
        for c in self.data.colnames:
            if self.data[c].unit is not None:
                self.data_reb[c] = self.data_reb[c]*self.data[c].unit

        self.data=self.data_reb
        
        self.set_fake_error(self.fake_error)

        self.data.sort('nu_data')
       
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

                msk=np.ones(len(self.data),dtype=np.bool)

            else:
                msk=get_freq_range_msk(self.data['nu_data'], nu_range)
        else:
              msk=get_data_set_msk(self.data, data_set)
        
        if msk is not None:
            self.data['dnuFnu_data'][msk]=np.sqrt(self.data['dnuFnu_data'][msk] * self.data['dnuFnu_data'][msk] + (self.data['nuFnu_data'][msk] * self.data['nuFnu_data'][msk] * syst * syst))
            self.data['nuFnu_data_log'][msk], self.data['dnuFnu_data_log'][msk]=self.lin_to_log(val=self.data['nuFnu_data'][msk], err=self.data['dnuFnu_data'][msk])

        else:
            self.data['dnuFnu_data']=np.sqrt(self.data['dnuFnu_data'] * self.data['dnuFnu_data'] + (self.data['nuFnu_data'] * self.data['nuFnu_data'] * syst * syst))
            self.data['nuFnu_data_log'], self.data['dnuFnu_data_log']=self.lin_to_log(val=self.data['nuFnu_data'], err=self.data['dnuFnu_data'])

        self.data.sort('nu_data')
        
    
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
                msk=get_freq_range_msk(self.data['nu_data'], nu_range)
            
        else:
            msk=get_data_set_msk(self.data, data_set)
                
        if data_msk is not None:
            
            if msk is not None:
                msk=msk*data_msk
            else:
                msk=data_msk
        
        if msk is not None:
            self.data['dnuFnu_data'][msk]= error_value * self.data['nuFnu_data'][msk]
            self.data['nuFnu_data_log'][msk], self.data['dnuFnu_data_log'][msk]=self.lin_to_log(val=self.data['nuFnu_data'][msk], err=self.data['dnuFnu_data'][msk])
        
        else:
            self.data['dnuFnu_data']= error_value * self.data['nuFnu_data']
            self.data['nuFnu_data_log'], self.data['dnuFnu_data_log']=self.lin_to_log(val=self.data['nuFnu_data'], err=self.data['dnuFnu_data'])

       

    def set_fake_error(self,val):
        """
        Sets the value for the fake error
        """
        self.fake_error=val
        self.data['dnuFnu_fake_log']= np.ones(self.data['nu_data_log'].size) * self.fake_error
        self.data['dnuFnu_fake']= self.data['nuFnu_data'] * self.fake_error
        
    
    def get_data_points(self,log_log=False,skip_UL=False,frame='obs',density=False):
        """
        Gives data point
        """

        if   skip_UL==True:
            msk= self.data['UL'] == False
        else:
            msk=np.ones(self.data['nu_data_log'].size, dtype=bool)

        check_frame(frame)
        if frame == 'obs':
            if   log_log  is True:
                _x = self.data['nu_data_log'][msk]
                _dx = self.data['dnu_data_log'][msk]
                _y = self.data['nuFnu_data_log'][msk]
                _dy =  self.data['dnuFnu_data_log'][msk]

            else:
                _x = self.data['nu_data'][msk]
                _dx = self.data['dnu_data'][msk]
                _y = self.data['nuFnu_data'][msk]
                _dy = self.data['dnuFnu_data'][msk]

        elif frame == 'src':

            dl = self.cosmo.get_DL_cm(self.z)

            _x = convert_nu_to_src(self.data['nu_data'][msk], self.z, 'obs')
            _dx = convert_nu_to_src(self.data['dnu_data'][msk], self.z, 'obs')
            _y = convert_nuFnu_to_nuLnu_src(self.data['nuFnu_data'][msk], self.z, 'obs', dl)
            _dy = convert_nuFnu_to_nuLnu_src(self.data['dnuFnu_data'][msk], self.z, 'obs', dl)

            if log_log is True:
                _x, _dx = self.lin_to_log(_x, _dx)
                _y, _dy = self.lin_to_log(_y, _dy)

        else:
            unexpected_behaviour()

        if density is True:

            if log_log is True:
                _y = _y - _x
            else:
                _y = _y / _x
                _dy=_dy/_x

        return _x ,_y, _dx, _dy
    
    def show_data_sets(self):
        shown=[]
        print('current datasets')
        for entry in self.data['data_set']:
            if entry not in shown:
                shown.append(entry)
                print ('dataset', entry)
                
    
    def get_data_sets(self):
        shown=[]
        for entry in np.unique(self.data['data_set']):
            if entry not in shown:
                shown.append(entry)
        
        return shown



    def plot_sed(self,plot_obj=None,frame='obs',color=None,fmt='o',ms=4,mew=0.5,figsize=None,show_dataset=False, density=False):
        if plot_obj is None:
            plot_obj = PlotSED(frame=frame, figsize=figsize,density=density)

        if show_dataset is False:

            plot_obj.add_data_plot(self, color=color, fmt=fmt, ms=ms, mew=mew, density=density)
        else:
            for ds in self.get_data_sets():
                self.filter_data_set(filters=ds,silent=True,exclude=False)
                plot_obj.add_data_plot(self, color=color, fmt=fmt, ms=ms, mew=mew ,label='dataset %s'%ds, density=density)
                self.reset_data()
            self.reset_data()

        return plot_obj

    
    def plot_time_spans(self,save_as=None):


        fig=plt.figure(figsize=(12,9))
        ax1 = fig.add_subplot(111)
        ax1.set_xlabel('MJD')
        data_sets=self.get_data_sets()
        y=0
        line_style='-'
        for data_set in data_sets:
            print (data_set)



            T1,T2,dT,n=self.get_time_span(data_set=data_set)
            print(T1,T2,dT,n)
            if T1!=-1:
                ax1.plot([T1,T2],[y,y],label=data_set,lw=3,marker='o')
                x=(T1+T2)/2
                ax1.text(x,y+0.3,data_set+' (%d)'%n)
                y=y+1

        ax1.set_ylim(-0.5,y+0.5)

        if save_as is not None:
            fig.savefig(save_as)

        return fig
       
        
    
    def get_time_span(self,data_set=None):
        return self.find_time_span(data_set=data_set,get_values=True)
    
    def show_time_span(self,data_set=None):
        self.find_time_span(data_set=data_set,silent=False)
    
    def find_time_span(self,data_set=None,silent=True,get_values=False):
        """
        returns Tstart, Tstop, and Delta T for the full data set (if no dat_set
        is provided), or for a specific data_set
        """

        time_span_found = False

        T1 = -1
        T2 = -1
        DT = -1
        nT = 0


        if data_set is None:
            #m1= self.data['T_start'] != 0

            T1 = self.data['T_start'].min()
            T2 = self.data['T_stop'].max()
            DT = T2 - T1
            nT = self.data['T_start'].size
            time_span_found=True

            
        elif  data_set  in self.data['data_set'].astype(str):
        
            #m1= self.data['T_start'] != 0
            m2= self.data['data_set'].astype(str) == data_set
            if self.data['data_set'][m2].size>0:
                try:

                    T1=self.data['T_start'][m2].min()
                    T2=self.data['T_stop'][m2].max()
                    DT=T2-T1
                    nT=self.data['T_start'][m2].size
                    time_span_found = True
                except:
                    time_span_found = False
                    T1=-1
                    T2=-1
                    DT=-1
                    nT=0
                    print ('something wrong with T_start and T_stop columns, check the values please, for data_set=',data_set)
        else:
            time_span_found = False
            print ("no data found for this selection, data_set= ",data_set)
            if data_set not in self.data['data_set']:
                print ("the data_set %s is not present in the data"%data_set)
                print ("possible data_set: ")
                print (self.show_data_sets())
                    

        
        if silent==False and time_span_found is True:
            print ("T_start=%f T_stop=%f DT=%f, number of points=%d"%(T1, T2, DT, nT))

        if get_values==True:
            return T1,T2,DT,nT
            

    def lin_to_log(self,val=None,err=None):
        
       return lin_to_log(val=val,err=err)
       
    
    def log_to_lin(self,log_val=None,log_err=None):
        
        return log_to_lin(log_val=log_val,log_err=log_err)
        

    @property
    def gammapy_table(self):
        return self.get_gammapy_table()


    def get_gammapy_table(self):
        table = Table()
        for c in self.metadata:
            table.meta[c]=self.metadata[c]
                
        table['e_ref']=self.data['nu_data'].to("eV", equivalencies=u.spectral())
      
        table["e2dnde"] = self.data['nuFnu_data']
        table["e2dnde_err"] = self.data['dnuFnu_data']
        
        table.meta["SED_TYPE"] = "e2dnde"
        
        return table


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