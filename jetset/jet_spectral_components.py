
__author__ = "Andrea Tramacere"


import numpy as np
import os
from astropy.table import Table
from numpy.core._multiarray_umath import zeros, log10
from scipy import interpolate

# on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

# if on_rtd == True:
#     try:
#         from .jetkernel import jetkernel as BlazarSED
#     except ImportError:
#         from .mock import jetkernel as BlazarSED
# else:

from .jetkernel import jetkernel as BlazarSED

from . import spectral_shapes
from .jetkernel_models_dic import nuFnu_obs_dict, n_seed_dic
from .plot_sedfit import PlotSpecComp,PlotSeedPhotons
from .utils import check_frame, unexpected_behaviour

__all__=['JetSeedPhotons','JetSpecComponent','SpecCompList']

class JetSeedPhotons(object):
    """

    """
    def __init__(self,name,blob_object,var_name=None):
        self.name = name

        self._blob_object = blob_object
        self._n_name, self._nu_name = n_seed_dic[self.name]

        self.n_ptr = getattr(blob_object, self._n_name)

        self.nu_ptr = getattr(blob_object, self._nu_name)
        #self.SED = spectral_shapes.SED(name=self.name)
        if var_name is not None:
            self._var_name=var_name

        self.fill(emiss_lim=self._blob_object.emiss_lim)

    def fill(self,log_log=False,emiss_lim=0):
        self.nu,self.n=self.get_spectral_points(log_log=log_log,emiss_lim=emiss_lim)

    def get_spectral_points(self,log_log=False,emiss_lim=0):

        #try:

        size=self._blob_object.nu_grid_size
        x=zeros(size)
        y=zeros(size)

        for i in range(size):
            x[i]=BlazarSED.get_spectral_array(self.nu_ptr,self._blob_object,i)
            y[i]=BlazarSED.get_spectral_array(self.n_ptr,self._blob_object,i)

            #print("->%e %e"%(x[i],y[i]))

        msk_nan=np.isnan(x)
        msk_nan+=np.isnan(y)
        #print('emiss lim',self.get_emiss_lim())
        x[msk_nan]=0.
        y[msk_nan]=emiss_lim

        msk=y<emiss_lim


        y[msk]=emiss_lim



        if log_log==True:
            msk = y <= 0.
            y[msk] = emiss_lim

            #x=x[msk]
            #    y=y[msk]

            x=log10(x)
            y=log10(y)



        return x,y

        #except:
        #    raise RuntimeError ('model evaluation failed in get_spectral_points')


    def plot(self, y_min=None,y_max=None):
        self.fill(emiss_lim=self._blob_object.emiss_lim)
        p=PlotSeedPhotons()
        p.plot(nu=self.nu,nuFnu=self.n,y_min=y_min,y_max=y_max)

        return p



class JetSpecComponent(object):
    """
    """

    def __repr__(self):
        return str(self.show())

    #def __str__(self):
    #    return str(self.show())


    def __init__(self,jet_obj,name,blob_object,var_name=None,state_dict=None,state=None):

        self.name=name
        self.jet_obj=jet_obj

        self._blob_object=blob_object
        self._nuFnu_name, self._nu_name=nuFnu_obs_dict[self.name]

        self.nuFnu_ptr=getattr(blob_object,self._nuFnu_name)

        self.nu_ptr=getattr(blob_object,self._nu_name)

        self.SED=spectral_shapes.SED(name=self.name,beaming=jet_obj.get_beaming())
        self.seed_field=None

        # self._nu_start_src_name, self._nu_stop_src_name = nu_src_start_stop_dict[self.name]
        #
        # self.nu_ptr_start = getattr(blob_object, self._nu_name)
        # self.nu_ptr_stop = getattr(blob_object, self._nu_name)
        #
        # self._nu_start_src = 'auto'
        # self._nu_stop_src = 'auto'
        # self._nu_start_obs = 'auto'
        # self._nu_stop_obs = 'auto'

        if name in n_seed_dic.keys():
            self.seed_field=JetSeedPhotons(name,blob_object)


        if var_name is not None:
            self._var_name=var_name

            if state_dict is None:
                self._state_dict = dict()
                self._state_dict['on'] = 1
                self._state_dict['off'] = 0
            else:
                self._state_dict=state_dict
            self.state='on'
        else:
            self._state_dict = {}
            self._var_name=None
            self._state='on'

        if state is not None and self._state_dict != {}:
            self.state=state

    # @property
    # def nu_boundaries(self,frame='obs'):
    #     check_frame(frame)
    #     if frame == 'src':
    #         return self._nu_start_src, self._nu_stop_src
    #     else:
    #         return self._nu_start_obs, self._nu_stop_obs
    #
    # @nu_boundaries.setter
    # def nu_boundaries(self,nu_start=None,nu_stop=None, frame='obs'):
    #     check_frame(frame)
    #     if frame == 'obs':
    #         self._nu_start_obs = nu
    #         self._nu_start_src = convert_nu_to_src(nu,self.jet_obj.get_par_by_type('redshift').val,'obs')
    #     else:
    #         self._nu_start_src = nu
    #         self._nu_start_obs = convert_nu_to_src(nu, self.jet_obj.get_par_by_type('redshift').val, 'obs')

    def get_emiss_lim(self,seed=False):
        return self._blob_object.emiss_lim


    def fill_SED(self,log_log=False,lin_nu=None,skip_zeros=False):

        x,y=self.get_SED_points( log_log=log_log,lin_nu=lin_nu,skip_zeros=skip_zeros)
        self.SED.beaming=self.jet_obj.get_beaming()
        self.SED.fill(nu=x,nuFnu=y,log_log=log_log)
        self.SED.fill_nuLnu(z=self.jet_obj.get_par_by_type('redshift').val,dl=self.jet_obj.get_DL_cm())

        if self.seed_field is not None:
            self.seed_field.fill(log_log=log_log)



    def get_SED_points(self, log_log=False, lin_nu=None,interp='linear',skip_zeros=False):

        size = self._blob_object.nu_grid_size
        x = zeros(size)
        y = zeros(size)

        for i in range(size):
            x[i] = BlazarSED.get_spectral_array(self.nu_ptr, self._blob_object, i)
            y[i] = BlazarSED.get_spectral_array(self.nuFnu_ptr, self._blob_object, i)


        msk_nan = np.isnan(x)
        msk_nan += np.isnan(y)

        x[msk_nan] = 0.
        y[msk_nan] = self.get_emiss_lim()

        msk = y < self.get_emiss_lim()
        y[msk] = self.get_emiss_lim()

        msk_zeros = y > self.get_emiss_lim()

        if lin_nu is not None:
            #f_interp=interpolate.Akima1DInterpolator(log10(x), log10(y))
            f_interp = interpolate.interp1d(log10(x), log10(y), bounds_error=False, kind=interp)
            y = np.power(10., f_interp(log10(lin_nu)))
            x=lin_nu
            msk_nan = np.isnan(y)
            y[msk_nan] = 0.
            msk_zeros = y > self.get_emiss_lim()
            y[~msk_zeros] = 0.

        if log_log == True:
            msk = y <= 0.
            y[msk] = -1.0E10
            msk_zeros = y >  self.get_emiss_lim()
            x = log10(x)
            y = log10(y)


        if skip_zeros is True:
            _x = x[msk_zeros]
            _y = y[msk_zeros]
        else:
            _x = x
            _y = y

        return _x, _y




    def update(self):
        size = self._blob.nu_grid_size
        x = zeros(size)
        y = zeros(size)

        for i in range(size):
            x[i] = BlazarSED.get_spectral_array(self.nu_ptr, self._blob, i)
            y[i] = BlazarSED.get_spectral_array(self.nuFnu_ptr, self._blob, i)


    def show(self):
        print('name                :',self.name)
        print('var name            :',self._var_name)
        print('state               :',self._state)
        #print('nu_start (src frame):', self._nu_start_src)
        #print('nu_stop  (src frame):', self._nu_stop_src)
        if self._state_dict is not None:
            print('allowed states :',[k for k in self._state_dict.keys()])


    @property
    def state(self,):
        return self._state

    @state.setter
    def state(self, val):
        if self._state_dict!={}:
            if val not in self._state_dict.keys():
                raise RuntimeError('val', val, 'not in allowed', self._state_dict.keys())
            self._state = val
            if self._var_name is not None:
                setattr(self._blob_object, self._var_name, self._state_dict[val])
        else:
            raise Warning('the state of the spectral component',self.name,' can not be changed')

    def get_var_state(self,):
        if self._var_name is not None:
            return  getattr(self._blob_object,self._var_name)
        else:
            return None



    def plot(self, y_min=None,y_max=None):
        p=PlotSpecComp()
        p.plot(nu=self.SED.nu.value,nuFnu=self.SED.nuFnu.value,y_min=y_min,y_max=y_max)

        return p


class SpecCompList(object):

    def __init__(self,sc_list):
        self._sc_list=sc_list
        self._table=None


    def show(self):
        for sc in self._sc_list:
            sc.show()

    def __repr__(self):
        return str(self.show())

    def build_table(self, restframe='obs'):

        _names = ['nu']
        _cols=[]

        check_frame(restframe)
        if restframe=='obs':
           _cols.append(self._sc_list[0].SED.nu)
        elif restframe=='src':
            _cols.append(self._sc_list[0].SED.nu_src)
        else:
            unexpected_behaviour()

        for ID,sc in enumerate(self._sc_list):
            _names.append(sc.name)
            if restframe == 'obs':
                _cols.append(sc.SED.nuFnu)
            else:
                _cols.append(sc.SED.nuLnu_src)



        _meta=dict(src_name=sc.jet_obj.name)
        _meta['redshift']=sc.jet_obj.get_par_by_type('redshift').val
        _meta['restframe']= restframe
        self._table = Table(_cols, names=_names,meta=_meta)

    def get_spectral_component_by_name(self,name,verbose=True):
        for i in range(len(self._sc_list)):
            if self._sc_list[i].name==name:
                return self._sc_list[i]
        else:
            if verbose==True:
                print ("no spectral components with name %s found"%name)


    @property
    def table(self):
        return self._table