
__author__ = "Andrea Tramacere"




# on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

# if on_rtd == True:
#     try:
#         from .jetkernel import jetkernel as BlazarSED
#     except ImportError:
#         from .mock import jetkernel as BlazarSED
# else:


import  warnings
import numpy as np
import time as system_time
import  ast
import copy
import dill as pickle


from tqdm.auto import tqdm

import threading

from astropy import constants as const

from .model_parameters import _show_table

from .jetkernel import jetkernel as BlazarSED
from .utils import   get_info, unexpected_behaviour

from astropy.table import Table
from astropy.units import Unit


from .jet_paramters import JetModelDictionaryPar,JetParameter,JetModelParameterArray
from .jet_emitters import  EmittersArrayDistribution
from .model_parameters import ModelParameter
from .model_manager import FitModel

from .plot_sedfit import  BasePlot,PlotPdistr,PlotTempEvDiagram,PlotTempEvEmitters, PlotSED


__all__=['JetTimeEvol','TimeEmittersDistribution', 'TimeEvolvingRegion']


class ProgressBarTempEV(object):

    def __init__(self, target_class,N):
        self.target_class = target_class
        self.N = N
        self.pbar = tqdm(total=self.N)


        self.stop = False
        self.pbar.n=0

    def update(self):
        step_tqdm = self.target_class.temp_ev.T_COUNTER +1 - self.pbar.n
        self.pbar.update(step_tqdm)

    def finalzie(self,max_try=10):
        n_try=1
        while (self.pbar.n<self.N and n_try<max_try):
            system_time.sleep(.1)
            n_try += 1
            self.update()
        self.pbar.display()

    def run(self):
        self.pbar.n = 0
        while (self.stop is False):
            self.update()
            system_time.sleep(.1)

# def _get_time_slice_from_time(time_samples_array, time):
#     ID = None
#     if time > time_samples_array[-1] or time < time_samples_array[0]:
#         warnings.warn('time must be within time start/stop range', time_samples_array[0], time_samples_array[-1])
#     else:
#         ID = (np.abs(time_samples_array - time)).argmin()
#     return ID
#
# def _get_time_from_time_slice(time_samples_array,time_slice):
#     ID = None
#     if time_slice > time_samples_array.size-1 or time_slice<0:
#         warnings.warn('time_slice must be within time 0/%d', time_samples_array.size)
#     else:
#         ID = time_slice
#     return ID
#
# def _get_time_slice_samples_from_time(time_samples_array, time,time_bin=None):
#     if time > time_samples_array[-1] or time < time_samples_array[0]:
#         warnings.warn('time must be within time start/stop range', time_samples_array[0], time_samples_array[-1])
#
#     ID = [np.argmin(np.fabs(time_samples_array - time))]
#     if time_bin is not None:
#         ID_l = np.argwhere(np.logical_and(time_samples_array >= time, time_samples_array <= time + time_bin))
#         if len(ID_l) > 0:
#             ID = ID_l.ravel()
#     return ID


class TimeEmittersDistribution(object):

    def __init__(self, jet, time_size, gamma_size):
        self.jet=jet
        self.time_blob=np.zeros(time_size)
        self.gamma=np.zeros(gamma_size)
        self.n_gamma = np.zeros((time_size, (gamma_size)) )

    #@property
    #def time_blob(self):
    #    return self.time_blob

    @property
    def time_src(self):
        return self.time_blob / self.jet.get_beaming()

    @property
    def time_obs(self):
        return self.time_blob * (1 + self.jet.parameters.z_cosm.val) / self.jet.get_beaming()

    def _time_obs_to_blob(self,time_obs):
        return time_obs / (1 + self.jet.parameters.z_cosm.val) * self.jet.get_beaming()

    def _time_src_to_blob(self,time_src):
        return time_src  * self.jet.get_beaming()

    def _time_blob_to_src(self,time_blob):
        return time_blob  / self.jet.get_beaming()
    
    def _time_blob_to_obs(self,time_obs):
        return time_obs / (1 + self.jet.parameters.z_cosm.val) * self.jet.get_beaming()

    def _get_time_slice_samples(self, time, frame, time_bin=None):
        """
        Returns the Time slice/s corresponding to a given time for *Sampled* time

        Parameters
        ----------
        time
        frame
        time_bin

        Returns
        -------

        """
        if frame == 'blob':
            t_ev_time = self.time_blob
        elif frame == 'src':
            #time=self._time_src_to_blob(time)
            t_ev_time = self.time_src
        elif frame == 'obs':
            #time=self._time_obs_to_blob(time)
            t_ev_time = self.time_obs
        else:
            raise  RuntimeError("frame must be 'blob' or 'src' or 'obs'")

        if time > t_ev_time[-1] or time < t_ev_time[0]:
            raise  RuntimeError('time in %s frame, must be within time start/stop range=[%e,%e]' % (frame, t_ev_time[0], t_ev_time[-1]))

        ID = np.array(np.argmin(np.fabs(t_ev_time - time))).ravel()
        if time_bin is not None:
            ID_l = np.argwhere(np.logical_and(t_ev_time >= time, t_ev_time <= time + time_bin))
            if len(ID_l) > 0:
                ID = ID_l.ravel()

        if ID.size == 0:
            raise  RuntimeError('time in %s frame, must be within time start/stop range=[%e,%e]' % (frame, t_ev_time.min(), t_ev_time.max()))
        return ID

    def _get_time_samples(self, time_slice,time_slice_bin=None):
        """
        Returns the blob time/s for a given time slice for *Sampled* time

        Parameters
        ----------
        time_slice
        time_slice_bin

        Returns
        -------

        """
        if time_slice == -1:
            time = np.array([self.time_blob[time_slice]])
            time_slice_out = np.array([time_slice])

        elif time_slice > self.time_blob.size - 1 or time_slice < 0:
            raise RuntimeError('time_slice must be -1, i.e. last slice, or within  [0, %d]' % self.time_blob.size)
            #warnings.warn('time_slice must be -1, i.e. last slice, or within  [0, %d]' % self.time_blob.size)
            #time = np.array([None])
            #time_slice_out =  np.array([None])
        else:
            time = np.array([self.time_blob[time_slice]])
            time_slice_out =  np.array([time_slice])


        if time_slice_bin is not None:
            time_slice_out=np.arange(time_slice, self.time_blob.size, time_slice_bin)
            time = self.time_blob[time_slice_out]

        if None in time or None in time_slice_out :
            warnings.warn('some time intervals are our of boundaries')

        return time,time_slice_out


class TimeEvolvingRegion(object):

    def __init__(self,temp_ev,jet,region_type,build_cached=True,):

        self._set_region_type(region_type)
        self._set_up(temp_ev=temp_ev,
                     jet=jet,
                     build_cached=build_cached)
        self.time_sampled_emitters=None

    def _set_region_type(self,region_type):
        self._chek_region_type(region_type)
        self._region_type = region_type

    def _chek_region_type(self,region_type):
        if region_type=='rad':
            pass
        elif region_type=='acc':
            pass
        else:
            raise  RuntimeError("region must be 'acc' or 'rad' ")

    @property
    def num_seds(self):
        if hasattr(self.temp_ev,'parameters'):
            N = self.temp_ev.parameters.num_samples.val
        else:
            N = None

        return N

    @property
    def region_type(self):
        return self._region_type

    @property
    def mult_factor(self):
        if self._region_type == 'acc':
            R, _mult_factor = self.temp_ev._get_R_acc_sphere()
        else:
            _mult_factor = 1
        return _mult_factor

    def _set_up(self,temp_ev,jet,build_cached=True,):
        self.temp_ev = temp_ev
        self.seds_array = None
        if self.num_seds is not None:
            self.time_array = np.zeros(self.num_seds, dtype=np.double)
        else:
            self.time_array = None
        self.jet = jet
        self.cached = False
        if build_cached is True:
            self.build_cached_SEDs()

    def _fill_temp_ev_array_post_run(self,):
        gamma_size = self.temp_ev._temp_ev.gamma_grid_size
        time_size = self.temp_ev.parameters.num_samples.val
        _t=np.zeros(time_size)
        for i in range(time_size):
            _t[i]=BlazarSED.get_temp_ev_N_time_array(self.temp_ev._temp_ev.N_time, self.temp_ev._temp_ev, i)            
        
        _t=_t[_t>=0.]
        time_size=_t.size
        self.temp_ev.parameters.num_samples.val=time_size
        
        self.time_sampled_emitters = TimeEmittersDistribution(jet=self.jet, time_size=time_size, gamma_size=gamma_size)

        if self.region_type=='rad':
            _ptr=self.temp_ev._temp_ev.N_rad_gamma
        elif  self.region_type=='acc':
            _ptr = self.temp_ev._temp_ev.N_acc_gamma
        else:
            unexpected_behaviour()

        for i in range(time_size):
            self.time_sampled_emitters.time_blob[i]=BlazarSED.get_temp_ev_N_time_array(self.temp_ev._temp_ev.N_time, self.temp_ev._temp_ev, i)            
            for j in range(gamma_size):
                self.time_sampled_emitters.gamma[j] = BlazarSED.get_temp_ev_gamma_array(self.temp_ev._temp_ev.gamma, self.temp_ev._temp_ev, j)
                self.time_sampled_emitters.n_gamma[i,j] = BlazarSED.get_temp_ev_N_gamma_array(_ptr,  self.temp_ev._temp_ev,i,j)



    def _update(self,temp_ev=None,jet=None,build_cached=True):
        if jet is None:
            jet=self.jet
        if temp_ev is None:
            temp_ev=self.temp_ev

        self._set_up(temp_ev=temp_ev,
                     jet=jet,
                     build_cached=build_cached)

    def _set_jet_post_run(self):
        _jet_emitters_distr = EmittersArrayDistribution(name='time_dep',
                                                             gamma_array=np.copy(self.time_sampled_emitters.gamma),
                                                             n_gamma_array=np.copy(self.time_sampled_emitters.n_gamma[-1]),
                                                             gamma_grid_size=self.temp_ev.jet_gamma_grid_size,
                                                             normalize=False)
        _jet_emitters_distr._fill()
        _jet_emitters_distr._fill()
        self.jet.set_emitters_distribution(_jet_emitters_distr)

    def _post_run_update(self,build_cached=True):
        self._set_jet_post_run()
        self._update(build_cached=build_cached)

    def set_time(self, time_slice=None, time=None, frame='blob'):
        if (time_slice is None and time is None) or (time_slice is not None and time is not None):
            raise RuntimeError('you can use either the N-th time slice, or the time in seconds')

        if time is not None:
            time_slice = self.time_sampled_emitters._get_time_slice_samples(time=time, frame=frame)[0]
            t = time
            if time_slice is None:
                raise RuntimeError('time=%e out of boundaries'%time)

            if frame == 'src':
                t = self.time_sampled_emitters._time_src_to_blob(time)
            if frame == 'obs':
                t = self.time_sampled_emitters._time_obs_to_blob(time)

        if time_slice is not None:
            time, time_ids = self.time_sampled_emitters._get_time_samples(time_slice=time_slice)
            t = time[0]
            if time is None:
                raise RuntimeError('time_slice=%d out of boundaries' % time_slice)
        if self._region_type == 'rad':
            self.jet.parameters.R.val = self.temp_ev._get_R_rad_sphere(t)
            self.jet.parameters.B.val = self.temp_ev._get_B_rad(t)
        elif self._region_type == 'acc':
            self.jet.parameters.R.val,_ = self.temp_ev._get_R_acc_sphere()
            self.jet.parameters.B.val = self.temp_ev._get_B_acc()

        self.jet.emitters_distribution._set_arrays(self.time_sampled_emitters.gamma,
                                                   self.time_sampled_emitters.n_gamma[time_slice].flatten())
        self.jet.emitters_distribution._fill()

    def get_SED(self,
                comp,
                frame,
                time_slice=None,
                time_slice_bin=None,
                time=None,
                time_bin=None,
                use_cached=False,
                average =False):

        if (time_slice is not None and time is not None):
            raise RuntimeError(
                'you can to pass either the N-th time slice "time_slice", or the blob time in seconds "time" ')

        if time_slice is None:
            time_samples = self.time_sampled_emitters._get_time_slice_samples(time, frame=frame, time_bin=time_bin)
        else:
            time_samples = [time_slice]
            if time_slice_bin is not None:
                time_samples=np.arange(time_slice, time_slice + time_slice_bin, 1)

        selected_seds_list=[]

        for ts in time_samples:
            if ts is not None:
                if use_cached is True and self.seds_array is not None:
                    selected_seds_list.append([sed for sed in self.seds_array[ts] if sed.name == comp][0])
                else:
                    self.set_time(time_slice=ts)
                    self.jet.eval()
                    for s in self.jet.spectral_components_list:
                        if s.name == comp:
                            selected_seds_list.append(s.SED)

        sed=None
        if len(selected_seds_list)>0:
            if len(selected_seds_list) == 1:
                sed = selected_seds_list[0]
            else:
                sed=copy.deepcopy(selected_seds_list[0])

                if average is True:
                    for ID in range(len(selected_seds_list)):
                        sed._nuFnu += selected_seds_list[ID]._nuFnu
                        sed._nuLnu_src += selected_seds_list[ID]._nuLnu_src


                    sed._nuFnu *= self.mult_factor/len(selected_seds_list)
                    sed._nuLnu_src *= self.mult_factor / len(selected_seds_list)

        return sed

    def build_cached_SEDs(self):
        print('caching SED for each saved distribution: start')
        self.seds_array = [None ] * self.num_seds
        pbar = tqdm(total=self.num_seds)
        pbar.n = 0
        for T in range(self.num_seds):
            self.set_time(time_slice=T)
            self.jet.eval()
            self.seds_array[T] = [copy.deepcopy(s.SED) for s in self.jet.spectral_components_list]
            pbar.update(1)
        pbar.display()
        print('caching SED for each saved distribution: done')
        self.cached=True

    def make_lc(self,
                nu1,
                nu2=None,
                comp='Sum',
                #region='rad',
                t1=None,
                t2=None,
                delta_t_out=None,
                cross_time_slices=1000,
                frame='obs',
                eval_cross_time=True,
                use_cached=False,
                R=None,
                name=None,
                density=False):

        beaming = self.jet.get_beaming()

        if frame == 'obs':
            _u = Unit('erg s-1')
            _t = self.time_sampled_emitters.time_obs
        elif frame == 'src':
            _u = Unit('erg s-1')
            _t = self.time_sampled_emitters.time_src
        elif frame == 'blob':
            _u = Unit('erg s-1 cm-2')
            _t = self.time_sampled_emitters.time_blob
        else:
            raise RuntimeError('rest_frame must be src or blob or obs')

        if t1 is None:
            t1 = _t[0]
        if t2 is None:
            t2 = _t[-1]

        selected_time_slices = np.logical_and(_t >= t1, _t <= t2)

        if delta_t_out is None:
            time_array_out = _t[selected_time_slices]
        else:
            time_array_out = np.arange(t1, t2, delta_t_out)

        time_array = _t[selected_time_slices]
        flux_array = np.zeros(time_array.shape)

        for ID, time_slice in enumerate(np.argwhere(selected_time_slices).ravel()):
            s = self.get_SED(comp, frame=frame, time_slice=time_slice, use_cached=use_cached)
            if frame == 'src':
                msk = s.nu_src.value >= nu1
                if nu2 is not None:
                    msk *= s.nu_src.value <= nu2
                x = s.nu_src[msk]
                y = s.nuLnu_src[msk] / x
            elif frame == 'obs':
                msk = s.nu.value >= nu1
                if nu2 is not None:
                    msk *= s.nu.value <= nu2
                x = s.nu[msk]
                y = s.nuFnu[msk] / x
            elif frame == 'blob':
                msk = s.nu_src.value >= nu1
                if nu2 is not None:
                    msk *= s.nu_src.value <= nu2
                x = s.nu_src[msk] / (beaming)
                y = s.nuLnu_src[msk] / s.nu_src[msk] / (beaming * beaming * beaming)
            else:
                raise RuntimeError('rest frame must be src or obs or blob')

            if density is True:
                y= y/x
            
            if nu2 is None:
                flux_array[ID] = y.value[0]
            else:
                flux_array[ID] = np.trapz(y.value, x.value)

        if density is True:
            un=Unit('erg s-1 cm-2 Hz-1')
        else:
             un=Unit('erg s-1 cm-2')

        flux_array_out = np.interp(time_array_out, time_array, flux_array, left=0, right=0)

        if eval_cross_time is True:
            t_cr = np.copy(time_array_out)
            if frame == 'src':
                t_cr_blob = self.time_sampled_emitters._time_src_to_blob(t_cr)
            if frame == 'obs':
                t_cr_blob = self.time_sampled_emitters._time_obs_to_blob(t_cr)
            flux_array_out,t_cr_blob = self.eval_cross_time(t_cr_blob, flux_array_out, n_slices=cross_time_slices, R=R)

            if frame == 'src':
                time_array_out = self.time_sampled_emitters._time_blob_to_src(t_cr_blob)
            if frame == 'obs':
                time_array_out = self.time_sampled_emitters._time_blob_to_src(t_cr_blob)

        return Table([time_array_out * Unit('s'), flux_array_out * un], names=('time', 'flux'),
                     meta={'name': name})

    def _lc_weight(self, delay_r, R, delta_R):
        return 3 * (R - delay_r) * (R + delay_r) * delta_R / (4 * R * R * R)

    def eval_cross_time(self, t, lc, n_slices=1000, R=None):
        if R is None:
            R = np.zeros(np.size(t))
            for ID,_t in enumerate(t):
                R[ID] = self.temp_ev._get_R_rad_sphere(t[ID])

        c = const.c.cgs.value
        delay_max=2 * R / c
        delta_R = (2 * R) / n_slices
        delay_t_delay = np.linspace(0, delay_max, n_slices)
        delay_r = R - delay_t_delay * c
        weight_array = self._lc_weight(delay_r, R, delta_R)       
        lc_out = np.zeros(lc.shape)
        t_out=np.zeros(t.shape)  
        for ID, lc_in in enumerate(lc_out):
            t_interp = np.linspace(t[ID] -delay_max[ID], t[ID], n_slices)
            m = np.atleast_1d( t_interp > t[0])            
            if m.sum()==0:
                lc_out[ID]=0
            else:
                lc_interp = np.interp(t_interp[m], t, lc, left=0, right=0)
                lc_out[ID] = np.sum(lc_interp * weight_array[m,ID])
            t_out[ID]=t[ID]
        return lc_out,t_out


class JetTimeEvol(object):
    """


        Parameters
        ----------
        jet
        Q_inj
        name
        inplace
        jet_gamma_grid_size
        log_sampling




        The model parameters can be set after object creation, eg:

        .. code-block:: python

            from jetset.jet_emitters_factory import EmittersFactory
            q_inj=InjEmittersFactory().create_inj_emitters('pl',emitters_type='electrons',normalize=True)
            q_inj.parameters.gmax.val=5


            from jetset.jet_timedep import  JetTimeEvol
            from jetset.jet_model import Jet

            my_jet=Jet()

            temp_ev_acc=JetTimeEvol(jet=my_jet_acc,Q_inj=q_inj,inplace=True)

            temp_ev_acc.parameters.duration.val=1E4


        model parameters
        ----------------
        duration: float (s)
            duration of the evolution in seconds

        gmin_grid: float
            lower bound of the gamma grid for the temporal evolution

        gmax_grid: float
            upper bound of the gamma grid for the temporal evolution

        gamma_grid_size: int
            size of the gamma grid

        TStart_Acc: float (s)
            starting time of the acceleration process

        TStop_Acc: float (s)
            stop time of the acceleration process

        TStart_Inj: float (s)
            starting time of the injection process

        TStop_Inj: float (s)
            stop time of the injection process

        T_esc: float (R/c)
            escape time scale in  (R/c)

        Esc_Index: float
            power-law index for the escapce term t_esc(gamma)=T_esc*gamma^(Esc_Index)

        t_D0: float (s)
            acceleration time scale for diffusive acceleration term

        t_A0: float (s)
            acceleration time scale for the systematic acceleration term

        Diff_Index: float (s)
            power-law index for the diffusive acceleration term

        Acc_Index: float (s)
            power-law index for the systematic acceleration term

        Lambda_max_Turb: float (cm)
            max scale of the turbulence

        Lambda_choer_Turb_factor: float (cm)
           max scale   for the quasi-linear resonant regime

        t_size: int
            number of steps for the time grid

        num_samples: int
            number of the time steps saved in the output

        L_inj: float (erg/s)
            if not None, the inj function is rescaled to match L_inj
    """

    def __init__(self,
                 jet_rad,
                 only_radiation=False, 
                 Q_inj=None,
                 name='jet_time_ev',
                 inplace=True,
                 log_sampling=False,
                 jet_gamma_grid_size=200,
                 setup=True):


        self._temp_ev = BlazarSED.MakeTempEv()
        if setup is True:
            self._setup_(jet_rad,Q_inj,name,log_sampling,jet_gamma_grid_size,inplace,only_radiation)

    def _setup_(self,
                jet_rad,
                Q_inj,
                name,
                log_sampling,
                jet_gamma_grid_size,
                inplace,
                only_radiation):

        if inplace is True:
            _jet_rad = jet_rad.clone()
        else:
            _jet_rad = jet_rad

        self._only_radiation=only_radiation
        self._bkp_acc_region=None
        self.rad_region = TimeEvolvingRegion(temp_ev=self, jet=_jet_rad, build_cached=False, region_type='rad')
        self.acc_region=None
        if self._only_radiation is False:
            _jet_acc = copy.deepcopy(self.rad_region.jet)
            _jet_acc.name = self.rad_region.jet.name + 'acc_region'
            self.acc_region = TimeEvolvingRegion(temp_ev=self, jet=_jet_acc, build_cached=False, region_type='acc')
        else:
            _jet_acc = copy.deepcopy(self.rad_region.jet)
            self.acc_region=None
            self._bkp_acc_region=TimeEvolvingRegion(temp_ev=self, jet=_jet_acc, build_cached=False, region_type='acc')

        self._m = FitModel()
        self._m.add_component(self.rad_region.jet)
        if self.acc_region is not None:
            self._m.add_component(self.acc_region.jet)
            self._m.parameters.link_par('z_cosm', [self.acc_region.jet.name], self.rad_region.jet.name)
        
        
        self.Q_inj=Q_inj
        self.name=name

        self._custom_q_jnj_profile = None
        self._custom_acc_profile = None
        self.time_steps_array = None
        self.IC_cooling = 'off'
        self.Sync_cooling = 'on'
        self.Adiabatic_cooling='on'
        self.region_expansion = 'off'
        self.log_sampling = log_sampling
        #self.time_sampled_emitters = None

        self.parameters = JetModelParameterArray(model=self)
        temp_ev_dict = self._build_par_dict()
        self.parameters.add_par_from_dict(temp_ev_dict,self,'_temp_ev',JetParameter)

        self.parameters.R_rad_start.val=self.rad_region.jet.parameters.R.val
        if self.acc_region is not None:
            self.parameters.B_acc.val = self.acc_region.jet.parameters.B.val
        
        self.parameters.B_rad.val = self.rad_region.jet.parameters.B.val

        self.jet_gamma_grid_size = jet_gamma_grid_size
        
    def __getstate__(self):
        return self._serialize_model()

    def __setstate__(self, state):
        self.__init__(jet_rad=None,inplace=False,setup=False)
        self._decode_model(state)

    def _serialize_model(self):

        _model = {}
        _model['version']=get_info()['version']
        _model['name'] = self.name
        _model['internals']={}
        #_model['internals']['_jet_rad']=self._jet_rad
        #_model['internals']['_jet_acc'] = self._jet_acc
        _model['internals']['Q_inj'] = self.Q_inj
        _model['internals']['_custom_q_jnj_profile'] = self._custom_q_jnj_profile
        _model['internals']['_custom_acc_profile'] = self._custom_acc_profile
        _model['internals']['_only_radiation'] = self._only_radiation
        _model['internals']['_bkp_acc_region'] = self._bkp_acc_region
        _model['internals']['time_steps_array'] = self.time_steps_array
        _model['internals']['IC_cooling'] = self.IC_cooling
        _model['internals']['Sync_cooling'] = self.Sync_cooling
        _model['internals']['Adiabatic_cooling'] = self.Adiabatic_cooling
        _model['internals']['region_expansion'] = self.region_expansion
        _model['internals']['log_sampling'] = self.log_sampling
        #_model['internals']['time_sampled_emitters'] = self.time_sampled_emitters
        _model['internals']['acc_region'] = self.acc_region
    
        _model['internals']['rad_region'] = self.rad_region
        #_model['internals']['_acc_seds_mult_factor'] = self._acc_seds_mult_factor
        _model['internals']['jet_gamma_grid_size'] = self.jet_gamma_grid_size

        _model['pars'] = {}
        _model['pars']=self.parameters._serialize_pars()

        return _model

    def _decode_model(self,_model):

        #if 'version' in _model.keys():
        #    self._set_version(_model['version'])
        #else:
        #    self._set_version(v='unknown(<1.1.2)')
        
        self.name = _model['name']
        self.acc_region=_model['internals']['acc_region']
        self.rad_region=_model['internals']['rad_region']
        self._only_radiation=_model['internals']['_only_radiation']
        self._bkp_acc_region=_model['internals']['_bkp_acc_region'] 
        self.parameters = JetModelParameterArray(model=self)
        temp_ev_dict = self._build_par_dict()
        self.parameters.add_par_from_dict(temp_ev_dict, self, '_temp_ev', JetParameter)

        _par_dict = _model['pars']
        _non_user_dict={}
        _user_dict={}
        for k, v in _model['pars'].items():
            if v['par_type'] == 'user_defined':
                _user_dict[k]=v
            else:
                _non_user_dict[k]=v
        self.parameters._decode_pars(_non_user_dict)

        for k, v in _user_dict.items():
            v['name']=k
            self.parameters.add_par(ModelParameter(**v))

        _par_dict = _model['internals']
        for k in _par_dict.keys():
            #print('k,v',k,v)
            setattr(self,k,_par_dict[str(k)])

    def save_model(self, file_name):

        pickle.dump(self, open(file_name, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)


    @classmethod
    def load_model(cls, file_name):
        try:
            c = pickle.load(open(file_name, "rb"))
            c.init_TempEv()
            return c
        except Exception as e:
            raise RuntimeError(e)


    def get_region(self, region):
        if region == 'rad':
            _reg = self.rad_region
        elif region == 'acc':
            _reg = self.acc_region
        else:
            raise RuntimeError('region must be rad or acc')
        return _reg

    def init_TempEv(self):
        BlazarSED.Init(self.rad_region.jet._blob, self.rad_region.jet.get_DL_cm())
        if self.acc_region is not None:
            BlazarSED.Init(self.acc_region.jet._blob, self.acc_region.jet.get_DL_cm())
        else:
             BlazarSED.Init(self._bkp_acc_region.jet._blob, self._bkp_acc_region.jet.get_DL_cm())

        self.temp_ev.R_H_rad_start = self.parameters.R_H_rad_start.val
        self._init_temp_ev()
        self._custom_q_jnj_profile=None
        self._custom_acc_profile=None
        self._set_inj_time_profile(user_defined_array=self._custom_q_jnj_profile)
        self._set_acc_time_profile(user_defined_array=self._custom_acc_profile)
        self._fill_temp_ev_array_pre_run()
        self.tempev_table
       
    def run(self,
            only_injection=True, 
            do_injection=True,
            cache_SEDs_rad=True, 
            cache_SEDs_acc=False):
        """

        Parameters
        ----------
        only_injection : if TRUE radiative regions is considered empty and is filled by particles escaped from acceleration region,
        otherwise particle escaped from acc regions are injected and mixed with pre-existing particles in the radiative region.
        if Q_inj is not provided in the object, then  only_injection will be considered FALSE
        do_injection : if TRUE, particles are injected in the acceleration region, otherwise acc region is skipped, and only particles
        in the radiative region are cooled.  if Q_inj is not provided in the object, then  do_injection will be considered FALSE
        cache_SEDs_rad : if TRUE, all the SEDs of the radiative region will be cached
        cache_SEDs_acc : if TRUE, all the SEDs of the acceleration region will be cached

        Returns
        -------

        """
        self.cache_SEDs_rad = cache_SEDs_rad
        if self.acc_region is not None:
            self.cache_SEDs_acc = cache_SEDs_acc
        print('temporal evolution running')
        self.init_TempEv()
        pbar=ProgressBarTempEV(target_class=self, N=self.parameters.t_size.val)
        t1 = threading.Thread(target=pbar.run)
        t1.start()
        do_injection = all((do_injection ,(self.Q_inj!=None)))
        only_injection = all((only_injection ,(self.Q_inj!=None)))
        if self._only_radiation is True and self._only_radiation is True:
            do_injection=2
        elif self._only_radiation is False and do_injection is True:
            do_injection=1
        if self.acc_region is not None:
            BlazarSED.Run_temp_evolution(self.rad_region.jet._blob, self.acc_region.jet._blob, self._temp_ev, int(only_injection), int(do_injection))
        else:
            BlazarSED.Run_temp_evolution(self.rad_region.jet._blob, self._bkp_acc_region.jet._blob, self._temp_ev, int(only_injection), int(do_injection))
        pbar.stop = True
        pbar.finalzie()
        #self._fill_temp_ev_array_post_run()
        self.rad_region._fill_temp_ev_array_post_run()
        if self.acc_region is not None:
            self.acc_region._fill_temp_ev_array_post_run()
        print('temporal evolution completed')

        #self._set_jet_post_run()
        self.rad_region._post_run_update(build_cached=cache_SEDs_rad)
        if self.acc_region is not None:
            self.acc_region._post_run_update(build_cached=cache_SEDs_acc)


    def _init_temp_ev(self):
        self.time_steps_array = np.linspace(0., self._temp_ev.duration, self.parameters.t_size.val, dtype=np.double)
        if self.Q_inj is not None:
            setattr(self._temp_ev,'Q_inj_jetset_gamma_grid_size',self.Q_inj._gamma_grid_size)
        BlazarSED.Init_Q_inj(self._temp_ev)
        if self.acc_region is not None:
            BlazarSED.Init_temp_evolution(self.rad_region.jet._blob, self.acc_region.jet._blob, self._temp_ev, self.rad_region.jet.get_DL_cm())
        else:
            BlazarSED.Init_temp_evolution(self.rad_region.jet._blob, self._bkp_acc_region.jet._blob, self._temp_ev, self.rad_region.jet.get_DL_cm())

        if self.Q_inj is not None:
            if self.acc_region is not None:
                self.Q_inj._set_L_inj(self.parameters.L_inj.val,self.V_acc())
            else:
                self.Q_inj._set_L_inj(self.parameters.L_inj.val,self.V_rad())
            Ne_custom_ptr = getattr(self._temp_ev, 'Q_inj_jetset')
            gamma_custom_ptr = getattr(self._temp_ev, 'gamma_inj_jetset')
            for ID in range(self.Q_inj._gamma_grid_size):
                BlazarSED.set_q_inj_user_array(gamma_custom_ptr, self._temp_ev, self.Q_inj.gamma_e[ID], ID)
                BlazarSED.set_q_inj_user_array(Ne_custom_ptr, self._temp_ev, self.Q_inj.n_gamma_e[ID], ID)

    def _fill_temp_ev_array_pre_run(self):
        size = BlazarSED.static_ev_arr_grid_size
        self.gamma_pre_run=np.zeros(size)
        self.t_Sync_cool_pre_run = np.zeros(size)
        self.t_D_pre_run = np.zeros(size)
        self.t_DA_pre_run = np.zeros(size)
        self.t_A_pre_run = np.zeros(size)
        self.t_Esc_rad_pre_run = np.zeros(size)
        self.t_Esc_acc_pre_run = np.zeros(size)
        self.R_H_t_pre_run = np.zeros(size)
        self.R_t_pre_run = np.zeros(size)
        self.B_t_pre_run = np.zeros(size)
        for i in range(size):
            self.gamma_pre_run[i] = BlazarSED.get_temp_ev_array_static(self._temp_ev.g, i)
            self.t_Sync_cool_pre_run[i] = BlazarSED.get_temp_ev_array_static(self._temp_ev.t_Sync_cool, i)
            self.t_D_pre_run[i] = BlazarSED.get_temp_ev_array_static(self._temp_ev.t_D, i)
            self.t_DA_pre_run[i] = BlazarSED.get_temp_ev_array_static(self._temp_ev.t_DA, i)
            self.t_A_pre_run[i] = BlazarSED.get_temp_ev_array_static(self._temp_ev.t_A, i)
            self.t_Esc_acc_pre_run[i] = BlazarSED.get_temp_ev_array_static(self._temp_ev.t_Esc_acc, i)
            self.t_Esc_rad_pre_run[i] = BlazarSED.get_temp_ev_array_static(self._temp_ev.t_Esc_rad, i)
            self.R_H_t_pre_run[i] =  BlazarSED.get_temp_ev_array_static(self._temp_ev.R_H_t_pre, i)
            self.R_t_pre_run[i] =  BlazarSED.get_temp_ev_array_static(self._temp_ev.R_t_pre, i)
            self.B_t_pre_run[i] =  BlazarSED.get_temp_ev_array_static(self._temp_ev.B_t_pre, i)

        for i in range(self.parameters.t_size.val):
            T_inj_profile_ptr = getattr(self._temp_ev, 'T_inj_profile')
            BlazarSED.set_temp_ev_Time_array(T_inj_profile_ptr, self._temp_ev, self.custom_q_jnj_profile[i], i)
            Acc_profile_ptr = getattr(self._temp_ev, 'T_acc_profile')
            BlazarSED.set_temp_ev_Time_array(Acc_profile_ptr, self._temp_ev, self.custom_acc_profile[i], i)

    # def _fill_temp_ev_array_post_run(self):
    #     if self.acc_region is not None:
    #         self.acc_region._set_time_sampled_emitters('acc')
    #     self.rad_region._set_time_sampled_emitters('rad')
    #
    #     # gamma_size = self._temp_ev.gamma_grid_size
    #     # time_size = self.parameters.num_samples.val
    #     # self.time_sampled_emitters=TimeEmittersDistribution(time_size=time_size, gamma_size=gamma_size)
    #     #
    #     # self.time_sampled_emitters.time = np.zeros(time_size)
    #     # #self.gamma = np.zeros(gamma_size)
    #     # #self.N_gamma = np.zeros((time_size,(gamma_size)))
    #     # for i in range(time_size):
    #     #     self.time_sampled_emitters.time[i]=BlazarSED.get_temp_ev_N_time_array(self._temp_ev.N_time, self._temp_ev, i)
    #     #     for j in range(gamma_size):
    #     #         self.time_sampled_emitters.gamma[j] = BlazarSED.get_temp_ev_gamma_array(self._temp_ev.gamma, self._temp_ev, j)
    #     #         self.time_sampled_emitters.n_gamma_rad[i,j] = BlazarSED.get_temp_ev_N_gamma_array(self._temp_ev.N_rad_gamma, self._temp_ev,i,j)
    #     #         self.time_sampled_emitters.n_gamma_acc[i, j] = BlazarSED.get_temp_ev_N_gamma_array(self._temp_ev.N_acc_gamma,self._temp_ev, i, j)

    # def _set_jet_post_run(self):
    #     self._jet_emitters_distr = EmittersArrayDistribution(name='time_dep',
    #                                                          gamma_array=np.copy(self.time_sampled_emitters.gamma),
    #                                                          n_gamma_array=np.copy(self.time_sampled_emitters.n_gamma_rad[-1]),
    #                                                          gamma_grid_size=self.jet_gamma_grid_size,
    #                                                          normalize=False)
    #     self._jet_emitters_distr._fill()
    #     self._jet_emitters_distr._fill()
    #     self._jet_rad.set_emitters_distribution(self._jet_emitters_distr)
    #     self._jet_emitters_distr = EmittersArrayDistribution(name='time_dep',
    #                                                          gamma_array=np.copy(self.time_sampled_emitters.gamma),
    #                                                          n_gamma_array=np.copy(
    #                                                              self.time_sampled_emitters.n_gamma_acc[-1]),
    #                                                          gamma_grid_size=self.jet_gamma_grid_size,
    #                                                          normalize=False)
    #     self.acc_region.jet.set_emitters_distribution(self._jet_emitters_distr)

    @property
    def Delta_R_acc(self):
        return self.parameters.Delta_R_acc.val

    @property
    def time_steps_array(self):
        return self._time_steps_array

    @time_steps_array.setter
    def time_steps_array(self, v):
        self._time_steps_array = v

    @property
    def temp_ev(self):
        return self._temp_ev

    def _get_R_rad_sphere(self, time):
        R = self.parameters.R_rad_start.val
        if self.region_expansion == 'on':
            R = BlazarSED.eval_R_jet_t(self.rad_region.jet._blob, self.temp_ev, time)

        return R

    def _get_B_rad(self, time):
        if self.region_expansion == 'on':
            R = self._get_R_rad_sphere(time)
            return BlazarSED.eval_B_jet_t(self.rad_region.jet._blob, self.temp_ev, R, time)
        else:
            return self.parameters.B_rad.val

    def _get_B_acc(self):
        return self.parameters.B_acc.val

    def _get_adiab_cooling_time_from_R(self, R):
        # R = BlazarSED.eval_R_jet_t(self._jet_rad._blob, self.temp_ev, R_H)
        return BlazarSED.Adiabatic_Cooling_time(self.temp_ev, self.rad_region.jet._blob, R)

    def _get_R_acc_sphere(self):
        R = self.parameters.R_rad_start.val
        V_shell = R * R * self.Delta_R_acc
        R1 = (4 / 3) * np.pi * (self.Delta_R_acc / 0.5) ** 3
        N = V_shell / R1
        return self.Delta_R_acc * 0.5, N

    def V_acc(self):
        R = self.parameters.R_rad_start.val
        return R * R * self.Delta_R_acc*np.pi

    def V_rad(self):
        R = self.parameters.R_rad_start.val
        return (4/3)*np.pi*R*R*R

    def show_pars(self, sort_key='par type'):
        self.parameters.show_pars(sort_key=sort_key)

    @property
    def log_sampling(self):
        """
        logarithmically spaced bool
        Returns
        -------

        """
        return self._log_sampling

    @log_sampling.setter
    def log_sampling(self, v):
        """
        if True, the time grid is logarithmically spaced
        Returns
        -------

        """
        if v is False or v is True:
            pass
        else:
            raise RuntimeError('this parameter must be bolean')

        self._log_sampling = v
        self.temp_ev.LOG_SET= np.int(v)

    @property
    def t_unit_rad(self):
        """
        blob time in R/c
        Returns
        -------

        """
        return self._temp_ev.t_unit_rad

    @property
    def t_unit_acc(self):
        """
        blob time in R/c
        Returns
        -------

        """
        return self._temp_ev.t_unit_acc

    @property
    def delta_t(self):
        """
        blob delta t in s
        Returns
        -------

        """
        return self._temp_ev.deltat

    @property
    def IC_cooling(self):
        return self._IC_cooling

    @IC_cooling.setter
    def IC_cooling(self, val):
        state_dict = dict((('on', 1), ('off', 0)))
        if val not in state_dict.keys():
            raise  RuntimeError('allowed values are on/off')
        self._IC_cooling=val
        self.temp_ev.do_Compton_cooling=state_dict[val]

    @property
    def region_expansion(self):
        return self._region_expansion

    @region_expansion.setter
    def region_expansion(self, val):
        state_dict = dict((('on', 1), ('off', 0)))
        if val not in state_dict.keys():
            raise RuntimeError('allowed values are on/off')
        self._region_expansion = val
        self.temp_ev.do_Expansion = state_dict[val]

    @property
    def Sync_cooling(self):
        return self._Sync_cooling

    @Sync_cooling.setter
    def Sync_cooling(self, val):
        state_dict = dict((('on', 1), ('off', 0)))
        if val not in state_dict.keys():
            raise RuntimeError('allowed values are on/off')
        self._Sync_cooling = val
        self.temp_ev.do_Sync_cooling = state_dict[val]

    @property
    def Adiabatic_cooling(self):
        return self._Adiabatic_cooling

    @Sync_cooling.setter
    def Adiabatic_cooling(self, val):
        state_dict = dict((('on', 1), ('off', 0)))
        if val not in state_dict.keys():
            raise RuntimeError('allowed values are on/off')
        self._Adiabatic_cooling = val
        self.temp_ev.do_Adiabatic_cooling = state_dict[val]


    @property
    def custom_q_jnj_profile(self):
        return self._custom_q_jnj_profile

    @custom_q_jnj_profile.setter
    def custom_q_jnj_profile(self,user_defined_array):
        self._set_inj_time_profile(user_defined_array)

    @property
    def custom_acc_profile(self):
        return self._custom_acc_profile

    @custom_acc_profile.setter
    def custom_acc_profile(self, user_defined_array):
        self._set_acc_time_profile(np.double(user_defined_array>0))

    def _set_inj_time_profile(self,user_defined_array=None):
        if user_defined_array is None:
            self._custom_q_jnj_profile = np.zeros(self.parameters.t_size.val, dtype=np.double)
            msk= self.time_steps_array >= self._temp_ev.TStart_Inj
            msk*= self.time_steps_array <= self._temp_ev.TStop_Inj
            self._custom_q_jnj_profile[msk] = 1.0
        else:
            if np.shape(user_defined_array)!=(self.parameters.t_size.val,):
                raise  RuntimeError('user_defined_array must be 1d array with size =',self._temp_ev.T_SIZE)
            self._custom_q_jnj_profile = np.double(user_defined_array)

    def _set_acc_time_profile(self,user_defined_array=None):
        self._custom_acc_profile = np.zeros(self._temp_ev.T_SIZE, dtype=np.double)

        if user_defined_array is None and self._only_radiation is False:
            self._custom_acc_profile = np.zeros(self._temp_ev.T_SIZE, dtype=np.double)
            msk= self.time_steps_array >= self._temp_ev.TStart_Acc
            msk*= self.time_steps_array <= self._temp_ev.TStop_Acc
            self._custom_acc_profile[msk] = 1.0
        elif user_defined_array is not None and self._only_radiation is False:
            if np.shape(user_defined_array)!=(self.parameters.t_size.val,):
                raise  RuntimeError('user_defined_array must be 1d array with size =',self._temp_ev.T_SIZE)
            self._custom_acc_profile = np.double(user_defined_array)

    # def get_SED(self,
    #             comp,
    #             region='rad',
    #             time_slice=None,
    #             time_slice_bin=None,
    #             time=None,
    #             time_bin=None,
    #             use_cached=False,
    #             average=False,
    #             frame='obs'):
    #
    #     if region == 'rad':
    #         sed = self.rad_region.get_SED(comp=comp,
    #                                       time_slice=time_slice,
    #                                       time_slice_bin=time_slice_bin,
    #                                       time=time,
    #                                       time_bin=time_bin,
    #                                       use_cached=use_cached,
    #                                       average=average,
    #                                       frame=frame)
    #     elif region == 'acc':
    #         sed = self.acc_region.get_SED(comp=comp,
    #                                       time_slice=time_slice,
    #                                       time_slice_bin=time_slice_bin,
    #                                       time=time,
    #                                       time_bin=time_bin,
    #                                       use_cached=use_cached,
    #                                       average=average,
    #                                       frame=frame)
    #     else:
    #         raise RuntimeError("region must be 'acc' or 'rad'")
    #     return sed

    # def make_lc(self,
    #             nu1,
    #             nu2=None,
    #             comp='Sum',
    #             region='rad',
    #             t1=None,
    #             t2=None,
    #             delta_t_out=None,
    #             cross_time_slices=1000,
    #             rest_frame='obs',
    #             eval_cross_time=True,
    #             use_cached=False,
    #             R=None,
    #             name=None):
    #
    #     return self.get_region(region).make_lc(nu1,
    #                                     nu2 = nu2,
    #                                     comp = comp,
    #                                     region = region,
    #                                     t1 = t1,
    #                                     t2 = t2,
    #                                     delta_t_out = delta_t_out,
    #                                     cross_time_slices = cross_time_slices,
    #                                     rest_frame = rest_frame,
    #                                     eval_cross_time = eval_cross_time,
    #                                     use_cached = use_cached,
    #                                     R = R,
    #                                     name = name)

    def eval_L_tot_inj(self):
        if self.Q_inj is not None and self.acc_region is not None:
            return self.Q_inj.eval_U_q() * self.V_acc()
        elif self.Q_inj is not None and self._only_radiation is True:
            return self.Q_inj.eval_U_q() * self.V_rad()
        else:
            return None

    def _get_time_slice_T_array(self, time_blob):
        """
        Returns the time slice for a given time for *Non-sampled* times

        Parameters
        ----------
        time_blob

        Returns
        -------

        """
        r = min(0, np.int(np.log10(self.delta_t))-1)
        if r<0:
            _time = np.round(time_blob, -r)
        else:
            _time = time_blob

        if _time== -1:
            _id = -1

        if _time > self.time_steps_array[-1] or _time < self.time_steps_array[0]:
            raise RuntimeError('time=%e must be within time start/stop range =[%e,%e]'%(_time,self.time_steps_array[0], self.time_steps_array[-1]))

        dt = (self.time_steps_array[-1] - self.time_steps_array[0]) / self.time_steps_array.size
        _id = np.int(_time/dt)

        if _id<0:
            _id =  0

        if _id > self.time_steps_array.size - 1:
            _id = self.time_steps_array.size - 1

        return _id

    def set_time(self, time_slice=None, time=None, frame='blob'):
        self.rad_region.set_time(time_slice=time_slice,time=time,frame=frame)
        if self.acc_region is not None:
            self.acc_region.set_time(time_slice=time_slice, time=time, frame=frame)

    # def set_time(self,time_slice=None,time=None,frame='blob',region='rad'):
    #     if (time_slice is None and time is None) or (time_slice is not None and time is not None):
    #         raise RuntimeError('you can use either the N-th time slice, or the time in seconds')
    #
    #
    #     region = self.get_region(region)
    #
    #     if time is not None:
    #         if frame == 'obs':
    #             time = time
    #         elif frame == 'src':
    #             _t = self.temp_ev.time_sampled_emitters.time_src
    #
    #         elif frame == 'blob':
    #             pass
    #         else:
    #             raise RuntimeError('rest_frame must be src or blob or obs')
    #
    #
    #         time_slice=region.time_sampled_emitters._get_time_slice_samples(time)[0]
    #
    #     if time_slice is not  None:
    #         time,time_ids = region.time_sampled_emitters._get_time_samples(time_slice=time_slice)
    #         time = time[0]
    #
    #     self.rad_region.jet.parameters.R.val=self._get_R_rad_sphere(time)
    #     self.rad_region.jet.parameters.B.val = self._get_B_rad(time)
    #     self.rad_region.jet.emitters_distribution._set_arrays(self.time_sampled_emitters.gamma, self.time_sampled_emitters.n_gamma_rad[time_slice].flatten())
    #     self.rad_region.jet.emitters_distribution._fill()
    #
    #     R,N_spheres = self._get_R_acc_sphere()
    #     self._acc_seds_mult_factor = N_spheres
    #     self.acc_region.jet.parameters.R.val = R
    #     self.rad_region.jet.parameters.B.val = self._get_B_acc()
    #     self.acc_region.jet.emitters_distribution._set_arrays(self.time_sampled_emitters.gamma,
    #                                                    self.time_sampled_emitters.n_gamma_acc[time_slice].flatten())
    #     self.acc_region.jet.emitters_distribution._fill()


    def plot_time_profile(self,figsize=(8,8),dpi=120):
        p=PlotTempEvDiagram(figsize=figsize,dpi=dpi,expanding_region=self.region_expansion=='on')
        p.plot(self.time_steps_array,
               self.custom_q_jnj_profile,
               self.custom_acc_profile,
               self.R_t_pre_run,
               self.B_t_pre_run,
               self.R_H_t_pre_run)
        return p

    def plot_tempev_emitters(self,region='rad',figsize=(8,8),dpi=120,energy_unit='gamma',loglog=True,plot_Q_inj=True,pow=None):
        region=self.get_region(region)

        p=PlotTempEvEmitters(figsize=figsize,dpi=dpi,loglog=loglog)
        p.plot_distr(temp_ev=self,region=region,energy_unit=energy_unit,plot_Q_inj=plot_Q_inj,pow=pow)

        return p

    def plot_tempev_model(self,
                          comp='Sum',
                          region='rad',
                          frame='obs',
                          t1=None,
                          t2=None,
                          time_slice=None,
                          time_slice_bin=None,
                          time=None,
                          time_bin=None,
                          density=False,
                          use_cached=False,
                          sed_data=None,
                          plot_obj=None,
                          average=False):


        if plot_obj is None:
            plot_obj=PlotSED(frame=frame)

        region = self.get_region(region)
        plot_obj.plot_tempev_model(temp_ev=self,
                                   region=region,
                                   comp=comp,
                                   t1=t1,
                                   t2=t2,
                                   time_slice=time_slice,
                                   time_slice_bin=time_slice_bin,
                                   time=time,
                                   time_bin=time_bin,
                                   density=density,
                                   sed_data=sed_data,
                                   use_cached=use_cached,
                                   average=average,)

        return plot_obj

    def plot_pre_run_plot(self,figsize=(8,6),dpi=120):
        p=BasePlot(figsize=figsize,dpi=dpi)
        p.ax.loglog(self.gamma_pre_run, self.t_Sync_cool_pre_run, label='t coool, synch.')
        if self._only_radiation is False:
            p.ax.loglog(self.gamma_pre_run, self.t_D_pre_run, label='t acc. diffusive')
            p.ax.loglog(self.gamma_pre_run, self.t_DA_pre_run, label='t acc. diffusive syst.')
            p.ax.loglog(self.gamma_pre_run, self.t_A_pre_run, label='t acc. systematic')
            p.ax.loglog(self.gamma_pre_run, self.t_Esc_acc_pre_run, label='t esc R acc.')
            p.ax.axvline(self._temp_ev.gamma_eq_t_A, ls='--')
            p.ax.axvline(self._temp_ev.gamma_eq_t_D, ls='--')
        
        p.ax.loglog(self.gamma_pre_run, self.t_Esc_rad_pre_run, label='t esc R rad.')
        if self.region_expansion=='on':
            p.ax.axhline(self._get_adiab_cooling_time_from_R(R=self.parameters.R_rad_start.val), label='t adiab (exp. start)', ls='-', color='black', lw=0.5)
            p.ax.axhline(self._get_adiab_cooling_time_from_R(R=self.R_t_pre_run.max()),
                         label='t adiab (exp. stop)', ls='-.', color='black', lw=0.5)
        p.ax.axhline(self.delta_t, label='delta t', ls=':',color='magenta',lw=0.5)
       

        p.ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), ncol=1, prop={'size':10})
        p.ax.set_xlabel('$gamma$')
        p.ax.set_ylabel('time (s)')
        return p

    def show_model(self, getstring=False, names_list=None, sort_key=None):
        print('-'*80)
        print("JetTimeEvol model description")
        print('-'*80)
       
        self.tempev_table
        print(" ")
        print("physical setup: ")
        print("")
        print('-'*80)
        _show_table(self._tempev_table)
        print("")
        print("model parameters: ")
        print("")
        print('-'*80)
        self.show_pars()

    def _build_row_dict(self, name, type='', unit='', val='', val_by=None, islog=False, unit1=''):
        row_dict = {}
        row_dict['name'] = name
        row_dict['par type'] = type
        row_dict['units'] = unit
        row_dict['val'] = val
        if val_by is not None:
            val_by = val / val_by
        row_dict['val*'] = val_by
        row_dict['units*'] = unit1
        row_dict['log'] = islog
        return row_dict

    def _build_tempev_table(self, jet_rad,jet_acc):

        _names = ['name', 'par type', 'val', 'units', 'val*', 'units*', 'log', ]

        rows = []
        rows.append(self._build_row_dict('delta t', 'time', 's', val=self.delta_t, val_by=self.t_unit_rad,
                                         unit1='R/c', islog=False))
        rows.append(self._build_row_dict('log. sampling', 'time', '', val=self.log_sampling, islog=False))
        rows.append(
            self._build_row_dict('R/c', 'time', 's', val=self.t_unit_rad, val_by=self.t_unit_rad, unit1='R/c',
                                 islog=False))
        
        rows.append(self._build_row_dict('IC cooling', '', '', val=self.IC_cooling, islog=False))
        rows.append(self._build_row_dict('Sync cooling', '', '', val=self.Sync_cooling, islog=False))
        rows.append(self._build_row_dict('Adiab. cooling', '', '', val=self.Adiabatic_cooling, islog=False))
        rows.append(self._build_row_dict('Reg. expansion', '', '', val=self.region_expansion, islog=False))

        if self.acc_region is not None:
            rows.append(self._build_row_dict('Diff coeff', '', 's-1', val=self._temp_ev.Diff_Coeff, islog=False))
            rows.append(self._build_row_dict('Acc coeff', '', 's-1', val=self._temp_ev.Acc_Coeff, islog=False))
            rows.append(self._build_row_dict('Diff index', '', '', val=self._temp_ev.Diff_Index, islog=False))
            rows.append(self._build_row_dict('Acc index', '', 's-1', val=self._temp_ev.Acc_Index, islog=False))
            rows.append(
                self._build_row_dict('Tesc acc', 'time', 's', val=self._temp_ev.T_esc_Coeff_acc, val_by=self.t_unit_acc,
                                     unit1='R_acc/c', islog=False))
            rows.append(
            self._build_row_dict('Eacc max', 'energy', 'erg', val=self._temp_ev.E_acc_max, islog=False))

        rows.append(
            self._build_row_dict('Tesc rad', 'time', 's', val=self._temp_ev.T_esc_Coeff_rad, val_by=self.t_unit_rad,
                                 unit1='R/c', islog=False))

        

        if jet_acc is not None:
            rows.append(
                self._build_row_dict('Delta R acc', 'accelerator_width', 'cm', val=self.parameters.Delta_R_acc.val))

        if jet_acc is not None:
            rows.append(
                self._build_row_dict('B acc', 'magnetic field', 'cm', val=self.parameters.B_acc.val))

        
        rows.append(
            self._build_row_dict('R_rad rad start', 'region_position', 'cm', val=self.parameters.R_rad_start.val))

        rows.append(
            self._build_row_dict('R_H rad start', 'region_position', 'cm',
                                 val=self.parameters.R_H_rad_start.val))

        if self.region_expansion=='on':
            rows.append(
                self._build_row_dict('beta exp.', 'region_position', 'v/c',
                                    val=self.parameters.beta_exp_R.val, unit1='cm/s',val_by=1.0/const.c.cgs ))

        if jet_acc is not None:
            rows.append(
                self._build_row_dict('T_A0=1/ACC_COEFF', 'time', 's', val=self._temp_ev.t_A0, val_by=self.t_unit_rad,
                                    unit1='R/c', islog=False))
            rows.append(
                self._build_row_dict('T_D0=1/DIFF_COEFF', 'time', 's', val=self._temp_ev.t_D0, val_by=self.t_unit_rad,
                                    unit1='R/c', islog=False))
            rows.append(self._build_row_dict('T_DA0=1/(2*DIFF_COEFF)', 'time', 's', val=self._temp_ev.t_DA0,
                                            val_by=self.t_unit_rad, unit1='R/c', islog=False))

            rows.append(self._build_row_dict('gamma Lambda Turb.  max', '', '', val=self._temp_ev.Gamma_Max_Turb_L_max,
                                            islog=False))
            rows.append(self._build_row_dict('gamma Lambda Coher. max', '', '', val=self._temp_ev.Gamma_Max_Turb_L_coher,
                                            islog=False))

            rows.append(self._build_row_dict('gamma eq Syst. Acc (synch. cool)', '', '', val=self._temp_ev.gamma_eq_t_A,
                                            islog=False))
            rows.append(self._build_row_dict('gamma eq Diff. Acc (synch. cool)', '', '', val=self._temp_ev.gamma_eq_t_D,
                                            islog=False))

        if jet_acc is not None:
            t = BlazarSED.Sync_tcool(jet_acc._blob, self._temp_ev.gamma_eq_t_D)
            rows.append(self._build_row_dict('T cooling(gamma_eq=gamma_eq_Diff)', '', 's', val=t, islog=False))
            t = BlazarSED.Sync_tcool(jet_acc._blob, self._temp_ev.gamma_eq_t_A)
            rows.append(self._build_row_dict('T cooling(gamma_eq=gamma_eq_Sys)', '', 's', val=t, islog=False))
            t = BlazarSED.Sync_tcool(jet_acc._blob, self._temp_ev.gamma_eq_t_A)
        rows.append(
            self._build_row_dict('T min. synch. cooling', '', 's', val=self.t_Sync_cool_pre_run.min(), islog=False))

        if self.Q_inj is not None:
            rows.append(self._build_row_dict('L inj (electrons)', 'injected lum.', 'erg/s', val=self.eval_L_tot_inj(),
                                             islog=False))
            #rows.append(self._build_row_dict('E_inj (electrons)', '', 'erg', val=self.eval_E_tot_inj(), islog=False))

        self._tempev_table = Table(rows=rows, names=_names, masked=False)

        self._fromat_column_entry(self._tempev_table)

    @property
    def tempev_table(self):
        if self.acc_region is not None:
            self._build_tempev_table(self.rad_region.jet, self.acc_region.jet)
        else:
            self._build_tempev_table(self.rad_region.jet, None)

    def _fromat_column_entry(self, t):

        for n in t.colnames:

            for ID, v in enumerate(t[n].data):
                try:
                    c = ast.literal_eval(t[n].data[ID])
                    if type(c) == int:
                        t[n].data[ID] = '%d' % c
                    else:
                        t[n].data[ID] = '%e' % c
                except:
                    pass

    def _build_par_dict(self):
        """
        These are the model parameters to set after object creation:
        eg:
            temp_ev_acc.parameters.duration.val=1E4

        model parameters
        ----------------
        duration: float (s)
            duration of the evolution in seconds

        gmin_grid: float
            lower bound of the gamma grid for the temporal evolution

        gmax_grid: float
            upper bound of the gamma grid for the temporal evolution

        gamma_grid_size: int
            size of the gamma grid

        TStart_Acc: float (s)
            starting time of the acceleration process

        TStop_Acc: float (s)
            stop time of the acceleration process

        TStart_Inj: float (s)
            starting time of the injection process

        TStop_Inj: float (s)
            stop time of the injection process

        T_esc: float (R/c)
            escape time scale in  (R/c)

        Esc_Index: float
            power-law index for the escapce term t_esc(gamma)=T_esc*gamma^(Esc_Index)

        t_D0: float (s)
            acceleration time scale for diffusive acceleration term

        t_A0: float (s)
            acceleration time scale for the systematic acceleration term

        Diff_Index: float (s)
            power-law index for the diffusive acceleration term

        Acc_Index: float (s)
            power-law index for the systematic acceleration term

        Lambda_max_Turb: float (cm)
            max scale of the turbulence

        Lambda_choer_Turb_factor: float (cm)
           max scale   for the quasi-linear resonant regime

        t_size: int
            number of steps for the time grid

        num_samples: int
            number of the time steps saved in the output

        log_sampling: bool
           if True, the time grid is logarithmically spaced

        L_inj: float (erg/s)
            if not None, the inj function is rescaled to match L_inj
        """

        model_dict_temp_ev = {}

        model_dict_temp_ev['duration'] = JetModelDictionaryPar(ptype='time_grid', vmin=0, vmax=None, punit='s', froz=True, log=False,val=1E5)

        model_dict_temp_ev['gmin_grid'] = JetModelDictionaryPar(ptype='gamma_grid', vmin=0, vmax=None, punit='', froz=True, log=False, val=1E1,jetkernel_par_name='gmin_griglia')

        model_dict_temp_ev['gmax_grid'] = JetModelDictionaryPar(ptype='gamma_grid', vmin=0, vmax=None, punit='', froz=True, log=False, val=1E8,jetkernel_par_name='gmax_griglia')

        model_dict_temp_ev['gamma_grid_size'] = JetModelDictionaryPar(ptype='gamma_grid', vmin=0, vmax=None, punit='', froz=True, log=False, val=5000,jetkernel_par_name='gamma_grid_size')


        if self.acc_region is not None:
            model_dict_temp_ev['TStart_Acc'] = JetModelDictionaryPar(ptype='time_grid', vmin=0, vmax=None, punit='s', froz=True,
                                                            log=False,val=0)

            model_dict_temp_ev['TStop_Acc'] = JetModelDictionaryPar(ptype='time_grid', vmin=0, vmax=None, punit='s', froz=True,
                                                            log=False,val=1E5)

        model_dict_temp_ev['TStart_Inj'] = JetModelDictionaryPar(ptype='time_grid', vmin=0, vmax=None, punit='s', froz=True,
                                                         log=False,val=0)

        model_dict_temp_ev['TStop_Inj'] = JetModelDictionaryPar(ptype='time_grid', vmin=0, vmax=None, punit='s', froz=True,
                                                         log=False,val=1E5)

        if self.acc_region is not None:
            model_dict_temp_ev['T_esc_acc'] = JetModelDictionaryPar(ptype='escape_time', vmin=None, vmax=None, punit='(R_acc/c)', froz=True,
                                                         log=False,val=2.0,jetkernel_par_name='T_esc_Coeff_R_by_c_acc')
            
            model_dict_temp_ev['Esc_Index_acc'] = JetModelDictionaryPar(ptype='fp_coeff_index', vmin=None, vmax=None, punit='', froz=True,
                                                            log=False,val=0.0)                                 

            model_dict_temp_ev['Esc_Index_acc'] = JetModelDictionaryPar(ptype='fp_coeff_index', vmin=None, vmax=None, punit='', froz=True,
                                                            log=False,val=0.0) 
            
            model_dict_temp_ev['t_D0'] = JetModelDictionaryPar(ptype='acceleration_time', vmin=0, vmax=None, punit='s', froz=True,
                                                            log=False,val=1E4)

            model_dict_temp_ev['t_A0'] = JetModelDictionaryPar(ptype='acceleration_time', vmin=0, vmax=None, punit='s', froz=True,
                                                                log=False,val=1E3)

            model_dict_temp_ev['Diff_Index'] = JetModelDictionaryPar(ptype='fp_coeff_index', vmin=0, vmax=None, punit='s', froz=True,
                                                        log=False,val=2)

            model_dict_temp_ev['Acc_Index'] = JetModelDictionaryPar(ptype='fp_coeff_index', vmin=None, vmax=None, punit='', froz=True,
                                                        log=False,val=1)

            model_dict_temp_ev['Delta_R_acc'] = JetModelDictionaryPar(ptype='accelerator_width', vmin=0, vmax=None, punit='cm',
                                                                jetkernel_par_name='Delta_R_acc',
                                                                froz=True,
                                                                log=False, val=1E13)

            model_dict_temp_ev['B_acc'] = JetModelDictionaryPar(ptype='magnetic_field', vmin=0, vmax=None,
                                                              punit='G',
                                                              jetkernel_par_name='B_acc',
                                                              froz=True,
                                                              log=False, val=0.1)
            
            model_dict_temp_ev['E_acc_max'] = JetModelDictionaryPar(ptype='acc_energy', vmin=0, vmax=None, punit='erg',
                                                       froz=True,
                                                       log=False, val=1E60)
            
            model_dict_temp_ev['Lambda_max_Turb'] = JetModelDictionaryPar(ptype='turbulence_scale', vmin=0, vmax=None, punit='cm', froz=True,
                                                          log=False,val=1E15)

            model_dict_temp_ev['Lambda_choer_Turb_factor'] = JetModelDictionaryPar(ptype='turbulence_scale', vmin=0, vmax=None, punit='cm',
                                                                    froz=True,log=False, val=0.1)



        model_dict_temp_ev['T_esc_rad'] = JetModelDictionaryPar(ptype='escape_time', vmin=None, vmax=None, punit='(R/c)', froz=True,
                                                   log=False, val=2.0, jetkernel_par_name='T_esc_Coeff_R_by_c_rad')

       

        model_dict_temp_ev['Esc_Index_rad'] = JetModelDictionaryPar(ptype='fp_coeff_index', vmin=None, vmax=None, punit='',
                                                                froz=True,
                                                                log=False, val=0.0)

      

        model_dict_temp_ev['R_rad_start'] = JetModelDictionaryPar(ptype='region_size', vmin=0, vmax=None,
                                                            punit='cm',
                                                            jetkernel_par_name='R_rad_start',
                                                            froz=True,
                                                            log=False, val=1E15)

        model_dict_temp_ev['R_H_rad_start'] = JetModelDictionaryPar(ptype='region_position', vmin=0, vmax=None,
                                                                  punit='cm',
                                                                  jetkernel_par_name='R_H_rad_start',
                                                                  froz=True,
                                                                  log=False, val=1E17)

        model_dict_temp_ev['m_B'] = JetModelDictionaryPar(ptype='magnetic_field_index', vmin=1, vmax=2, punit='',
                                                                froz=True,
                                                                log=False, val=1)

       
        model_dict_temp_ev['t_jet_exp'] = JetModelDictionaryPar(ptype='exp_start_time', vmin=0, vmax=None,
                                                            punit='s',
                                                            jetkernel_par_name='t_jet_exp',
                                                            froz=True,
                                                            log=False, val=1E14)

        model_dict_temp_ev['beta_exp_R'] = JetModelDictionaryPar(ptype='beta_expansion', vmin=0, vmax=1,
                                                                punit='v/c',
                                                                jetkernel_par_name='v_exp_by_c',
                                                                froz=True,
                                                                log=False, val=0)

       

        model_dict_temp_ev['B_rad'] = JetModelDictionaryPar(ptype='magnetic_field', vmin=0, vmax=None,
                                                            punit='G',
                                                            jetkernel_par_name='B_rad',
                                                            froz=True,
                                                            log=False, val=0.1)


        

     
        model_dict_temp_ev['t_size'] = JetModelDictionaryPar(ptype='time_grid', vmin=0, vmax=None,
                                                                         punit='',
                                                                         froz=True,
                                                                         log=False,val=1000,
                                                                        jetkernel_par_name='T_SIZE')

        model_dict_temp_ev['num_samples'] = JetModelDictionaryPar(ptype='time_ev_output', vmin=0, vmax=None, punit='', froz=True, log=False, val=50,jetkernel_par_name='NUM_SET')

        #model_dic['log_sampling'] = JetModelDictionaryPar(ptype='time_ev_output', vmin=0, vmax=None, punit='', froz=True, log=False, val=0,allowed_values=[0,1],jetkernel_par_name='LOG_SET')

        #model_dic['IC_cooling'] = JetModelDictionaryPar(ptype='time_ev_output', vmin=0, vmax=None, punit='', froz=True, log=False, val=0,allowed_values=[0,1],jetkernel_par_name='do_Compton_cooling')

        model_dict_temp_ev['L_inj'] = JetModelDictionaryPar(ptype='inj_luminosity', vmin=0, vmax=None,
                                                        punit='erg/s',
                                                        froz=True,
                                                        log=False,
                                                        val=1E39)

        return model_dict_temp_ev
