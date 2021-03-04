
__author__ = "Andrea Tramacere"




# on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

# if on_rtd == True:
#     try:
#         from .jetkernel import jetkernel as BlazarSED
#     except ImportError:
#         from .mock import jetkernel as BlazarSED
# else:

from .jetkernel import jetkernel as BlazarSED


import numpy as np
import time
import  ast


from tqdm.autonotebook import tqdm

import threading

from astropy import constants as const

from .model_parameters import _show_table


from astropy.table import  Table


from .jet_paramters import JetModelDictionaryPar,JetParameter,JetModelParameterArray
from .jet_emitters import  ArrayDistribution, EmittersArrayDistribution

from .plot_sedfit import  BasePlot,PlotPdistr,PlotTempEvDiagram,PlotTempEvEmitters


__all__=['JetTimeEvol']


class ProgressBar(object):

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
            time.sleep(.1)
            n_try += 1
            self.update()


    def run(self):
        self.pbar.n = 0
        while (self.stop is False):
            self.update()
            time.sleep(.1)



class JetTimeEvol(object):
    """


        Parameters
        ----------
        jet
        Q_inj
        name
        inplace
        jet_gamma_grid_size






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

        T_SIZE: int
            number of steps for the time grid

        NUM_SET: int
            number of the time steps saved in the output

        LOG_SET: bool
           if True, the time grid is logarithmically spaced

        L_inj: float (erg/s)
            if not None, the inj function is rescaled to match L_inj
    """

    def __init__(self,
                 jet,
                 Q_inj=None,
                 #flag='tests',
                 name='jet_time_ev',
                 inplace=True,
                 jet_gamma_grid_size=200):


        self._temp_ev = BlazarSED.MakeTempEv()
        if inplace is True:
            self.jet=jet.clone()
        else:
            self.jet=jet

        self.Q_inj=Q_inj
        self.name=name

        self.parameters = JetModelParameterArray(model=self)
        self.parameters.add_par_from_dict(self._build_par_dict(),self,'_temp_ev',JetParameter)
        self.jet_gamma_grid_size = jet_gamma_grid_size





    @property
    def temp_ev(self):
        return self._temp_ev

    def show_pars(self, sort_key='par type'):
        self.parameters.show_pars(sort_key=sort_key)

    def init_TempEv(self,custom_inj_profile=None,custom_acc_profile=None):
        BlazarSED.Init(self.jet._blob, self.jet.get_DL_cm())
        self._init_temp_ev()
        self._fill_temp_ev_array_pre_run()
        #if custom_ini_profile is not None:
        self._build_tempev_table(self.jet)


    def run(self,only_injection=True):
        print('temporal evolution running')
        self.init_TempEv()
        pbar=ProgressBar(target_class=self,N=self.parameters.T_SIZE.val)
        t1 = threading.Thread(target=pbar.run)
        t1.start()
        BlazarSED.Run_temp_evolution(self.jet._blob, self._temp_ev, int(only_injection))
        pbar.stop = True
        pbar.finalzie()
        self._fill_temp_ev_array_post_run()
        print('temporal evolution completed')

        self._set_jet_post_run()

    def _set_jet_post_run(self):
        self._jet_emitters_distr = EmittersArrayDistribution(name='time_dep',
                                                             gamma_array=np.copy(self.gamma),
                                                             n_gamma_array=np.copy(self.N_gamma[-1]),
                                                             gamma_grid_size=self.jet_gamma_grid_size,
                                                             normalize=False)
        self._jet_emitters_distr._fill()
        self._jet_emitters_distr._fill()
        self.jet.set_emitters_distribution(self._jet_emitters_distr)


    def _init_temp_ev(self):
        if self.Q_inj is not None:
            setattr(self._temp_ev,'Q_inj_jetset_gamma_grid_size',self.Q_inj._gamma_grid_size)
        BlazarSED.Init_Q_inj(self._temp_ev)
        BlazarSED.Init_temp_evolution(self.jet._blob, self._temp_ev, self.jet.get_DL_cm())
        if self.Q_inj is not None:
            self.Q_inj._set_L_inj(self.parameters.L_inj.val,BlazarSED.V_sphere(self.jet.parameters.R.val))
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
        self.t_Esc_pre_run = np.zeros(size)


        for i in range(size):
            self.gamma_pre_run[i] = BlazarSED.get_temp_ev_array_static(self._temp_ev.g, i)
            self.t_Sync_cool_pre_run[i] = BlazarSED.get_temp_ev_array_static(self._temp_ev.t_Sync_cool, i)
            self.t_D_pre_run[i] = BlazarSED.get_temp_ev_array_static(self._temp_ev.t_D, i)
            self.t_DA_pre_run[i] = BlazarSED.get_temp_ev_array_static(self._temp_ev.t_DA, i)
            self.t_A_pre_run[i] = BlazarSED.get_temp_ev_array_static(self._temp_ev.t_A, i)
            self.t_Esc_pre_run[i] = BlazarSED.get_temp_ev_array_static(self._temp_ev.t_Esc, i)



    def eval_L_tot_inj(self):
        if self.Q_inj is not None and self.jet is not None:
            return self.Q_inj.eval_U_q() * BlazarSED.V_sphere(self.jet.parameters.R.val)

    def set_time(self,T_step=None):
        self.jet.emitters_distribution._set_arrays(self.gamma, self.N_gamma[T_step].flatten())
        self.jet.emitters_distribution._fill()

    def inj_time_profile(self,user_defined_arry=None):
        if user_defined_arry is None:
            self._inj_t_prof=np.zeros(self._temp_ev.T_SIZE, dtype=np.float)
            _t=np.linspace(0.,self._temp_ev.duration,self._temp_ev.T_SIZE)
            msk=_t>=self._temp_ev.TStart_Inj_Inj
            msk*=_t<=self._temp_ev.TStop_Inj
            self._inj_t_prof[msk]=1.0
        else:
            if np.size(user_defined_arry)!=self._temp_ev.T_SIZE:
                raise  RuntimeError('inj_time_profile must match T_SIZE',self._temp_ev.T_SIZE)
            self._inj_t_prof=user_defined_arry

    def _fill_temp_ev_array_post_run(self):
        gamma_size = self._temp_ev.gamma_grid_size
        time_size = self._temp_ev.NUM_SET
        self.N_time = np.zeros(time_size)
        self.gamma = np.zeros(gamma_size)
        self.N_gamma = np.zeros((time_size,(gamma_size)))
        for i in range(time_size):
            self.N_time[i]=BlazarSED.get_temp_ev_N_time_array(self._temp_ev.N_time, self._temp_ev, i)
            for j in range(gamma_size):
                self.gamma[j] = BlazarSED.get_temp_ev_gamma_array(self._temp_ev.gamma, self._temp_ev, j)
                self.N_gamma[i,j] = BlazarSED.get_temp_ev_N_gamma_array(self._temp_ev.N_gamma, self._temp_ev,i,j)

    def plot_time_profile(self,figsize=(8,8),dpi=120):
        p=PlotTempEvDiagram(figsize=figsize,dpi=dpi)
        p.plot(self.parameters.duration.val,
               self.parameters.TStart_Inj.val,
               self.parameters.TStop_Inj.val,
               self.parameters.TStart_Acc.val,
               self.parameters.TStop_Acc.val,)
        return p

    def plot_TempEv_emitters(self,figsize=(8,8),dpi=120,energy_unit='gamma',loglog=True,plot_Q_inj=True,pow=None):
        p=PlotTempEvEmitters(figsize=figsize,dpi=dpi,loglog=loglog)
        p.plot_distr(self,energy_unit=energy_unit,plot_Q_inj=plot_Q_inj,pow=pow)

        return p


    def plot_pre_run_plot(self,figsize=(8,6),dpi=120):
        p=BasePlot(figsize=figsize,dpi=dpi)
        p.ax.loglog(self.gamma_pre_run, self.t_Sync_cool_pre_run, label='t coool, synch.')
        p.ax.loglog(self.gamma_pre_run, self.t_D_pre_run, label='t acc. diffusive')
        p.ax.loglog(self.gamma_pre_run, self.t_DA_pre_run, label='t acc. diffusive syst.')
        p.ax.loglog(self.gamma_pre_run, self.t_A_pre_run, label='t acc. systematic')
        p.ax.loglog(self.gamma_pre_run, self.t_Esc_pre_run, label='t esc')
        p.ax.axvline(self._temp_ev.gamma_eq_t_A, ls='--')
        p.ax.axvline(self._temp_ev.gamma_eq_t_D, ls='--')

        p.ax.legend()
        return p


    def make_lc(self, t1, t2, delta_t_out, nu1, nu2, comp='Sum', delta_T=1, n_slices=1000,rest_frame='obs',eval_cross_time=True):
        _T_array = np.arange(t1, t2, delta_T, dtype=np.int)
        _flux_array=np.zeros(_T_array.shape)
        for ID, T in enumerate(_T_array):
            self.set_time(T)
            self.jet.eval()
            s = getattr(self.jet.spectral_components, comp)
            if rest_frame=='src':
                msk = s.SED.nu_src.value >= nu1
                msk *= s.SED.nu_src.value <= nu2
                x = s.SED.nu_src[msk]
                y = s.SED.nuLnu_src[msk] / x
            elif rest_frame == 'obs':
                msk = s.SED.nu.value >= nu1
                msk *= s.SED.nu.value <= nu2
                x = s.SED.nu[msk]
                y = s.SED.nuFnu[msk] / x
            else:
                raise  RuntimeError('rest frame must be src or obs')
            _flux_array[ID] = np.trapz(y.value, x.value)


        _time_array = self.N_time[_T_array]
        self.time_array = np.arange(_time_array[0], _time_array[-1], delta_t_out)
        self.flux_array=np.interp(self.time_array, _time_array, _flux_array)

        if eval_cross_time is True:
            if rest_frame=='src':
                beaming = 1
                z=0
            else:
                beaming = None
                z=None
            self.time_array, self.flux_array=self.lc_t_obs(self.time_array,self.flux_array,n_slices=n_slices, beaming=beaming, z=z)

        return np.copy(self.time_array), np.copy(self.flux_array)

    def _lc_weight(self, delay_r, R, delta_R):
        return 3*(R-delay_r)*(R+delay_r)*delta_R/(4*R*R*R)


    def lc_t_obs(self,t,lc,n_slices=1000,R=None,beaming=None, z=None):
        if R is None:
            R=self.jet.parameters.R.val

        if beaming is  None:
            beaming=self.jet.get_beaming()

        if z is None:
            z=self.jet.parameters.z_cosm.val

        c = const.c.cgs.value

        delta_R = (2 * R) / n_slices
        delay_t = np.linspace(0, 2 * R / c, n_slices)
        delay_r = R - delay_t * c
        weight_array = self._lc_weight(delay_r, R, delta_R)
        lc_out = np.zeros(lc.shape)
        for ID, lc_in in enumerate(lc):
            t_interp = np.linspace(t[ID] - 2 * R / c, t[ID], n_slices)
            m = t_interp > t[0]
            lc_interp = np.interp(t_interp[m], t, lc)
            lc_out[ID] = np.sum(lc_interp * weight_array[m])

        return t*(1+z)/beaming, lc_out

    def plot_model(self, T1, T2, step,p=None):
        t_array = np.linspace(T1, T2, step, dtype=np.int)
        for ID, T in enumerate(t_array):
            #ed = EmittersArrayDistribution(name='time_dep',
            #                               gamma_array=np.copy(self.temp_ev.gamma),
            #                               n_gamma_array=np.copy(self.temp_ev.N_gamma[T]),
            #                               gamma_grid_size=self.gamma_grid_size)
            #ed._fill()
            #self.jet.set_emitters_distribution(ed)
            self.set_time(T)
            self.jet.eval()
            label = None
            ls = '-'
            if self.N_time[T] > self.parameters.TStop_Acc.val:
                color = 'g'
                ls = '--'
            else:
                color = 'b'
                ls = '-'

            if ID == 0:
                ls = '.'
                label = 'start'
                color = 'black'
            if ID == t_array.size - 1:
                ls = '*-'
                label = 'stop'

            p = self.jet.plot_model(plot_obj=p, comp='Sum', auto_label=False, label=label, line_style=ls, color=color)

    def show_model(self, jet=None, getstring=False, names_list=None, sort_key=None):
        if jet is None:
            jet = self.jet
        print('-'*80)
        print("JetTimeEvol model description")
        print('-'*80)
        self._build_tempev_table(jet)
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

    def _build_tempev_table(self, jet):

        _names = ['name', 'par type', 'val', 'units', 'val*', 'units*', 'log', ]

        rows = []
        rows.append(self._build_row_dict('delta t', 'time', 's', val=self._temp_ev.deltat, val_by=self._temp_ev.t_unit,
                                         unit1='R/c', islog=False))
        rows.append(
            self._build_row_dict('R/c', 'time', 's', val=self._temp_ev.t_unit, val_by=self._temp_ev.t_unit, unit1='R/c',
                                 islog=False))
        rows.append(self._build_row_dict('Diff coeff', '', 's-1', val=self._temp_ev.Diff_Coeff, islog=False))
        rows.append(self._build_row_dict('Acc coeff', '', 's-1', val=self._temp_ev.Acc_Coeff, islog=False))

        rows.append(self._build_row_dict('Diff index', '', '', val=self._temp_ev.Diff_Index, islog=False))
        rows.append(self._build_row_dict('Acc index', '', 's-1', val=self._temp_ev.Acc_Index, islog=False))

        rows.append(
            self._build_row_dict('Tesc', 'time', 's', val=self._temp_ev.T_esc_Coeff, val_by=self._temp_ev.t_unit,
                                 unit1='R/c', islog=False))
        rows.append(
            self._build_row_dict('T_A0=1/ACC_COEFF', 'time', 's', val=self._temp_ev.t_A0, val_by=self._temp_ev.t_unit,
                                 unit1='R/c', islog=False))
        rows.append(
            self._build_row_dict('T_D0=1/DIFF_COEFF', 'time', 's', val=self._temp_ev.t_D0, val_by=self._temp_ev.t_unit,
                                 unit1='R/c', islog=False))
        rows.append(self._build_row_dict('T_DA0=1/(2*DIFF_COEFF)', 'time', 's', val=self._temp_ev.t_DA0,
                                         val_by=self._temp_ev.t_unit, unit1='R/c', islog=False))

        rows.append(self._build_row_dict('gamma Lambda Turb.  max', '', '', val=self._temp_ev.Gamma_Max_Turb_L_max,
                                         islog=False))
        rows.append(self._build_row_dict('gamma Lambda Coher. max', '', '', val=self._temp_ev.Gamma_Max_Turb_L_coher,
                                         islog=False))

        rows.append(self._build_row_dict('gamma eq Syst. Acc (synch. cool)', '', '', val=self._temp_ev.gamma_eq_t_A,
                                         islog=False))
        rows.append(self._build_row_dict('gamma eq Diff. Acc (synch. cool)', '', '', val=self._temp_ev.gamma_eq_t_D,
                                         islog=False))

        t = BlazarSED.Sync_tcool(jet._blob, self._temp_ev.gamma_eq_t_D)
        rows.append(self._build_row_dict('T cooling(gamma_eq=gamma_eq_Diff)', '', 's', val=t, islog=False))
        t = BlazarSED.Sync_tcool(jet._blob, self._temp_ev.gamma_eq_t_A)
        rows.append(self._build_row_dict('T cooling(gamma_eq=gamma_eq_Sys)', '', 's', val=t, islog=False))
        t = BlazarSED.Sync_tcool(jet._blob, self._temp_ev.gamma_eq_t_A)
        rows.append(
            self._build_row_dict('T min. synch. cooling', '', 's', val=self.t_Sync_cool_pre_run.min(), islog=False))

        if self.Q_inj is not None:
            rows.append(self._build_row_dict('L inj (electrons)', 'injected lum.', 'erg/s', val=self.eval_L_tot_inj(),
                                             islog=False))
            #rows.append(self._build_row_dict('E_inj (electrons)', '', 'erg', val=self.eval_E_tot_inj(), islog=False))

        self._tempev_table = Table(rows=rows, names=_names, masked=False)

        self._fromat_column_entry(self._tempev_table)

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

        T_SIZE: int
            number of steps for the time grid

        NUM_SET: int
            number of the time steps saved in the output

        LOG_SET: bool
           if True, the time grid is logarithmically spaced

        L_inj: float (erg/s)
            if not None, the inj function is rescaled to match L_inj
        """

        model_dic = {}

        model_dic['duration'] = JetModelDictionaryPar(ptype='time_grid', vmin=0, vmax=None, punit='s', froz=True, log=False,val=1E5)

        model_dic['gmin_grid'] = JetModelDictionaryPar(ptype='gamma_grid', vmin=0, vmax=None, punit='', froz=True, log=False, val=1E1,jetkernel_par_name='gmin_griglia')

        model_dic['gmax_grid'] = JetModelDictionaryPar(ptype='gamma_grid', vmin=0, vmax=None, punit='', froz=True, log=False, val=1E8,jetkernel_par_name='gmax_griglia')

        model_dic['gamma_grid_size'] = JetModelDictionaryPar(ptype='gamma_grid', vmin=0, vmax=None, punit='', froz=True, log=False, val=5000,jetkernel_par_name='gamma_grid_size')

        model_dic['TStart_Acc'] = JetModelDictionaryPar(ptype='time_grid', vmin=0, vmax=None, punit='s', froz=True,
                                                         log=False,val=0)

        model_dic['TStop_Acc'] = JetModelDictionaryPar(ptype='time_grid', vmin=0, vmax=None, punit='s', froz=True,
                                                         log=False,val=1E5)

        model_dic['TStart_Inj'] = JetModelDictionaryPar(ptype='time_grid', vmin=0, vmax=None, punit='s', froz=True,
                                                         log=False,val=0)

        model_dic['TStop_Inj'] = JetModelDictionaryPar(ptype='time_grid', vmin=0, vmax=None, punit='s', froz=True,
                                                         log=False,val=1E5)

        model_dic['T_esc'] = JetModelDictionaryPar(ptype='escape_time', vmin=None, vmax=None, punit='(R/c)', froz=True,
                                                         log=False,val=2.0,jetkernel_par_name='T_esc_Coeff_R_by_c')

        model_dic['Esc_Index'] = JetModelDictionaryPar(ptype='fp_coeff_index', vmin=None, vmax=None, punit='', froz=True,
                                                            log=False,val=0.0)

        model_dic['t_D0'] = JetModelDictionaryPar(ptype='acceleration_time', vmin=0, vmax=None, punit='s', froz=True,
                                                            log=False,val=1E4)

        model_dic['t_A0'] = JetModelDictionaryPar(ptype='acceleration_time', vmin=0, vmax=None, punit='s', froz=True,
                                                            log=False,val=1E3)

        model_dic['Diff_Index'] = JetModelDictionaryPar(ptype='fp_coeff_index', vmin=0, vmax=None, punit='s', froz=True,
                                                     log=False,val=2)

        model_dic['Acc_Index'] = JetModelDictionaryPar(ptype='fp_coeff_index', vmin=None, vmax=None, punit='', froz=True,
                                                     log=False,val=1)

        model_dic['Lambda_max_Turb'] = JetModelDictionaryPar(ptype='turbulence_scale', vmin=0, vmax=None, punit='cm', froz=True,
                                                          log=False,val=1E15)

        model_dic['Lambda_choer_Turb_factor'] = JetModelDictionaryPar(ptype='turbulence_scale', vmin=0, vmax=None, punit='cm',
                                                                froz=True,log=False, val=0.1)

        model_dic['T_SIZE'] = JetModelDictionaryPar(ptype='time_grid', vmin=0, vmax=None,
                                                                         punit='',
                                                                         froz=True,
                                                                         log=False,val=1000)

        model_dic['NUM_SET'] = JetModelDictionaryPar(ptype='time_ev_output', vmin=0, vmax=None, punit='', froz=True, log=False, val=50)

        model_dic['LOG_SET'] = JetModelDictionaryPar(ptype='time_ev_output', vmin=0, vmax=None, punit='', froz=True, log=False, val=50,allowed_values=[0,1])

        model_dic['L_inj'] = JetModelDictionaryPar(ptype='inj_luminosity', vmin=0, vmax=None,
                                                        punit='erg/s',
                                                        froz=True,
                                                        log=False,
                                                        val=1E39)

        return model_dic
