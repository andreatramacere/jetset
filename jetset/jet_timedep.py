#from __future__ import absolute_import, division, print_function

#from builtins import (str, open, super, range,
#                      object, map)




__author__ = "Andrea Tramacere"



import os

on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if on_rtd == True:
    try:
        from .jetkernel import jetkernel as BlazarSED
    except ImportError:
        from .mock import jetkernel as BlazarSED
else:
    from .jetkernel import jetkernel as BlazarSED


import numpy as np

from .model_parameters import ModelParameterArray, ModelParameter

from .output import makedir,WorkPlace

from astropy.table import  Table

from .utils import safe_run,set_str_attr, old_model_warning

from .jet_paramters import JetModelDictionaryPar,JetParameter,JetModelParameterArray
from .jet_emitters import  ArrayDistribution, EmittersArrayDistribution

from .plot_sedfit import  BasePlot,PlotPdistr,PlotTempEvDiagram,PlotTempEvEmitters

import  ast

import copy


__all__=['JetTimeEvol']






class JetTimeEvol(object):

    def __init__(self,
                 jet,
                 #flag='tests',
                 name='jet_time_ev'):
                 #inplace=True):

        self._temp_ev = BlazarSED.MakeTempEv()
        #if inplace is True:
        #self.jet=jet.clone(jet)
        #else:
        self.jet=jet


        self.name=name

        self.parameters = JetModelParameterArray(model=self)
        self.parameters.add_par_from_dict(self._build_par_dict(),self,'_temp_ev',JetParameter)
        self.Q_inj=None
        #self.init_TempEv()

    def _build_par_dict(self):

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


    def show_pars(self, sort_key='par type'):
        self.parameters.show_pars(sort_key=sort_key)

    def init_TempEv(self, skip_jet=False):
        #if jet is None:
        #    jet=self.jet
        self.jet.set_blob()
        print("--> self.jet._blob.E_tot_e ",self.jet._blob.E_tot_e,self._temp_ev.Q_scaling_factor)
        BlazarSED.Init_temp_evolution(self.jet._blob, self._temp_ev, self.jet.get_DL_cm())
        print("--> self.jet._blob.E_tot_e ", self.jet._blob.E_tot_e,self._temp_ev.Q_scaling_factor)
        self._fill_temp_ev_array_pre_run()
        print("--> self.jet._blob.E_tot_e ", self.jet._blob.E_tot_e,self._temp_ev.Q_scaling_factor)
        self._build_tempev_table(self.jet)
        print("--> self.jet._blob.E_tot_e ", self.jet._blob.E_tot_e,self._temp_ev.Q_scaling_factor,self._temp_ev.deltat)
        print("--> self._Q_inj.max() ",self._Q_inj.max())
        #self.Q_inj = ArrayDistribution(np.copy(self.jet.emitters_distribution.gamma_e),
        #                               np.copy(self.jet.emitters_distribution.n_gamma_e))
        self.Q_inj = ArrayDistribution(self.gamma, self._Q_inj)
        if skip_jet is False:
            ed = EmittersArrayDistribution(name='time_dep',
                                           gamma_array=np.copy(self.Q_inj.e_array),
                                           n_gamma_array=np.copy(self.Q_inj.n_array))
            ed._fill()
            self.jet.set_emitters_distribution(ed)
            self.jet.eval()
        print("--> self._Q_inj.max() ", self._Q_inj.max())
    def run(self):
        #self.init_TempEv(skip_jet=True)
        #self.set_path(path,clean_work_dir=clean_work_dir)
        #self.set_flag(flag)
        BlazarSED.Run_temp_evolution(self.jet._blob, self._temp_ev)
        self._fill_temp_ev_array_post_run()

    def _fill_temp_ev_array_pre_run(self):
        size = BlazarSED.static_ev_arr_grid_size
        self.gamma_pre_run=np.zeros(size)
        self.t_Sync_cool_pre_run = np.zeros(size)
        self.t_D_pre_run = np.zeros(size)
        self.t_DA_pre_run = np.zeros(size)
        self.t_A_pre_run = np.zeros(size)
        self.t_Esc_pre_run = np.zeros(size)
        gamma_size = self._temp_ev.gamma_grid_size
        self.gamma = np.zeros(gamma_size)
        self._Q_inj = np.zeros(gamma_size)

        for i in range(size):
            self.gamma_pre_run[i] = BlazarSED.get_temp_ev_array_static(self._temp_ev.g, i)
            self.t_Sync_cool_pre_run[i] = BlazarSED.get_temp_ev_array_static(self._temp_ev.t_Sync_cool, i)
            self.t_D_pre_run[i] = BlazarSED.get_temp_ev_array_static(self._temp_ev.t_D, i)
            self.t_DA_pre_run[i] = BlazarSED.get_temp_ev_array_static(self._temp_ev.t_DA, i)
            self.t_A_pre_run[i] = BlazarSED.get_temp_ev_array_static(self._temp_ev.t_A, i)
            self.t_Esc_pre_run[i] = BlazarSED.get_temp_ev_array_static(self._temp_ev.t_Esc, i)

        for j in range(gamma_size):
            self.gamma[j] = BlazarSED.get_Q_inj_array(self._temp_ev.gamma, self._temp_ev, j)
            self._Q_inj[j] = BlazarSED.get_Q_inj_array(self._temp_ev.Q_inj, self._temp_ev, j)

    def show_model(self,jet=None,getstring=False,names_list=None,sort_key=None):
        if jet is None:
            jet=self.jet
        print("-------------------------------------------------------------------------------------------------------------------")
        print("JetTimeEvol model description")
        print("-------------------------------------------------------------------------------------------------------------------")
        self._build_tempev_table(jet)
        print(" ")
        print("physical  setup: ")
        print("")
        self._tempev_table.pprint_all()


        print("")
        print("model parameters: ")
        print("")
        self.show_pars()

    def _build_row_dict(self, name, type='', unit='', val='', val_by=None, islog=False,unit1=''):
        row_dict = {}
        row_dict['name'] = name
        row_dict['par type'] = type
        row_dict['units'] = unit
        row_dict['val'] =  val
        if val_by is not None:
            val_by= val/val_by
        row_dict['val*'] = val_by
        row_dict['units*'] = unit1
        row_dict['log'] = islog
        return  row_dict

    def _build_tempev_table(self, jet):

        _names = ['name', 'par type', 'val','units', 'val*', 'units*' ,'log',]

        rows=[]
        rows.append(self._build_row_dict('delta t', 'time', 's', val=self._temp_ev.deltat, val_by=self._temp_ev.t_unit,unit1='R/c', islog=False))
        rows.append(self._build_row_dict('R/c', 'time', 's', val=self._temp_ev.t_unit, val_by=self._temp_ev.t_unit,unit1='R/c', islog=False))
        rows.append(self._build_row_dict('Diff coeff', '', 's-1', val=self._temp_ev.Diff_Coeff, islog=False))
        rows.append(self._build_row_dict('Acc coeff', '', 's-1', val=self._temp_ev.Acc_Coeff, islog=False))

        rows.append(self._build_row_dict('Diff index', '', '', val=self._temp_ev.Diff_Index, islog=False))
        rows.append(self._build_row_dict('Acc index', '', 's-1', val=self._temp_ev.Acc_Index, islog=False))

        rows.append(self._build_row_dict('Tesc', 'time', 's', val=self._temp_ev.T_esc_Coeff, val_by=self._temp_ev.t_unit,unit1='R/c', islog=False))
        rows.append(self._build_row_dict('T_A0=1/ACC_COEFF', 'time', 's', val=self._temp_ev.t_A0, val_by=self._temp_ev.t_unit,unit1='R/c', islog=False))
        rows.append(self._build_row_dict('T_D0=1/DIFF_COEFF', 'time', 's', val=self._temp_ev.t_D0, val_by=self._temp_ev.t_unit,unit1='R/c', islog=False))
        rows.append(self._build_row_dict('T_DA0=1/(2*DIFF_COEFF)', 'time', 's', val=self._temp_ev.t_DA0, val_by=self._temp_ev.t_unit,unit1='R/c', islog=False))

        rows.append(self._build_row_dict('gamma Lambda Turb.  max', '', '', val=self._temp_ev.Gamma_Max_Turb_L_max, islog=False))
        rows.append(self._build_row_dict('gamma Lambda Coher. max', '', '', val=self._temp_ev.Gamma_Max_Turb_L_coher, islog=False))

        rows.append(self._build_row_dict('gamma eq Syst. Acc (synch. cool)', '', '', val=self._temp_ev.gamma_eq_t_A, islog=False))
        rows.append(self._build_row_dict('gamma eq Diff. Acc (synch. cool)', '', '', val=self._temp_ev.gamma_eq_t_D, islog=False))

        t = BlazarSED.Sync_tcool(jet._blob, self._temp_ev.gamma_eq_t_D)
        rows.append(self._build_row_dict('T cooling(gamma_eq=gamma_eq_Diff)', '', 's', val=t, islog=False))
        t = BlazarSED.Sync_tcool(jet._blob, self._temp_ev.gamma_eq_t_A)
        rows.append(self._build_row_dict('T cooling(gamma_eq=gamma_eq_Sys)', '', 's', val=t, islog=False))
        t = BlazarSED.Sync_tcool(jet._blob, self._temp_ev.gamma_eq_t_A)
        rows.append(self._build_row_dict('T min. synch. cooling', '', 's', val=self.t_Sync_cool_pre_run.min(), islog=False))


        rows.append(self._build_row_dict('L inj', 'injected lum.', 'erg/s', val=self._temp_ev.L_inj, islog=False))
        rows.append(self._build_row_dict('E_tot (electrons)/(delta t)', '', 'erg/s', val=jet._blob.E_tot_e*self._temp_ev.Q_scaling_factor / self._temp_ev.deltat, islog=False))


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





    def set_time(self,t=None,T_step=None):
        self.jet.emitters_distribution._array_gamma = self.gamma
        self.jet.emitters_distribution._array_n_gamma = self.N_gamma[T_step].flatten()




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
        self.N_time  = np.zeros(time_size)
        self.N_gamma = np.zeros((time_size,(gamma_size)))
        for i in range(time_size):
            self.N_time[i]=BlazarSED.get_temp_ev_N_time_array(self._temp_ev.N_time, self._temp_ev, i)
            for j in range(gamma_size):
                self.gamma[j] = BlazarSED.get_temp_ev_gamma_array(self._temp_ev.gamma, self._temp_ev, j)
                self.N_gamma[i,j] = BlazarSED.get_temp_ev_N_gamma_array(self._temp_ev.N_gamma, self._temp_ev,i,j)


    def plot_Inj_profile(self,figsize=(8,8),dpi=120):
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

    def plot_Q_inj(self,p=None,figsize=(8,6),dpi=120, energy_unit='gamma',loglog=True):
        if p is None:
            p = PlotPdistr(figsize=figsize,dpi=dpi,injection=True,loglog=loglog)

        p.plot_distr(self.Q_inj.e_array,
                     self.Q_inj.n_array,
                     particle='electrons',
                     energy_unit=energy_unit)


        return p

    # def set_path(self, path, clean_work_dir=True):
    #     if path is None:
    #         path=self.path
    #
    #     if path.endswith('/'):
    #         pass
    #     else:
    #         path += '/'
    #
    #     set_str_attr(self._temp_ev, 'path', path)
    #     makedir(path, clean_work_dir=clean_work_dir)



    # def set_flag(self,flag=None):
    #     if flag is None:
    #         flag=self.flag
    #
    #     self._temp_ev.STEM=flag


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


