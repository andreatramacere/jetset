#from __future__ import absolute_import, division, print_function

#from builtins import (bytes, str, open, super, range,
#                      zip, round, input, int, pow, object, map, zip)




__author__ = "Andrea Tramacere"



import os
import numpy as np
from numpy import log10,array,zeros,power,shape
from inspect import signature

from .jetkernel_models_dic import gamma_dic_e ,gamma_dic_p, gamma_dic_pp_e_second, available_N_distr, N_distr_descr, available_emitters_type
from .plot_sedfit import PlotPdistr
from .jet_paramters import *
from .utils import set_str_attr
from .model_parameters import ModelParameterArray, ModelParameter


on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
if on_rtd == True:
    try:
        from jetkernel import jetkernel as BlazarSED
    except ImportError:
        from .mock import jetkernel as BlazarSED
else:
    from jetkernel import jetkernel as BlazarSED


__all__=['JetkernelEmittersDistribution', 'EmittersDistribution', 'ArrayDistribution']




class ArrayDistribution(object):

    def __init__(self, e_array, n_array, gamma_grid_size=None):
        self.e_array = e_array
        self.n_array = n_array
        _size = e_array.size

        if n_array.size != _size:
            raise RuntimeError('e_array and n_array must have same size')
        self.size = _size
        if gamma_grid_size is None:
                self.gamma_grid_size = _size * 2 + 1


class EmittersDistribution(object):

    def __init__(self,
                 name,
                 spectral_type,
                 jet=None,
                 gamma_grid_size=1000,
                 log_values=False,
                 emitters_type='electrons',
                 normalize=False,
                 skip_build=False):

        self._spectral_type = None
        self._allowed_spectral_types=['bkn','plc','lp','lppl','pl']
        self._set_emitters_type(emitters_type)
        self._set_spectral_type(spectral_type)
        self._Norm=1.0
        if skip_build is False:
            self._build(jet,name,log_values, gamma_grid_size,normalize)

    def _copy_from_jet(self, jet):
        self._name = jet.emitters_distribution._name
        self._log_values = jet.emitters_distribution._log_values
        self._gamma_grid_size = jet.emitters_distribution._log_values
        self._gamma_grid_size = jet.emitters_distribution._gamma_grid_size
        self.normalize = jet.emitters_distribution.normalize

        for par in self.parameters.par_array:
            p = jet.emitters_distribution.parameters.get_par_by_name(par.name)
            par.val = p.val

    def _build(self,jet,name,log_values, gamma_grid_size,normalize):
        self._user_defined=True
        self._name = name
        self._log_values = log_values
        self._gamma_grid_size=gamma_grid_size
        self.parameters = ModelParameterArray()
        self.set_parameters_dict()
        self.normalize = normalize

        self.set_jet(jet)

    def _set_spectral_type(self,spectral_type):
        if spectral_type not in self._allowed_spectral_types:
            raise RuntimeError('spectral_type=',spectral_type,' not in allowed',self._allowed_spectral_types)
        else:
            self._spectral_type = spectral_type

    @property
    def spectral_type(self):
        return self._spectral_type

    def _update_parameters_dict(self):
        for par in self.parameters.par_array:
            self._parameters_dict[par.name].val = par.val
            self._parameters_dict[par.name].vmin = par.val_min
            self._parameters_dict[par.name].vmax = par.val_max
            self._parameters_dict[par.name].log = par.islog
            self._parameters_dict[par.name].punit = par.units
            #self._parameters_dict[par.name].is_in_jetkernel = par.is_in_jetkernel
            #self._parameters_dict[par.name].allowed_values = par.allowed_values

            self._parameters_dict[par.name]._is_dependent = par._is_dependent
            self._parameters_dict[par.name]._depending_par = par._depending_par
            self._parameters_dict[par.name]._func = par._func
            self._parameters_dict[par.name]._master_par = par._master_par

    def add_par(self, name, par_type, val, vmax, vmin, unit='', log=False, frozen=False):
        if log is True:
            val = np.log10(val)

        if name not in self._parameters_dict.keys():
            self._parameters_dict[name]=JetModelDictionaryPar(ptype=par_type,
                                                              vmin=vmin,
                                                              vmax=vmax,
                                                              log=log,
                                                              val=val,
                                                              punit=unit,
                                                              froz=frozen,
                                                              is_in_jetkernel=False)
        else:
            raise ValueError('par',name,'already assigned')

        self.parameters.add_par(ModelParameter(name=name,
                                               par_type=par_type,
                                               val=val,
                                               val_min=vmin,
                                               val_max=vmax,
                                               units=unit,
                                               frozen=frozen,
                                               log=log,))

    def set_parameters_dict(self,):
        model_dict, a_h, b_h, a_l, b_l, a_t, b_t = self._get_base_dict_and_bounds()

        model_dict['gmin'].val = 2.0
        model_dict['gmin'].is_in_jetkernel = True

        model_dict['gmax'].val = 1E6
        model_dict['gmax'].is_in_jetkernel = True

        model_dict['N'].val = 100.
        model_dict['N'].is_in_jetkernel = True

        if 'NH_pp' in model_dict.keys():
            model_dict['NH_pp'].val = 1.
            model_dict['NH_pp'].is_in_jetkernel = True


        self._parameters_dict=model_dict

        for k,par in model_dict.items():
            #print('==> k',k,'par',par)
            if self._log_values is True:
                val = np.log10(par.val)
            else:
                val = par.val

            self.parameters.add_par(ModelParameter(name=k,
                                                   par_type=par.ptype,
                                                   val=val,
                                                   val_min=par.vmin,
                                                   val_max=par.vmax,
                                                   units=par.punit,
                                                   frozen=par.froz,
                                                   log=par.log))


    def set_distr_func(self,distr_func):
        if distr_func is not None:
            self._validate_func(distr_func)

        self.distr_func = distr_func

    def _validate_func(self,distr_func):
        s=signature(distr_func)
        for param in s.parameters.values():
            if param.name not in self._parameters_dict.keys() and param.name!='gamma':
                raise RuntimeError('distr_func parmaters',param.name,'not in parameters_dict.keys()',self._parameters_dict.keys())
        if 'gamma' not in s.parameters.keys():
            raise RuntimeError('gamma not present in distr_func signature')


    def _set_emitters_type(self,emitters_type):
        self.available_emitters_type = available_emitters_type
        self._check_emitters_type(emitters_type)
        self.emitters_type = emitters_type

    def _get_base_dict_and_bounds(self):
        model_dic={}
        a_h, b_h = self.set_bounds(1, 1E15, log_val=self._log_values)
        a_l, b_l = self.set_bounds(1, 1E9, log_val=self._log_values)
        a_t, b_t = self.set_bounds(1, 1E9, log_val=self._log_values)
        model_dic['gmin'] = JetModelDictionaryPar(ptype='low-energy-cut-off', vmin=a_l, vmax=b_l,
                                                  punit='lorentz-factor', log=self._log_values)
        model_dic['gmax'] = JetModelDictionaryPar(ptype='high-energy-cut-off', vmin=a_h, vmax=b_h,
                                                  punit='lorentz-factor', log=self._log_values)

        if self.emitters_type =='electrons':
            model_dic['N'] = JetModelDictionaryPar(ptype='emitters_density',vmin=0,vmax=None,punit='cm^-3')

        elif self.emitters_type =='protons':
            model_dic['N'] = JetModelDictionaryPar(ptype='emitters_density', vmin=0, vmax=None, punit='cm^-3')
            model_dic['NH_pp'] = JetModelDictionaryPar(ptype='target_density', vmin=0, vmax=None, punit='cm^-3',froz =False)

        return model_dic, a_h, b_h, a_l, b_l, a_t, b_t

    def set_jet(self,jet):
        if jet is not None:
            name = 'from_array'
            #if emitters_type != 'electrons':
            #    raise RuntimeError('not implemented yet for', emitters_type, 'only for electrons')

            self._jet = jet
            set_str_attr(jet._blob, 'DISTR', name)
            set_str_attr(jet._blob, 'PARTICLE', self.emitters_type)


            p = self._jet.get_par_by_name('gmin')
            if p is not None:
                p.set(value=self.parameters.get_par_by_name('gmin').val_lin)
            else:
                p=self.parameters.get_par_by_name('gmin')

            setattr(jet._blob, 'gmin', p.val_lin)

            p = self._jet.get_par_by_name('gmax')
            if p is not None:
                p.set(value=self.parameters.get_par_by_name('gmax').val_lin)
            else:
                p = self.parameters.get_par_by_name('gmax')
            setattr(jet._blob, 'gmax', p.val_lin)

            if self.emitters_type == 'electrons':
                BlazarSED.build_Ne_custom(jet._blob, self._gamma_grid_size)
            elif self.emitters_type == 'protons':
                BlazarSED.build_Ne_custom(jet._blob, self._gamma_grid_size)
                BlazarSED.build_Np_custom(jet._blob, self._gamma_grid_size)
            else:
                raise RuntimeError('emitters type',self.emitters_type,'not among valid: electrons,protons')
        else:
            self._jet=jet

    def _eval_func(self,gamma):
        p_dict = {}
        _skip=['gmin','gmax','N','NH_pp']
        for par in self.parameters.par_array:
            if par.name not in _skip:
                p_dict[par.name] = par.val
        p_dict['gamma']=gamma

        return self.distr_func(**p_dict)

    @property
    def name(self):
        return self._name

    def _check_emitters_type(self, emitters_type):
        if emitters_type not in available_emitters_type:
            raise RuntimeError("type of emitters :%s, not allowed" % emitters_type, "please choose among: ",
                               self.available_emitters_type)
        else:
            pass

    def update(self):
        self._fill()

    def set_grid_size(self,gamma_grid_size):
        if gamma_grid_size is not None:
            if self._jet is not  None:
              setattr(self._jet._blob,'gamma_grid_size' ,gamma_grid_size)

        self._fill()

    def set_grid(self):

        gmin = self.parameters.get_par_by_name('gmin').val_lin
        gmax = self.parameters.get_par_by_name('gmax').val_lin

        self._gamma_grid = np.logspace(np.log10(gmin), np.log10(gmax), self._gamma_grid_size)

    def eval_N_tot(self):
        if self.emitters_type=='electrons':
            return np.trapz(self.n_gamma_e,self.gamma_e)
        elif self.emitters_type=='protons':
            return np.trapz(self.n_gamma_p, self.gamma_p)
        else:
            raise  RuntimeError('emitters type',self.emitters_type, 'not valid')

    def _fill(self,normalize=False):
        self.set_grid()



        self.f=self._eval_func(gamma=self._gamma_grid)
        self._Norm=1.0

        if self.normalize is True:
            self._Norm=1.0/np.trapz(self.f,self._gamma_grid)

        self.f = self.f*self._Norm*self.parameters.get_par_by_name('N').val

        if self.emitters_type == 'electrons':
            self.gamma_e = zeros(self._gamma_grid_size)
            self.n_gamma_e = zeros(self._gamma_grid_size)
            self.gamma_e = self._gamma_grid
            self.n_gamma_e = self.f

            self.n_gamma_e[np.isnan(self.n_gamma_e)] = 0
            self.n_gamma_e[np.isinf(self.n_gamma_e)] = 0

        if self.emitters_type == 'protons':
            self.gamma_p = zeros(self._gamma_grid_size)
            self.n_gamma_p = zeros(self._gamma_grid_size)
            self.gamma_e_second_inj = zeros(self._gamma_grid_size)
            self.n_gamma_e_second_inj = zeros(self._gamma_grid_size)

            self.gamma_p = self._gamma_grid
            self.n_gamma_p = self.f

            self.n_gamma_p[np.isnan(self.n_gamma_p)] = 0
            self.n_gamma_p[np.isinf(self.n_gamma_p)] = 0

            self.n_gamma_e_second_inj[np.isnan(self.n_gamma_e_second_inj)] = 0
            self.n_gamma_e_second_inj[np.isinf(self.n_gamma_e_second_inj)] = 0

        if self._jet is not None:

            if self.emitters_type == 'electrons':
                # this is not the size used in jetkernel, because it is the python array without the simpson grid
                size = self._gamma_grid_size
                Ne_custom_ptr = getattr(self._jet._blob, 'Ne_custom')
                gamma_custom_ptr = getattr(self._jet._blob, 'gamma_e_custom')

                for ID in range(size):
                    BlazarSED.set_elec_custom_array(gamma_custom_ptr, self._jet._blob, self._gamma_grid[ID], ID)
                    BlazarSED.set_elec_custom_array(Ne_custom_ptr, self._jet._blob, self.f[ID], ID)


            if self.emitters_type == 'protons':

                #this is not the size used in jetkernel, because it is the python array without the simpson grid
                size = self._gamma_grid_size

                Np_custom_ptr = getattr(self._jet._blob, 'Np_custom')
                gamma_custom_ptr = getattr(self._jet._blob, 'gamma_p_custom')
                for ID in range(size):
                    BlazarSED.set_elec_custom_array(gamma_custom_ptr, self._jet._blob, self._gamma_grid[ID], ID)
                    BlazarSED.set_elec_custom_array(Np_custom_ptr, self._jet._blob, self.f[ID], ID)




    def _set_blob(self):
        if self.emitters_type == 'electrons':

            self._Ne_name, self._gammae_name = gamma_dic_e['electron_distr']
            self.Ne_ptr = getattr(self._jet._blob, self._Ne_name)
            self.e_gamma_ptr = getattr(self._jet._blob, self._gammae_name)

        elif self.emitters_type == 'protons':
            self._Ne_name, self._gammae_name = gamma_dic_e['electron_distr']
            self._Np_name, self._gammap_name = gamma_dic_p['proton_distr']
            self._Q_inj_e_second_name, self._gammae_inj_sec_name = gamma_dic_pp_e_second['e_second_inj']

            self.Np_ptr = getattr(self._jet._blob, self._Np_name)
            self.p_gamma_ptr = getattr(self._jet._blob, self._gammap_name)
            self.Ne_ptr = getattr(self._jet._blob, self._Ne_name)
            self.e_gamma_ptr = getattr(self._jet._blob, self._gammae_name)

            self._Q_inj_e_second_ptr = getattr(self._jet._blob, self._Q_inj_e_second_name)
            self.e_inj_second_gamma_ptr = getattr(self._jet._blob, self._gammae_inj_sec_name)

        size = self._jet._blob.gamma_grid_size
        self.gamma_e = zeros(size)
        self.n_gamma_e = zeros(size)
        self.gamma_p = zeros(size)
        self.n_gamma_p = zeros(size)
        self.gamma_e_second_inj = zeros(size)
        self.n_gamma_e_second_inj = zeros(size)
        self.gamma_e_second_inj = zeros(size)
        self.q_gamma_e_second_inj = zeros(size)

        for ID in range(size):
            self.gamma_e[ID] = BlazarSED.get_elec_array(self.e_gamma_ptr, self._jet._blob, ID)
            self.n_gamma_e[ID] = BlazarSED.get_elec_array(self.Ne_ptr, self._jet._blob, ID)

        if self.emitters_type == 'protons':
            for ID in range(size):
                self.gamma_p[ID] = BlazarSED.get_elec_array(self.p_gamma_ptr, self._jet._blob, ID)
                self.n_gamma_p[ID] = BlazarSED.get_elec_array(self.Np_ptr, self._jet._blob, ID)

            for ID in range(size):
                self.gamma_e_second_inj[ID] = BlazarSED.get_elec_array(self.e_inj_second_gamma_ptr, self._jet._blob, ID)
                self.n_gamma_e_second_inj[ID] = BlazarSED.get_elec_array(self._Q_inj_e_second_ptr, self._jet._blob, ID)



    def set_bounds(self,a,b,log_val=False):
        if log_val == False:
            return [a,b]

        else:
            return np.log10([a,b])

    def plot(self, p=None, y_min=None, y_max=None, x_min=None, x_max=None, energy_unit='gamma'):

        self.update()
        if self._jet is not None:
            self._set_blob()

        if p is None:
            p = PlotPdistr()
        if self.emitters_type == 'electrons':
            p.plot_distr(self.gamma_e,
                         self.n_gamma_e,
                         y_min=y_min,
                         y_max=y_max,
                         x_min=x_min,
                         x_max=x_max,
                         particle='electrons',
                         energy_unit=energy_unit)

        if self.emitters_type == 'protons':
            p.plot_distr(self.gamma_p,
                         self.n_gamma_p,
                         y_min=y_min,
                         y_max=y_max,
                         x_min=x_min,
                         x_max=x_max,
                         particle='protons',
                         energy_unit=energy_unit)

            p.plot_distr(self.gamma_e_second_inj,
                         self.n_gamma_e_second_inj,
                         y_min=y_min,
                         y_max=y_max,
                         x_min=x_min,
                         x_max=x_max,
                         particle='electrons',
                         energy_unit=energy_unit)

        return p

    def plot2p(self, p=None, y_min=None, y_max=None, x_min=None, x_max=None, energy_unit='gamma'):
        self.update()
        if self._jet is not None:
            self._set_blob()

        if p is None:
            p = PlotPdistr()

        p.plot_distr2p(self.gamma_e,
                       self.n_gamma_e,
                       y_min=y_min,
                       y_max=y_max,
                       x_min=x_min,
                       x_max=x_max,
                       particle='electrons',
                       energy_unit=energy_unit)

        if self.emitters_type == 'protons':
            p.plot_distr(self.gamma_p,
                         self.n_gamma_p,
                         y_min=y_min,
                         y_max=y_max,
                         x_min=x_min,
                         x_max=x_max,
                         particle='protons',
                         energy_unit=energy_unit)

            p.plot_distr(self.gamma_e_second_inj,
                         self.n_gamma_e_second_inj,
                         y_min=y_min,
                         y_max=y_max,
                         x_min=x_min,
                         x_max=x_max,
                         particle='electrons',
                         energy_unit=energy_unit)

        return p

    def plot3p(self, p=None, y_min=None, y_max=None, x_min=None, x_max=None, energy_unit='gamma'):
        self.update()
        if self._jet is not None:
            self._set_blob()

        if p is None:
            p = PlotPdistr()

        p.plot_distr3p(self.gamma_e,
                       self.n_gamma_e,
                       y_min=y_min,
                       y_max=y_max,
                       x_min=x_min,
                       x_max=x_max,
                       particle='electrons',
                       energy_unit=energy_unit)

        if self.emitters_type == 'protons':
            p.plot_distr(self.gamma_p,
                         self.n_gamma_p,
                         y_min=y_min,
                         y_max=y_max,
                         x_min=x_min,
                         x_max=x_max,
                         particle='protons',
                         energy_unit=energy_unit)

            p.plot_distr(self.gamma_e_second_inj,
                         self.n_gamma_e_second_inj,
                         y_min=y_min,
                         y_max=y_max,
                         x_min=x_min,
                         x_max=x_max,
                         particle='electrons',
                         energy_unit=energy_unit)

        return p


class JetkernelEmittersDistribution(EmittersDistribution):

    def __init__(self,name,jet,gamma_grid_size=None,log_values=False,emitters_type='electrons'):

        super(JetkernelEmittersDistribution, self).__init__(name,
                                                            spectral_type=name,
                                                            skip_build=True,
                                                            emitters_type=emitters_type)
        self._build( jet, name, log_values, gamma_grid_size,emitters_type)

    def _build(self, jet, name, log_values, gamma_grid_size,emitters_type):
        self._set_N_distr_name(name)
        self._user_defined = False
        self._name = name
        self._log_values = log_values
        self._set_jet(jet, name,log_values=log_values, emitters_type=emitters_type)
        self.set_grid_size(gamma_grid_size)

        self._parameters_dict= self._build_emitters_distribution_dict(name, emitters_type=emitters_type)

    def _set_jet(self,jet, name, log_values=False, emitters_type='electrons'):
        self._jet = jet
        #jet.set_emitters_distribution(self, log_values=log_values, emitters_type=emitters_type)
        set_str_attr(jet._blob, 'DISTR', name)
        set_str_attr(jet._blob, 'PARTICLE', emitters_type)

    def _set_N_distr_name(self,name):
        self.available_N_distr = available_N_distr
        if name == 'from_array':
            pass
        elif name not in self.available_N_distr:
            raise RuntimeError ("electron distribution model %s not allowed" % name, "please choose among: ", self.available_N_distr)
        else:
            pass

    def set_grid_size(self,gamma_grid_size):
        if gamma_grid_size is not None:
            setattr(self._jet._blob,'gamma_grid_size' ,gamma_grid_size)

        self._set_blob()
        self._fill()

    def update(self):
        self._set_blob()
        self._fill()

    @staticmethod
    def available_distributions():
        for k in N_distr_descr.keys():
            print('%s: %s'%(k,N_distr_descr[k]))

    @classmethod
    def from_array(cls, jet, custom_Ne, emitters_type='electrons'):
        # TODO update for protons

        if emitters_type != 'electrons':
            raise RuntimeError('not implemented yet for', emitters_type, 'only for electrons')

        name = 'from_array'

        BlazarSED.build_Ne_custom(jet._blob, custom_Ne.size)
        Ne_custom_ptr = getattr(jet._blob, 'Ne_custom')
        gamma_custom_ptr = getattr(jet._blob, 'gamma_e_custom')
        jet._blob.N=1.0
        jet._blob.N_0 = 1.0
        jet.Norm_distr=0
        for ID in range(custom_Ne.size):
            BlazarSED.set_elec_custom_array(gamma_custom_ptr, jet._blob, custom_Ne.e_array[ID], ID)
            BlazarSED.set_elec_custom_array(Ne_custom_ptr, jet._blob, custom_Ne.n_array[ID], ID)

        setattr(jet._blob, 'gmin', custom_Ne.e_array[0])
        setattr(jet._blob, 'gmax', custom_Ne.e_array[-1])

        return cls(name, jet, custom_Ne.gamma_grid_size)

    def _set_blob(self):

        if self.emitters_type=='electrons':
            BlazarSED.InitNe(self._jet._blob)
            self._Ne_name, self._gammae_name = gamma_dic_e['electron_distr']
            self.Ne_ptr = getattr(self._jet._blob, self._Ne_name)
            self.e_gamma_ptr = getattr(self._jet._blob, self._gammae_name)

        elif self.emitters_type=='protons':
            BlazarSED.Init_Np_Ne_pp(self._jet._blob)
            self._Ne_name, self._gammae_name = gamma_dic_e['electron_distr']
            self._Np_name, self._gammap_name = gamma_dic_p['proton_distr']
            self._Q_inj_e_second_name, self._gammae_inj_sec_name = gamma_dic_pp_e_second['e_second_inj']

            self.Np_ptr = getattr(self._jet._blob, self._Np_name)
            self.p_gamma_ptr = getattr(self._jet._blob, self._gammap_name)
            self.Ne_ptr = getattr(self._jet._blob, self._Ne_name)
            self.e_gamma_ptr = getattr(self._jet._blob, self._gammae_name)

            self._Q_inj_e_second_ptr = getattr(self._jet._blob, self._Q_inj_e_second_name)
            self.e_inj_second_gamma_ptr = getattr(self._jet._blob, self._gammae_inj_sec_name)
        else:
            raise RuntimeError('')

    def _fill(self):
        size = self._jet._blob.gamma_grid_size
        self.gamma_e = zeros(size)
        self.n_gamma_e = zeros(size)
        self.gamma_p = zeros(size)
        self.n_gamma_p = zeros(size)
        self.gamma_e_second_inj = zeros(size)
        self.n_gamma_e_second_inj = zeros(size)

        for ID in range(size):
            self.gamma_e[ID] = BlazarSED.get_elec_array(self.e_gamma_ptr, self._jet._blob, ID)
            self.n_gamma_e[ID] = BlazarSED.get_elec_array(self.Ne_ptr, self._jet._blob, ID)

        if self.emitters_type=='protons':
            for ID in range(size):
                self.gamma_p[ID] = BlazarSED.get_elec_array(self.p_gamma_ptr, self._jet._blob, ID)
                self.n_gamma_p[ID] = BlazarSED.get_elec_array(self.Np_ptr, self._jet._blob, ID)

            for ID in range(size):
                self.gamma_e_second_inj[ID] = BlazarSED.get_elec_array(self.e_inj_second_gamma_ptr, self._jet._blob, ID)
                self.n_gamma_e_second_inj[ID] = BlazarSED.get_elec_array(self._Q_inj_e_second_ptr, self._jet._blob, ID)

        self.n_gamma_p[np.isnan(self.n_gamma_p)]=0
        self.n_gamma_p[np.isinf(self.n_gamma_p)]=0

        self.n_gamma_e[np.isnan(self.n_gamma_e)] = 0
        self.n_gamma_e[np.isinf(self.n_gamma_e)] = 0


    def _build_emitters_distribution_dict(self, distribution_name, emitters_type='electrons'):
        """
        Builds the dictionary to init the :class:`.JetParameter`
        objects   for the electron  distribution:

        The following :class:`.JetParameter`: objects the *do not* depend
        on the type of electron distribution

                - N, particle density in cm^-3

                - gmin, the minimum value of the electron Lorentz factor

                - gmax, the maximum value of the electron Lorentz factor

        The following :class:`.JetParameter`: objects *depend* on the type of electron distribution:

            - **power law**, electron_distribution='pl'

               - p

            - **broken power-law**, electron_distribution= **'bkn'**

                - p
                - p_1
                - gamma_break

            - **log-parabola**, electron_distribution= **'lp'**

                - r
                - s
                - gamma0_log_parab (fixed)

            - **log-parabola** with a low energy power-law tail, electron_distribution= **'lppl'**

                - r
                - s
                - gamma0_log_parab

            - **log-parabola** defined by peak energy, electron_distribution= **'lpep'**

                - r
                - s
                - gammap_log_parab,

            - **power-law cut-off**, lectron_distribution= **'plc'**

                - p
                - gamma_cut

        """

        model_dic, a_h, b_h, a_l, b_l, a_t, b_t= self._get_base_dict_and_bounds()

        a_h,b_h=self.set_bounds(1,1E15,log_val=self._log_values)
        a_l, b_l = self.set_bounds(1, 1E9, log_val=self._log_values)
        a_t, b_t = self.set_bounds(1, 1E9, log_val=self._log_values)

        if distribution_name == 'pl':
            model_dic['p'] = JetModelDictionaryPar(ptype='HE_spectral_slope',vmin=-10.,vmax=10.,punit='')

        if distribution_name == 'bkn':
            model_dic['p'] = JetModelDictionaryPar(ptype='LE_spectral_slope',vmin=-10.,vmax=10.,punit='')
            model_dic['p_1'] = JetModelDictionaryPar(ptype='HE_spectral_slope',vmin=-10.,vmax=10.,punit='')
            model_dic['gamma_break'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)

        if distribution_name == 'lp':
            model_dic['s'] = JetModelDictionaryPar(ptype='LE_spectral_slope', vmin=-10., vmax=10., punit='')
            model_dic['r'] = JetModelDictionaryPar(ptype='spectral_curvature', vmin=-15., vmax=15., punit='')
            model_dic['gamma0_log_parab'] = JetModelDictionaryPar(ptype='turn-over-energy', vmin=a_t, vmax=b_t,punit='lorentz-factor', log=self._log_values,froz=True)

        if distribution_name == 'lppl':
            model_dic['s'] = JetModelDictionaryPar(ptype='LE_spectral_slope',vmin=-10.,vmax=10.,punit='')
            model_dic['r'] = JetModelDictionaryPar(ptype='spectral_curvature',vmin=-15.,vmax=15.,punit='')
            model_dic['gamma0_log_parab'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)

        if distribution_name == 'lpep':
            model_dic['r'] =  JetModelDictionaryPar(ptype='spectral_curvature',vmin=-15.,vmax=15.,punit='')
            model_dic['gammap_log_parab'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)

        if distribution_name == 'plc':
            model_dic['p'] = JetModelDictionaryPar(ptype='LE_spectral_slope', vmin=-10., vmax=10., punit='')
            model_dic['gamma_cut'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)

        if distribution_name == 'spitkov':
            model_dic['spit_index'] = JetModelDictionaryPar(ptype='LE_spectral_slope',vmin=-10.,vmax=10.,punit='')
            model_dic['spit_temp'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)
            model_dic['spit_gamma_th'] =JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)

        if distribution_name == 'lppl_pile_up':
            model_dic['s'] =JetModelDictionaryPar(ptype='LE_spectral_slope',vmin=-10.,vmax=10.,punit='')
            model_dic['r'] =JetModelDictionaryPar(ptype='spectral_curvature',vmin=-15.,vmax=15.,punit='')
            model_dic['gamma0_log_parab'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)
            model_dic['gamma_inj'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)
            model_dic['gamma_pile_up'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)
            model_dic['ratio_pile_up'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=0.0,vmax=None,punit='')
            model_dic['alpha_pile_up'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=0.0,vmax=10.0,punit='')

        if distribution_name == 'bkn_pile_up':
            model_dic['p'] = JetModelDictionaryPar(ptype='LE_spectral_slope',vmin=-10.,vmax=10.,punit='')
            model_dic['p_1'] = JetModelDictionaryPar(ptype='LE_spectral_slope',vmin=-10.,vmax=10.,punit='')
            model_dic['gamma_break'] =  JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)
            model_dic['gamma_pile_up'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)
            model_dic['gamma_pile_up_cut'] =JetModelDictionaryPar(ptype='turn-over-energy',vmin=a_t,vmax=b_t,punit='lorentz-factor',log=self._log_values)
            model_dic['alpha_pile_up'] = JetModelDictionaryPar(ptype='turn-over-energy',vmin=0.0,vmax=10.0,punit='')

        if distribution_name == 'from_array':
            model_dic.pop('gmin')
            model_dic.pop('gmax')

        return model_dic



