__author__ = "Andrea Tramacere"

import os
import numpy as np
from inspect import signature
from numba import njit
from numba.extending import is_jitted
import copy
import warnings
from astropy.constants import m_e,m_p,c
from scipy import interpolate
from .jetkernel_models_dic import gamma_dic_e, gamma_dic_p, gamma_dic_pp_e_second, gamma_dic_e_equilibrium, available_emitters_type
from .plot_sedfit import PlotPdistr
from .jet_paramters import *
from .utils import set_str_attr, get_nested_attr, set_nested_attr
from .model_parameters import ModelParameterArray, ModelParameter
from .jet_kernel_tools import get_emitters_c_array1d_fast as get_emitters
from .jet_kernel_tools import set_emitters_c_array1d_fast as set_emitters

on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
if on_rtd is True:
    try:
        from .jetkernel import jetkernel as BlazarSED
    except ImportError:
        from .mock import jetkernel as BlazarSED
else:
    from .jetkernel import jetkernel as BlazarSED

__all__ = ['EmittersDistribution',
           'BaseEmittersDistribution',
           'ArrayDistribution',
           'EmittersArrayDistribution',
           'InjEmittersDistribution',
           'InjEmittersArrayDistribution']


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


def check_par_name(method):
    def inner(ref,  *args, **kwargs):
        #print(ref,args,kwargs)
        if args[0]  in ref._skip:
            raise RuntimeError('par name',args[0],'can not be in proteced names',ref._skip)
        else:
            return method(ref,*args, **kwargs)

    return inner


class BaseEmittersDistribution(object):

    def __init__(self,
                 name,
                 spectral_type,
                 gamma_grid_size=200,
                 log_values=False,
                 emitters_type='electrons',
                 skip_build=False,
                 normalize=False):

        self._spectral_type = None
        self._allowed_spectral_types = ['bkn', 'plc', 'lp', 'lppl', 'pl', 'lpep', 'array','user_defined']
        self._set_emitters_type(emitters_type)
        self._set_spectral_type(spectral_type)
        self._Norm = 1.0
        self._skip = ['gmin', 'gmax', 'N', 'NH_pp', 'Q']
        if skip_build is False:
            self._build( name, log_values, gamma_grid_size, normalize)

    def _build(self, name, log_values, gamma_grid_size, normalize):
        self._user_defined = True
        self._name = name
        self._log_values = log_values
        self._gamma_grid_size = gamma_grid_size
        self.parameters = ModelParameterArray()
        self._set_base_pars()
        self.normalize = normalize
        self._reset_equilibrium_state()

    def _reset_equilibrium_state(self):
        # Shared state used by plotting/readback when cooling-equilibrium overlays are available.
        self.gamma_cooling_eq = None
        self.gamma_cooling_eq_second = None
        self._primaries_done = False
        self._secondaries_done = False


    def _set_spectral_type(self,spectral_type):
        if spectral_type not in self._allowed_spectral_types:
            raise RuntimeError('spectral_type=',spectral_type,' not in allowed',self._allowed_spectral_types)
        else:
            self._spectral_type = spectral_type


    @staticmethod
    def spectral_types_obs_constrain():
        return ['bkn', 'plc', 'lp', 'lppl', 'pl', 'lpep', 'array']

    @property
    def spectral_type(self):
        return self._spectral_type

    @check_par_name
    def add_par(self, name, par_type, val, vmax, vmin, unit='', log=False, frozen=False):
        

        self.parameters.add_par(ModelParameter(name=name,
                                               par_type=par_type,
                                               val=val,
                                               val_min=vmin,
                                               val_max=vmax,
                                               units=unit,
                                               frozen=frozen,
                                               log=log,))

    def set_distr_func(self,distr_func):
        
        if distr_func is not None:
            self._validate_func(distr_func)
        self.distr_func=distr_func

    def _validate_func(self,distr_func):
        s=signature(distr_func)
        for param in s.parameters.values():
            if self.parameters.get_par_by_name(param.name) is None and param.name!='gamma':
                raise RuntimeError('distr_func parameter',param.name,'not among parameters names',self.parameters.names)
        if 'gamma' not in s.parameters.keys():
            raise RuntimeError('gamma not present in distr_func signature')

    def _set_emitters_type(self,emitters_type):
        self.available_emitters_type = available_emitters_type
        self._check_emitters_type(emitters_type)
        self.emitters_type = emitters_type


    def _set_base_pars(self):
        #model_dic = {}
        a_h, b_h = self.set_bounds(1, 1E15, log_val=self._log_values)
        a_l, b_l = self.set_bounds(1, 1E9, log_val=self._log_values)
        #a_t, b_t = self.set_bounds(1, 1E9, log_val=self._log_values)

        self.parameters.add_par( ModelParameter(name='gmin', par_type='low-energy-cut-off',val=2, val_min=a_l, val_max=b_l,
                                                  units='lorentz-factor', log=self._log_values))

        self.parameters.add_par( ModelParameter(name='gmax', par_type='high-energy-cut-off',val=1E5, val_min=a_h, val_max=b_h,
                                                  units='lorentz-factor', log=self._log_values))

        self.parameters.add_par( ModelParameter(name='N', par_type='emitters_density', val=1, val_min=0, val_max=None, units='cm^-3'))


    def set_bounds(self,a,b,log_val=False):
        if log_val == False:
            return [a,b]

        else:
            return np.log10([a,b])

    def _set_log_val(self,a,log_val=False):
        if log_val == False:
            return a

        else:
            return np.log10(a)


    def _eval_func(self,gamma):
        p_dict = {}
        for par in self.parameters.par_array:
            if par.name not in self._skip:
                p_dict[par.name] = par.val_lin
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

    def set_grid(self):
        gmin = self.parameters.get_par_by_name('gmin').val_lin
        gmax = self.parameters.get_par_by_name('gmax').val_lin
        self._gamma_grid = np.logspace(np.log10(gmin), np.log10(gmax), self._gamma_grid_size)

    def eval_N(self):
        #NOTE: commented until is fixed the setting to zero of Ne_jetset for protons
        #self.update()
        if self.emitters_type=='electrons':
            return np.trapezoid(self.n_gamma_e,self.gamma_e)
        elif self.emitters_type=='protons':
            return np.trapezoid(self.n_gamma_p, self.gamma_p)
        #NOTE: to be added for leptonic-equilibrium
        #elif self.emitters_type=='electrons-equilibrium':
        #    return np.trapezoid(self.n_gamma_e,self.gamma_e)
        else:
            raise  RuntimeError('emitters type',self.emitters_type, 'not valid')

    def eval_U(self, gmin=None, gmax=None):
        #NOTE: commented until is fixed the setting to zero of Ne_jetset for protons
        #self.update()
        if self.emitters_type == 'electrons':# or self.emitters_type == 'electrons-equilibrium':
            x=self.gamma_e
            y=self.n_gamma_e * self.gamma_e
            cost=m_e.cgs.value * (c.cgs ** 2).value

        elif self.emitters_type == 'protons':
            x = self.gamma_p
            y = self.n_gamma_p * self.gamma_p
            cost =  m_p.cgs.value * (c.cgs ** 2).value

        else:
            raise  RuntimeError('emitters type',self.emitters_type, 'not valid')

        msk = np.ones(x.shape, dtype=bool)
        if gmin is not None:
            msk = x>=gmin
        if gmax is not None:
            msk = np.logical_and(msk, x<=gmax)

        return np.trapezoid(y[msk],x[msk]) * cost

    def _fill(self,):
        self.set_grid()

        self.f=self._eval_func(gamma=self._gamma_grid)
        self._Norm=1.0

        if self.normalize is True:
            self._Norm=1.0/np.trapezoid(self.f,self._gamma_grid)

        self.f = self.f*self._Norm*self.parameters.get_par_by_name('N').val

        #NOTE: to be added for leptonic-equilibrium
        # if self.emitters_type == 'electrons-equilibrium':
        #     #TODO implement here
        #     if self._primaries_done is False:
        #         self.gamma_e_inj = np.zeros(self._gamma_grid_size)
        #         self.n_gamma_e_inj = np.zeros(self._gamma_grid_size)

        if self.emitters_type == 'electrons':
            self.gamma_e = np.zeros(self._gamma_grid_size)
            self.n_gamma_e = np.zeros(self._gamma_grid_size)
            self.gamma_e = self._gamma_grid
            self.n_gamma_e = self.f

            self.n_gamma_e[np.isnan(self.n_gamma_e)] = 0
            self.n_gamma_e[np.isinf(self.n_gamma_e)] = 0

        if self.emitters_type == 'protons':
            self.gamma_p = np.zeros(self._gamma_grid_size)
            self.n_gamma_p = np.zeros(self._gamma_grid_size)
            if self._secondaries_done is False:
                self.gamma_e_second_inj = np.zeros(self._gamma_grid_size)
                self.n_gamma_e_second_inj = np.zeros(self._gamma_grid_size)

            self.gamma_p = self._gamma_grid
            self.n_gamma_p = self.f

            self.n_gamma_p[np.isnan(self.n_gamma_p)] = 0
            self.n_gamma_p[np.isinf(self.n_gamma_p)] = 0

            self.n_gamma_e_second_inj[np.isnan(self.n_gamma_e_second_inj)] = 0
            self.n_gamma_e_second_inj[np.isinf(self.n_gamma_e_second_inj)] = 0


    def _plot(self,m, p, y_min=None, y_max=None, x_min=None, x_max=None, energy_unit='gamma',label=None, loglog=False):

        if hasattr(self,'_jet'):
            if self._jet is not None:
                self._set_blob()
                self._jet.eval()
        #NOTE: commented until is fixed the setting to zero of Ne_jetset for protons
        #self.update()
        if self.emitters_type == 'electrons':
            if label is None:
                label = 'electrons'
            m(self.gamma_e,
                         self.n_gamma_e,
                         y_min=y_min,
                         y_max=y_max,
                         x_min=x_min,
                         x_max=x_max,
                         particle='electrons',
                         energy_unit=energy_unit,
                         label=label)
            if getattr(self, '_primaries_done', False) is True:
                if self.gamma_cooling_eq is not None:
                    if energy_unit != 'gamma':
                        eq = self.gamma_cooling_eq * (m_e * c * c).to(energy_unit).value
                    else:
                        eq = self.gamma_cooling_eq
                    if loglog is True:
                        eq = np.log10(eq)
                    p.ax.axvline(eq, ls='--', label='cooling eq. primary', lw=0.5, c='r')

                if hasattr(self, 'gamma_e_inj') and hasattr(self, 'n_gamma_e_inj'):
                    m(self.gamma_e_inj,
                      self.n_gamma_e_inj,
                      y_min=y_min,
                      y_max=y_max,
                      x_min=x_min,
                      x_max=x_max,
                      particle='electrons',
                      energy_unit=energy_unit,
                      label='electrons (inj)')

        if self.emitters_type == 'protons':
            if label is None:
                label = 'protons'
            m(self.gamma_p,
                         self.n_gamma_p,
                         y_min=y_min,
                         y_max=y_max,
                         x_min=x_min,
                         x_max=x_max,
                         particle='protons',
                         energy_unit=energy_unit,
                         label=label)
            if getattr(self, '_secondaries_done', False) is True:
                if hasattr(self,'n_gamma_e'):
                    m(self.gamma_e,
                                 self.n_gamma_e,
                                 y_min=y_min,
                                 y_max=y_max,
                                 x_min=x_min,
                                 x_max=x_max,
                                 particle='electrons',
                                 energy_unit=energy_unit,
                                 label='electrons sec.')
                    if energy_unit != 'gamma':

                        eq= self.gamma_cooling_eq_second* (m_e * c * c).to(energy_unit).value
                    else:
                        eq = self.gamma_cooling_eq_second
                    if loglog is True:
                        eq = np.log10(eq)
                    p.ax.axvline(eq, ls='--', label='cooling. eq. second.', lw=0.5, c='r')
                
                if hasattr(self, 'gamma_e_second_inj'):
                    m(self.gamma_e_second_inj,
                      self.n_gamma_e_second_inj,
                      y_min=y_min,
                      y_max=y_max,
                      x_min=x_min,
                      x_max=x_max,
                      particle='electrons',
                      energy_unit=energy_unit,
                      label='electrons sec. (inj)')
            

        return p

    def plot(self, p=None, y_min=None, y_max=None, x_min=None, x_max=None, energy_unit='gamma',label=None,loglog=False):
        if p is None:
            p = PlotPdistr(loglog=loglog)
        m=getattr(p,'plot_distr')
        self._plot(m,p,y_min=y_min,y_max=y_max,x_min=x_min,x_max=x_max,energy_unit=energy_unit,label=label,loglog=loglog)
        return p

    def plot2p(self, p=None, y_min=None, y_max=None, x_min=None, x_max=None, energy_unit='gamma',label=None,loglog=False):
        if p is None:
            p = PlotPdistr(loglog=loglog)
        m=getattr(p,'plot_distr2p')
        self._plot(m,p,y_min=y_min,y_max=y_max,x_min=x_min,x_max=x_max,energy_unit=energy_unit,label=label,loglog=loglog)

        return p

    def plot3p(self, p=None, y_min=None, y_max=None, x_min=None, x_max=None, energy_unit='gamma',label=None,loglog=False):
        if p is None:
            p = PlotPdistr(loglog=loglog)
        m=getattr(p,'plot_distr3p')
        self._plot(m,p,y_min=y_min,y_max=y_max,x_min=x_min,x_max=x_max,energy_unit=energy_unit,label=label,loglog=loglog)
        return p


class EmittersDistribution(BaseEmittersDistribution):

    def __init__(self,
                 name,
                 spectral_type,
                 jet=None,
                 gamma_grid_size=200,
                 log_values=False,
                 emitters_type='electrons',
                 normalize=False,
                 skip_build=False):

        super(EmittersDistribution, self).__init__(name,
                                                   spectral_type=spectral_type,
                                                   emitters_type=emitters_type,
                                                   gamma_grid_size=gamma_grid_size,
                                                   skip_build=True,
                                                   normalize=normalize,
                                                   log_values=log_values)


        if isinstance(self,EmittersArrayDistribution):
            self._array_gamma = None
            self._array_n_gamma = None
        else:
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
            par.set(val=p.val, skip_dep_par_warning=True)

    def update(self):
        if self.emitters_type == 'electrons' and hasattr(self, '_jet') and self._jet is not None:
            if getattr(self._jet._blob.emitters, 'do_equilibrium', 0) == 1:
                # In equilibrium mode the physical electron distribution comes from C readback.
                self._set_blob()
                return
        self._fill()

    def _build(self, jet, name, log_values, gamma_grid_size, normalize):
        self._user_defined=True
        self._name = name
        self._log_values = log_values
        self._gamma_grid_size=gamma_grid_size
        self.parameters = ModelParameterArray()
        self.set_parameters_dict()
        self.normalize = normalize
        self.gamma_cooling_eq = None
        self.gamma_cooling_eq_second=None
        self._primaries_done = False
        self._secondaries_done=False
        self.set_jet(jet)

    def _update_parameters_dict(self):
        for par in self.parameters.par_array:
            self._parameters_dict[par.name].val = par.val
            self._parameters_dict[par.name].vmin = par.val_min
            self._parameters_dict[par.name].vmax = par.val_max
            self._parameters_dict[par.name].log = par.islog
            self._parameters_dict[par.name].punit = par.units
            self._parameters_dict[par.name]._is_dependent = par._is_dependent
            self._parameters_dict[par.name]._depending_pars = par._depending_pars
            self._parameters_dict[par.name]._func = par._func
            self._parameters_dict[par.name]._master_pars = par._master_pars

    def add_par(self, name, par_type, val, vmax, vmin, unit='', log=False, frozen=False):
        #if log is True:
        #    val = np.log10(val)

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

        #print('==> add par',name,val,log)
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
        if isinstance(self,EmittersArrayDistribution):
            model_dict['N'].val = 1.0
        else:
            model_dict['N'].val = 100.
        model_dict['N'].is_in_jetkernel = True

        if 'NH_pp' in model_dict.keys():
            model_dict['NH_pp'].val = 1.
            model_dict['NH_pp'].is_in_jetkernel = True

        self._parameters_dict=model_dict
        for k,par in model_dict.items():
            if par.log is True:
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

    def _set_base_pars(self):
        raise RuntimeError('allowed only in parent class')

    def _get_base_dict_and_bounds(self):
        model_dic = {}
        a_h, b_h = self.set_bounds(1, 1E15, log_val=self._log_values)
        a_l, b_l = self.set_bounds(1, 1E9, log_val=self._log_values)
        a_t, b_t = self.set_bounds(1, 1E9, log_val=self._log_values)
        model_dic['gmin'] = JetModelDictionaryPar(ptype='low-energy-cut-off', vmin=a_l, vmax=b_l,
                                                  punit='lorentz-factor', log=self._log_values,
                                                  jetkernel_par_name='emitters.gmin')
        model_dic['gmax'] = JetModelDictionaryPar(ptype='high-energy-cut-off', vmin=a_h, vmax=b_h,
                                                  punit='lorentz-factor', log=self._log_values,
                                                  jetkernel_par_name='emitters.gmax')

        model_dic['N'] = JetModelDictionaryPar(ptype='emitters_density', vmin=0, vmax=None, punit='cm^-3',
                                               jetkernel_par_name='emitters.N')

        if self.emitters_type =='protons':
            model_dic['NH_pp'] = JetModelDictionaryPar(ptype='target_density', vmin=0, vmax=None, punit='cm^-3',
                                                       froz=False, jetkernel_par_name='PP_gamma.NH_pp')

        if isinstance(self, EmittersArrayDistribution):
            model_dic['N'] = JetModelDictionaryPar(ptype='scaling_factor', vmin=0, vmax=None, punit='',
                                                   jetkernel_par_name='emitters.N')
        return model_dic, a_h, b_h, a_l, b_l, a_t, b_t

    def set_jet(self, jet):
        if jet is not None:

            #name passed to the C code
            #do not change
            name = 'jetset'

            self._jet = jet
            set_str_attr(jet._blob, 'core.DISTR', name)
            set_str_attr(jet._blob, 'core.PARTICLE', self.emitters_type)

            p = self._jet.get_par_by_name('gmin')
            if p is not None:
                p.set(val=self.parameters.get_par_by_name('gmin').val_lin)
            else:
                p=self.parameters.get_par_by_name('gmin')
            set_nested_attr(jet._blob, 'emitters.gmin', p.val_lin)

            p = self._jet.get_par_by_name('gmax')
            if p is not None:
                p.set(val=self.parameters.get_par_by_name('gmax').val_lin)
            else:
                p = self.parameters.get_par_by_name('gmax')
            set_nested_attr(jet._blob, 'emitters.gmax', p.val_lin)

            self._jet._blob.emitters.gamma_grid_size = self._gamma_grid_size
        else:
            self._jet=jet


    def set_grid_size(self,gamma_grid_size):
        if gamma_grid_size is not None:
            if self._jet is not  None:
              set_nested_attr(self._jet._blob, 'emitters.gamma_grid_size', gamma_grid_size)
        self._fill()

    def set_grid(self):
        if self._jet is None:
            gmin = self.parameters.get_par_by_name('gmin').val_lin
            gmax = self.parameters.get_par_by_name('gmax').val_lin
            self._gamma_grid = np.logspace(np.log10(gmin), np.log10(gmax), self._gamma_grid_size)
        else:
            BlazarSED.setNgrid(self._jet._blob)
            if self.emitters_type == 'electrons':
                #print('==> Build Ne start')
                BlazarSED.build_Ne_jetset(self._jet._blob)
                gamma_ptr = get_nested_attr(self._jet._blob, 'emitters.griglia_gamma_jetset_Ne_log')
                #print('==> Build Ne done')
            elif self.emitters_type == 'protons':
                BlazarSED.build_Np_jetset(self._jet._blob)
                gamma_ptr = get_nested_attr(self._jet._blob, 'emitters.griglia_gamma_jetset_Np_log')
            #NOTE: to be added for leptonic-equilibrium
            # elif self.emitters_type == 'electrons-equilibrium':
            #     BlazarSED.build_Ne_jetset(self._jet._blob)
            #     gamma_ptr = get_nested_attr(self._jet._blob, 'emitters.griglia_gamma_jetset_Ne_log')
            #     #print('==> Build Ne done')
            else:
                raise RuntimeError('emitters type', self.emitters_type, 'not valid')

            size = self._jet._blob.emitters.gamma_grid_size
            self._gamma_grid = np.zeros(size)
            for ID in range(size):
                self._gamma_grid[ID] = BlazarSED.get_elec_array(gamma_ptr, self._jet._blob, ID)
            self._gamma_grid_size=self._jet.gamma_grid_size
           



    def _fill(self,):
        super(EmittersDistribution, self)._fill()
        if self._jet is not None:

            if self.emitters_type == 'electrons':
                size = self._gamma_grid_size
                Ne_ptr = get_nested_attr(self._jet._blob, 'emitters.Ne_jetset')
                set_emitters(Ne_ptr,self._jet._blob,size,self.f)

            if self.emitters_type == 'protons':
                size = self._gamma_grid_size
                Np_ptr = get_nested_attr(self._jet._blob, 'emitters.Np_jetset')
                set_emitters(Np_ptr,self._jet._blob,size,self.f)
 


    def _set_blob(self):
        size = self._jet._blob.emitters.gamma_grid_size
        if self.emitters_type == 'electrons':

            self._Ne_name, self._gammae_name = gamma_dic_e['electron_distr']
            self.Ne_ptr = get_nested_attr(self._jet._blob, self._Ne_name)
            self.e_gamma_ptr = get_nested_attr(self._jet._blob, self._gammae_name)
            self.gamma_e,self.n_gamma_e=get_emitters(self.e_gamma_ptr,self.Ne_ptr, self._jet._blob,size)
            self._primaries_done = False
            self.gamma_cooling_eq = None
            if getattr(self._jet._blob.emitters, 'do_equilibrium', 0) == 1:
                self._Q_inj_e_name, self._gammae_inj_name = gamma_dic_e_equilibrium['e_inj']
                self._Q_inj_e_ptr = get_nested_attr(self._jet._blob, self._Q_inj_e_name)
                self.e_inj_gamma_ptr = get_nested_attr(self._jet._blob, self._gammae_inj_name)
                self.gamma_e_inj, self.n_gamma_e_inj = get_emitters(self.e_inj_gamma_ptr,
                                                                    self._Q_inj_e_ptr,
                                                                    self._jet._blob,
                                                                    size)
                self.gamma_cooling_eq = self._jet._blob.emitters.gamma_cooling_eq
                self._primaries_done = True
        
        elif self.emitters_type == 'protons':
            self._Ne_name, self._gammae_name = gamma_dic_e['electron_distr']
            self._Np_name, self._gammap_name = gamma_dic_p['proton_distr']
            self._Q_inj_e_second_name, self._gammae_inj_sec_name = gamma_dic_pp_e_second['e_second_inj']

            self.Np_ptr = get_nested_attr(self._jet._blob, self._Np_name)
            self.p_gamma_ptr = get_nested_attr(self._jet._blob, self._gammap_name)
            self.Ne_ptr = get_nested_attr(self._jet._blob, self._Ne_name)
            self.e_gamma_ptr = get_nested_attr(self._jet._blob, self._gammae_name)
            self._Q_inj_e_second_ptr = get_nested_attr(self._jet._blob, self._Q_inj_e_second_name)
            self.e_inj_second_gamma_ptr = get_nested_attr(self._jet._blob, self._gammae_inj_sec_name)

        else:
            raise RuntimeError(f"emitters type {self.emitters_type} not valid ", available_emitters_type)
       
        
        #NOTE: this is needed to get the pointers to Ne also in the case of protons
        self.gamma_e,self.n_gamma_e=get_emitters(self.e_gamma_ptr,self.Ne_ptr, self._jet._blob,size)
        if self.emitters_type == 'protons':
            self.gamma_p,self.n_gamma_p=get_emitters(self.p_gamma_ptr,self.Np_ptr, self._jet._blob,size)
            self.gamma_e_second_inj,self.n_gamma_e_second_inj=get_emitters(self.e_inj_second_gamma_ptr,self._Q_inj_e_second_ptr, self._jet._blob,size)
            self.gamma_cooling_eq_second= self._jet._blob.emitters.gamma_cooling_eq
            self._secondaries_done = True
       
        
        
        
              
     
    
    def _activate_numba(self):
        self._py_distr_func=copy.deepcopy(self.distr_func)
        try:
            if is_jitted(self.distr_func) is False:
                self.distr_func = njit(self.distr_func,fastmath=True, cache=False)
        except Exception as e:
           print("failed to activate numba for",str(e))


class EmittersArrayDistribution(EmittersDistribution):
    def __init__(self,
                 name,
                 jet=None,
                 emitters_type='electrons',
                 normalize=False,
                 skip_build=False,
                 gamma_array=None,
                 n_gamma_array=None,
                 gamma_grid_size=None):

        super(EmittersArrayDistribution, self).__init__(name,
                                                        spectral_type='array',
                                                        emitters_type=emitters_type)

        if gamma_array is not None and n_gamma_array is not None:
            if self._spectral_type != 'array':
                raise RuntimeError('you can pass  gamma_array and n_gamma_array only for array distribution')
            self._set_arrays(gamma_array,n_gamma_array)
        elif (gamma_array is not None) and (n_gamma_array is not None) is False:
            raise RuntimeError('you have to pass both gamma_array and n_gamma_array for array distribution')

        if gamma_grid_size is None:
            gamma_grid_size=gamma_array.size
        if skip_build is False:
            self._build(jet, name, False, gamma_grid_size, normalize)
        self.parameters.gmin.val=self._array_gamma[0]
        self.parameters.gmax.val = self._array_gamma[-1]
        self.parameters.N.val=1
        self.set_distr_func(self._array_func)

    def _set_arrays(self,gamma_array,n_gamma_array):
        _ids = np.argsort(gamma_array)
        self._array_gamma = gamma_array[_ids]
        self._array_n_gamma = n_gamma_array[_ids]

    def _array_func(self,gamma):
        msk_nan=self._array_n_gamma>0
        if msk_nan.sum()>1:
            f_interp = interpolate.interp1d(np.log10(self._array_gamma[msk_nan]), np.log10(self._array_n_gamma[msk_nan]), bounds_error=False, kind='linear')
            y = np.power(10., f_interp(np.log10(gamma)))
            msk_nan = np.isnan(y)
            y[msk_nan] = 0
        else:
            y = np.zeros(gamma.size)

        return y

    def _activate_numba(self):
        pass

class InjEmittersDistribution(BaseEmittersDistribution):

    def __init__(self,
                 name,
                 spectral_type,
                 gamma_grid_size=200,
                 log_values=False,
                 skip_build=False,
                 normalize=True,
                 emitters_type='electrons'):

        super(InjEmittersDistribution, self).__init__(name=name,
                                                      spectral_type=spectral_type,
                                                      emitters_type='electrons',
                                                      gamma_grid_size=gamma_grid_size,
                                                      skip_build=True,
                                                      normalize=True,
                                                      log_values=log_values)

        if skip_build is False:
            self._build(name, log_values, gamma_grid_size, normalize=True)


    def _build(self, name, log_values, gamma_grid_size, normalize):
        self._user_defined = True
        self._name = name
        self._log_values = log_values
  
        self.parameters = ModelParameterArray()
        self.set_parameters_dict()
        self.parameters.add_par( ModelParameter(name='Q', par_type='emitters_density', val=1E-3, val_min=0, val_max=None, units='cm-3 s-1'))
        self._parameters_dict['Q'] = JetModelDictionaryPar(ptype='emitters_density',
                                                           vmin=0,
                                                           vmax=None,
                                                           punit='cm-3 s-1',
                                                           val=1E-3,
                                                           log=False,
                                                           froz=False,
                                                           is_in_jetkernel=False)
        self.normalize = normalize
        self.gamma_cooling_eq = None
        self.gamma_cooling_eq_second=None
        self._primaries_done = False
        self._secondaries_done=False
        self._gamma_grid_size = gamma_grid_size
        self._reset_equilibrium_state()
    
    #def _set_base_pars(self):
    #
    #
    #    self.parameters.add_par( JetModelDictionaryPar(name='Q', par_type='emitters_density', val=1E-3, val_min=0, val_max=None, units='cm-3 s-1'))
        
    def set_parameters_dict(self,):
        model_dict, a_h, b_h, a_l, b_l, a_t, b_t = self._get_base_dict_and_bounds()
        model_dict['gmin'].val = 2.0
        model_dict['gmin'].is_in_jetkernel = True

        model_dict['gmax'].val = 1E6
        model_dict['gmax'].is_in_jetkernel = True

        self._parameters_dict=model_dict
        for k,par in model_dict.items():
            if par.log is True:
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

    def _get_base_dict_and_bounds(self):
        model_dic = {}
        a_h, b_h = self.set_bounds(1, 1E15, log_val=self._log_values)
        a_l, b_l = self.set_bounds(1, 1E9, log_val=self._log_values)
        a_t, b_t = self.set_bounds(1, 1E9, log_val=self._log_values)
        model_dic['gmin'] = JetModelDictionaryPar(ptype='low-energy-cut-off', vmin=a_l, vmax=b_l,
                                                  punit='lorentz-factor', log=self._log_values,
                                                  jetkernel_par_name='emitters.gmin')
        model_dic['gmax'] = JetModelDictionaryPar(ptype='high-energy-cut-off', vmin=a_h, vmax=b_h,
                                                  punit='lorentz-factor', log=self._log_values,
                                                  jetkernel_par_name='emitters.gmax')
        return model_dic, a_h, b_h, a_l, b_l, a_t, b_t

    def _fill(self,):
        self.set_grid()
        self.f=self._eval_func(gamma=self._gamma_grid)
        self._Norm=1.0
        if self.normalize is True:
            _integ = np.trapezoid(self.f,self._gamma_grid)
            if _integ > 0:
                self._Norm=1.0/_integ
        self.f = self.f*self._Norm*self.parameters.get_par_by_name('Q').val
        self.gamma_e = self._gamma_grid.copy()
        self.n_gamma_e = np.asarray(self.f, dtype=np.float64)
        self.n_gamma_e[np.isnan(self.n_gamma_e)] = 0
        self.n_gamma_e[np.isinf(self.n_gamma_e)] = 0

        if hasattr(self, '_jet') and self._jet is not None:
            if getattr(self._jet._blob.emitters, 'do_equilibrium', 0) == 1:
                size = self._gamma_grid_size
                q_ptr = get_nested_attr(self._jet._blob, 'emitters.Q_inj_e_primaries')
                set_emitters(q_ptr,self._jet._blob,size,self.n_gamma_e)

    def eval_U_q(self):
        self._fill()
        x=self.gamma_e
        y=self.f*x
        cost=m_e.cgs.value * (c.cgs ** 2).value
        return np.trapezoid(y,x) * cost
 

    def _set_L_inj(self, L_inj_target_erg, volume):
        if L_inj_target_erg >0:
             self.parameters.Q.val *=L_inj_target_erg/(self.eval_U_q() * volume )
        #NOTE: commented until is fixed the setting to zero of Ne_jetset for protons
        #self.update()
    
    #NOTE: not used, to be removed
    def set_temp_ev(self):
        self.e_gamma_ptr = getattr(self._temp_ev, self._gammae_name)
        self._Q_inj_e_second_ptr = get_nested_attr(self._temp_ev._blob, self._Q_inj_e_second_name)

    def add_par(self, name, par_type, val, vmax, vmin, unit='', log=False, frozen=False):
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



class InjEmittersArrayDistribution(InjEmittersDistribution):
    def __init__(self,
                 name,
                 emitters_type='electrons',
                 normalize=False,
                 skip_build=False,
                 gamma_array=None,
                 n_gamma_array=None,
                 gamma_grid_size=None):

        super(InjEmittersArrayDistribution, self).__init__(name,
                                                        spectral_type='array',
                                                        emitters_type=emitters_type)

        if gamma_array is not None and n_gamma_array is not None:
            if self._spectral_type != 'array':
                raise RuntimeError('you can pass  gamma_array and n_gamma_array only for array distribution')
            self._set_arrays(gamma_array,n_gamma_array)
        elif (gamma_array is not None) and (n_gamma_array is not None) is False:
            raise RuntimeError('you have to pass both gamma_array and n_gamma_array for array distribution')

        if gamma_grid_size is None:
            gamma_grid_size=gamma_array.size
        if skip_build is False:
            self._build( name, False, gamma_grid_size, normalize)
        self.parameters.gmin.val=self._array_gamma[0]
        self.parameters.gmax.val = self._array_gamma[-1]
        self.parameters.Q.val =1
        self.set_distr_func(self._array_func)

    def _set_arrays(self,gamma_array,n_gamma_array):
        _ids = np.argsort(gamma_array)
        self._array_gamma = gamma_array[_ids]
        self._array_n_gamma = n_gamma_array[_ids]

    def _array_func(self,gamma):
        msk_nan=self._array_n_gamma>0
        f_interp = interpolate.interp1d(np.log10(self._array_gamma[msk_nan]), np.log10(self._array_n_gamma[msk_nan]), bounds_error=False, kind='linear')
        y = np.power(10., f_interp(np.log10(gamma)))
        msk_nan = np.isnan(y)
        y[msk_nan] = 0

        return y
    
    def _activate_numba(self):
        pass
