__author__ = "Andrea Tramacere"

import numpy as np
import numba as nb
import pprint
from jetset.jet_emitters import EmittersDistribution, InjEmittersDistribution
__all__=['EmittersFactory']

_available_dict = {'lp': 'log-parabola',
                   'pl': 'powerlaw',
                   'lppl': 'log-parabola with low-energy powerlaw branch',
                   'lpep': 'log-parabola defined by peak energy',
                   'plc': 'powerlaw with cut-off',
                   'bkn': 'broken powerlaw',
                   'superexp': 'powerlaw with super-exp cut-off'}


@nb.njit(fastmath=True, cache=True)
def distr_func_bkn(gamma_break, gamma, p, p_1):
    f = np.zeros(gamma.shape)
    m = gamma < gamma_break
    f[m] = np.power(gamma[m], -p)
    f[~m] = np.power(gamma_break, -(p - p_1)) * np.power(gamma[~m], -p_1)
    return f


@nb.njit(fastmath=True, cache=True)
def distr_func_pl(gamma, p, ):
    return np.power(gamma, -p)


@nb.njit(fastmath=True, cache=True)
def distr_func_plc(gamma, gamma_cut, p,):
    return np.power(gamma, -p) * np.exp(-(gamma / gamma_cut) )


@nb.njit(fastmath=True, cache=True)
def distr_func_super_exp(gamma, gamma_cut, p, a):
    return np.power(gamma, -p) * np.exp(-(1 / a) * (gamma / gamma_cut) ** a)


@nb.njit(fastmath=True, cache=True)
def distr_func_lp(gamma, gamma0_log_parab, r, s):
    return np.power((gamma / gamma0_log_parab), (-s - r * np.log10(gamma / gamma0_log_parab)))


@nb.njit(fastmath=True, cache=True)
def distr_func_lep(gamma, gamma_p, r):
    return np.power(10., (-r * np.power(np.log10(gamma / gamma_p), 2)))


@nb.njit(fastmath=True, cache=True)
def distr_func_lppl(gamma, gamma0_log_parab, r, s):
    f = np.zeros(gamma.shape)
    m = gamma < gamma0_log_parab
    f[m] = np.power(gamma[m]/gamma0_log_parab, -s)
    f[~m] = np.power(gamma[~m]/gamma0_log_parab, (-s -r*np.log10(gamma[~m] / gamma0_log_parab )))
    return f


class EmittersFactory:
    def __repr__(self):
        return str(pprint.pprint(self._available_dict))

    def __init__(self):
        self._available_dict=_available_dict

        self._func_dict = {'pl': self._create_pl,
                           'bkn': self._create_bkn,
                           'plc': self._create_plc,
                           'lp': self._create_lp,
                           'lppl': self._create_lppl,
                           'lpep': self._create_lpep,
                           'superexp': self._create_super_exp}

        self._set_emitters_class()

    def _set_emitters_class(self):
        if type(self) == EmittersFactory:
            self._emitters_class=EmittersDistribution
        elif type(self) ==InjEmittersFactory:
            self._emitters_class = InjEmittersDistribution
        else:
            raise RuntimeError('class instance', type(self), 'not valid')

    @staticmethod
    def available_distributions():
        for k in _available_dict.keys():
            print('%s: %s' % (k, _available_dict[k]))

    @staticmethod
    def available_distributions_list():
        return  _available_dict.keys()

    def create_emitters(self,
                        name,
                        gamma_grid_size=200,
                        log_values=False,
                        emitters_type='electrons',
                        normalize=True,
                        skip_build=False):

        if name not in self._available_dict.keys():
            raise RuntimeError('name', name, 'not among available', self._available_dict.keys())

        return self._func_dict[name](gamma_grid_size=gamma_grid_size,
                                     log_values=log_values,
                                     normalize=normalize,
                                     skip_build=skip_build,
                                     emitters_type=emitters_type)

    def _create_bkn(self,gamma_grid_size,log_values,normalize,skip_build,emitters_type):

        n_e_bkn = self._emitters_class(name='bkn',
                                       spectral_type='bkn',
                                       normalize=normalize,
                                       emitters_type=emitters_type,
                                       log_values=log_values,
                                       skip_build=skip_build,
                                       gamma_grid_size=gamma_grid_size)

        a_t, b_t = n_e_bkn.set_bounds(1, 1E9, log_val=n_e_bkn._log_values)
        gamma_break_val = n_e_bkn._set_log_val(1E4,log_val=log_values)
        n_e_bkn.add_par('gamma_break', par_type='turn-over-energy', val=gamma_break_val, vmin=a_t, vmax=b_t, unit='lorentz-factor',log=log_values)
        n_e_bkn.add_par('p', par_type='LE_spectral_slope', val=2.5, vmin=-10., vmax=10, unit='')
        n_e_bkn.add_par('p_1', par_type='HE_spectral_slope', val=3.5, vmin=-10., vmax=10, unit='')
        n_e_bkn.set_distr_func(distr_func_bkn)
        return n_e_bkn

    def _create_pl(self,gamma_grid_size, log_values, normalize, skip_build, emitters_type):

        n_e_pl = self._emitters_class(name='pl',
                                       spectral_type='pl',
                                       normalize=normalize,
                                       emitters_type=emitters_type,
                                       log_values=log_values,
                                       skip_build=skip_build,
                                       gamma_grid_size=gamma_grid_size)

        a_t, b_t = n_e_pl.set_bounds(1, 1E9, log_val=n_e_pl._log_values)

        n_e_pl.add_par('p', par_type='LE_spectral_slope', val=2.0, vmin=-10., vmax=10, unit='')
        n_e_pl.set_distr_func(distr_func_pl)

        return n_e_pl


    def _create_plc(self, gamma_grid_size,log_values,normalize,skip_build,emitters_type):

        n_e_plc = n_e_bkn = self._emitters_class(name='plc',
                                             spectral_type='plc',
                                             normalize=normalize,
                                             emitters_type=emitters_type,
                                             log_values=log_values,
                                             skip_build=skip_build,
                                             gamma_grid_size=gamma_grid_size)

        a_t, b_t = n_e_plc.set_bounds(1, 1E9, log_val=n_e_plc._log_values)
        gamma_cut_val = n_e_plc._set_log_val(1E4,log_val=log_values)
        n_e_plc.add_par('gamma_cut', par_type='turn-over-energy', val=gamma_cut_val, vmin=a_t, vmax=b_t,
                              unit='lorentz-factor',log=log_values)
        n_e_plc.add_par('p', par_type='LE_spectral_slope', val=2.0, vmin=-10., vmax=10, unit='')
        n_e_plc.set_distr_func(distr_func_plc)
        return n_e_plc


    def _create_super_exp(self, gamma_grid_size, log_values, normalize, skip_build, emitters_type):

        n_e_super_exp  = self._emitters_class(name='super_exp',
                                             spectral_type='plc',
                                             normalize=normalize,
                                             emitters_type=emitters_type,
                                             log_values=log_values,
                                             skip_build=skip_build,
                                             gamma_grid_size=gamma_grid_size)

        a_t, b_t = n_e_super_exp.set_bounds(1, 1E9, log_val=n_e_super_exp._log_values)
        gamma_cut_val = n_e_super_exp._set_log_val(1E4,log_val=log_values)
        n_e_super_exp.add_par('gamma_cut', par_type='turn-over-energy', val=gamma_cut_val, vmin=a_t, vmax=b_t,
                              unit='lorentz-factor',log=log_values)
        n_e_super_exp.add_par('p', par_type='LE_spectral_slope', val=2.0, vmin=-10., vmax=10, unit='')
        n_e_super_exp.add_par('a', par_type='spectral_curvature', val=1.0, vmin=0., vmax=100., unit='')
        n_e_super_exp.set_distr_func(distr_func_super_exp)

        return n_e_super_exp

    def _create_lp(self, gamma_grid_size,log_values,normalize,skip_build,emitters_type):

        n_lp  = self._emitters_class(name='lp',
                                    spectral_type='lp',
                                    normalize=normalize,
                                    emitters_type=emitters_type,
                                    log_values=log_values,
                                    skip_build=skip_build,
                                    gamma_grid_size=gamma_grid_size)

        a_t, b_t = n_lp.set_bounds(1, 1E9, log_val=n_lp._log_values)
        gamma0_log_parab_val = n_lp._set_log_val(1E4,log_val=log_values)
        n_lp.add_par('gamma0_log_parab', par_type='turn-over-energy', val=gamma0_log_parab_val, vmin=a_t, vmax=b_t,
                              unit='lorentz-factor',log=log_values)
        n_lp.add_par('s', par_type='LE_spectral_slope', val=2.0, vmin=-10., vmax=10, unit='')
        n_lp.add_par('r', par_type='spectral_curvature', val=0.4, vmin=-15., vmax=15., unit='')

        n_lp.set_distr_func(distr_func_lp)

        return n_lp

    def _create_lpep(self,gamma_grid_size, log_values, normalize, skip_build, emitters_type):

        n_lep  = self._emitters_class(name='lpep',
                                    spectral_type='lp',
                                    normalize=normalize,
                                    emitters_type=emitters_type,
                                    log_values=log_values,
                                    skip_build=skip_build,
                                    gamma_grid_size=gamma_grid_size)

        a_t, b_t = n_lep.set_bounds(1, 1E9, log_val=n_lep._log_values)
        gamma_p_val = n_lep._set_log_val(1E4,log_val=log_values)
        n_lep.add_par('gamma_p', par_type='turn-over-energy', val=gamma_p_val, vmin=a_t, vmax=b_t,
                     unit='lorentz-factor',log=log_values)
        n_lep.add_par('r', par_type='spectral_curvature', val=0.4, vmin=-15., vmax=15., unit='')
        n_lep.set_distr_func(distr_func_lep)

        return n_lep

    def _create_lppl(self, gamma_grid_size, log_values, normalize, skip_build, emitters_type):

        n_lppl  = self._emitters_class(name='lppl',
                                    spectral_type='lp',
                                    normalize=normalize,
                                    emitters_type=emitters_type,
                                    log_values=log_values,
                                    skip_build=skip_build,
                                    gamma_grid_size=gamma_grid_size)

        a_t, b_t = n_lppl.set_bounds(1, 1E9, log_val=n_lppl._log_values)
        gamma0_log_parab_val = n_lppl._set_log_val(1E4,log_val=log_values)

        n_lppl.add_par('gamma0_log_parab', par_type='turn-over-energy', val=gamma0_log_parab_val, vmin=a_t, vmax=b_t,
                     unit='lorentz-factor',log=log_values)
        n_lppl.add_par('s', par_type='LE_spectral_slope', val=2.0, vmin=-10., vmax=10, unit='')
        n_lppl.add_par('r', par_type='spectral_curvature', val=0.4, vmin=-15., vmax=15., unit='')

        n_lppl.set_distr_func(distr_func_lppl)

        return n_lppl


class InjEmittersFactory(EmittersFactory):

    def __init__(self):
        super(InjEmittersFactory, self).__init__()



    def create_inj_emitters(self,
                        name,
                        gamma_grid_size=200,
                        log_values=False,
                        emitters_type='electrons',
                        normalize=True,
                        skip_build=False):

        if name not in self._available_dict.keys():
            raise RuntimeError('name', name, 'not among available', self._available_dict.keys())
        return self._func_dict[name](gamma_grid_size=gamma_grid_size,
                                     log_values=log_values,
                                     normalize=normalize,
                                     skip_build=skip_build,
                                     emitters_type=emitters_type)