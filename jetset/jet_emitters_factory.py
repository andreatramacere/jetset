__author__ = "Andrea Tramacere"

import numpy as np
import numba as nb
import pprint
from jetset.jet_emitters import EmittersDistribution
__all__=['EmittersFactory']

_available_dict = {'lp': 'log-parabola',
                                'pl': 'powerlaw',
                                'lppl': 'log-parabola with low-energy powerlaw branch',
                                'lpep': 'log-parabola defined by peak energy',
                                'plc': 'powerlaw with cut-off',
                                'bkn': 'broken powerlaw',
                                'superexp':'powerlaw with super-exp cut-off'}
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
                           'sperexp': self._create_super_exp}

    @staticmethod
    def available_distributions():
        for k in _available_dict.keys():
            print('%s: %s' % (k, _available_dict[k]))

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

    @staticmethod
    def _create_bkn(gamma_grid_size,log_values,normalize,skip_build,emitters_type):
        @nb.njit(fastmath=True, cache=True)
        def distr_func_bkn(gamma_break, gamma, s1, s2):
            f = np.zeros(gamma.shape)
            m = gamma < gamma_break
            f[m] = np.power(gamma[m], -s1)
            f[~m] = np.power(gamma_break, -(s1 - s2)) * np.power(gamma[~m], -s2)
            return f

        n_e_bkn = EmittersDistribution(name='bkn',
                                       spectral_type='bkn',
                                       normalize=normalize,
                                       emitters_type=emitters_type,
                                       log_values=log_values,
                                       skip_build=skip_build,
                                       gamma_grid_size=gamma_grid_size)

        a_t, b_t = n_e_bkn.set_bounds(1, 1E9, log_val=n_e_bkn._log_values)

        n_e_bkn.add_par('gamma_break', par_type='turn-over-energy', val=1E4, vmin=a_t, vmax=b_t, unit='lorentz-factor')
        n_e_bkn.add_par('s1', par_type='LE_spectral_slope', val=2.5, vmin=-10., vmax=10, unit='')
        n_e_bkn.add_par('s2', par_type='HE_spectral_slope', val=3.5, vmin=-10., vmax=10, unit='')
        n_e_bkn.set_distr_func(distr_func_bkn)
        return n_e_bkn

    @staticmethod
    def _create_pl(gamma_grid_size, log_values, normalize, skip_build, emitters_type):
        @nb.njit(fastmath=True, cache=True)
        def distr_func_super_exp(gamma, p, ):
            return np.power(gamma, -p)

        n_e_pl = EmittersDistribution(name='pl',
                                       spectral_type='pl',
                                       normalize=normalize,
                                       emitters_type=emitters_type,
                                       log_values=log_values,
                                       skip_build=skip_build,
                                       gamma_grid_size=gamma_grid_size)

        a_t, b_t = n_e_pl.set_bounds(1, 1E9, log_val=n_e_pl._log_values)

        n_e_pl.add_par('p', par_type='LE_spectral_slope', val=2.0, vmin=-10., vmax=10, unit='')
        n_e_pl.set_distr_func(distr_func_super_exp)

        return n_e_pl

    @staticmethod
    def _create_plc(gamma_grid_size,log_values,normalize,skip_build,emitters_type):
        @nb.njit(fastmath=True, cache=True)
        def distr_func_plc(gamma, gamma_cut, p,):
            return np.power(gamma, -p) * np.exp(-(gamma / gamma_cut) )

        n_e_plc = EmittersDistribution(name='plc',
                                             spectral_type='plc',
                                             normalize=normalize,
                                             emitters_type=emitters_type,
                                             log_values=log_values,
                                             skip_build=skip_build,
                                             gamma_grid_size=gamma_grid_size)

        a_t, b_t = n_e_plc.set_bounds(1, 1E9, log_val=n_e_plc._log_values)
        n_e_plc.add_par('gamma_cut', par_type='turn-over-energy', val=1E4, vmin=a_t, vmax=b_t,
                              unit='lorentz-factor',log=log_values)
        n_e_plc.add_par('p', par_type='LE_spectral_slope', val=2.0, vmin=-10., vmax=10, unit='')
        n_e_plc.set_distr_func(distr_func_plc)
        return n_e_plc

    @staticmethod
    def _create_super_exp(gamma_grid_size, log_values, normalize, skip_build, emitters_type):
        @nb.njit(fastmath=True, cache=True)
        def distr_func_super_exp(gamma, gamma_cut, p, a):
            return np.power(gamma, -p) * np.exp(-(1 / a) * (gamma / gamma_cut) ** a)

        n_e_super_exp = EmittersDistribution(name='plc',
                                             spectral_type='plc',
                                             normalize=normalize,
                                             emitters_type=emitters_type,
                                             log_values=log_values,
                                             skip_build=skip_build,
                                             gamma_grid_size=gamma_grid_size)

        a_t, b_t = n_e_super_exp.set_bounds(1, 1E9, log_val=n_e_super_exp._log_values)
        n_e_super_exp.add_par('gamma_cut', par_type='turn-over-energy', val=1E4, vmin=a_t, vmax=b_t,
                              unit='lorentz-factor',log=log_values)
        n_e_super_exp.add_par('p', par_type='LE_spectral_slope', val=2.0, vmin=-10., vmax=10, unit='')
        n_e_super_exp.add_par('a', par_type='spectral_curvature', val=1.0, vmin=0., vmax=100., unit='')
        n_e_super_exp.set_distr_func(distr_func_super_exp)

        return n_e_super_exp

    @staticmethod
    def _create_lp(gamma_grid_size,log_values,normalize,skip_build,emitters_type):
        @nb.njit(fastmath=True, cache=True)
        def distr_func_lp(gamma, gamma0_log_parab, r, s):
            return np.power((gamma / gamma0_log_parab), (-s - r * np.log10(gamma / gamma0_log_parab)))

        n_lp = EmittersDistribution(name='lp',
                                    spectral_type='lp',
                                    normalize=normalize,
                                    emitters_type=emitters_type,
                                    log_values=log_values,
                                    skip_build=skip_build,
                                    gamma_grid_size=gamma_grid_size)

        a_t, b_t = n_lp.set_bounds(1, 1E9, log_val=n_lp._log_values)
        n_lp.add_par('gamma0_log_parab', par_type='turn-over-energy', val=1E4, vmin=a_t, vmax=b_t,
                              unit='lorentz-factor',log=log_values)
        n_lp.add_par('s', par_type='LE_spectral_slope', val=2.0, vmin=-10., vmax=10, unit='')
        n_lp.add_par('r', par_type='spectral_curvature', val=1.0, vmin=-15., vmax=15., unit='')

        n_lp.set_distr_func(distr_func_lp)

        return n_lp

    @staticmethod
    def _create_lpep(gamma_grid_size, log_values, normalize, skip_build, emitters_type):
        @nb.njit(fastmath=True, cache=True)
        def distr_func_lp_ep(gamma, gamma_p, r):
            return np.power(10., (-r * np.power(np.log10(gamma /gamma_p), 2)))

        n_lpep = EmittersDistribution(name='lpep',
                                    spectral_type='lp',
                                    normalize=normalize,
                                    emitters_type=emitters_type,
                                    log_values=log_values,
                                    skip_build=skip_build,
                                    gamma_grid_size=gamma_grid_size)

        a_t, b_t = n_lpep.set_bounds(1, 1E9, log_val=n_lpep._log_values)
        n_lpep.add_par('gamma_p', par_type='turn-over-energy', val=1E4, vmin=a_t, vmax=b_t,
                     unit='lorentz-factor',log=log_values)
        n_lpep.add_par('r', par_type='spectral_curvature', val=1.0, vmin=-15., vmax=15., unit='')
        n_lpep.set_distr_func(distr_func_lp_ep)

        return n_lpep

    @staticmethod
    def _create_lppl(gamma_grid_size, log_values, normalize, skip_build, emitters_type):
        @nb.njit(fastmath=True, cache=True)
        def distr_func_lp(gamma, gamma0_log_parab, r, s):
            f = np.zeros(gamma.shape)
            m = gamma < gamma0_log_parab
            f[m] = np.power(gamma[m]/gamma0_log_parab, -s)
            f[~m] = np.power(gamma[~m]/gamma0_log_parab, (-s -r*np.log10(gamma[~m] / gamma0_log_parab )))
            return f

        n_lppl = EmittersDistribution(name='lppl',
                                    spectral_type='lp',
                                    normalize=normalize,
                                    emitters_type=emitters_type,
                                    log_values=log_values,
                                    skip_build=skip_build,
                                    gamma_grid_size=gamma_grid_size)

        a_t, b_t = n_lppl.set_bounds(1, 1E9, log_val=n_lppl._log_values)

        n_lppl.add_par('gamma0_log_parab', par_type='turn-over-energy', val=1E4, vmin=a_t, vmax=b_t,
                     unit='lorentz-factor',log=log_values)
        n_lppl.add_par('s', par_type='LE_spectral_slope', val=2.0, vmin=-10., vmax=10, unit='')
        n_lppl.add_par('r', par_type='spectral_curvature', val=1.0, vmin=-15., vmax=15., unit='')

        n_lppl.set_distr_func(distr_func_lp)

        return n_lppl

    @staticmethod
    def _create_lppl(gamma_grid_size, log_values, normalize, skip_build, emitters_type):
        @nb.njit(fastmath=True, cache=True)
        def distr_func_lp(gamma, gamma0_log_parab, r, s):
            f = np.zeros(gamma.shape)
            m = gamma < gamma0_log_parab
            f[m] = np.power(gamma[m] / gamma0_log_parab, -s)
            f[~m] = np.power(gamma[~m] / gamma0_log_parab, (-s - r * np.log10(gamma[~m] / gamma0_log_parab)))
            return f

        n_lppl = EmittersDistribution(name='lppl',
                                      spectral_type='lp',
                                      normalize=normalize,
                                      emitters_type=emitters_type,
                                      log_values=log_values,
                                      skip_build=skip_build,
                                      gamma_grid_size=gamma_grid_size)

        a_t, b_t = n_lppl.set_bounds(1, 1E9, log_val=n_lppl._log_values)

        n_lppl.add_par('gamma0_log_parab', par_type='turn-over-energy', val=1E4, vmin=a_t, vmax=b_t,
                       unit='lorentz-factor',log=log_values)
        n_lppl.add_par('s', par_type='LE_spectral_slope', val=2.0, vmin=-10., vmax=10, unit='')
        n_lppl.add_par('r', par_type='spectral_curvature', val=1.0, vmin=-15., vmax=15., unit='')

        n_lppl.set_distr_func(distr_func_lp)

        return n_lppl
