"""Internal absorption modeling utilities for jet spectral components."""

import os
import ctypes

on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if on_rtd:
    from .mock import jetkernel as BlazarSED
else:
    from .jetkernel import jetkernel as BlazarSED


import numpy as np


def _get_c_array_read_only(ptr, size):
    size = int(size)
    if size <= 0 or int(ptr) == 0:
        return np.zeros(0, dtype=np.float64)
    arr = (ctypes.c_double * size).from_address(int(ptr))
    return np.ctypeslib.as_array(arr).copy()


class InternalAbsorption(object):


    """Compute internal gamma-gamma absorption from jet seed photon fields.

    Notes
    -----
    Evaluates optical depth ``tau(nu)`` using BLR or DT photon distributions
    through the C backend and provides attenuation factors for integration into
    jet spectral component calculations.
    """
    def __init__(self,
                 jet,
                 nu_min=None,
                 seed_photons_name='BLR',
                 N_soft=50,
                 N_hard=50,
                 N_R_H=20,
                 N_theta=20,
                 use_R_H_profile_extrapolation=False,
                 ):
        """Create a new `InternalAbsorption` instance.
        
        Parameters
        ----------
        jet : object
            Jet model instance.
        nu_min : object, optional
            Minimum frequency in Hz.
        seed_photons_name : str, optional
            Seed-photon field identifier (for example ``BLR`` or ``DT``).
        N_soft : int, optional
            Number of soft-photon energy samples.
        N_hard : int, optional
            Number of hard-photon energy samples.
        N_R_H : int, optional
            Number of distance samples along ``R_H``.
        N_theta : int, optional
            Number of angular samples.
        use_R_H_profile_extrapolation : bool, optional
            If ``True``, enable r h profile extrapolation.
        """
        if seed_photons_name in ['BLR','DT']:
            self._seed_photons_name=seed_photons_name
        else:
            raise RuntimeError('seed_photons_name %s not valid'%seed_photons_name)

        self._N_soft=N_soft
        self._N_hard=N_hard

        self._N_R_H=N_R_H
        self._N_theta=N_theta
        self._nu_min=nu_min
        self._jet=jet
        self._use_R_H_profile_extrapolation=use_R_H_profile_extrapolation

    def eval_tau_photons(self,
                         nu_src,
                         R_H,
                         peak=False,
                         use_R_H_profile_extrapolation=False):

        """Evaluate tau photons.
        
        Parameters
        ----------
        nu_src : object
            Source-frame frequency array in Hz.
        R_H : object
            Distance from black hole in cm.
        peak : bool, optional
            If ``True``, use peak-optimized seed-photon sampling.
        use_R_H_profile_extrapolation : bool, optional
            If ``True``, enable r h profile extrapolation.

        Notes
        -----
        Internal absorption is evaluated via
        ``eval_internal_abs_tau_isolated`` only.
        
        Returns
        -------
        object
            Computed value.
        """
        nu_src = np.atleast_1d(np.asarray(nu_src, dtype=np.float64))
        if nu_src.size == 0:
            return np.zeros(0, dtype=np.float64), np.zeros(0, dtype=np.float64)

        nu_src_max = float(np.max(nu_src))
        r_h_override = -1.0 if R_H is None else float(R_H)

        nu_min = -1.0 if self._nu_min is None else float(self._nu_min)

        if not hasattr(BlazarSED, 'eval_internal_abs_tau_isolated'):
            raise RuntimeError('jetkernel extension is missing eval_internal_abs_tau_isolated; rebuild the C extension.')

        tau_size = BlazarSED.eval_internal_abs_tau_isolated(
            self._jet._blob,
            self._seed_photons_name,
            nu_min,
            int(self._N_soft),
            int(self._N_hard),
            int(self._N_R_H),
            int(self._N_theta),
            int(bool(use_R_H_profile_extrapolation)),
            int(bool(peak)),
            nu_src_max,
            r_h_override,
        )

        if tau_size < 0:
            raise RuntimeError('internal absorption C evaluation failed for component %s' % self._seed_photons_name)

        if self._seed_photons_name == "BLR":
            comp = self._jet._blob.core.internal_abs.BLR
        else:
            comp = self._jet._blob.core.internal_abs.DT

        nu_tau = _get_c_array_read_only(comp.nu_tau, comp.tau_size)
        tau = _get_c_array_read_only(comp.tau, comp.tau_size)

        return tau,nu_tau

    def _get_default_nu_src(self):
        sed_nu_src = None
        if hasattr(self._jet, 'spectral_components'):
            try:
                sed_nu_src = self._jet.spectral_components.Sum.SED.nu_src
            except Exception:
                sed_nu_src = None

        if sed_nu_src is not None:
            try:
                sed_nu_src = np.atleast_1d(np.asarray(sed_nu_src.value, dtype=np.float64))
            except Exception:
                sed_nu_src = None

        if sed_nu_src is not None and sed_nu_src.size > 0 and np.all(np.isfinite(sed_nu_src)):
            return sed_nu_src

        nu_src, _ = self._jet._prepare_nu_model(nu=None, loglog=False)
        nu_src = np.atleast_1d(np.asarray(nu_src, dtype=np.float64))
        if nu_src.size == 0 or np.any(~np.isfinite(nu_src)):
            raise RuntimeError('invalid jet frequency grid for internal absorption evaluation')
        return nu_src


    def eval(self, get_tau=False, lin_nu=None, peak=False):
        """Evaluate model output.
        
        Parameters
        ----------
        get_tau : bool, optional
            If ``True``, return optical depth instead of attenuation.
        lin_nu : object, optional
            Linear-frequency array in Hz.
        peak : bool, optional
            If ``True``, use peak-optimized seed-photon sampling.
        
        Returns
        -------
        object
            Computed value.
        """

        if lin_nu is None:
            nu_src = self._get_default_nu_src()
        else:
            nu_src = np.atleast_1d(np.asarray(lin_nu, dtype=np.float64))
        
       
        tau,nu_tau=self.eval_tau_photons(nu_src,
                                        R_H=self._jet.parameters.R_H.val,
                                        peak=peak,
                                        use_R_H_profile_extrapolation=self._use_R_H_profile_extrapolation)
        EPS = 1e-300   
        nu_src_pos = np.maximum(nu_src, EPS)
        nu_tau_pos = np.maximum(nu_tau, EPS)
        tau_pos    = np.maximum(tau,    EPS)

        tau_interp = 10**np.interp(
            np.log10(nu_src_pos),
            np.log10(nu_tau_pos),
            np.log10(tau_pos),
            left=np.log10(EPS),
            right=None
        )
        
        if get_tau:
            return tau_interp,nu_src_pos

        return np.exp(-tau_interp),nu_src_pos
        
