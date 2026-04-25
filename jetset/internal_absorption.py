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
                         skip_check=True,
                         peak=False,
                         use_R_H_profile_extrapolation=False,
                         use_isolated_eval=None):

        """Evaluate tau photons.
        
        Parameters
        ----------
        nu_src : object
            Source-frame frequency array in Hz.
        R_H : object
            Distance from black hole in cm.
        skip_check : bool, optional
            If ``True``, skip check.
        peak : bool, optional
            If ``True``, use peak-optimized seed-photon sampling.
        use_R_H_profile_extrapolation : bool, optional
            If ``True``, enable r h profile extrapolation.
        use_isolated_eval : bool or None, optional
            Controls which C backend path is used:
            ``None`` (default) uses isolated when available, otherwise fallback;
            ``True`` forces isolated path; ``False`` forces non-isolated path.
        
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

        # Kept for API compatibility. C-side IA cache/state is authoritative.
        _ = skip_check

        if not hasattr(BlazarSED, 'eval_internal_abs_tau'):
            raise RuntimeError('jetkernel extension is missing eval_internal_abs_tau; rebuild the C extension.')

        nu_min = -1.0 if self._nu_min is None else float(self._nu_min)

        if use_isolated_eval not in (None, True, False):
            raise RuntimeError('use_isolated_eval must be None, True, or False')

        isolated_available = hasattr(BlazarSED, 'eval_internal_abs_tau_isolated')
        if use_isolated_eval is None:
            use_isolated = isolated_available
        else:
            use_isolated = bool(use_isolated_eval)

        if use_isolated:
            if not isolated_available:
                raise RuntimeError('isolated IA evaluation requested but not available; rebuild the C extension.')
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
        else:
            if R_H is not None:
                self._jet.set_par('R_H', val=R_H)
            tau_size = BlazarSED.eval_internal_abs_tau(
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
            )
            # Non-isolated fallback evaluates IA on the live blob and can
            # perturb seed-field buffers; rebuild them at the model R_H.
            if hasattr(BlazarSED, 'spectra_External_Fields'):
                BlazarSED.spectra_External_Fields(1, self._jet._blob, 0)

        if tau_size < 0:
            raise RuntimeError('internal absorption C evaluation failed for component %s' % self._seed_photons_name)

        if self._seed_photons_name == "BLR":
            comp = self._jet._blob.core.internal_abs.BLR
        else:
            comp = self._jet._blob.core.internal_abs.DT

        nu_tau = _get_c_array_read_only(comp.nu_tau, comp.tau_size)
        tau = _get_c_array_read_only(comp.tau, comp.tau_size)

        return tau,nu_tau



    def eval(self, get_tau=False,skip_check=True,lin_nu=None,peak=False,use_isolated_eval=None):
        """Evaluate model output.
        
        Parameters
        ----------
        get_tau : bool, optional
            If ``True``, return optical depth instead of attenuation.
        skip_check : bool, optional
            If ``True``, skip check.
        lin_nu : object, optional
            Linear-frequency array in Hz.
        peak : bool, optional
            If ``True``, use peak-optimized seed-photon sampling.
        use_isolated_eval : bool or None, optional
            Passed to ``eval_tau_photons`` to select isolated/non-isolated C path.
        
        Returns
        -------
        object
            Computed value.
        """

        if lin_nu is None:
            nu_src=self._jet.spectral_components.Sum.SED.nu_src.value
        else:
            nu_src=lin_nu
            nu_src=np.atleast_1d(nu_src)
        
       
        tau,nu_tau=self.eval_tau_photons(nu_src,
                                        R_H=self._jet.parameters.R_H.val,
                                        skip_check=skip_check,
                                        peak=peak,
                                        use_R_H_profile_extrapolation=self._use_R_H_profile_extrapolation,
                                        use_isolated_eval=use_isolated_eval)
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
        
