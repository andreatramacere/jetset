import os

# Prefer the portable numba backend for this kernel path to avoid OpenMP runtime
# deprecation chatter (e.g. omp_set_nested) in environments using Intel OMP.
os.environ.setdefault("NUMBA_THREADING_LAYER", "workqueue")

import os

on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if on_rtd:
    from .mock import jetkernel as BlazarSED
else:
    from .jetkernel.jetkernel import HPLANCK as h
    from .jetkernel.jetkernel import MEC2 as mec2
    from .jetkernel.jetkernel import SIGTH 
    from .jetkernel import jetkernel as BlazarSED
    H_OVER_MEC2 = h / mec2
    SIGMA_PREF = 0.75 * SIGTH * 0.5




from .jet_kernel_tools import get_spectral_c_array_read_only
from .utils import get_nested_attr
from numba import njit, prange
import numpy as np
import warnings



@njit(fastmath=True)
def _sigma_numba(s, pref):
    beta = np.sqrt(1.0 - 1.0 / s)
    term = (3.0 - beta ** 4) * np.log((1.0 + beta) / (1.0 - beta)) - 2.0 * beta * (2.0 - beta ** 2)
    return pref * (1.0 - beta ** 2) * term


@njit(parallel=True, fastmath=True)
def _compute_tau_numba(nu_src,
                       nu_grid,
                       n_grid,
                       mu_grid,
                       R_H_grid,
                       h_over_mec2,
                       sigma_pref):

    n_gamma = nu_src.shape[0]
    n_rh, n_theta, n_soft = nu_grid.shape
    two_pi = 2.0 * np.pi
    tau = np.zeros(n_gamma)

    for gamma_idx in prange(n_gamma):
        eps_gamma = nu_src[gamma_idx] * h_over_mec2
        tau_gamma = 0.0
        prev_rh_integral = 0.0
        for rh_idx in range(n_rh):
            mu_vals = mu_grid[rh_idx]
            mu_integral = 0.0
            prev_mu_integral = 0.0
            for theta_idx in range(n_theta):
                mu_val = mu_vals[theta_idx]
                one_minus_mu = 1.0 - mu_val
                if one_minus_mu < 1e-20:
                    one_minus_mu = 1e-20

                nu_integral = 0.0
                prev_nu_val = nu_grid[rh_idx, theta_idx, 0]
                prev_integrand = 0.0
                for soft_idx in range(n_soft):
                    nu_val = nu_grid[rh_idx, theta_idx, soft_idx]
                    eps_soft = nu_val * h_over_mec2
                    s_val = eps_gamma * eps_soft * one_minus_mu * 0.5

                    integrand = 0.0
                    if s_val >= 1.0:
                        integrand = _sigma_numba(s_val, sigma_pref) * n_grid[rh_idx, theta_idx, soft_idx] * one_minus_mu

                    if soft_idx > 0:
                        dnu = nu_val - prev_nu_val
                        nu_integral += 0.5 * (prev_integrand + integrand) * dnu

                    prev_nu_val = nu_val
                    prev_integrand = integrand

                if theta_idx > 0:
                    dmu = mu_val - mu_vals[theta_idx - 1]
                    mu_integral += 0.5 * (prev_mu_integral + nu_integral) * dmu

                prev_mu_integral = nu_integral

            if rh_idx > 0:
                d_rh = R_H_grid[rh_idx] - R_H_grid[rh_idx - 1]
                tau_gamma += 0.5 * (prev_rh_integral + mu_integral) * d_rh

            prev_rh_integral = mu_integral

        tau[gamma_idx] = two_pi * tau_gamma

    return tau


class InternalAbsorption(object):


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
        if seed_photons_name in ['BLR','DT']:
            self._seed_photons_name=seed_photons_name
        else:
            raise RuntimeError('seed_photons_name %s not valid'%seed_photons_name)

        self._N_soft=N_soft
        self._N_hard=N_hard

        self._N_R_H=N_R_H
        self._N_theta=N_theta
        self._nu_min=nu_min
        jet.skip_internal_absorption_serial=True
        self._jet=jet.clone()
        jet.skip_internal_absorption_serial=False
        self._jet_orig=jet
        self._old_tau=None
        self._old_nu_tau=None
        self._parameters_old=None
        self._use_R_H_profile_extrapolation=use_R_H_profile_extrapolation
        
        #self._DT_pars=['R_H','tau_DT','R_DT','T_DT']
    
    
    def _update_parameters_old(self):
        self._parameters_old={}
        for p_orig in self._jet_orig.parameters.par_array:
            if not p_orig.frozen and not p_orig._is_dependent:
                self._parameters_old[p_orig.name]=p_orig.val

    def _check_eval_needed(self,ptype):
        changed=False
        if self._parameters_old is None:
            return True
        if self._old_tau is None:
            return True
        pars_new=self._jet_orig.parameters.get_pars_by_type(ptype)
        
        pars_new.extend(self._jet_orig.parameters.get_pars_by_type('Disk'))
        pars_new.extend([self._jet_orig.parameters.get_par_by_name('R_H')])
        for p_new in pars_new:
            if p_new.name in self._parameters_old.keys():   
                if p_new.val != self._parameters_old[p_new.name]:
                    changed=True
            else:
                changed=True

        return changed


    def _get_R_seed(self):
        if self._seed_photons_name not in ["BLR","DT"]:
             raise RuntimeError('seed_photons_name %s not valid'%self._seed_photons_name)

        if self._seed_photons_name == "BLR":
            R_seed=self._jet.parameters.R_BLR_out.val
        elif self._seed_photons_name == "DT":
            R_seed=self._jet.parameters.R_DT.val

        return R_seed
    
    def _get_R_seed_field(self,R_seed,N_soft,peak):
        nu_soft_0,n_soft_0=self.get_n(R_H=R_seed/1000,
                        seed_photons_name=self._seed_photons_name,
                        N_soft=N_soft,
                        peak=peak)
        
        return nu_soft_0,n_soft_0


    
    def eval_tau_photons(self,
                         nu_src,
                         R_H,
                         skip_check=True,
                         peak=False,
                         use_R_H_profile_extrapolation=False):

        for p_orig in self._jet_orig.parameters.par_array:
            if not p_orig.frozen and not p_orig._is_dependent:
                p=self._jet.get_par_by_name(p_orig.name)
                p.val=p_orig.val
        
        R_H_range=np.logspace(0,3,self._N_R_H)*R_H
        
        if peak is True:
            N_soft=1
        else:
            N_soft=self._N_soft
        R_seed=self._get_R_seed()

        nu_soft_0,n_soft_0=self._get_R_seed_field(R_seed,N_soft,peak)
        
            
        if self._nu_min is None:
            nu_min=1E40/(nu_soft_0.max())
        else:
            nu_min = self._nu_min

        nu_tau=np.logspace(np.log10(nu_min),np.log10(nu_src.max()),self._N_hard)

        if  skip_check:
            tau_changed=True
        else:
            tau_changed=self._check_eval_needed(ptype=self._seed_photons_name)
    
        self._update_parameters_old()

        if self._old_tau is not None:
            if not tau_changed and np.array_equal(self._old_nu_tau, nu_tau):
                return self._old_tau,nu_tau
    

        shape=(self._N_R_H,self._N_theta,N_soft)
        nu_soft=np.zeros(shape)
        n_soft=np.zeros(shape)
        mu_range=np.zeros(shape)

       
        for ID_RH,R_H in enumerate(R_H_range):
            self._jet.set_par('R_H',val=R_H)
           
            if use_R_H_profile_extrapolation:
                if R_H<=R_seed:
                    R_x=1
                else:
                    R_x=R_seed/R_H

                n_soft[ID_RH]=n_soft_0*(R_x)**2
                nu_soft[ID_RH]=nu_soft_0

            else:
                nu_soft[ID_RH],n_soft[ID_RH]=self.get_n(R_H=R_H,
                                            seed_photons_name=self._seed_photons_name,
                                            N_soft=N_soft,
                                            peak=peak)
            
            if R_H<R_seed:
                mu_max=1.0
                mu_min=-1
            else:
                mu_max=1.0
                mu_min = np.sqrt(1.0 - (R_seed / R_H)**2)
            

            mu_range[ID_RH]=(np.ones((self._N_theta,N_soft)).T*np.linspace(mu_min,mu_max,self._N_theta).T).T

        nu_tau = np.atleast_1d(nu_tau)
        #if nu_min is not None:
        nu_tau = nu_tau[nu_tau >= nu_min]

        if nu_tau.size == 0:
            return np.zeros(0, dtype=np.float64),nu_tau

        nu_soft_grid = np.ascontiguousarray(nu_soft, dtype=np.float64)
        n_soft_grid = np.ascontiguousarray(n_soft, dtype=np.float64)
        mu_grid = np.ascontiguousarray(mu_range[:, :, 0], dtype=np.float64)
        R_H_grid = np.ascontiguousarray(R_H_range, dtype=np.float64)
        nu_tau_grid = np.ascontiguousarray(nu_tau, dtype=np.float64)

        try:
            tau = _compute_tau_numba(
                nu_tau_grid,
                nu_soft_grid,
                n_soft_grid,
                mu_grid,
                R_H_grid,
                H_OVER_MEC2,
                SIGMA_PREF,
            )
        except Exception as exc:
            warnings.warn(
                f"Falling back to NumPy tau computation because the numba implementation failed: {exc}",
                RuntimeWarning,
                stacklevel=2,
            )
            one_minus_mu = np.clip(1.0 - mu_range, 1e-20, None)
            eps_soft = nu_soft_grid * H_OVER_MEC2
            eps_gamma = (nu_tau_grid * H_OVER_MEC2)[:, None, None, None]
            s = eps_gamma * eps_soft[None, :, :, :] * one_minus_mu[None, :, :, :] / 2.0
            sigma_vals = self.sigma(s)
            integrand = sigma_vals * n_soft_grid[None, :, :, :] * one_minus_mu[None, :, :, :]
            int_over_nu = np.trapezoid(integrand, nu_soft_grid[None, :, :, :], axis=-1)
            int_over_mu = np.trapezoid(int_over_nu, mu_range[None, :, :, 0], axis=-1)
            tau = 2.0 * np.pi * np.trapezoid(int_over_mu, R_H_grid, axis=-1)

        return tau,nu_tau



    def sigma(self, s):
        """
        Pair-production cross section [cm^2], s = dimensionless CM energy squared
        s must satisfy s >= 1
        """
        out = np.zeros_like(s)
        mask = s >= 1.0
        if not np.any(mask):
            return out
        sm = s[mask]
        beta = np.sqrt(1.0 - 1.0/sm)
        pref = 0.75 * SIGTH * 0.5  # 3/16 * σ_T = 0.1875 σ_T
        term = (3 - beta**4) * np.log((1 + beta)/(1 - beta)) - 2*beta*(2 - beta**2)
        out[mask] = pref * (1 - beta**2) * term
        return out

    def get_n(self,
            seed_photons_name,
            N_soft,
            R_H,
            peak=False,
            rescale=True):
                    
        self._jet.set_par('R_H',val=R_H)
        BlazarSED.Build_I_nu_Disk(self._jet._blob)
        if seed_photons_name == "BLR":
            n_name='BLR.spec.n_nu_DRF'
            nu_name='BLR.spec.nu_DRF'
            BlazarSED.Build_I_nu_BLR(self._jet._blob)
            nu_start=self._jet._blob.BLR.spec.nu_min_DRF
            nu_stop=self._jet._blob.BLR.spec.nu_max_DRF
        elif seed_photons_name== "DT":
            n_name='DT.spec.n_nu_DRF'
            nu_name='DT.spec.nu_DRF'
            BlazarSED.Build_I_nu_DT(self._jet._blob)
            nu_start=self._jet._blob.DT.spec.nu_min_DRF
            nu_stop=self._jet._blob.DT.spec.nu_max_DRF
        else:
            raise RuntimeError('seed_photons_name %s not valid'%seed_photons_name)
        
        n_ptr = get_nested_attr(self._jet._blob, n_name)
        nu_ptr = get_nested_attr(self._jet._blob, nu_name)
        
        size=self._jet._blob.core.nu_grid_size
       
        x,y=get_spectral_c_array_read_only(nu_ptr,n_ptr,size)
        msk=np.logical_and(x>=nu_start,x<=nu_stop)
        x=x[msk]
        y=y[msk]
        if peak is True:
            idx=np.argmax(y)
            scale_factor=np.trapezoid(y,x)
            x=np.atleast_1d(x[idx])
            y=np.atleast_1d(y[idx])
            if rescale is True:
                scale_factor=y*x/scale_factor
                y=y/scale_factor
            
            
        else:
            x=x[y>y.max()/100]
            y=y[y>y.max()/100]
            x_int=np.logspace(np.log10(x[0]),np.log10(x[-1]),N_soft)
            y_int = np.interp(np.log10(x_int),np.log10(x), np.log10(y),left=1E-200,right=1E-200)
            y_int = np.power(10., y_int)
            x=x_int
            y=y_int
    
        return x,y



    def eval(self, get_tau=False,skip_check=True,lin_nu=None,peak=False):
        """
        """

        if lin_nu is None:
            nu_src=self._jet_orig.spectral_components.Sum.SED.nu_src.value
        else:
            nu_src=lin_nu
            nu_src=np.atleast_1d(nu_src)
        
       
        tau,nu_tau=self.eval_tau_photons(nu_src,
                                        R_H=self._jet_orig.parameters.R_H.val,
                                        skip_check=skip_check,
                                        peak=peak,
                                        use_R_H_profile_extrapolation=self._use_R_H_profile_extrapolation)

        self._old_tau=tau
        self._old_nu_tau=nu_tau
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
        
