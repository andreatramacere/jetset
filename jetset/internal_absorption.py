from .base_model import MultiplicativeModel
from jetset.jetkernel.jetkernel import HPLANCK as h
from jetset.jetkernel.jetkernel import MEC2 as mec2
from jetset.jetkernel.jetkernel import SIGTH 
from jetset.jetkernel import jetkernel as BlazarSED
from .jet_kernel_tools import get_spectral_c_array_read_only
from numba import njit, prange
import numpy as np
import warnings

H_OVER_MEC2 = h / mec2
SIGMA_PREF = 0.75 * SIGTH * 0.5


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
                 nu_min,
                 seed_photons_name='BLR',
                 N_soft=50,
                 N_hard=50,
                 N_R_H=20,
                 N_theta=20,
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
        self._parameters_old=None
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


    def eval_tau_photons(self,
                        R_H,
                        nu_src,
                        nu_min,
                        skip_check=True,
                        peak=False):

        for p_orig in self._jet_orig.parameters.par_array:
            if not p_orig.frozen and not p_orig._is_dependent:
                p=self._jet.get_par_by_name(p_orig.name)
                p.val=p_orig.val
        
        if  skip_check:
            tau_changed=True
        else:
            tau_changed=self._check_eval_needed(ptype=self._seed_photons_name)
    
        self._update_parameters_old()

        if not tau_changed:
            return self._old_tau
    
        R_H_range=np.logspace(0,5,self._N_R_H)*R_H
        if peak is True:
            N_soft=1
        else:
            N_soft=self._N_soft

        
        shape=(self._N_R_H,self._N_theta,N_soft)
        nu=np.zeros(shape)
        n=np.zeros(shape)
        mu_range=np.zeros(shape)
        for ID_RH,R_H in enumerate(R_H_range):
            self._jet.set_par('R_H',val=R_H)
            nu[ID_RH],n[ID_RH]=self.get_n(R_H=R_H,
                                        seed_photons_name=self._seed_photons_name,
                                        N_soft=N_soft,
                                        peak=peak)
            
            if self._seed_photons_name == "BLR":
                if R_H<self._jet.parameters.R_BLR_in.val:
                    mu_max=1.0
                    mu_min=-1
                else:
                    mu_max=1.0
                    mu_min = np.sqrt(1.0 - (self._jet.parameters.R_BLR_out.val / R_H)**2)
            
            elif self._seed_photons_name == "DT":
                if R_H<self._jet.parameters.R_DT.val:
                    mu_max=1.0
                    mu_min=-1
                else:
                    mu_max=1.0
                    mu_min = np.sqrt(1.0 - (self._jet.parameters.R_DT.val / R_H)**2)
            else:
                raise RuntimeError('seed_photons_name %s not valid'%self._seed_photons_name)
     
            mu_range[ID_RH]=(np.ones((self._N_theta,N_soft)).T*np.linspace(mu_min,mu_max,self._N_theta).T).T

        nu_src = np.atleast_1d(nu_src)
        if nu_min is not None:
            nu_src = nu_src[nu_src >= nu_min]

        if nu_src.size == 0:
            return np.zeros(0, dtype=np.float64)

        nu_grid = np.ascontiguousarray(nu, dtype=np.float64)
        n_grid = np.ascontiguousarray(n, dtype=np.float64)
        mu_grid = np.ascontiguousarray(mu_range[:, :, 0], dtype=np.float64)
        R_H_grid = np.ascontiguousarray(R_H_range, dtype=np.float64)
        nu_src_grid = np.ascontiguousarray(nu_src, dtype=np.float64)

        try:
            tau = _compute_tau_numba(
                nu_src_grid,
                nu_grid,
                n_grid,
                mu_grid,
                R_H_grid,
                H_OVER_MEC2,
                SIGMA_PREF,
            )
        except Exception as exc:
            warnings.warn(
                f"Falling back to NumPy tau computation because the optimized implementation failed: {exc}",
                RuntimeWarning,
                stacklevel=2,
            )
            one_minus_mu = np.clip(1.0 - mu_range, 1e-20, None)
            eps_soft = nu_grid * H_OVER_MEC2
            eps_gamma = (nu_src_grid * H_OVER_MEC2)[:, None, None, None]
            s = eps_gamma * eps_soft[None, :, :, :] * one_minus_mu[None, :, :, :] / 2.0
            sigma_vals = self.sigma(s)
            integrand = sigma_vals * n_grid[None, :, :, :] * one_minus_mu[None, :, :, :]
            int_over_nu = np.trapz(integrand, nu_grid[None, :, :, :], axis=-1)
            int_over_mu = np.trapz(int_over_nu, mu_range[None, :, :, 0], axis=-1)
            tau = 2.0 * np.pi * np.trapz(int_over_mu, R_H_grid, axis=-1)

        return tau



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
            n_name='n_%s_DRF'%self._seed_photons_name
            nu_name='nu_%s_disk_RF'%self._seed_photons_name
            BlazarSED.Build_I_nu_BLR(self._jet._blob)
            nu_start=self._jet._blob.nu_start_BLR_disk_RF
            nu_stop=self._jet._blob.nu_stop_BLR_disk_RF
        elif seed_photons_name== "DT":
            n_name='n_%s_DRF'%self._seed_photons_name
            nu_name='nu_%s_disk_RF'%self._seed_photons_name
            BlazarSED.Build_I_nu_DT(self._jet._blob)
            nu_start=self._jet._blob.nu_start_DT_DRF
            nu_stop=self._jet._blob.nu_stop_DT_DRF
        else:
            raise RuntimeError('seed_photons_name %s not valid'%seed_photons_name)
        
        n_ptr = getattr(self._jet._blob, n_name)
        nu_ptr = getattr(self._jet._blob, nu_name)
        
        size=self._jet._blob.nu_grid_size
        #x=np.zeros(size)
        #y=np.zeros(size)
        
    
        #for i in range(size):
        #x=BlazarSED.get_spectral_array_np(nu_ptr,self._jet._blob)
        #y=BlazarSED.get_spectral_array_np(n_ptr,self._jet._blob)
        x,y=get_spectral_c_array_read_only(nu_ptr,n_ptr,size)
        msk=np.logical_and(x>=nu_start,x<=nu_stop)
        x=x[msk]
        y=y[msk]
        if peak is True:
            id=np.argmax(y)
            scale_factor=np.trapz(y,x)
            x=np.atleast_1d(x[id])
            y=np.atleast_1d(y[id])
            if rescale is True:
                scale_factor=y*x/scale_factor
                y=y/scale_factor
            
            
        else:
            x=x[y>y.max()/100000]
            y=y[y>y.max()/100000]
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
        nu_abs=np.logspace(np.log10(self._nu_min),np.log10(nu_src.max()),self._N_hard)
        tau=np.zeros(nu_src.shape)
        tau=self.eval_tau_photons(nu_src=nu_abs,
                                  nu_min=self._nu_min,
                                  R_H=self._jet_orig.parameters.R_H.val,
                                  skip_check=skip_check,
                                  peak=peak,
                                  )
        self._old_tau=tau
        tau_interp = np.interp(nu_src, nu_abs, tau,left=0,right=0)
       
     
        
        if get_tau:

            return tau_interp,nu_src
        
