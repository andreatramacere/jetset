from .base_model import MultiplicativeModel
from jetset.jetkernel.jetkernel import HPLANCK as h
from jetset.jetkernel.jetkernel import MEC2 as mec2
from jetset.jetkernel.jetkernel import SIGTH 
from jetset.jetkernel import jetkernel as BlazarSED

import numpy as np


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
            p_old = self._parameters_old.get_par_by_name(p_new.name)
            if p_old is None:
                changed=True
            if p_new.val != p_old.val:
                changed=True
            #print('p_old',p_old.name,p_old.val,'p_new',p_new.name,p_new.val)
        #print('-> testing',changed)
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
    
        self._parameters_old=self._jet_orig.parameters

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
        tau=np.zeros(nu_src.size)
        integrand_nu=np.zeros(n.shape)
        
        for ID_RH,R_H in enumerate(R_H_range):
            self._jet.set_par('R_H',val=R_H)
            nu[ID_RH],n[ID_RH]=self.get_n(R_H=R_H,
                                        seed_photons_name=self._seed_photons_name,
                                        N_soft=N_soft,
                                        peak=peak)
            
            if self._seed_photons_name == "BLR":
                if R_H<self._jet.parameters.R_BLR_in.val:
                    mu_max=1.0
                else:
                    mu_max=self._jet.parameters.R_BLR_in.val/R_H 
            
            elif self._seed_photons_name == "DT":
                if R_H<self._jet.parameters.R_DT.val:
                    mu_max=1.0
                else:
                    mu_max=self._jet.parameters.R_DT.val/R_H
            else:
                raise RuntimeError('seed_photons_name %s not valid'%self._seed_photons_name)

            mu_range[ID_RH]=(np.ones((self._N_theta,N_soft)).T*np.linspace(-mu_max,mu_max,self._N_theta).T).T
            
        for ID,nu_gamma in enumerate(nu_src):
            if nu_gamma<nu_min:
                tau[ID]=0
            else:
                eps_gamma=np.ones(shape)*nu_gamma*h/mec2
                eps_soft=nu*h/mec2
                E_th=mec2**2/(eps_soft*mec2*(1-mu_range))
                s =eps_soft * eps_gamma* (1 - mu_range)/ 2   
                s[eps_gamma*mec2<E_th]=0
                integrand_nu=self.sigma(s)*n*(1+mu_range)
                if peak is False:
                    int_d_nu=np.trapz(integrand_nu,nu,axis=2)
                    int_d_mu=np.trapz(int_d_nu,mu_range[::,::,0],axis=1)
                else:
                    int_d_nu=integrand_nu*nu
                    int_d_mu=np.trapz(int_d_nu[::,::,0],mu_range[::,::,0],axis=1) 
                tau[ID]=np.pi*2*np.trapz(int_d_mu,R_H_range)
        
        return tau



    def sigma(self,s=np.array([[[]]])):
        a = 3 / 8 * SIGTH /(s*s*s)
        b = (s**2+s-0.5)
        s_sqr=np.sqrt(s)
        s_sqr_1=np.sqrt(s-1)
        c=np.log((s_sqr+s_sqr_1)/(s_sqr-s_sqr_1))
        d=-np.sqrt((s*(s-1)))*(s+1)
        
        e = a * (b*c+d)
        e[s<1] = 0

        return e

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
        x=np.zeros(size)
        y=np.zeros(size)
        
    
        for i in range(size):
            x[i]=BlazarSED.get_spectral_array(nu_ptr,self._jet._blob,i)
            y[i]=BlazarSED.get_spectral_array(n_ptr,self._jet._blob,i)
        msk=np.logical_and(x>=nu_start,y<=nu_stop)
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
            x=x[y>y.max()/100]
            y=y[y>y.max()/100]
            x_int=np.logspace(np.log10(x[0]),np.log10(x[-1]),N_soft)
            y_int = np.interp(np.log10(x_int),np.log10(x), np.log10(y),left=1E-200,right=1E-200)
            y_int = np.power(10., y_int)
            x=x_int
            y=y_int
    
        return x,y



    def eval(self, get_tau=False,skip_check=True,lin_nu=None):
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
                                  )
        self._old_tau=tau
        tau_interp = np.interp(nu_src, nu_abs, tau,left=0,right=0)
       
     
        
        if get_tau:

            return tau_interp,nu_src
        