#!/usr/bin/env python

"""
Module: obs_constrain
===================================================================

This module contains all the classes necessary to constrain 
the model parmateres starting from the SED shape




Classes and Inheritance Structure
-------------------------------------------------------------------

.. inheritance-diagram:: BlazarSEDFit.obs_constrain
   


Classes relations
----------------------------------------------

.. figure::  classes_obs_constrain.png
   :align:   center     


  
.. autosummary::
   ObsConstrain
   
    
Module API
-------------------------------------------------------------------

"""


from __future__ import absolute_import, division, print_function

from builtins import (bytes, str, open, super, range,
                      zip, round, input, int, pow, object, map, zip)

__author__ = "Andrea Tramacere"

import math as m

import scipy as sp

from scipy import double,logspace

from numpy import polyfit,polyval


import numpy as n

from .frame_converter import convert_nu_to_blob

from . import sed_models_dic as Model_dic

#import jet_wrapper
#BlazarSED = jet_wrapper.importer()

try:
    from .jetkernel import jetkernel as BlazarSED
except ImportError:
    from .mock import jetkernel as BlazarSED

#from model_manager import FitModel

from .jet_model import Jet

from .sed_shaper import index_array,peak_values

from .output import section_separator,WorkPlace,makedir


__all__=['ObsConstrain','check_boundaries','check_gamma_tansp','check_t_var',
         'constr_B_from_nu_peaks','constr_R_from_CD','find_B_from_nu_p_S','find_gamma0',
         'find_gamma_3p_SSC','find_gamma_Synch','find_HE_cut_off','find_s','find_s1',
         'find_turn_over','get_Comp_factor','get_R_tvar','get_U_Sync_from_Ph',
         'rescale_Ne','set_gmin_from_nu_cut_IR']

class ObsConstrain(object):
    """
    doc
    """
   
    
    
    def __init__(self,B_range=None,
                 distr_e=None,
                 t_var_sec=None,
                 nu_cut_IR=None,
                 beaming=None,
                 theta=None,
                 bulk_factor=None,
                 obj_class=None,
                 restframe=None,
                 z=None,
                 obspar_workplace=None,
                 **keywords):

        if  'SEDShape' not in keywords :  

            self.indices=index_array()
            
            keys = sorted(keywords.keys())
            for kw in keys:
                print ("kw",kw)
                if kw=='indices':
                    self.indices=keywords[kw]
                    
            for index in self.indices.idx_array:
                index.show_val()
            
            self.beta_S=keywords['beta_S']
            self.nu_p_S_obs=keywords['nu_p_S_obs']
            self.nu_p_IC_obs=keywords['nu_p_IC_obs']
            self.nuFnu_p_S_obs=keywords['nuFnu_p_S_obs']
            self.nuFnu_p_IC_obs=keywords['nuFnu_p_IC_obs']
            self.S_nu_max=keywords['S_nu_max']
            
            self.z=z
            self.class_obj=obj_class

            
        else:
            
            self.SEDShape=keywords['SEDShape']
            self.indices=self.SEDShape.indices
           
            self.beta_S = -self.SEDShape.S_peak.curvature
            self.nu_p_S_obs = n.power(10.,self.SEDShape.S_peak.nu_p_val)
            self.nuFnu_p_S_obs = n.power(10.,self.SEDShape.S_peak.nuFnu_p_val)
            self.S_nu_max=n.power(10.,self.SEDShape.S_nu_max)
            self.S_LE_slope=self.SEDShape.S_LE_slope
            self.nu_p_IC_obs  = n.power(10.,self.SEDShape.IC_peak.nu_p_val)
            self.nuFnu_p_IC_obs = n.power(10.,self.SEDShape.IC_peak.nuFnu_p_val)
            
           
            
            self.class_obj=self.SEDShape.obj_class
            self.z=self.SEDShape.sed_data.z
            self.rest_frame='obs'
            
           

        
        
            
        self.beaming=beaming
        self.theta=theta
        self.bulk_factor=bulk_factor

        if self.theta is None and self.bulk_factor is None and self.beaming is not None:
            self.beaming_expr='delta'
        elif self.theta is not None and self.bulk_factor is not None and self.beaming is None:
            self.beaming_expr='bulk_theta'
        else:
            raise RuntimeError('''either you provide the beaming value, or both theta and the bulk factor values ''')
                
        
        self.B_min=B_range[0]
        self.B_max=B_range[1]
        self.B_start = (self.B_min + self.B_max)/2.0

        self.t_var_sec=t_var_sec
        
        self.distr_e=distr_e
        
        self.nu_cut_IR=nu_cut_IR
                        
            
        if obspar_workplace is None:
            obspar_workplace=WorkPlace()
            self.out_dir=obspar_workplace.out_dir
            self.flag=obspar_workplace.flag
        else:
            self.out_dir=obspar_workplace.out_dir
            self.flag=obspar_workplace.flag
    
    def constrain_SSC_model(self,name=None,jet_model=None,params_grid_size=10,electron_distribution_log_values=False):
        """
        constarin SSC model paramters
        """
        
         
        
        print  (section_separator)
        print ("***  constrains parameters from observable ***")
        print
        
        
        model=self.get_model_constraint(name=name,
                                        jet_model=jet_model,
                                        params_grid_size=params_grid_size,
                                        electron_distribution_log_values=electron_distribution_log_values)
         
        
        print
        print  (section_separator)

        return model
    
    
    
    
    def constrain_SSC_EC_model(self,name=None,jet_model=None,EC_componets_list=['BLR'],params_grid_size=10):
        """
        constarin SSC model paramters
        """
        
         
        
        print  (section_separator)
        print  ("***  constrains parameters from observable ***")
        print
        
        
        model=self.get_model_constraint(name=name,jet_model=jet_model,EC_componets_list=EC_componets_list,params_grid_size=params_grid_size)
        
         

        print
        print  (section_separator)

        return model
     
     
     
    
    def get_model_constraint(self,name=None,
                             jet_model=None,
                             EC_componets_list=None,
                             params_grid_size=10,
                             electron_distribution_log_values=False):
        
        if name is None:
            name=self.distr_e
            
        out_dir='%s/obs_constrain_%s/'%(self.out_dir,name)
        makedir(out_dir)
        
        if jet_model is None:
            if self.beaming_expr=='delta':
                jet_model=Jet(name=name, electron_distribution=self.distr_e,electron_distribution_log_values=electron_distribution_log_values)
            elif self.beaming_expr=='bulk_theta': 
                jet_model=Jet(name=name, electron_distribution=self.distr_e,beaming_expr='bulk_theta',electron_distribution_log_values=electron_distribution_log_values)
            
            else:
                raise RuntimeError('''wrong beaming_expr value=%s, allowed 'delta' or 'bulk_theta' '''%self.beaming_expr)
        
        else:
            self.distr_e=jet_model.get_electron_distribution_name()

    
        nu_p_EC_seed_field=None
        if EC_componets_list is not None:
            jet_model.add_EC_component(EC_componets_list)
            jet_model.set_par('L_disk',val=self.SEDShape.L_Disk)
            jet_model.set_par('T_disk_max',val=self.SEDShape.T_Disk)
            nu_p_EC_seed_field=self.SEDShape.nu_p_Disk
     
        #setting the Jet object 
        path_initial=jet_model.get_path()
        flag_initial=jet_model.get_flag()
        
        jet_model.set_path(out_dir)
        jet_model.set_flag(name)
        
        
        IC_nu_size_initial=jet_model.get_IC_nu_size()
        nu_seed_size_initial=jet_model.get_seed_nu_size()
        
        jet_model.set_IC_nu_size(50)
        jet_model.set_seed_nu_size(50)

        jet_model.show_pars()
        
        #beaming
        print ("---> ***  emitting region parameters  ***")
       
       
        #SEtting emetting_region parameters
        
        #!!Controlla beaming max
        #beaming
        
        if self.beaming_expr=='delta':
            beaming_par=jet_model.get_par_by_type('beaming')
            beaming_par.set(val=self.beaming)
            print ("--->",beaming_par.get_description())
            
           
        elif self.beaming_expr=='bulk_theta': 
            bulk_factor_par=jet_model.get_par_by_type('jet-bulk-factor')
            theta_par=jet_model.get_par_by_type('jet-viewing-angle')
            theta_par.set(val=self.theta)
            bulk_factor_par.set(val=self.bulk_factor)
            print ("--->",theta_par.get_description())
            print ("--->",bulk_factor_par.get_description())
            self.beaming=jet_model.get_beaming(theta=self.theta,bulk_factor=self.bulk_factor,beaming=-1.0)
            print ("---> beaming set to",self.beaming)
        else:
            raise RuntimeError('''wrong beaming_expr value=%s, allowed 'delta' or 'bulk_theta' '''%self.beaming_expr)
        
            
        #z
        z_par=jet_model.get_par_by_type('redshift')
        if z_par is not None:
            print ("---> setting par type redshift, corresponding to par %s"%(z_par.name))
            z_par.set(val=self.z)
            print ("---> ",z_par.get_description())
            print
        #B
        B_par=jet_model.get_par_by_type('magnetic_field')
        if B_par is not None:
            print ("---> setting par type magnetic_field, corresponding to par %s"%(B_par.name))
            B_par.set(val=self.B_start)
            print ("---> ",B_par.get_description())
            print


        #R
        R_tvar=get_R_tvar(self.beaming,self.t_var_sec,+self.z)
        R_par=jet_model.get_par_by_type('region_size')
        if R_par is not None:
            print ("---> setting par type region_size, corresponding to par %s"%(R_par.name))

            R_par.set(val=set_lin_log_val(R_par,R_tvar))

            print ("---> ",R_par.get_description())
            print
        
        
        
        print
        
        print ("---> *** electron distribution parameters ***")
        
        #elec distr law
        print ("---> distribution type: ",jet_model.get_electron_distribution_name())
        
       
        
        #r
        curvature_par=jet_model.get_par_by_type('spectral_curvature')
        if curvature_par is not None:
        
            if self.beta_S is not None:
                curvature_par.set(val=self.beta_S*5.0)
                print ("---> r elec. spec. curvature =%e"%curvature_par.val)
            else:
                curvature_par.set(val=0.4*5.0)
                print ("---> beta Sync not provided, using 0.4, corresponding to r=2.0")
 
            print ("---> setting par type curvature, corresponding to par %s"%(curvature_par.name))
            print ("---> ",curvature_par.get_description())
            print
        
        #s
        LE_spectral_slope_par=jet_model.get_par_by_type('LE_spectral_slope')
        if LE_spectral_slope_par is not None:
            s_Plank,s_X,s_Fermi,index_s=find_s(self.class_obj,self.nu_p_S_obs,self.S_LE_slope,self.indices) 
            LE_spectral_slope_par.set(val=index_s)
            print ("---> setting par type LE_spectral_slope, corresponding to par %s"%(LE_spectral_slope_par.name))
            print ("---> ",LE_spectral_slope_par.get_description())
            print


        #s1
        HE_spectral_slope_par=jet_model.get_par_by_type('HE_spectral_slope')
        if HE_spectral_slope_par is not None:
            index_s1=find_s1(self.class_obj,self.indices) 
            HE_spectral_slope_par.set(val=index_s1)
            print ("---> setting par type LE_spectral_slope, corresponding to par %s"%(HE_spectral_slope_par.name))
            print ("---> ",HE_spectral_slope_par.get_description())
            print()

        gamma_3p_Sync=find_gamma_Synch(self.nu_p_S_obs,self.rest_frame,B_par.val,self.beaming,z_par.val)
        print ("---> gamma_3p_Sync= %e, assuming B=%e"%(gamma_3p_Sync,B_par.val))

        
        
        #gmax start
        gmax=find_HE_cut_off(self.distr_e,self.S_nu_max,self.rest_frame,B_par.val,self.beaming,z_par.val)
        print ("---> gamma_max=%e from nu_max_Sync= %e, using B=%e"%(gmax,self.S_nu_max,B_par.val))
        HE_cut_off_par=jet_model.get_par_by_type('high-energy-cut-off')
        if HE_cut_off_par is not None:
            HE_cut_off_par.set(val=set_lin_log_val(HE_cut_off_par,gmax))

            print ("---> setting par type high-energy-cut-off, corresponding to par %s"%(HE_cut_off_par.name))
            print ("---> ",HE_cut_off_par.get_description())
            print ()
        
        #gmin start
        gmin= set_gmin_from_nu_cut_IR(self.nu_cut_IR,self.rest_frame,B_par.val,self.beaming,z_par.val)  
        par_LE_cut_off_par=jet_model.get_par_by_type('low-energy-cut-off')
        if par_LE_cut_off_par is not None:
            par_LE_cut_off_par.set(val=set_lin_log_val(par_LE_cut_off_par,gmin))

            print ("---> setting par type low-energy-cut-off, corresponding to par %s"%(par_LE_cut_off_par.name))
            print ("---> ",par_LE_cut_off_par.get_description())
            print ()
        
        
        #turn-over-energy from gamma_3p_Sync
        turn_over_par=jet_model.get_par_by_type('turn-over-energy')
        if turn_over_par is not None:
            _t=find_turn_over(jet_model,self.distr_e,gamma_3p_Sync)
            turn_over_par.set(val=set_lin_log_val(turn_over_par,_t))

            print ("---> setting par type turn-over energy, corresponding to par %s"%(turn_over_par.name))
            print ("---> using gamma_3p_Sync=",gamma_3p_Sync)
            print ("---> ",turn_over_par.get_description())
            print ()
         
   
        
        #turn-over-energy from gamma_3p_SSC
        gamma_3p_SSC= find_gamma_3p_SSC(self.nu_p_S_obs,self.nu_p_IC_obs,self.rest_frame,gamma_3p_Sync,self.beaming,z_par.val,nu_p_EC_seed_field=nu_p_EC_seed_field)
        print ("---> gamma_3p_SSCc= %e",gamma_3p_SSC)
       
        if turn_over_par is not None:
            _t=find_turn_over(jet_model,self.distr_e,gamma_3p_SSC)
            turn_over_par.set(val=set_lin_log_val(turn_over_par,_t))

            print ("---> setting par type turn-over energy, corresponding to par %s"%(turn_over_par.name))
            print ("---> using gamma_3p_SSC=",gamma_3p_SSC)
            print ("---> ",turn_over_par.get_description())
            print ()
        
        print ()
       
       
       
        #N
        N_par=jet_model.get_par_by_type('electron_density')
        if N_par is not None:
            N,ratio=rescale_Ne(jet_model,self.nuFnu_p_S_obs,self.rest_frame)
            print ("---> setting par type electron_density, corresponding to par %s"%(N_par.name))
            N_par.set(val=N)
            print('--->',N_par.get_description())
        #find B
        #print "estimate B from nu_p_S, and gamma_3p_SSC"
        B=find_B_from_nu_p_S(self.nu_p_S_obs,gamma_3p_SSC,self.rest_frame,self.beaming,z_par.val)
        print ("---> B from nu_p_S=%e"%B)
       

        print ("---> get B from best matching of nu_p_IC")
        if B_par is not None:
            if EC_componets_list is None:
                B_from_nu_peaks,failed=constr_B_from_nu_peaks (jet_model,self.nu_p_S_obs,self.nu_p_IC_obs,self.rest_frame,self.B_min,self.B_max,self.beaming,params_grid_size)
            else:
                B_from_nu_peaks,failed=constr_B_from_nu_peaks (jet_model,self.nu_p_S_obs,self.nu_p_IC_obs,self.rest_frame,self.B_min,self.B_max,self.beaming,params_grid_size,EC=True)

            B_par.set(val=B_from_nu_peaks)
            print ("---> setting par type magnetic_field, corresponding to par %s"%(B_par.name))
            print ("---> ",B_par.get_description())
            print ()

            if failed==False:
                print ("---> best B found: ",B_par.get_description())
            else:
                print ("---> constrain failed, B set to: ",B_par.get_description())
                print ()
                
    
                
            print ()
            print ("---> update pars for new B ")
            
            #update gmin from new B
            gmin= set_gmin_from_nu_cut_IR(self.nu_cut_IR,self.rest_frame,B_par.val,self.beaming,z_par.val)  
            if par_LE_cut_off_par is not None:
                par_LE_cut_off_par.set(val=set_lin_log_val(par_LE_cut_off_par,gmin))

                print ("---> setting par type low-energy-cut-off, corresponding to par %s"%(par_LE_cut_off_par.name))
                print ("---> ",par_LE_cut_off_par.get_description())
                print ()
         
            #update gamma_3p and gamma_cut from new B
            gamma_3p_Sync=find_gamma_Synch(self.nu_p_S_obs,self.rest_frame,B_par.val,self.beaming,z_par.val)
            if turn_over_par is not None:
                _t=find_turn_over(jet_model,self.distr_e,gamma_3p_Sync)
                turn_over_par.set(val=set_lin_log_val(turn_over_par,_t))
                print ("---> setting par type low-energy-cut-off, corresponding to par %s"%(turn_over_par.name))
                print ("---> using gamma_3p_Sync=",gamma_3p_Sync)
                print ("---> ",turn_over_par.get_description())
                print
                
            
            #update gmax for New B
            gmax=find_HE_cut_off(self.distr_e,self.S_nu_max,self.rest_frame,B_par.val,self.beaming,z_par.val)
            print ("---> gamma_max=%e from nu_max_Sync= %e, using B=%e"%(gmax,self.S_nu_max,B_par.val))
            if HE_cut_off_par is not None:
                HE_cut_off_par.set(val=set_lin_log_val(HE_cut_off_par,gmax))

                print ("---> setting par type high-energy-cut-off, corresponding to par %s"%(HE_cut_off_par.name))
                print ("---> ",HE_cut_off_par.get_description())
                print ()
            
            #update N for New B
            if N_par is not None:
                N,ratio=rescale_Ne(jet_model,self.nuFnu_p_S_obs,self.rest_frame)
                print ("---> setting par type electron_density, corresponding to par %s"%(N_par.name))
                N_par.set(val=N)
                
            
                
        
        
        #Improve R according to CD
        print ("---> get R from Compoton Dominance (CD)")
        R_start=R_tvar
        if EC_componets_list  is None:
            R_from_CD,failed=constr_R_from_CD(jet_model,self.nuFnu_p_S_obs,self.nuFnu_p_IC_obs,self.nu_p_IC_obs,self.rest_frame,R_tvar,params_grid_size)
        else:
            R_from_CD,failed=constr_R_from_CD(jet_model,self.nuFnu_p_S_obs,self.nuFnu_p_IC_obs,self.nu_p_IC_obs,self.rest_frame,R_tvar,params_grid_size,EC=True)

        if failed==False:
            

            if R_par is not None:
                print ("---> setting par type region_size, corresponding to par %s"%(R_par.name))
                R_par.set(val=set_lin_log_val(R_par,R_from_CD))
                print ("---> ",R_par.get_description())
                print ()


                if N_par is not None:
                    N,ratio=rescale_Ne(jet_model,self.nuFnu_p_S_obs,self.rest_frame)
                    print ("---> setting par type electron_density, corresponding to par %s"%(N_par.name))
                    N_par.set(val=N)
        else:
            R_par.set(val=set_lin_log_val(R_par,R_start))
            print ("---> constrain failed, R unchanged: ")
            print ("---> ",R_par.get_description())
            print ()
            
         
        #Check tau_gamma_gamma and t_var
        #check_gamma_tansp(jet_model,self.beaming,sp.logspace(21,28,10),self.rest_frame)
        print ("---> t_var (days)", check_t_var(R_par.val_lin,self.beaming,z_par.val)/86400.)
    
                
        jet_model.set_flag('obs_constrain_final')
       
        print ()
        
        print ("show pars")
        jet_model.show_pars()
        
        #sets to initial values
        
        #if self.nu_p_IC_obs is not None:
        jet_model.set_IC_nu_size(IC_nu_size_initial)
        jet_model.set_seed_nu_size(nu_seed_size_initial)
        
        
        print ("eval_model")
        jet_model.eval(fill_SED=True)
        
        jet_model.set_flag(flag_initial)
        jet_model.set_path(path_initial)
        
        
#        if self.SEDShape.L_host!=None:
#           
#            SED_model=FitModel(jet=jet_model,template=self.SEDShape.host)
#       
#        else:
#        
#            SED_model=FitModel(jet=jet_model,template=self.SEDShape.host)

            
        return jet_model
        
    





def check_t_var(R,beaming,z):
    return R*(1+z)/(BlazarSED.vluce_cm*beaming)




def check_gamma_tansp(jet,beaming_val,nu_IC_data,rest_frame):
    """retrun tau_gamma_gamma for a given  IC frequency 
    
    Args:
        blob:
        temap_ev:
        nu_IC_data: data freq of the IC component to check fro tau_gamma_gamma
        
    Returns:
        tau_gamma_gamma
        
    """
    #################################
    # This function must not        #
    # change jet.blob attributes    #
    # all the changed values must   #
    # be set back to their original #
    # values                        #
    #################################
    
    flag_initial=jet.get_flag()
    
    jet.set_flag('tau_gamma_gamma')

    IC_initial=jet.get_IC_mode()

    jet.set_IC_mode('off')
    
    
    #nu_IC_blob=convert_nu_to_blob(nu_IC_data,rest_frame,blob.beam_obj,blob.z_cosm)
    #target_nu_blob = SED.MEC2 * SED.MEC2 / (SED.HPLANCK * SED.HPLANCK * nu_IC_blob)
    target_nu_obs=BlazarSED.MEC2 * BlazarSED.MEC2 / (BlazarSED.HPLANCK * BlazarSED.HPLANCK * nu_IC_data)
    #print nu_IC_blob,target_nu,target_nu*nu_IC_blob

   
    jet.eval()
    
    
    nu_obs=sp.array([])
    nuFnu_obs=sp.array([])
    nu_obs,nuFnu_obs=jet.get_SED_points(log_log=True,name='Sync')
   
    #for i in range(BlazarSED.GetNuIntMaxSynch()):
    #    
    #     if BlazarSED.GetSEDSynch(i)>0:
    #         nu_obs=sp.append(nu_obs,sp.log10(BlazarSED.GetNuObsSynch(i)))
    #        nuFnu_obs = sp.append(nuFnu_obs, sp.log10(BlazarSED.GetSEDSynch(i)))

     
    nuFnu_target = sp.power(10.0,sp.interp(sp.log10(target_nu_obs),nu_obs,nuFnu_obs))

    if sp.shape(target_nu_obs)==():
        target_nu_obs=sp.array([target_nu_obs])
        
        nuFnu_target=sp.array([nuFnu_target])
        
        nu_IC_data=sp.array([nu_IC_data ])
        
    for i in range(len(target_nu_obs)):
        #print nuFnu_target,blob.beam_obj,blob.z_cosm,blob.dist
        
        
        
        z_val=jet.get_par_by_type('redshift').val
        
        DL_cm_val=jet.get_DL_cm()
        
        nuLnu_target_blob= BlazarSED.nuFnu_obs_to_nuLnu_blob(nuFnu_target[i],beaming_val,z_val,DL_cm_val)
        
        R_val=jet.get_par_by_type('region_size').val_lin
        
        tau_gamma_gamma=BlazarSED.SIGTH/5*nuLnu_target_blob/(4*sp.pi * BlazarSED.vluce_cm * R_val * BlazarSED.MEC2 )
        
        print ("     -->target comoving freq=%e"%target_nu_obs[i], "observed freq=%e"%nu_IC_data[i], "tau_gamma=%e"%(tau_gamma_gamma))
    
    
    jet.set_IC_mode(IC_initial)
    jet.set_flag(flag_initial)



def set_gmin_from_nu_cut_IR(nu_cut_IR,rest_frame,B,beaming,z):
        if nu_cut_IR is not None:
            return find_gamma_Synch (nu_cut_IR,rest_frame,B,beaming,z)
            
        else:
            print ("--> !! No nu_cut_IR provided set gmin to 1")
            return 1
           
        
        
        

def find_turn_over(jet,distr_e,gamma_3p):
    
    if distr_e=='lppl' or distr_e=='lp':
        r=jet.get_par_by_type('spectral_curvature').val
        s=jet.get_par_by_type('LE_spectral_slope').val
        return find_gamma0(r,s,gamma_3p)
        

    elif distr_e=='lpep':
        r=jet.get_par_by_type('spectral_curvature').val
        return gamma_3p/10**(1.5/r)

    #!! TO FIX
    elif distr_e=='plc':
        return gamma_3p*2

    
    elif distr_e=='bkn':
        "!!!!!!!!! MUST FIX ACCORDING P1 and P2!!!!!!"
        return gamma_3p



def find_HE_cut_off(distr_e,nu_S_max,rest_frame,B,beaming,z):
    gamma_max_Sync=find_gamma_Synch(nu_S_max,rest_frame,B ,beaming,z)
    
    if distr_e=='lppl' or distr_e=='lp' or distr_e=='lpep':
        return gamma_max_Sync
    
    elif  distr_e=='bkn' or distr_e=='plc' or distr_e=='pl':
        return gamma_max_Sync/1.5
    

def find_gamma0(r,s,gamma_3p):
    """returns the value of gamma_0  for
    a log_par+pl distribution
    
    Args:
        r: curvature
        s: spectral index in the PL branch
        gamma_3p: peak value requested for n(gamma)gamma^3

    Returns:
        
        
    """
    if (r!=0.0):
        c=(3-s)/(2*r)
        gamma_0= gamma_3p/(pow(10,c))
    else:
        gamma_0= gamma_3p

    if gamma_0<1.0:
        gamma_0=1.0
    
    return gamma_0




def find_B_from_nu_p_S(nu_p_S,gamma_3p,rest_frame,beaming,z):
    """returns B according to Ep_S and gamma_3p


        Args:
            nu_p_S:  Synchrotron peack frequency
            gamma_3p: peak value of n(gamma)gamma^3
            

            re_eval: def==True, set the flag to re_evaluate find_gamma_3p_Synch, after updating B

        Returns:
            B
    """
    nu_p_blob=convert_nu_to_blob(nu_p_S,rest_frame,beaming,z)

    #first tryal for B
    B= (nu_p_blob)/(gamma_3p*gamma_3p*3.7E6 )
    
    return B    
    

def  find_gamma_Synch (nu_S,rest_frame,B,beaming,z):
    """returns the value of gamma corresponding to the  Synch freq
    
        
        Args:
            B: magnetic field
            nu_S:  Synchrotron  frequency
            rest_frame: rest frame cooresponding to nu_p_S
            z: redshift
            beamign: beaming factor
            
        Returns:
            gamma
            
    """
    nu_S_blob=convert_nu_to_blob(nu_S,rest_frame,beaming,z)
     
    return m.sqrt(nu_S_blob/(3.7E6*B))



def find_gamma_3p_SSC(nu_p_S,nu_p_IC,rest_frame,gamma_3p_Sync,beaming,z,nu_p_EC_seed_field=None):
    """returns the value of gamma_3p from nu_p_S/nu_p_IC


        Args:
            nu_p_S :  Synchrotron peack frequency
            nu_p_IC:  IC peack frequency
            rest_frame: rest frame cooresponding to peak frequencies

                    
        Returns:
            gamma_3p_SSC
    """
    #print "*** find_gamma_3p_SSC ***"
    #print "nu_p_S",nu_p_S
    #print "nu+p_IC",nu_p_IC
    #print "beaming",beaming
    #print "z",z
    
    

    
    #print "--> gamma_TH, gamma_KN_Ep_C",g_TH,gamma_KN
    if nu_p_EC_seed_field is None:
        g_TH=(4./3)*m.sqrt(nu_p_IC/nu_p_S)
        nu_p_seed_blob=convert_nu_to_blob(nu_p_S,rest_frame,beaming,z)
       
    else:
        nu_p_seed_blob=nu_p_EC_seed_field/(1+z)*beaming
        g_TH=(4./3)*m.sqrt(nu_p_IC/nu_p_seed_blob)
    #print"--> comp_fac,gamma_TH", comp_factor,g_TH,
    comp_factor=get_Comp_factor(gamma_3p_Sync,nu_p_seed_blob)
    print ("nu_p_seed_blob",nu_p_seed_blob)

    gamma_3p_SSC_TH=g_TH
    #print "--> gamma_3p_SSC_TH",gamma_3p_SSC_TH,g_TH
    print ("COMP FACTOR",comp_factor,gamma_3p_SSC_TH)
    
    if comp_factor>0.01:
        #print "--> correction", 0.2*pow(comp_factor,-0.45)
        gamma_3p_SSC=g_TH/(0.2*pow(comp_factor,-0.45))
    else:
        gamma_3p_SSC=gamma_3p_SSC_TH
    
    #   print "--> gamma_3p_SSC=",gamma_3p_SSC
    return gamma_3p_SSC



def get_Comp_factor(gamma,nu_p_S_blob):
    """returns the compton factor = nu_blob_seed_Synch*hplanck/(mec2)
        
        Args:
            gamma: gamma of the up-scattering electrons
            u_p_S_blob: Synchrotron peack frequency in the blob rest frame
           
        Returns:
         compton factor
    """
    # print "--> Ep_S_blob ",nu_p_S_blob
    return gamma*nu_p_S_blob*BlazarSED.HPLANCK/(BlazarSED.MEC2)





def find_s(class_obj,nu_p_S_obs,S_LE_slope,indices):
    """Find the index of the low-energy PL branch of n(gamma), from PL fit over various instrument bands

    Args:
        class_obj: object  class
        indices_array

    Returns:
        s_Planck, s_X, s_Fermi, s

    """
    
    #1 get S index from nu_p_S_obs
    #2 get IC index from nu_p_IC_obs

    
    if (class_obj=='LSP' or class_obj=='ISP' or class_obj=='HSP'):
        
        s_radio_mm=None
        if indices.get_by_name('radio_mm')is not None:
            if indices.get_by_name('radio_mm').val is not None:
                radio_mm_index=indices.get_by_name('radio_mm').val.spectral
                s_radio_mm=(1.0-2*radio_mm_index)
                print ("---> s_radio_mm",radio_mm_index,s_radio_mm)

                
        
        s_X=None
        if indices.get_by_name('X')is not None:
            if indices.get_by_name('X').val is not None:
                X_index=indices.get_by_name('X').val.spectral
                s_X=(1.0-2*X_index)
                print ("---> s_X",s_X)
    


      
        s_Fermi=None
        if indices.get_by_name('Fermi') is not None:
            if indices.get_by_name('Fermi').val is not None:
                Fermi_index=indices.get_by_name('Fermi').val.spectral
                s_Fermi=(0.4+ Fermi_index*-1*1.60)
                print ("---> s_Fermi",s_Fermi)
    
        
        s_UV_X=None
        if indices.get_by_name('UV_X') is not None:
            if indices.get_by_name('UV_X').val is not None:
                UV_X_index=indices.get_by_name('UV_X').val.spectral
                s_UV_X=(1+ UV_X_index*-1*2.0)
                print ("---> s_UV_X",s_UV_X)
        
        s_Opt_UV=None
        if indices.get_by_name('Opt_UV') is not None:
            if indices.get_by_name('Opt_UV').val is not None:
                Opt_UV_index=indices.get_by_name('Opt_UV').val.spectral
                s_Opt_UV=(1 -1.0*Opt_UV_index*2.0)
                print ("---> s_Opt_UV",Opt_UV_index,s_Opt_UV)
        
        if class_obj=='LSP' :
            #print s_X
            if s_X is not None and s_X-3.0<0:
                #if X_index+1.0>0.2:
                s=s_X
                print ("---> s from X_index",s)
            elif s_radio_mm is not None and s_radio_mm+2>0.0:
                s=s_radio_mm
                print ("---> s from radio_mm_index")
            else:   
                s=2.0
        
        if class_obj=='ISP':
            if s_radio_mm is not None and s_radio_mm+2>0.2:
                s=s_radio_mm
                print ("---> s from radio_mm_index")
            else:   
                s=2.0
        
#        if class_obj=='ISP':
#            if s_Planck!=None:
#                s=s_Planck
#            elif s_X!=None:
#                s=s_X
#            else:
#                s=2.0   
        
        s_SED_shape_photon=S_LE_slope-2.0
        s_SED_shape=-2.0*s_SED_shape_photon-1.0
        print ("---> s from synch log-log fit",s_SED_shape)

        if class_obj=='HSP':
           
            s_UV=None
            
            if s_Opt_UV is None and s_UV_X is not None:
                s_UV=s_UV_X
            elif s_Opt_UV is not None and s_UV_X is None:
                s_UV=s_Opt_UV
            elif s_Opt_UV is not None and s_UV_X is not None:
                s_UV=s_UV_X
                
            if s_Fermi is not None and s_UV is not None:
                s=(s_Fermi +s_UV )/(2.0)
                print ("---> s from (s_Fermi + s_UV)/2")
            elif s_UV is not None :
                s=s_UV
                print ("---> s from  s_UV")

            elif s_Fermi is not None :
                s=s_Fermi
                print ("---> s from Fermi")

            
            elif s_radio_mm is not None :
                s=s_radio_mm
                print ("---> s from radio_mm")

            else:
                s=2.0 
   
    print ("---> power-law index s, class obj=%s s chosen is %f"%(class_obj,s))
    #print s_mm, s_X, s_Fermi, s
    return s_radio_mm, s_X, s_Fermi, s



def find_s1(class_obj,indices):
    #print "---> !!! fake function, still to develop"
    val=3.5
    print ("---> set s1 to %f"%val)
    return val

def rescale_Ne(jet,Lp_S,rest_frame):
    """Rescales N.blob to get the wanted Lp_S_blob


        Args:
            blob: SED module class
            temp_ev: SED module class
            Lp_S_obs_blob: peak luminosity of the Synch component 
            rest_frame: rest_frame of the observed data
            


        Returns:
            Lp_S_blob/Lp_S_model_blob
    """
    #################################
    # This function must not        #
    # change jet.blob attributes    #
    # all the changed values must   #
    # be set back to their original #
    # values                        #
    #################################

    #blob.STEM='constr_rescale_Ne-pre'
    flag_initial=jet.get_flag()
    
    N_initial=jet.get_par_by_name('N').val
    #print('-->3',N_initial)
    jet.set_flag('constr_R')
    
    
    IC_initial=jet.get_IC_mode()

    jet.set_IC_mode('off')
    
    jet.set_flag('constr_rescale_Ne')

    jet.eval()
    
    Lp_S_model=jet.get_SED_peak( Model_dic.Sync_nuFnu_p_dic[rest_frame])

    #print('-->3', Lp_S_model)
    N_new=N_initial*Lp_S/Lp_S_model

    jet.set_par('N',val=N_new)

    
    
    jet.eval()
    
    #reset old values
    jet.set_IC_mode(IC_initial)
    jet.set_par('N',val=N_initial)
    jet.set_flag(flag_initial)
    
    Lp_S_model=jet.get_SED_peak( Model_dic.Sync_nuFnu_p_dic[rest_frame])
    
    return N_new, Lp_S/Lp_S_model



def get_U_Sync_from_Ph(jet,re_eval_sync=False):
    """returns U_synch by integrating Synch photon spectrum

    Args:
        blob: SED module class
        temp_ev: SED module class

        re_eval_sync: def==False, re_eval Synch spectrum and then get U_synch
    Returns:
        U_synch

    """
    if re_eval_sync==True:
        comp_old=jet.blob.do_SSC
        jet.blob.do_SSC=0
        jet.eval()
        jet.blob.do_SSC=comp_old
    
    return BlazarSED.Uph_Sync(jet.blob)


def check_boundaries(val,val_min,val_max,val_name):
    failed=False
    if val<val_min or val>val_max:
        failed=True
        print ("---> %s=%e, out of boundaries %e %e, rejected"%(val_name,val,val_min,val_max))
        if val<val_min:
            val=val_min
        elif val>val_max:
            val=val_max
        print ("     Best %s not found, (temporary set to %e)"%(val_name,val))
    else:        
        print ("     Best %s=%e"%(val_name,val))

    return val,failed

def get_R_tvar(beaming,t_var_sec,z):
    return BlazarSED.vluce_cm*beaming*t_var_sec/(1+z)

def set_lin_log_val(p,v):
    if p.islog is True:
        v = m.log10(v)

    return v

def constr_R_from_CD(jet,nuFnu_p_S,nuFnu_p_IC,nu_p_IC,rest_frame,R_tvar,params_grid_size,EC=False):
    #################################
    # This function must not        #
    # change jet.blob attributes    #
    # all the changed values must   #
    # be set back to their original #
    # values                        #
    #################################
    
    flag_initial=jet.get_flag()

    
    R_initial=jet.get_par_by_name('R').val_lin
    N_initial=jet.get_par_by_name('N').val

    jet.set_flag('constr_R')
    
    CD_model_log=[]
    CD_obs=nuFnu_p_IC/nuFnu_p_S
    R_min=R_tvar/1000
    R_min=1E13
    R_max=R_tvar*2
    R_grid=logspace(n.log10(R_min),n.log10(R_max),params_grid_size)
    f=open('%s/R_vs_CD.dat'%jet.get_path(),'w')
    for R in R_grid:

        jet.set_par('R',val=set_lin_log_val(jet.get_par_by_name('R'),R))
        #print('--> 1',jet.get_par_by_name('R').val)

        N_res,ratio=rescale_Ne(jet,nuFnu_p_S,rest_frame)
        
        jet.set_par('N',val=N_res)
        
        jet.eval()

        Lp_S=jet.get_SED_peak(Model_dic.Sync_nuFnu_p_dic[rest_frame])
        if EC==False:
            Lp_IC=jet.get_SED_peak(Model_dic.SSC_nuFnu_p_dic[rest_frame])
        else:
            nu_p_IC,Lp_IC=jet.get_SED_peak(freq_range=[0.01*nu_p_IC,10*nu_p_IC])

        #print('--> 2', N_res, ratio, Lp_S,Lp_IC,nuFnu_p_S)

        CD=Lp_IC/Lp_S
        CD_model_log.append(n.log10(CD))

        #print "     R=%e, CD_obs=%e, CD_model=%e"%(R,CD_obs,CD)
        #print(n.log10(R),n.log10(CD/CD_obs),file=f)

    f.close()

    R_grid_log=n.log10(R_grid)
    p=polyfit(CD_model_log,R_grid_log,2)
    #print "--> lll",CD_model_log,R_grid_log
    best_R=polyval(p,n.log10(CD_obs))
    
    Best_R=n.power(10.,best_R)
    Best_R,failed=check_boundaries(Best_R,R_min,R_max,'R')
    
    jet.set_flag(flag_initial)

    jet.set_par('R',val=set_lin_log_val(jet.get_par_by_name('R'),R_initial))
    jet.set_par('N',val=N_initial)
    
    return  Best_R,failed



def constr_B_from_nu_peaks(jet,nu_p_S,nu_p_IC,rest_frame,B_min,B_max,beaming,params_grid_size,EC=False):
    #################################
    # This function must not        #
    # change jet.blob attributes    #
    # all the changed values must   #
    # be set back to their original #
    # values                        #
    #################################

    flag_initial=jet.get_flag()
    
    B_initial=jet.get_par_by_name('B').val
    
    jet.set_flag('constr_B')
 
    nu_p_IC_model_log=[]
    
    B_grid=logspace(sp.log10(B_min),sp.log10(B_max),params_grid_size)
    
    f=open('%s/B_vs_nu_p_IC.dat'%jet.get_path(),'w')
    
    #fin the turn-over variable to update with B
    turn_over_energy=jet.get_par_by_type('turn-over-energy')
    
    if turn_over_energy is not None :
        turn_over_energy_initial=turn_over_energy.val_lin
    
    
    z=jet.get_par_by_type('redshift').val
    
    for B in B_grid:
        
        jet.set_par('B',val=B)
        
        
        gamma_3p=find_gamma_Synch(nu_p_S,rest_frame,B,beaming,z)

        #upda-turn-over variable
        if turn_over_energy is not None :
            turn_over_energy.set(val=set_lin_log_val(turn_over_energy,gamma_3p))

        
        jet.eval()   
        
        if EC==False:
            nu_p_model=jet.get_SED_peak(Model_dic.SSC_nu_p_dic[rest_frame])
        else:
             nu_p_model,nuFnu_p_model=jet.get_SED_peak(freq_range=[0.01*nu_p_IC,10*nu_p_IC])
        
        nu_p_IC_model_log.append(n.log10(nu_p_model))

        print(n.log10(B),n.log10(nu_p_model),file=f)
        #print n.log10(B),n.log10(nu_p_model)

    f.close()
    
    B_grid_log=n.log10(B_grid)
    
     
    p=polyfit(nu_p_IC_model_log,B_grid_log,2)    
    
    Best_B=polyval(p,n.log10(nu_p_IC))
    Best_B = n.power(10., Best_B)
    Best_B,failed = check_boundaries(Best_B, B_min, B_max, 'B')
    
    
    #REST CHANGED VALUES
    jet.set_flag(flag_initial)
    jet.set_par('B',val=B_initial)

    if turn_over_energy is not None :
        turn_over_energy.set(val=set_lin_log_val(turn_over_energy,turn_over_energy_initial))

        
    return  Best_B,failed
