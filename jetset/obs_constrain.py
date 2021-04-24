
__author__ = "Andrea Tramacere"

import math as m


from numpy import polyfit,polyval


import numpy as np

from .frame_converter import convert_nu_to_blob

from . import jetkernel_models_dic as Model_dic
#
# on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
#
# if on_rtd == True:
#     try:
#         from .jetkernel import jetkernel as BlazarSED
#     except ImportError:
#         from .mock import jetkernel as BlazarSED
# else:

from .jetkernel import jetkernel as BlazarSED

from .jet_model import Jet

from .sed_shaper import index_array,peak_values

from .output import section_separator,WorkPlace,makedir

from .utils import *

__all__=['ObsConstrain','check_boundaries', 'check_t_var',
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
                 #restframe=None,
                 z=None,
                 obspar_workplace=None,
                 **keywords):

        if 'SEDShape' not in keywords :

            self.indices=index_array()
            
            keys = sorted(keywords.keys())
            for kw in keys:
                print("kw",kw)
                if kw =='indices':
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
            self.nu_p_S_obs = np.power(10., self.SEDShape.S_peak.nu_p_val)
            self.nuFnu_p_S_obs = np.power(10., self.SEDShape.S_peak.nuFnu_p_val)
            self.S_nu_max=np.power(10., self.SEDShape.S_nu_max)
            self.S_LE_slope=self.SEDShape.S_LE_slope
            self.nu_p_IC_obs  = np.power(10., self.SEDShape.IC_peak.nu_p_val)
            self.nuFnu_p_IC_obs = np.power(10., self.SEDShape.IC_peak.nuFnu_p_val)

            self.class_obj = self.SEDShape.obj_class
            self.z=self.SEDShape.sed_data.z
            self.rest_frame='obs'
            check_frame(self.rest_frame)

        if self.z<0.:
            raise  RuntimeError('redshift value must be>0, please check z in sed_data or z in class constructor parameter')
        
        
            
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

    def constrain_SSC_model(self,name=None,jet_model=None,params_grid_size=10,electron_distribution_log_values=False,silent=False):

        """
        constarin SSC model paramters
        """
        print(section_separator)
        print("***  constrains parameters from observable ***")
        print()

        model=self.get_model_constraint(name=name,
                                        jet_model=jet_model,
                                        params_grid_size=params_grid_size,
                                        electron_distribution_log_values=electron_distribution_log_values,
                                        silent=silent)
        print()
        print(section_separator)

        return model
    
    
    
    
    def constrain_SSC_EC_model(self,name=None,
                               jet_model=None,
                               EC_componets_list=['EC_BLR'],
                               params_grid_size=10,
                               electron_distribution_log_values=False,
                               R_H=None,
                               silent=False,
                               R_H_within_BLR=False,
                               R_H_within_DT=False,
                               disk_type='BB'):
        """
        constarin SSC model paramters
        """
        
         
        
        print(section_separator)
        print ("***  constrains parameters from observable ***")
        print()
        
        
        model=self.get_model_constraint(name=name,
                                        jet_model=jet_model,
                                        EC_componets_list=EC_componets_list,
                                        params_grid_size=params_grid_size,
                                        electron_distribution_log_values=electron_distribution_log_values,
                                        silent=silent,
                                        R_H=R_H,
                                        R_H_within_BLR=R_H_within_BLR,
                                        R_H_within_DT=R_H_within_DT,
                                        disk_type=disk_type)
        
         

        print()
        print (section_separator)

        return model
     
     
     
    
    def get_model_constraint(self,
                             name=None,
                             jet_model=None,
                             EC_componets_list=None,
                             params_grid_size=10,
                             electron_distribution_log_values=False,
                             silent=False,
                             R_H=None,
                             R_H_within_BLR=False,
                             R_H_within_DT=False,
                             disk_type='BB'):
        


        #out_dir='%s/obs_constrain_%s/'%(self.out_dir,name)
        #makedir(out_dir)

        if silent is False:
            print(section_separator)
            print("---> ***  emitting region parameters  ***")

        if jet_model is None:
            if self.beaming_expr=='delta':
                jet_model=Jet(name=name, electron_distribution=self.distr_e,electron_distribution_log_values=electron_distribution_log_values)
            elif self.beaming_expr=='bulk_theta':
                jet_model=Jet(name=name, electron_distribution=self.distr_e,beaming_expr='bulk_theta',electron_distribution_log_values=electron_distribution_log_values)

            else:
                raise RuntimeError('''wrong beaming_expr value=%s, allowed 'delta' or 'bulk_theta' '''%self.beaming_expr)
        if jet_model.emitters_distribution.spectral_type not in jet_model.emitters_distribution.spectral_types_obs_constrain():
            raise RuntimeError('''to use osb constrain the spectral type of emitters has to be ''' %jet_model.emitters_distribution.spectral_types_obs_constrain())

        if R_H is not None:
            jet_model.set_par('R_H',val=R_H)
        self.emitters_distr_spectral_type=jet_model.emitters_distribution.spectral_type

        nu_p_EC_seed_field=None
        if EC_componets_list is not None:
            jet_model.add_EC_component(EC_componets_list,disk_type=disk_type)
            if hasattr(jet_model.parameters, 'L_Disk'):
                jet_model.set_par('L_Disk',val=self.SEDShape.L_Disk)
            if hasattr(jet_model.parameters, 'T_Disk'):
                jet_model.set_par('T_Disk', val=self.SEDShape.T_Disk)
                if silent is False:
                    print('---> EC set L_D ',jet_model.parameters.L_Disk.val)
                    print('---> EC set T_D', jet_model.parameters.T_Disk.val)
                if hasattr(jet_model.parameters, 'R_BLR_in'):
                    R_BLR_in=1E17*np.sqrt(self.SEDShape.L_Disk/1E45)
                    jet_model.set_par('R_BLR_in', val=R_BLR_in)
                    jet_model.set_par('R_BLR_out', val=R_BLR_in*2)

                    if silent is False:
                        print('---> EC set R_BLR_in', jet_model.parameters.R_BLR_in.val)
                        print('---> EC set R_BLR_out', jet_model.parameters.R_BLR_out.val)
                    if jet_model.parameters.R_H.val > R_BLR_in and R_H_within_BLR is True:
                        jet_model.parameters.R_H.val = R_BLR_in * 0.8
                        if silent is False:
                            print('---> moved R_H within BLR to R_H=%e'%jet_model.parameters.R_H.val)
                if hasattr(jet_model.parameters, 'R_DT'):
                    R_DT=2.5E18*np.sqrt(self.SEDShape.L_Disk/1E45)
                    jet_model.set_par('R_DT', val = R_DT)
                    if silent is False:
                        print('---> EC set R_DT', jet_model.parameters.R_DT.val)
                    if jet_model.parameters.R_H.val > R_DT and R_H_within_DT is True:
                        jet_model.parameters.R_H.val = R_DT * 0.8
                        if silent is False:
                            print('---> moved R_H within DT to R_H=%e'%jet_model.parameters.R_H.val)


            nu_p_EC_seed_field=self.SEDShape.nu_p_Disk

        if silent is False:
            print()
        #setting the Jet object 
        #path_initial=jet_model.get_path()
        #flag_initial=jet_model.get_flag()
        
        #jet_model.set_path(out_dir)
        #jet_model.set_flag(name)

        IC_nu_size_initial=jet_model.get_IC_nu_size()
        nu_seed_size_initial=jet_model.get_seed_nu_size()


        jet_model.set_IC_nu_size(50)
        jet_model.set_seed_nu_size(50)

        #jet_model.show_pars()
        
        #beaming

       
       
        #SEtting emetting_region parameters
        
        #!!Controlla beaming max
        #beaming

        if self.beaming_expr=='delta':
            beaming_par=jet_model.get_par_by_type('beaming')
            beaming_par.set(val=self.beaming)
            #if silent is False:
            #    print("--->",beaming_par.get_description())
            
           
        elif self.beaming_expr=='bulk_theta': 
            bulk_factor_par=jet_model.get_par_by_type('jet-bulk-factor')
            theta_par=jet_model.get_par_by_type('jet-viewing-angle')
            theta_par.set(val=self.theta)
            bulk_factor_par.set(val=self.bulk_factor)
            #if silent is False:
            #    print("--->",theta_par.get_description())
            #    print("--->",bulk_factor_par.get_description())
            #self.beaming=jet_model.get_beaming(theta=self.theta,bulk_factor=self.bulk_factor,beaming=-1.0)
            self.beaming = jet_model.get_beaming()
            if silent is False:
                print("---> setting beaming  to",self.beaming)
                print()
        else:
            raise RuntimeError('''wrong beaming_expr value=%s, allowed 'delta' or 'bulk_theta' '''%self.beaming_expr)
        
            
        #z
        z_par=jet_model.get_par_by_type('redshift')
        if z_par is not None:
            if silent is False:
                print ("---> setting par type redshift, corresponding to par %s"%(z_par.name))
                print()
            z_par.set(val=self.z)
            #if silent is False:
            #    print("---> ",z_par.get_description())
            #    print()
        #B
        B_par=jet_model.get_par_by_type('magnetic_field')
        B_par.set(val=self.B_start)
        if B_par is not None:
            if silent is False:
                print("---> setting par type magnetic_field, corresponding to par %s=%e"%(B_par.name,B_par.val))
                print()



        #R
        R_tvar,completed=get_R_tvar(self.beaming,self.t_var_sec,+self.z)

        R_par=jet_model.get_par_by_type('region_size')
        R_par.set(val=set_lin_log_val(R_par, R_tvar))
        if R_par is not None and completed is True:
            if silent is False:

                print("---> setting par type region_size, corresponding to par %s=%e"%(R_par.name,R_par.val))
                print("---> completed", completed)
                print()



        if silent is False:
            print()
            print("---> *** electron distribution parameters ***")
            print('---> emitters distribution spectral type', jet_model.emitters_distribution.spectral_type)
            print('---> emitters distribution name', jet_model.emitters_distribution.name)
            print()

        #r
        curvature_par=jet_model.get_par_by_type('spectral_curvature')
        if curvature_par is not None:
        
            if self.beta_S is not None:
                curvature_par.set(val=self.beta_S*5.0)
                if silent is False:
                    print("---> r elec. spec. curvature =%e"%curvature_par.val)
            else:
                curvature_par.set(val=0.4*5.0)
                if silent is False:
                    print("---> beta Sync not provided, using 0.4, corresponding to r=2.0")
            if silent is False:
                print("---> setting par type curvature, corresponding to par %s"%(curvature_par.name))
                #print("---> ",curvature_par.get_description())
                print()
        
        #s
        LE_spectral_slope_par=jet_model.get_par_by_type('LE_spectral_slope')
        if LE_spectral_slope_par is not None:
            (s_Plank,s_X,s_Fermi,index_s),completed=find_s(self.class_obj,self.nu_p_S_obs,self.S_LE_slope,self.indices,silent=silent)
            LE_spectral_slope_par.set(val=index_s)
            if silent is False:
                print("---> setting par type LE_spectral_slope, corresponding to par %s"%(LE_spectral_slope_par.name))
                print("---> task completed",completed)
                #print("---> ",LE_spectral_slope_par.get_description())
                print()


        #s1
        HE_spectral_slope_par=jet_model.get_par_by_type('HE_spectral_slope')
        if HE_spectral_slope_par is not None:
            index_s1,completed=find_s1(self.class_obj,self.indices)
            HE_spectral_slope_par.set(val=index_s1)
            if silent is False:
                print("---> setting par type LE_spectral_slope, corresponding to par %s"%(HE_spectral_slope_par.name))
                print("---> task completed", completed)
                #print("---> ",HE_spectral_slope_par.get_description())
                print()

        gamma_3p_Sync,completed=find_gamma_Synch(self.nu_p_S_obs,self.rest_frame,B_par.val,self.beaming,z_par.val)
        if silent is False:
            print ("---> setting gamma_3p_Sync= %e, assuming B=%e"%(gamma_3p_Sync,B_par.val))
            print ("---> task completed",completed)
            print()
        
        #gmax start
        gmax, completed=find_HE_cut_off(self.emitters_distr_spectral_type,self.S_nu_max,self.rest_frame,B_par.val,self.beaming,z_par.val)
        if silent is False:
            print("---> gamma_max=%e from nu_max_Sync= %e, using B=%e"%(gmax,self.S_nu_max,B_par.val))
            print("---> task completed", completed)
        HE_cut_off_par=jet_model.get_par_by_type('high-energy-cut-off')
        if HE_cut_off_par is not None:
            HE_cut_off_par.set(val=set_lin_log_val(HE_cut_off_par,gmax))
            if silent is False:
                print("---> setting par type high-energy-cut-off, corresponding to par %s=%e"%(HE_cut_off_par.name,HE_cut_off_par.val))
                #print("---> ",HE_cut_off_par.get_description())
                print()
        
        #gmin start
        gmin,completed= set_gmin_from_nu_cut_IR(self.nu_cut_IR,self.rest_frame,B_par.val,self.beaming,z_par.val)
        par_LE_cut_off_par=jet_model.get_par_by_type('low-energy-cut-off')
        if par_LE_cut_off_par is not None:
            par_LE_cut_off_par.set(val=set_lin_log_val(par_LE_cut_off_par,gmin))
            if silent is False:
                print("---> setting par type low-energy-cut-off, corresponding to par %s=%e"%(par_LE_cut_off_par.name,par_LE_cut_off_par.val))
                print("---> task completed", completed)
                #print("---> ",par_LE_cut_off_par.get_description())
                print()

        #turn-over-energy from gamma_3p_Sync
        turn_over_par=jet_model.get_par_by_type('turn-over-energy')
        if turn_over_par is not None:
            _t,completed=find_turn_over(jet_model,self.emitters_distr_spectral_type,gamma_3p_Sync)
            turn_over_par.set(val=set_lin_log_val(turn_over_par,_t))

            if silent is False:
                print("---> setting par type turn-over energy, corresponding to par %s=%e"%(turn_over_par.name,turn_over_par.val))
                print("---> task completed", completed)
                print("---> using gamma_3p_Sync=",gamma_3p_Sync)
                #print("---> ",turn_over_par.get_description())
                print()

        #turn-over-energy from gamma_3p_SSC
        gamma_3p_SSC,completed= find_gamma_3p_SSC(self.nu_p_S_obs,self.nu_p_IC_obs,self.rest_frame,gamma_3p_Sync,self.beaming,z_par.val,nu_p_EC_seed_field=nu_p_EC_seed_field,silent=silent)
        if silent is False:
            print("---> determine gamma_3p_SSCc= %e"%gamma_3p_SSC)
            print("---> task completed", completed)
            print()
       
        if turn_over_par is not None:
            _t,completed=find_turn_over(jet_model,self.emitters_distr_spectral_type,gamma_3p_SSC)
            turn_over_par.set(val=set_lin_log_val(turn_over_par,_t))

            if silent is False:
                print("---> setting par type turn-over energy, corresponding to par %s=%e"%(turn_over_par.name,turn_over_par.val))
                print("---> task completed", completed)
                print("---> using gamma_3p_SSC=%e"%gamma_3p_SSC)
                #print("---> ",turn_over_par.get_description())
                print()

        if silent is False:
            print ()
       
       
       
        #N
        N_par=jet_model.get_par_by_type('emitters_density')
        if N_par is not None:
            completed=rescale_Ne(jet_model,self.nuFnu_p_S_obs,self.nu_p_S_obs,self.rest_frame)
            if silent is False:
                print("---> setting par type emitters_density, corresponding to par %s"%(N_par.name))
                print("---> to N=%e" % (N_par.val))
                print("---> task completed",completed)
            if silent is False:
                #print('--->',N_par.get_description())
                print()

        B,B=find_B_from_nu_p_S(self.nu_p_S_obs,gamma_3p_SSC,self.rest_frame,self.beaming,z_par.val)
        if silent is False:
            print("---> setting B from nu_p_S to B=%e"%B)
            print("---> to B=%e" %B)

        if silent is False:
            print("---> setting B from best matching of nu_p_IC")
            print()

        if B_par is not None:
            if EC_componets_list is None:
                (B_from_nu_peaks, failed), completed = constr_B_from_nu_peaks (jet_model,self.nu_p_S_obs,self.nu_p_IC_obs,self.rest_frame,self.B_min,self.B_max,self.beaming,params_grid_size,silent=silent)
            else:
                (B_from_nu_peaks, failed), completed = constr_B_from_nu_peaks (jet_model,self.nu_p_S_obs,self.nu_p_IC_obs,self.rest_frame,self.B_min,self.B_max,self.beaming,params_grid_size,EC=True,silent=silent)

            if completed is True and failed is False:
                B_par.set(val=B_from_nu_peaks)

            if silent is False:
                print("---> setting par type magnetic_field, corresponding to par %s"%(B_par.name))
                #print("---> ",B_par.get_description())
                print("---> task completed ",completed)

            if failed is False and completed is True:
                if silent is False:
                    print ("---> best B found: %e"%B_par.val)
            else:
                if silent is False:
                    print("---> constrain failed, B set to: %e"%B_par.val)
                    print()

            if silent is False:
                print()
                print("---> update pars for new B ")
            
            #update gmin from new B
            gmin,completed= set_gmin_from_nu_cut_IR(self.nu_cut_IR,self.rest_frame,B_par.val,self.beaming,z_par.val)
            if par_LE_cut_off_par is not None:
                par_LE_cut_off_par.set(val=set_lin_log_val(par_LE_cut_off_par,gmin))

                if silent is False:
                    print("---> setting par type low-energy-cut-off, corresponding to par %s"%(par_LE_cut_off_par.name))
                    print("---> task completed", completed)
                    print("---> set to %e"%par_LE_cut_off_par.val)
                    print()

            #update gamma_3p and gamma_cut from new B
            gamma_3p_Sync,completed=find_gamma_Synch(self.nu_p_S_obs,self.rest_frame,B_par.val,self.beaming,z_par.val)
            if turn_over_par is not None:
                _t,completed=find_turn_over(jet_model,self.emitters_distr_spectral_type,gamma_3p_Sync)
                turn_over_par.set(val=set_lin_log_val(turn_over_par,_t))

                if silent is False:
                    print("---> setting par type low-energy-cut-off, corresponding to par %s"%(turn_over_par.name))
                    print("---> task completed", completed)
                    print("---> task completed ", completed)
                    print("---> using gamma_3p_Sync=",gamma_3p_Sync)
                    print("---> to %e"%turn_over_par.val)
                    print()
                
            
            #update gmax for New B
            gmax,completed=find_HE_cut_off(self.emitters_distr_spectral_type,self.S_nu_max,self.rest_frame,B_par.val,self.beaming,z_par.val)
            if silent is False:
                print("---> gamma_max=%e from nu_max_Sync= %e, using B=%e"%(gmax,self.S_nu_max,B_par.val))
                print("---> task completed", completed)
            if HE_cut_off_par is not None:
                HE_cut_off_par.set(val=set_lin_log_val(HE_cut_off_par,gmax))

                if silent is False:
                    print("---> setting par type high-energy-cut-off, corresponding to par %s"%(HE_cut_off_par.name))
                    print("---> set to %e"%HE_cut_off_par.val)
                    print()
            
            #update N for New B
            if N_par is not None:
                completed = rescale_Ne(jet_model, self.nuFnu_p_S_obs, self.nu_p_S_obs, self.rest_frame)
                if silent is False:
                    print("---> setting par type emitters_density, corresponding to par %s"%(N_par.name))
                    print("---> to N=%e" % (N_par.val))
                    print("---> task completed",completed)
                    print()
                #N_par.set(val=N)

        #Improve R according to CD
        if silent is False:
            print ("---> setting R from Compton Dominance (CD)")

        R_start=R_tvar
        if EC_componets_list  is None:
            (R_from_CD,failed),completed=constr_R_from_CD(jet_model,self.nuFnu_p_S_obs,self.nu_p_S_obs,self.nuFnu_p_IC_obs,self.nu_p_IC_obs,self.rest_frame,R_tvar,params_grid_size,silent=silent)
        else:
            (R_from_CD,failed),completed=constr_R_from_CD(jet_model,self.nuFnu_p_S_obs,self.nu_p_S_obs,self.nuFnu_p_IC_obs,self.nu_p_IC_obs,self.rest_frame,R_tvar,params_grid_size,EC=True,silent=silent)

        if failed is False and completed is True:
            #
            if R_par is not None:
                R_par.set(val=set_lin_log_val(R_par,R_from_CD))
                if silent is False:
                    print("---> setting par type region_size, corresponding to par %s"%(R_par.name))
                    print("---> set to %e"%R_par.val)
                    print("---> task completed", completed)
                #if silent is False:
                #    print("---> ",R_par.get_description())
                #    print()

                if N_par is not None:
                    completed=rescale_Ne(jet_model,self.nuFnu_p_S_obs,self.nu_p_S_obs,self.rest_frame)
                    if silent is False:
                        print("---> updating setting par type emitters_density, corresponding to par %s"%(N_par.name))
                        print("---> set to %e"%N_par.val)
                        print("---> task completed", completed)
        else:
            R_par.set(val=set_lin_log_val(R_par,R_start))

            if silent is False:
                print("---> constrain failed, R unchanged: ")
                print("---> set to %e"%R_par.val)
                print()
            
         
        #Check tau_gamma_gamma and t_var
        #check_gamma_tansp(jet_model,self.beaming,sp.logspace(21,28,10),self.rest_frame)
        if silent is False:
            print("---> t_var (days)", check_t_var(R_par.val_lin,self.beaming,z_par.val)/86400.)
    
                
        #jet_model.set_flag('obs_constrain_final')

        if silent is False:
            print()
        
            print("show pars")

        jet_model.show_pars()
        
        #sets to initial values
        
        #if self.nu_p_IC_obs is not None:
        jet_model.set_IC_nu_size(IC_nu_size_initial)
        jet_model.set_seed_nu_size(nu_seed_size_initial)

        if silent is False:
           print("eval_model")

        jet_model.eval(fill_SED=True)
        
        #jet_model.set_flag(flag_initial)
        #jet_model.set_path(path_initial)
        
        
#        if self.SEDShape.L_host!=None:
#           
#            SED_model=FitModel(jet=jet_model,template=self.SEDShape.host)
#       
#        else:
#        
#            SED_model=FitModel(jet=jet_model,template=self.SEDShape.host)

            
        return jet_model
        

def run_task(func):

    def func_wrapper(*args, **kwargs):
        completed=True

        try:
            return func(*args, **kwargs),completed
        except Exception as e:
            print('Exception',e)
            return None, completed
        print("-" * 20)

    return func_wrapper


def check_t_var(R,beaming,z):
    return R*(1+z)/(BlazarSED.vluce_cm*beaming)

@run_task
def set_gmin_from_nu_cut_IR(nu_cut_IR,rest_frame,B,beaming,z):
        if nu_cut_IR is not None:
            _v, completed= find_gamma_Synch (nu_cut_IR,rest_frame,B,beaming,z)
            return _v
        else:
            print ("--> !! No nu_cut_IR provided set gmin to 1")
            return 1
           
        
        
        
@run_task
def find_turn_over(jet,distr_e,gamma_3p):
    
    if distr_e=='lppl' or distr_e=='lp':
        r=jet.get_par_by_type('spectral_curvature').val
        s=jet.get_par_by_type('LE_spectral_slope').val
        _v, completed =find_gamma0(r,s,gamma_3p)
        return _v

    elif distr_e=='lpep':
        r=jet.get_par_by_type('spectral_curvature').val
        return gamma_3p/10**(1.5/r)

    #!! TO FIX
    elif distr_e=='plc':
        return gamma_3p*2

    
    elif distr_e=='bkn':
        "!!!!!!!!! MUST FIX ACCORDING P1 and P2!!!!!!"
        return gamma_3p


@run_task
def find_HE_cut_off(distr_e,nu_S_max,rest_frame,B,beaming,z):
    gamma_max_Sync,completed=find_gamma_Synch(nu_S_max,rest_frame,B ,beaming,z)
    
    if distr_e=='lppl' or distr_e=='lp' or distr_e=='lpep':
        return gamma_max_Sync
    
    elif  distr_e=='bkn' or distr_e=='plc' or distr_e=='pl':
        return gamma_max_Sync
    raise RuntimeError ('no gmax found!',distr_e)

@run_task
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



@run_task
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
    
@run_task
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


@run_task
def find_gamma_3p_SSC(nu_p_S,nu_p_IC,rest_frame,gamma_3p_Sync,beaming,z,nu_p_EC_seed_field=None,silent=False):
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
    comp_factor,completed=get_Comp_factor(gamma_3p_Sync,nu_p_seed_blob)
    if silent is False:
        print ("---> nu_p_seed_blob=%e"%nu_p_seed_blob)

    gamma_3p_SSC_TH=g_TH
    #print "--> gamma_3p_SSC_TH",gamma_3p_SSC_TH,g_TH
    if silent is False:
        print ("---> COMPTON FACTOR=%e"%comp_factor,gamma_3p_SSC_TH)
    
    if comp_factor>0.01:
        #print "--> correction", 0.2*pow(comp_factor,-0.45)
        gamma_3p_SSC=g_TH/(0.2*pow(comp_factor,-0.45))
    else:
        gamma_3p_SSC=gamma_3p_SSC_TH
    
    #   print "--> gamma_3p_SSC=",gamma_3p_SSC
    return gamma_3p_SSC


@run_task
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




@run_task
def find_s(class_obj,nu_p_S_obs,S_LE_slope,indices,silent=False):
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
                if silent is False:
                    print ("---> s_radio_mm",radio_mm_index,s_radio_mm)

                
        
        s_X=None
        if indices.get_by_name('X')is not None:
            if indices.get_by_name('X').val is not None:
                X_index=indices.get_by_name('X').val.spectral
                s_X=(1.0-2*X_index)
                if silent is False:
                    print ("---> s_X",s_X)
    


      
        s_Fermi=None
        if indices.get_by_name('Fermi') is not None:
            if indices.get_by_name('Fermi').val is not None:
                Fermi_index=indices.get_by_name('Fermi').val.spectral
                s_Fermi=(0.4+ Fermi_index*-1*1.60)
                if silent is False:
                   print ("---> s_Fermi",s_Fermi)
    
        
        s_UV_X=None
        if indices.get_by_name('UV_X') is not None:
            if indices.get_by_name('UV_X').val is not None:
                UV_X_index=indices.get_by_name('UV_X').val.spectral
                s_UV_X=(1+ UV_X_index*-1*2.0)
                if silent is False:
                    print ("---> s_UV_X",s_UV_X)
        
        s_Opt_UV=None
        if indices.get_by_name('Opt_UV') is not None:
            if indices.get_by_name('Opt_UV').val is not None:
                Opt_UV_index=indices.get_by_name('Opt_UV').val.spectral
                s_Opt_UV=(1 -1.0*Opt_UV_index*2.0)
                if silent is False:
                    print ("---> s_Opt_UV",Opt_UV_index,s_Opt_UV)
        
        if class_obj=='LSP' :
            #print s_X
            if s_X is not None and s_X-3.0<0:
                #if X_index+1.0>0.2:
                s=s_X
                if silent is False:
                    print ("---> s from X_index",s)
            elif s_radio_mm is not None and s_radio_mm+2>0.0:
                s=s_radio_mm
                if silent is False:
                    print ("---> s from radio_mm_index")
            else:   
                s=2.0
        
        if class_obj=='ISP':
            if s_radio_mm is not None and s_radio_mm+2>0.2:
                s=s_radio_mm
                if silent is False:
                    print ("---> s from radio_mm_index")
            else:   
                s=2.0

        s_SED_shape_photon=S_LE_slope-2.0
        s_SED_shape=-2.0*s_SED_shape_photon-1.0
        if silent is False:
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
                if silent is False:
                    print ("---> s from (s_Fermi + s_UV)/2")
            elif s_UV is not None :
                s=s_UV
                if silent is False:
                    print ("---> s from  s_UV")

            elif s_Fermi is not None :
                s=s_Fermi
                if silent is False:
                    print ("---> s from Fermi")

            
            elif s_radio_mm is not None :
                s=s_radio_mm
                if silent is False:
                    print ("---> s from radio_mm")

            else:
                s=2.0
    if silent is False:
        print ("---> power-law index s, class obj=%s s chosen is %f"%(class_obj,s))
    res = (s_radio_mm, s_X, s_Fermi,s)
    return res


@run_task
def find_s1(class_obj,indices):
    #print "---> !!! fake function, still to develop"
    val=3.5
    #print ("---> set s1 to %f"%val)
    return val

@run_task
def rescale_Ne(jet,S_p,nu_p,rest_frame):
    """Rescales N.blob to get the wanted Lp_S_blob
    """

    if rest_frame=='obs':
        jet.set_N_from_nuFnu(nuFnu_obs=S_p,nu_obs=nu_p)
    else:
        jet.set_N_from_nuLnu(nuLnu_src=S_p, nu_src=nu_p)

    return


@run_task
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

def check_boundaries(val,val_min,val_max,val_name,silent=False):
    failed=False
    if val<val_min or val>val_max:
        failed=True
        if silent is False:
            print ("---> %s=%e, out of boundaries %e %e, rejected"%(val_name,val,val_min,val_max))

        if val<val_min:
            val=val_min
        elif val>val_max:
            val=val_max

        if silent is False:
            print ("     Best %s not found, (temporary set to %e)"%(val_name,val))
    else:
        if silent is False:
            print ("     Best %s=%e"%(val_name,val))

    return val,failed

@run_task
def get_R_tvar(beaming,t_var_sec,z):
    return BlazarSED.vluce_cm*beaming*t_var_sec/(1+z)

def set_lin_log_val(p,v):
    if p.islog is True:
        v = m.log10(v)
    return v

@run_task
def constr_R_from_CD(jet,nuFnu_p_S,nu_p_S,nuFnu_p_IC,nu_p_IC,rest_frame,R_tvar,params_grid_size,EC=False,silent=False):
    #################################
    # This function must not        #
    # change jet.blob attributes    #
    # all the changed values must   #
    # be set back to their original #
    # values                        #
    #################################
    R_initial=jet.get_par_by_name('R').val_lin
    N_initial=jet.get_par_by_name('N').val

    CD_model_log=[]
    CD_obs=nuFnu_p_IC/nuFnu_p_S
    R_min=R_tvar/1000
    R_min=1E13
    R_max=R_tvar*2
    R_grid=np.logspace(np.log10(R_min), np.log10(R_max), params_grid_size)
    for R in R_grid:

        jet.set_par('R',val=set_lin_log_val(jet.get_par_by_name('R'),R))

        completed = rescale_Ne(jet, nuFnu_p_S, nu_p_S, rest_frame)

        jet.eval()

        Lp_S=jet.get_SED_peak(Model_dic.Sync_nuFnu_p_dic[rest_frame])
        if EC==False:
            Lp_IC=jet.get_SED_peak(Model_dic.SSC_nuFnu_p_dic[rest_frame])
        else:
            nu_p_IC,Lp_IC=jet.get_SED_peak(freq_range=[0.01*nu_p_IC,10*nu_p_IC])
        CD=Lp_IC/Lp_S
        CD_model_log.append(np.log10(CD))


    R_grid_log=np.log10(R_grid)
    p=polyfit(CD_model_log,R_grid_log,2)
    best_R=polyval(p, np.log10(CD_obs))
    Best_R=np.power(10., best_R)
    Best_R,failed=check_boundaries(Best_R,R_min,R_max,'R',silent=silent)
    

    jet.set_par('R',val=set_lin_log_val(jet.get_par_by_name('R'),R_initial))
    jet.set_par('N',val=N_initial)
    res = (Best_R,failed)
    return res


@run_task
def constr_B_from_nu_peaks(jet,nu_p_S,nu_p_IC,rest_frame,B_min,B_max,beaming,params_grid_size,EC=False,silent=False):
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
    
    B_grid=np.logspace(np.log10(B_min),np.log10(B_max),params_grid_size)
    
    #f=open('%s/B_vs_nu_p_IC.dat'%jet.get_path(),'w')
    
    #fin the turn-over variable to update with B
    turn_over_energy=jet.get_par_by_type('turn-over-energy')
    
    if turn_over_energy is not None :
        turn_over_energy_initial=turn_over_energy.val_lin
    
    
    z=jet.get_par_by_type('redshift').val
    
    for B in B_grid:
        
        jet.set_par('B',val=B)
        
        
        gamma_3p,completed=find_gamma_Synch(nu_p_S,rest_frame,B,beaming,z)

        #upda-turn-over variable
        if turn_over_energy is not None :
            turn_over_energy.set(val=set_lin_log_val(turn_over_energy,gamma_3p))

        
        jet.eval()   
        
        if EC==False:
            nu_p_model=jet.get_SED_peak(Model_dic.SSC_nu_p_dic[rest_frame])
        else:
             nu_p_model,nuFnu_p_model=jet.get_SED_peak(freq_range=[0.01*nu_p_IC,10*nu_p_IC])
        
        nu_p_IC_model_log.append(np.log10(nu_p_model))

        #print(n.log10(B),n.log10(nu_p_model),file=f)
        #print n.log10(B),n.log10(nu_p_model)

    #f.close()
    
    B_grid_log=np.log10(B_grid)
    
     
    p=polyfit(nu_p_IC_model_log,B_grid_log,2)    
    
    Best_B=polyval(p, np.log10(nu_p_IC))
    Best_B = np.power(10., Best_B)
    Best_B,failed  = check_boundaries(Best_B, B_min, B_max, 'B',silent=silent)
    
    
    #REST CHANGED VALUES
    jet.set_flag(flag_initial)
    jet.set_par('B',val=B_initial)

    if turn_over_energy is not None :
        turn_over_energy.set(val=set_lin_log_val(turn_over_energy,turn_over_energy_initial))

    res= (Best_B,failed)
    return res

#def check_gamma_transp(jet, beaming_val, nu_IC_data, rest_frame):
#     """retrun tau_gamma_gamma for a given  IC frequency
#
#     Args:
#         blob:
#         temap_ev:
#         nu_IC_data: data freq of the IC component to check fro tau_gamma_gamma
#
#     Returns:
#         tau_gamma_gamma
#
#     """
#     #################################
#     # This function must not        #
#     # change jet.blob attributes    #
#     # all the changed values must   #
#     # be set back to their original #
#     # values                        #
#     #################################
#
#     flag_initial=jet.get_flag()
#
#     jet.set_flag('tau_gamma_gamma')
#
#     IC_initial=jet.get_IC_mode()
#
#     jet.set_IC_mode('off')
#
#
#     #nu_IC_blob=convert_nu_to_blob(nu_IC_data,rest_frame,blob.beam_obj,blob.z_cosm)
#     #target_nu_blob = SED.MEC2 * SED.MEC2 / (SED.HPLANCK * SED.HPLANCK * nu_IC_blob)
#     target_nu_obs=BlazarSED.MEC2 * BlazarSED.MEC2 / (BlazarSED.HPLANCK * BlazarSED.HPLANCK * nu_IC_data)
#     #print nu_IC_blob,target_nu,target_nu*nu_IC_blob
#
#
#     jet.eval()
#
#
#     #nu_obs=sp.array([])
#     #nuFnu_obs=sp.array([])
#     nu_obs,nuFnu_obs=jet.spectral_components.Sync.get_SED_points(log_log=True,name='Sync')
#
#     #for i in range(BlazarSED.GetNuIntMaxSynch()):
#     #
#     #     if BlazarSED.GetSEDSynch(i)>0:
#     #         nu_obs=sp.append(nu_obs,sp.log10(BlazarSED.GetNuObsSynch(i)))
#     #        nuFnu_obs = sp.append(nuFnu_obs, sp.log10(BlazarSED.GetSEDSynch(i)))
#
#
#     nuFnu_target = sp.power(10.0,sp.interp(sp.log10(target_nu_obs),nu_obs,nuFnu_obs))
#
#     if sp.shape(target_nu_obs)==():
#         target_nu_obs=sp.array([target_nu_obs])
#
#         nuFnu_target=sp.array([nuFnu_target])
#
#         nu_IC_data=sp.array([nu_IC_data ])
#
#     for i in range(len(target_nu_obs)):
#         #print nuFnu_target,blob.beam_obj,blob.z_cosm,blob.dist
#
#
#
#         z_val=jet.get_par_by_type('redshift').val
#
#         DL_cm_val=jet.get_DL_cm()
#
#         nuLnu_target_blob= BlazarSED.nuFnu_obs_to_nuLnu_blob(nuFnu_target[i],beaming_val,z_val,DL_cm_val)
#
#         R_val=jet.get_par_by_type('region_size').val_lin
#
#         tau_gamma_gamma=BlazarSED.SIGTH/5*nuLnu_target_blob/(4*sp.pi * BlazarSED.vluce_cm * R_val * BlazarSED.MEC2 )
#
#         print ("     -->target comoving freq=%e"%target_nu_obs[i], "observed freq=%e"%nu_IC_data[i], "tau_gamma=%e"%(tau_gamma_gamma))
#
#
#     jet.set_IC_mode(IC_initial)
#     jet.set_flag(flag_initial)