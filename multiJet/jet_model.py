"""
Module: jet_model
========================



Overview
--------
   
This module  provides an interface to call the BlazarSED code, setting the
physical parameters and running the code. The BlazarSED code is a numerical 
accurate C code, to evaluate SSC/EC emission processes in a relativistic jet. 
The python wrappper is  built using SWIG. 

The :class:`.Jet` is used to create a *jet* object, providing the high 
level interface to set the jet physical  paramters, and evaluater the SSC/EC 
emission processes.



Classes relations
---------------------------------------
.. figure::  classes_jet_model.png
   :align:   center  


Classes and Inheritance Structure
----------------------------------------------
.. inheritance-diagram:: BlazarSEDFit.jet_model


  


Module API
-----------

Summary
---------
.. autosummary::
   JetParameter
   Jet
    
Module API

"""
'''
Created on 2013 1 27

@author: orion
'''


#import jet_wrapper
from jetkernel import jetkernel as BlazarSED

import os

import  spectral_shapes  

import numpy as np
from numpy import log10,array,zeros,power,shape

from scipy.interpolate import interp1d

from model_parameters import ModelParameterArray, ModelParameter
from base_model import  Model

from output import makedir,workplace,clean_dir

from sed_models_dic import nuFnu_obs_dic

class JetParameter(ModelParameter):
    """
    
    This class is a subclass of the :class:`.ModelParameter` class,
    extending the base class to  handles SSC/EC parameters, 
    overriding the :meth:`.ModelParameter.set` in order to propagate the
    parameter value to the BlazarSED object instance
           
 
    """
    def __init__(self,blob,**keywords):
    
        self.__blob=blob
        
        self.par_type_list=['electron_density',
                            'Disk',
                            'BLR',
                            'DT',
                            'region_size',
                            'electron_energy',
                            'LE_spectral_slope',
                            'HE_spectral_slope',
                            'high-energy-cut-off',
                            'low-energy-cut-off',
                            'spectral_curvature',
                            'turn-over-energy',
                            'magnetic_field',
                            'beaming',
                            'jet-viewing-angle',
                            'jet-bulk-factor',
                            'redshift']
        
        if 'par_type' in keywords.keys() and keywords['par_type'] not in self.par_type_list:
            msg= "blob_parameter type %s not allowed"%keywords['par_type'] + "\n please choose among %s"%self.par_type_list
            raise ValueError("%s"%msg)
            
        
 
        super(JetParameter,self).__init__(  **keywords)
        
       
        if 'val' in keywords.keys():
            val=keywords['val']
            self.assign_val(self.name,val)
            
        
        
    #OVERRIDES Base Method
    def set(self,**keywords):
        """        
        overrides the  :meth:`.ModelParameter.set` method in order to propagate the
        parameter value to the BlazarSED object instance
        """
        super(JetParameter,self).set(**keywords )
        
        
        if 'val' in keywords.keys():
            
            self.assign_val(self.name,keywords['val']) 
        
    
    def assign_val(self,name,val):
        """
        sets the :class:`.JetParameter` value in the BlazarSED object
        """
    
    
        b=getattr(self.__blob,name)
        
        
           
        if type(b)==int:
            val=int(val)
        elif type(b)==float:
            val=float(val)
        elif type(b)==str:
            val=val
        
        setattr(self.__blob,name,val)
        
        
    


def build_emitting_region_dic(beaming_expr='delta'):
    """
    
    Builds a dictionary to init the :class:`.JetParameter` 
    objects   for the emitting region:
        
        - **R**, the radius of the emitting region in cm
        
        - **B**, the magnetic field in G
        
        - **beaming**, the beaming factor
        
        - **z**, the redshift

   
    """
    
    model_dic={}
    model_dic['R']=['region_size',0,None,'cm']
    model_dic['B']=['magnetic_field',0,None,'G']
    
    if beaming_expr=='bulk_theta':
        model_dic['theta']=['jet-viewing-angle',0.0,None,'deg']
        model_dic['BulkFactor']=['jet-bulk-factor',1.0,None,'Lorentz-factor']
    elif beaming_expr=='delta' or beaming_expr=='':
        model_dic['beam_obj']=['beaming',1,None,'']
    else:
        raise  RuntimeError('''wrong beaming_expr="%s" value, allowed 'delta' or 'bulk_theta' '''%(beaming_expr))
    
    model_dic['z_cosm']=['redshift',0,None,'']
    
    
    return model_dic
    




def build_ExtFields_dic(EC_model_list,allowed_EC_components_list):
    """
    """
    

        
    
    model_dic={}

    for EC_model in EC_model_list:
            
        if EC_model not in allowed_EC_components_list:
            raise RuntimeError("EC model %s not allowed"%EC_model,"please choose among ", allowed_EC_components_list)
     
     
        if EC_model=='Disk':
            model_dic['L_disk']=['Disk',0,None,'erg/s']
            model_dic['R_inner_Sw']=['Disk',0,None,'Sw. radii']
            model_dic['R_ext_Sw']=['Disk',0,None,'Sw. radii']
            model_dic['T_disk_max']=['Disk',0,None,'K']
            model_dic['accr_eff']=['Disk',0,None,'']
            model_dic['R_H']=['Disk',0,None,'cm']
    
                    
        if EC_model=='BLR':
            model_dic['tau_BLR']=['BLR',0.0,1.0,'']
            model_dic['R_BLR_in']=['BLR',0,None,'cm']
            model_dic['R_BLR_out']=['BLR',0,None,'cm']
    
    
            
        if EC_model=='DT':
            model_dic['T_DT']=['DT',0.0,None,'K']
            model_dic['R_DT']=['DT',0.0,None,'cm']
            model_dic['tau_DT']=['DT',0.0,1.0,'']
    
    
    return model_dic




def build_electron_distribution_dic(electron_distribution):
    """
    Builds the dictionary to init the :class:`.JetParameter` 
    objects   for the electron  distribution:
        
    The following :class:`.JetParameter`: objects the *do not* depend
    on the type of electron distribution
            
            - N, particle density in cm^-3
            
            - gmin, the minimum value of the electron Lorentz factor
            
            - gmax, the maximum value of the electron Lorentz factor
        
    The following :class:`.JetParameter`: objects *depend* on the type of electron distribution:
        
        - **power law**, electron_distribution='pl'
           
           - p
        
        - **broken power-law**, electron_distribution= **'bkn'**

            - p
            - p_1
            - gamma_break
        
        - **log-parabola**, electron_distribution= **'lp'**
            
            - r 
            - s
            - gamma0_log_parab (fixed)
            
        - **log-parabola** with a low energy power-law tail, electron_distribution= **'lppl'**
            
            - r 
            - s
            - gamma0_log_parab
        
        - **log-parabola** defined by peak energy, electron_distribution= **'lpep'**
            
            - r 
            - s
            - gammap_log_parab, 
        
        - **power-law cut-off**, lectron_distribution= **'plc'**
        
            - p
            - gamma_cut
    
    """
     
        
    model_dic={}
    model_dic['N']=['electron_density',0,None,'cm^-3']
    model_dic['gmin']=['low-energy-cut-off',1,None,'Lorentz-factor']
    model_dic['gmax']=['high-energy-cut-off',1,None,'Lorentz-factor']
    
    
 
        
    if electron_distribution=='pl':
        model_dic['p']=['HE_spectral_slope',-10,10,'']
            
            
    if electron_distribution=='bkn':
        model_dic['p']=['LE_spectral_slope',-10,10,'']
        model_dic['p_1']=['HE_spectral_slope',-10,10,'']
        model_dic['gamma_break']=['turn-over-energy',1,None,'Lorentz-factor']
          
          
    if electron_distribution=='lp':
        model_dic['s']=['LE_spectral_slope',-10,10,'']
        model_dic['r']=['spectral_curvature',-10,10,'']
        model_dic['gamma0_log_parab']=['turn-over-energy',1,None,'Lorentz-factor',True]
        
        
    if electron_distribution=='lppl':
        model_dic['s']=['LE_spectral_slope',-10,10,'']
        model_dic['r']=['spectral_curvature',-10,10,'']
        model_dic['gamma0_log_parab']=['turn-over-energy',1,None,'Lorentz-factor']
 
 
    if electron_distribution=='lpep':   
        model_dic['r']=['spectral_curvature',-10,10,'']
        model_dic['gammap_log_parab']=['turn-over-energy',1,None,'Lorentz-factor']   
    
        
    if electron_distribution=='plc':
        model_dic['p']=['LE_spectral_slope',-10,10,'']
        model_dic['gamma_cut']=['turn-over-energy',1,None,'Lorentz-factor']

    
    return model_dic
            







class JetSpecComponent(object):
    """
    """
    def __init__(self,name,blob_object):
        
        self.name=name
        
        
        self.__nuFnu_name, self.__nu_name=nuFnu_obs_dic[self.name]
        
        self.nuFnu_ptr=getattr(blob_object,self.__nuFnu_name)
        
        self.nu_ptr=getattr(blob_object,self.__nu_name)
    
        self.SED=spectral_shapes.SED(name=self.name)
    



    
class Jet(Model):
    """

    This class allows to build a ``jet`` model providing the interface to the
    BlazarSED code, giving  full access to the physical parameters and
    providing the methods to run the code.
    A :class:`Jet` object  will store the
    the physical parameters in  the :class:`.ModelParameterArray` class,
    that is a collection of :class:`JetParameter` objects.

    :param name: (str), name id for the model
    :param electron_distribution: (str), the type of electron distribution

    :ivar parameters: instance of the :class:`.ModelParameterArray`
        storing and handling the array of  :class:`JetParameters` instances

    **Examples**

    .. literalinclude:: ../../../examples/modules/jet_model/Jet_example.py


    """

    def __init__(self,name='test',electron_distribution='lp',beaming_expr='delta',jet_workplace=None,verbose=None,clean_work_dir=True, **keywords):
        """ Constructor
            :param name: (str), name id for the model
            :param electron_distribution: (str), the type of electron distribution

        """

        super(Jet,self).__init__(  **keywords)

        self.name = name

        self.model_type='jet'

        self.__scale='lin-lin'

        self.__blob = init_SED(verbose=verbose)

        #self.jet_wrapper_dir=os.path.dirname(__file__)+'/jet_wrapper'

        #os.environ['BLAZARSED']=self.jet_wrapper_dir
        #print ("BLAZARSED DIR",self.jet_wrapper_dir)

        if jet_workplace is None:
            out_dir=workplace.out_dir+'/'+self.name+'_BalzarSED_prod/'
        else:
            out_dir=jet_workplace.out_dir+'/'+self.name+'_BalzarSED_prod/'

        self.set_path(out_dir,clean_work_dir=clean_work_dir)

        self.set_flag(self.name)

        self.init_BlazarSED()

        self._allowed_EC_components_list=['BLR','DT','Star','CMB','Disk','All','CMB_stat']

        self.EC_components_list =[]
        self.spectral_components=[]


        self.add_spectral_component('Sum')
        self.add_spectral_component('Sync')
        self.add_spectral_component('SSC')


        self.SED=self.get_spectral_component_by_name('Sum').SED

        self.parameters = ModelParameterArray()

        self.__emitting_region_dic = None
        self.__electron_distribution_dic= None
        self.__external_photon_fields_dic= None



        self.set_electron_distribution(electron_distribution)

        self.set_emitting_region(beaming_expr)


        self.flux_plot_lim=1E-30

    def show_spectral_components(self):

        print "Spectral components for Jet model:%s"%(self.name)

        for comp in self.spectral_components:

            print "comp: %s "%(comp.name)

        print




    def get_spectral_component_by_name(self,name,verbose=True):
        for i in xrange(len(self.spectral_components)):
            if self.spectral_components[i].name==name:
                return self.spectral_components[i]
        else:
            if verbose==True:
                print "no spectral components with name %s found"%name
                print "names in array are:"
                self.show_spectral_components()

            return None



    def list_spectral_components(self):
        for i in xrange(len(self.spectral_components)):
            print self.spectral_components[i].name





    def del_spectral_component(self,name,verbose=True):

        comp=self.get_spectral_component_by_name(name,verbose=verbose)

        if comp is not None:
            self.spectral_components.remove(comp)






    def add_spectral_component(self,name):

        self.spectral_components.append(JetSpecComponent(name,self.__blob))







    def set_emitting_region(self,beaming_expr):

        if  self.__emitting_region_dic is not None:
            self.del_par_from_dic(self.__emitting_region_dic)

        setattr(self.__blob,'BEAMING_EXPR',beaming_expr)

        self.__emitting_region_dic=build_emitting_region_dic(beaming_expr=beaming_expr)

        self.add_par_from_dic(self.__emitting_region_dic)

    def set_N_from_L(self,L_0, nu_0):
        """
        returns the normalization to math luminosity L_0, given the
        """
        N = 1.0
        self.set_par('N', N)
        gamma_grid_size = self._Jet__blob.gamma_grid_size
        self._Jet__blob.gamma_grid_size = 100

        z = self._Jet__blob.z_cosm
        delta = self.get_beaming()
        nu_blob = nu_0 * (1 + z) / delta
        self.init_BlazarSED()

        L_out = BlazarSED.Lum_Sync_at_nu(self._Jet__blob, nu_blob) * delta ** 4
        if L_out <= 0.0:
            # N*=10
            # myJet.set_par('N', N)
            # myJet.init_BlazarSED()
            # L_out=BlazarSED.Lum_Sync_at_nu(myJet._Jet__blob,nu_blob)*delta**4
            # myJet.set_par('B',1.0)
            # myJet.init_BlazarSED()
            # myJet._Jet__blob.nu=nu_blob*1000000
            self.parameters.show_pars()
            print '-> j', BlazarSED.j_nu_Sync(self._Jet__blob)
            L_out = BlazarSED.Lum_Sync_at_nu(self._Jet__blob, nu_blob) * delta ** 4
            print '-> L', L_out
            # sys.exit()

        # print '-> N',BlazarSED.j_nu_Sync(myJet._Jet__blob)
        # if int(L_out)==0:
        #    print 'set N->',L_out,nu_0,L_0,nu_blob,delta,myJet.get_par_by_name('B').val,myJet.get_par_by_name('N').val,myJet._Jet__blob.nu_B

        ratio = (L_out / L_0)

        # myJet.set_par('N', 1.0/ratio)
        # myJet.init_BlazarSED()
        # L_out=BlazarSED.Lum_Sync_at_nu(myJet._Jet__blob,nu_blob)*delta**4
        # print 'check->',L_out,nu_blob,L_0,L_out/L_0

        self._Jet__blob.gamma_grid_size = gamma_grid_size
        print 'setting N to ', N / ratio
        self.set_par('N', val=N / ratio)

    def set_B_eq(self, L_0, nu_0, B_min):
        """
        returns equipartiont B
        """
        # print B_min

        N_pts = 20
        b_grid = np.logspace(np.log10(B_min), 1, N_pts)
        U_e = np.zeros(N_pts)
        U_B = np.zeros(N_pts)
        N = np.zeros(N_pts)
        self.set_par('B', b_grid[0])
        self.init_BlazarSED()

        for ID, b in enumerate(b_grid):
            self.set_par('B', b)
            self.set_par('N', 1.0)
            # print 'B_eq',ID
            self.set_N_from_L(L_0, nu_0)
            N[ID]=self.get_par_by_name('N').val
            self.init_BlazarSED()
            #
            U_e[ID] = self._Jet__blob.U_e
            U_B[ID] = self._Jet__blob.UB
            # delta=Jet.get_beaming()
            # print "check L_in=%4.4e L_out=%4.4e"%(L_0,(L_0/delta**4)/BlazarSED.Power_Sync_Electron(Jet._Jet__blob))
        # fig, ax = plt.subplots()
        # ax.plot(b_grid,U_e)
        # ax.plot(b_grid,U_B)
        # ax.plot(b_grid,U_B+U_e)
        # ax.semilogy()
        # ax.semilogx()
        # plt.show()
        # print b_grid[np.argmin(U_B+U_e)]

        ID_min = np.argmin(U_B + U_e)

        # print'-> ID_min' ,ID_min
        # if ID_min>5:
        #    fig, ax = plt.subplots()
        #    ax.plot(b_grid,U_e)
        #    ax.plot(b_grid,U_B)
        #    ax.plot(b_grid,U_B+U_e)
        #    ax.semilogy()
        #    ax.semilogx()
        #    plt.show()

        self.set_par('B', val=b_grid[ID_min])
        self.set_par('N', val=N[ID_min])
        print 'setting B to ',b_grid[ID_min]
        print 'setting N to ',N[ID_min]
        return


    def set_electron_distribution(self,electron_distribution):

        elec_models_list=['lp','pl','lppl','lpep','plc','bkn']


        if electron_distribution not in elec_models_list:
            print "electron distribution model %s not allowed"%electron_distribution
            print "please choose among: ", elec_models_list
            return

        if self.__electron_distribution_dic is not None:
            self.del_par_from_dic(self.__electron_distribution_dic)


        self.__electron_distribution=electron_distribution

        self.__blob.DISTR=electron_distribution

        self.__electron_distribution_dic = build_electron_distribution_dic( self.__electron_distribution)

        self.add_par_from_dic(self.__electron_distribution_dic)







    def del_EC_component(self,EC_components):

        _del_EC_components=EC_components.split(',')
        if len(EC_components)==1:
            _del_EC_components=[EC_components]



        if 'All' in EC_components:
            _del_EC_components=self._allowed_EC_components_list[::]


        for EC_component in _del_EC_components:


            if EC_component not in self._allowed_EC_components_list:
                raise RuntimeError("EC_component %s not allowed" % EC_component, "please choose among ",
                                   self._allowed_EC_components_list)

            if EC_component=='Disk':
                if self.get_spectral_component_by_name('Disk', verbose=False) is None:
                    self.__blob.do_EC_Disk=0
                    self.del_spectral_component('EC_Disk',verbose=False)
                    self.del_spectral_component('Disk',verbose=False)
                    self.EC_components_list.remove(EC_component)

            if EC_component=='BLR':
                if self.get_spectral_component_by_name('BLR', verbose=False) is None:
                    self.__blob.do_EC_BLR=0

                    self.del_spectral_component('EC_BLR',verbose=False)
                    self.EC_components_list.remove(EC_component)

            if EC_component=='DT':
                if self.get_spectral_component_by_name('DT', verbose=False) is None:
                    self.__blob.do_EC_DT=0
                    self.del_spectral_component('EC_DT',verbose=False)
                    self.del_spectral_component('DT',verbose=False)
                    self.EC_components_list.remove(EC_component)

            if EC_component=='CMB':
                if self.get_spectral_component_by_name('CMB', verbose=False) is None:
                    self.__blob.do_EC_CMB=0
                    self.del_spectral_component('EC_CMB',verbose=False)
                    self.EC_components_list.remove(EC_component)

            if EC_component=='CMB_stat':
                if self.get_spectral_component_by_name('CMB_stat', verbose=False) is None:
                    self.__blob.do_EC_CMB_stat=0
                    self.del_spectral_component('EC_CMB_stat',verbose=False)
                    self.EC_components_list.remove(EC_component)

        self.del_par_from_dic(build_ExtFields_dic(_del_EC_components,self._allowed_EC_components_list))






    def add_EC_component(self,EC_components):

        _add_EC_components = EC_components.split(',')
        if len(_add_EC_components) == 1:
            _add_EC_components = [EC_components]

        if 'All' in EC_components:
            _add_EC_components=self._allowed_EC_components_list[::]

        for EC_component in _add_EC_components:
            if EC_component not in self._allowed_EC_components_list:
                raise RuntimeError("EC_component %s not allowed" % EC_component, "please choose among ",
                                   self._allowed_EC_components_list)




            if EC_component=='Disk':
                self.__blob.do_EC_Disk=1
                if self.get_spectral_component_by_name('EC_Disk',verbose=False) is None:
                    self.add_spectral_component('EC_Disk')
                    self.EC_components_list.append(EC_component)

                if self.get_spectral_component_by_name('Disk',verbose=False) is None:
                    self.add_spectral_component('Disk')


            if EC_component=='BLR':
                self.__blob.do_EC_BLR=1
                if self.get_spectral_component_by_name('EC_BLR',verbose=False) is None:
                    self.add_spectral_component('EC_BLR')
                    self.EC_components_list.append(EC_component)

                if self.get_spectral_component_by_name('Disk',verbose=False) is None:
                    self.add_spectral_component('Disk')


            if EC_component=='DT':
                self.__blob.do_EC_DT=1
                if self.get_spectral_component_by_name('EC_DT',verbose=False) is None:
                    self.add_spectral_component('EC_DT')
                    self.EC_components_list.append(EC_component)
                if self.get_spectral_component_by_name('DT',verbose=False) is None:
                    self.add_spectral_component('DT')


            if EC_component=='CMB':
                self.__blob.do_EC_CMB=1
                if self.get_spectral_component_by_name('EC_CMB',verbose=False) is None:
                    self.add_spectral_component('EC_CMB')
                    self.EC_components_list.append(EC_component)

            if EC_component == 'CMB_stat':
                self.__blob.do_EC_CMB_stat = 1
                if self.get_spectral_component_by_name('EC_CMB_stat', verbose=False) is  None:
                    self.add_spectral_component('EC_CMB_stat')
                    self.EC_components_list.append(EC_component)



        self.add_par_from_dic(build_ExtFields_dic(self.EC_components_list,self._allowed_EC_components_list))




    def del_par_from_dic(self,model_dic):
        """
        """
        for key in model_dic.keys():

            par=self.parameters.get_par_by_name(key)

            if par is not None:
                self.parameters.del_par(par)




    def add_par_from_dic(self,model_dic):
        """
        add the :class:`.JetParameter` object to the :class:`.ModelParameterArray`
        usign the dictionaries built by the :func:`build_emitting_region_dic`
        and :func:`build_electron_distribution_dic`
        """
        for key in model_dic.keys():

            pname=key
            pname_test=self.parameters.get_par_by_name(pname)
            if pname_test is None:
                pval=getattr(self.__blob,pname)
                ptype=model_dic[key][0]
                vmin=model_dic[key][1]
                vmax=model_dic[key][2]
                punit=model_dic[key][3]

                froz=False

                if len(model_dic[key])>4:
                    froz=model_dic[key][4]


                self.parameters.add_par(JetParameter(self.__blob,name=pname,par_type=ptype,val=pval,val_min=vmin,val_max=vmax,units=punit,frozen=froz))



    #PROTECTED MEMBER ACCESS
    #!! controlla i paramteri ....
    def get_electron_distribution(self):
        return self.__electron_distribution



    def get_DL_cm(self):
        self.init_BlazarSED()
        return self.__blob.dist



    def get_beaming(self,theta=None,bulk_factor=None,beaming=None):

        if theta is None:
            theta=self.__blob.theta

        if bulk_factor is None:
            bulk_factor=self.__blob.BulkFactor

        if beaming is None:
            beaming=self.__blob.beam_obj

        BlazarSED.SetBeaming(self.__blob)

        return self.__blob.beam_obj

    def set_flag(self,flag):
        self.__blob.STEM=flag

    def get_flag(self):
        return self.__blob.STEM

    def get_path(self):
        return self.__blob.path

    def set_path(self,path,clean_work_dir=True):
        self.__blob.path=path
        makedir(path,clean_work_dir=clean_work_dir)

    def set_SSC_mode(self,val):
        self.__blob.do_SSC=val

    def get_SSC_mode(self):
        return self.__blob.do_SSC

    def set_IC_mode(self,val):
        self.__blob.do_IC=val

    def get_IC_mode(self):
        return self.__blob.do_IC

    def set_sync_mode(self,val):
        self.__blob.do_Sync=val

    def get_sync_mode(self):
        return self.__blob.do_Sync

    def set_IC_nu_size(self,val):
        self.__blob.nu_IC_size=val

    def get_IC_nu_size(self):
        return self.__blob.nu_IC_size

    def set_seed_nu_size(self,val):
        self.__blob.nu_seed_size=val

    def get_seed_nu_size(self):
        return self.__blob.nu_seed_size

    def set_gamma_grid_size(self,val):
        self.__blob.gamma_grid_size=val

    def get_gamma_grid_size(self):
        return self.__blob.gamma_grid_size


    def set_verbosity(self,val):
        self.__blob.verbose=val

    def get_verbosity(self):
        return  self.__blob.verbose

    def debug_synch(self):
        print "nu stop synch", self.__blob.nu_stop_Sync
        print "nu stop synch ssc", self.__blob.nu_stop_Sync_ssc
        print "ID MAX SYNCH", self.__blob.NU_INT_STOP_Sync_SSC

    def debug_SSC(self):
        print "nu start SSC", self.__blob.nu_start_SSC
        print "nu stop SSC", self.__blob.nu_stop_SSC
        print "ID MAX SSC", self.__blob.NU_INT_STOP_COMPTON_SSC

    def set_par(self,par_name,val):
        """
        shortcut to :class:`ModelParametersArray.set` method
        set a parameter value

        :param par_name: (srt), name of the parameter
        :param val: parameter value

        """

        self.parameters.set(par_name, val=val)



    def get_par_by_type(self,par_type):
        """

        get parameter by type

        """
        for param in self.parameters.par_array:
            if param.par_type==par_type:
                return param

        return None

    def get_par_by_name(self,par_name):
        """

        get parameter by type

        """
        for param in self.parameters.par_array:
            if param.name==par_name:
                return param

        return None


    def show_model(self):
        self.show_pars()


    def show_pars(self):
        """
        shortcut to :class:`ModelParametersArray.show_pars` method
        shows all the paramters in the model

        """

        print "-----------------------------------------------------------------------------------------"
        print "model parameters for jet model:"
        print


        print "electron distribution type = %s  "%(self.__electron_distribution)

        self.parameters.show_pars()

        print "-----------------------------------------------------------------------------------------"


    def PlotModel(self,Plot,clean=False,autoscale=False,label=None,comp=None):
        if Plot is not None:
            if clean==True:
                Plot.clean_model_lines()

            for i in xrange(len(self.spectral_components)):
                if (comp is  None) or ((comp is not None)  and (comp==self.spectral_components[i].name)):

                    if self.spectral_components[i].name=='Sum':
                        if label is None:
                            comp_label= 'Jet model'
                        else:
                            comp_label=label
                        line_style='-'
                    else:
                        comp_label=self.spectral_components[i].name
                        line_style='--'
                    Plot.add_model_plot(self.spectral_components[i].SED,autoscale=autoscale,line_style=line_style,label=comp_label)
                else:
                    pass
        else:
            print "the plot window is not defined"

    def init_BlazarSED(self):
        BlazarSED.Init(self.__blob)




    def eval(self,init=True,fill_SED=True,nu=None,get_model=False,loglog=False,plot=None,label=None,phys_output=False):
        """
        Runs the BlazarSED  code for the current `JetModel` instance.

        :param init: (boolean), "defualt=True" initializes the BlazarSED code
            for the current `Jet` instance parameters values.

        """




        if init==True:

            BlazarSED.Init(self.__blob)



        BlazarSED.Run_SED(self.__blob)


        if phys_output==True:
            BlazarSED.EnergeticOutput(self.__blob)

        nu_sed_sum,nuFnu_sed_sum= self.get_SED_points()

        if fill_SED==True:

            self.SED.fill(nu=nu_sed_sum,nuFnu=nuFnu_sed_sum)

            for i in xrange(len(self.spectral_components)):

                print ('fill name',self.spectral_components[i].name)
                nu_sed,nuFnu_sed= self.get_SED_points(name=self.spectral_components[i].name)

                self.spectral_components[i].SED.fill(nu=nu_sed,nuFnu=nuFnu_sed)

        if nu is None:



            if loglog==True:

                model=log10(nuFnu_sed_sum)
            else:
                model=nuFnu_sed_sum

        else :

            if shape(nu)==():

                nu=array([nu])

            if loglog==False:
                nu_log=log10(nu)
            else:
                nu_log=nu

            #nu_sed_log,nuFnu_sed_log= self.get_SED_points(log_log=True)

            nu_sed_log=np.log10(nu_sed_sum)
            nuFnu_sed_log=np.log10(nuFnu_sed_sum)

            #print "here: ",nu_log, nu_sed_log.min()
            msk= nu_log > nu_sed_log.min()


            msk*=  nu_log < nu_sed_log.max()



            if loglog==False:
                model=zeros(nu_log.size)

                nuFnu_log=np.log10(nuFnu_sed_sum)

                f_interp =interp1d(nu_sed_log,nuFnu_sed_log,bounds_error=True)


                model[msk] = power(10.0,f_interp(nu_log[msk]))

            else:
                model=zeros(nu_log.size) + np.log10(self.flux_plot_lim)

                f_interp =interp1d(nu_sed_log,nuFnu_sed_log,bounds_error=True)

                model[msk] = f_interp(nu_log[msk])


        if plot is not None:
            if label is None:
                label= self.name

            self.PlotModel(plot, clean=True, label=self.name)


        if get_model==True:
            return model
        else:
            return None





    def get_SED_points(self,log_log=False,name='Sum'):

        try:
            spec_comp=self.get_spectral_component_by_name(name)


            nuFnu_ptr=spec_comp.nuFnu_ptr
            nu_ptr=spec_comp.nu_ptr

            size=self.__blob.spec_array_size
            x=zeros(size)
            y=zeros(size)

            for i in xrange(size):
                x[i]=BlazarSED.get_freq_array(nu_ptr,self.__blob,i)
                y[i]=BlazarSED.get_freq_array(nuFnu_ptr,self.__blob,i)

                #print"%s %e %e"%(name,x[i],y[i])

            msk=y>self.flux_plot_lim

            x=x[msk]
            y=y[msk]



            if log_log==True:


                #x=x[msk]
                #    y=y[msk]

                x=log10(x)
                y=log10(y)


            return x,y

        except:

            print "no spectral model found"




    def get_SED_peak(self,peak_name=None,freq_range=None,log_log=False):

        if peak_name is not None and freq_range is not None:
            print "either you provide peak_name or freq_range"
            raise ValueError
        elif   peak_name is not None:
            try:
                if log_log==False:
                    return getattr(self.__blob,peak_name)
                else:
                     return np.log10(getattr(self.__blob,peak_name) )
            except:
                print "peak name %s, not found, check name"%peak_name
                raise ValueError

        else:
            x,y=self.get_SED_points(log_log=log_log)
            msk1=x>freq_range[0]
            msk2=x<freq_range[1]
            #print x,freq_range,x[msk1*msk2]
            y_m= y[msk1*msk2].max()
            x_id= np.argmax(y[msk1*msk2])
            return x[msk1*msk2][x_id],y_m

def init_SED(verbose=None):
    
   
    
    blob=BlazarSED.MakeBlob()
    #temp_ev=BlazarSED.MakeTempEv()

    if verbose is None:
        blob.verbose=0
    else:
        blob.verbose=verbose
    

    blob.path="./"
    
    
    blob.MODE='custom'
    blob.gamma_grid_size=1000
    
    blob.nu_IC_size =50
    blob.nu_seed_size =100
    
    blob.do_Sync=2
    
    blob.do_SSC=1
   
    
    blob.R=3.0e15
    
    blob.B=0.1
    
    blob.z_cosm=0.1

    blob.BulkFactor=10
    blob.theta=0.1
    
    
    blob.N=100
    
    blob.gmin=2

    blob.gmax=1e8
    
    
     
    
    blob.nu_start_Sync = 1e6
    blob.nu_stop_Sync = 1e20
    
    blob.nu_start_SSC=1e14
    blob.nu_stop_SSC=1e30
    
    blob.nu_start_grid=1e6
    blob.nu_stop_grid=1e30
    
    
    
    
    return blob




        
