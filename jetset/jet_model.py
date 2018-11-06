from __future__ import absolute_import, division, print_function

from builtins import (bytes, str, open, super, range,
                      zip, round, input, int, pow, object, map, zip)


"""
This module  provides an interface to call the BlazarSED code, setting the
physical parameters and running the code. The BlazarSED code is a numerical 
accurate C code, to evaluate SSC/EC emission processes in a relativistic jet. 
The python wrappper is  built using SWIG. 

The :class:`.Jet` is used to create a *jet* object, providing the high 
level interface to set the jet physical  paramters, and evaluater the SSC/EC 
emission processes.

"""



'''
Created on 2013 1 27

@author: andrea tramacere
'''



__author__ = "Andrea Tramacere"





import os
import json
import  sys
import pickle

import numpy as np
from numpy import log10,array,zeros,power,shape

from scipy.interpolate import interp1d

#import jet_wrapper
from .jetkernel import jetkernel as BlazarSED

from . import spectral_shapes

from .model_parameters import ModelParameterArray, ModelParameter
from .base_model import  Model

from .output import makedir,WorkPlace,clean_dir

from .sed_models_dic import nuFnu_obs_dic,gamma_dic

from  .plot_sedfit import Plot

__all__=['Jet','JetParameter','JetSpecComponent','ElectronDistribution','build_emitting_region_dic',
         'build_ExtFields_dic','init_SED']


def str_hook(pairs):
    new_pairs = []
    for key, value in pairs:
        if isinstance(value, unicode):
            value = value.encode('utf-8')
        if isinstance(key, unicode):
            key = key.encode('utf-8')
        new_pairs.append((key, value))
    return dict(new_pairs)

class JetParameter(ModelParameter):
    """
    
    This class is a subclass of the :class:`.ModelParameter` class,
    extending the base class to  handles SSC/EC parameters, 
    overriding the :meth:`.ModelParameter.set` in order to propagate the
    parameter value to the BlazarSED object instance
           
 
    """
    def __init__(self,blob,**keywords):
    
        self._blob=blob
        
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
    
    
        b=getattr(self._blob,name)
        
        
           
        if type(b)==int:
            val=int(val)
        elif type(b)==float:
            val=float(val)
        elif type(b)==str:
            val=val
        
        setattr(self._blob,name,val)
        
        
    


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





            






class ElectronDistribution(object):

    def __init__(self,name,jet,gamma_grid_size=None):

        self._name=name



        self._jet=jet
        set_str_attr(jet._blob,'DISTR',name)
        #set_str(jet._blob.DISTR,name)

        if gamma_grid_size is not None:
            self.set_grid_size(gamma_grid_size)

        else:
            self._set_blob()
            self._fill()

    def update(self):
        self._set_blob()
        self._fill()

    def set_grid_size(self,gamma_grid_size):
        setattr(self._jet._blob,'gamma_grid_size' ,gamma_grid_size)
        self._set_blob()
        self._fill()

    def _set_blob(self):
        BlazarSED.Init(self._jet._blob)
        self._N_name, self._gamma_name = gamma_dic['electron_distr']
        self.Ne_ptr = getattr(self._jet._blob, self._N_name)
        self.gamma_ptr = getattr(self._jet._blob, self._gamma_name)


    def _fill(self):
        size = self._jet._blob.gamma_grid_size
        self.gamma = zeros(size)
        self.n_gamma = zeros(size)

        for ID in range(size):
            self.gamma[ID] = BlazarSED.get_elec_array(self.gamma_ptr, self._jet._blob, ID)
            self.n_gamma[ID] = BlazarSED.get_elec_array(self.Ne_ptr, self._jet._blob, ID)


    def plot(self, ax=None, y_min=None,y_max=None):
        self.update()
        import  pylab as plt
        if ax is None:
            fig,ax= plt.subplots()

        ax.plot(np.log10(self.gamma),np.log10(self.n_gamma))
        ax.set_xlabel(r'log($\gamma$)')
        ax.set_ylabel(r'log(n($\gamma$))')
        ax.set_ylim(y_min,y_max)
        plt.show()

        return ax,fig


    def plot3p(self, ax=None,y_min=None,y_max=None):
        self.update()
        import  pylab as plt
        if ax is None:
            fig,ax= plt.subplots()

        ax.plot(np.log10(self.gamma),np.log10(self.n_gamma*self.gamma*self.gamma*self.gamma))
        ax.set_xlabel(r'log($\gamma$)')
        ax.set_ylabel(r'log(n($\gamma$))')
        ax.set_ylim(y_min, y_max)
        plt.show()

        return ax,fig


    def _build_electron_distribution_dic(self,electron_distribution_name):
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

        model_dic = {}
        model_dic['N'] = ['electron_density', 0, None, 'cm^-3']
        model_dic['gmin'] = ['low-energy-cut-off', 1, None, 'Lorentz-factor']
        model_dic['gmax'] = ['high-energy-cut-off', 1, None, 'Lorentz-factor']

        if electron_distribution_name == 'pl':
            model_dic['p'] = ['HE_spectral_slope', -10, 10, '']

        if electron_distribution_name == 'bkn':
            model_dic['p'] = ['LE_spectral_slope', -10, 10, '']
            model_dic['p_1'] = ['HE_spectral_slope', -10, 10, '']
            model_dic['gamma_break'] = ['turn-over-energy', 1, None, 'Lorentz-factor']

        if electron_distribution_name == 'lp':
            model_dic['s'] = ['LE_spectral_slope', -10, 10, '']
            model_dic['r'] = ['spectral_curvature', -15, 15, '']
            model_dic['gamma0_log_parab'] = ['turn-over-energy', 1, None, 'Lorentz-factor', True]

        if electron_distribution_name == 'lppl':
            model_dic['s'] = ['LE_spectral_slope', -10, 10, '']
            model_dic['r'] = ['spectral_curvature', -15, 15, '']
            model_dic['gamma0_log_parab'] = ['turn-over-energy', 1, None, 'Lorentz-factor']

        if electron_distribution_name == 'lpep':
            model_dic['r'] = ['spectral_curvature', -15, 15, '']
            model_dic['gammap_log_parab'] = ['turn-over-energy', 1, None, 'Lorentz-factor']

        if electron_distribution_name == 'plc':
            model_dic['p'] = ['LE_spectral_slope', -10, 10, '']
            model_dic['gamma_cut'] = ['turn-over-energy', 1, None, 'Lorentz-factor']

        if electron_distribution_name == 'spitkov':
            model_dic['spit_index'] = ['LE_spectral_slope', -10, 10, '']
            model_dic['spit_temp'] = ['turn-over-energy', 1, None, 'Lorentz-factor']
            model_dic['spit_gamma_th'] = ['turn-over-energy', 1, None, 'Lorentz-factor']

        if electron_distribution_name == 'lppl_pile_up':
            model_dic['s'] = ['LE_spectral_slope', -10, 10, '']
            model_dic['r'] = ['spectral_curvature', -15, 15, '']
            model_dic['gamma0_log_parab'] = ['turn-over-energy', 1, None, 'Lorentz-factor']
            model_dic['gamma_inj'] = ['turn-over-energy', 1, None, 'Lorentz-factor']

            model_dic['gamma_pile_up'] = ['turn-over-energy', 1, None, 'Lorentz-factor']
            model_dic['ratio_pile_up'] = ['turn-over-energy', 0, None, '']

            model_dic['alpha_pile_up'] = ['turn-over-energy', 0.0, 10, '']
           # model_dic['ratio_pile_up'] = ['turn-over-energy',0, None, '']

        if electron_distribution_name == 'bkn_pile_up':
            model_dic['p'] = ['LE_spectral_slope', -10, 10, '']
            model_dic['p_1'] = ['HE_spectral_slope', -10, 10, '']
            model_dic['gamma_break'] = ['turn-over-energy', 1, None, 'Lorentz-factor']

            model_dic['gamma_pile_up'] = ['turn-over-energy', 1, None, 'Lorentz-factor']
            model_dic['gamma_pile_up_cut'] = ['turn-over-energy', 1, None, 'Lorentz-factor']

            model_dic['alpha_pile_up'] = ['turn-over-energy', 0.0, 10, '']


        return model_dic






class JetSpecComponent(object):
    """
    """
    def __init__(self,name,blob_object):
        
        self.name=name
        
        
        self._nuFnu_name, self._nu_name=nuFnu_obs_dic[self.name]
        
        self.nuFnu_ptr=getattr(blob_object,self._nuFnu_name)
        
        self.nu_ptr=getattr(blob_object,self._nu_name)
    
        self.SED=spectral_shapes.SED(name=self.name)
    
    def plot(self):
        pass



class TempEvol(Model):

    def __init__(self,verbose=None,out_dir='temp_ev',flag='test',clean_work_dir=True):

        self._temp_ev = self.build(verbose=verbose)

        self.set_path(out_dir, clean_work_dir=clean_work_dir)

        self.set_flag(flag)

    def build_TempEv(self,duration=1E5,
                     TStart_Acc=0.,
                     TStop_Acc=0.,
                     TStart_Inj=0.,
                     TStop_Inj=0.,
                     T_esc_Coeff=1E60,
                     Esc_Index=0.,
                     Acc_Index=1.,
                     Diff_Index=2.,
                     T_SIZE=5000,
                     NUM_SET=50,
                     Lambda_max_Turb=1E30):


        self._temp_ev = BlazarSED.MakeTempEv()

        self._temp_ev.duration = duration
        self._temp_ev.TStart_Acc = TStart_Acc
        self._temp_ev.TStop_Acc = TStop_Acc
        self._temp_ev.TStart_Inj = TStart_Inj
        self._temp_ev.TStop_Inj = TStop_Inj
        self._temp_ev.T_esc_Coeff = T_esc_Coeff
        self._temp_ev.Esc_Index = Esc_Index
        self._temp_ev.Acc_Index = Acc_Index
        self._temp_ev.Diff_Index = Diff_Index
        self._temp_ev.T_SIZE = T_SIZE
        self._temp_ev.NUM_SET = NUM_SET
        self._temp_ev.Lambda_max_Turb=Lambda_max_Turb


    def set_path(self,path,clean_work_dir=True):
        set_str_attr(self._temp_ev,'path',path)
        #set_str(self._blob.path,path)
        makedir(path,clean_work_dir=clean_work_dir)

    def set_flag(self,flag):
        self._temp_ev.STEM=flag


    def run(self,jet):
        BlazarSED.temp_evolution(jet._blob, self._temp_ev)




class Jet(Model):
    """

    This class allows to build a ``jet`` model providing the interface to the
    BlazarSED code, giving  full access to the physical parameters and
    providing the methods to run the code.
    A :class:`Jet` object  will store the
    the physical parameters in  the :class:`.ModelParameterArray` class,
    that is a collection of :class:`JetParameter` objects.


    **Examples**




    """

    def __init__(self,name='test',electron_distribution='lp',beaming_expr='delta',jet_workplace=None,verbose=None,clean_work_dir=True, **keywords):


        super(Jet,self).__init__(  **keywords)

        self.name = name

        self.model_type='jet'

        self._scale='lin-lin'

        self._blob = self.build_blob(verbose=verbose)

        #self.jet_wrapper_dir=os.path.dirname(__file__)+'/jet_wrapper'

        #os.environ['BLAZARSED']=self.jet_wrapper_dir
        #print ("BLAZARSED DIR",self.jet_wrapper_dir)

        if jet_workplace is None:
            jet_workplace=WorkPlace()
            out_dir= jet_workplace.out_dir + '/' + self.name + '_jet_prod/'
        else:
            out_dir=jet_workplace.out_dir+'/'+self.name+'_jet_prod/'

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

        self._emitting_region_dic = None
        self._electron_distribution_dic= None
        self._external_photon_fields_dic= None




        self.set_electron_distribution(electron_distribution)

        self.set_emitting_region(beaming_expr)


        self.flux_plot_lim=1E-30

    def build_blob(self,verbose=None):

        blob = BlazarSED.MakeBlob()
        # temp_ev=BlazarSED.MakeTempEv()

        if verbose is None:
            blob.verbose = 0
        else:
            blob.verbose = verbose

        set_str_attr(blob, 'path', './')
        # set_str(blob.path,'./')

        set_str_attr(blob, 'MODE', 'custom')
        # blob.MODE='custom'
        blob.gamma_grid_size = 1000

        blob.nu_IC_size = 50
        blob.nu_seed_size = 100

        blob.do_Sync = 2

        blob.do_SSC = 1

        blob.R = 3.0e15

        blob.B = 0.1

        blob.z_cosm = 0.1

        blob.BulkFactor = 10
        blob.theta = 0.1

        blob.N = 100

        blob.gmin = 2

        blob.gmax = 1e8

        blob.nu_start_Sync = 1e6
        blob.nu_stop_Sync = 1e20

        blob.nu_start_SSC = 1e14
        blob.nu_stop_SSC = 1e30

        blob.nu_start_grid = 1e6
        blob.nu_stop_grid = 1e30

        return blob

    def save_model(self,file_name):
        _model={}
        _model['electron_distribution']=self._electron_distribution_name
        _model['beaming_expr']=self._beaming_expr
        _model['model_spectral_components']=self.get_spectral_component_names_list()
        _model['EC_components_list'] = self.EC_components_list
        _model['pars']={}

        for p in self.parameters.par_array:
            _model['pars'][p.name]=p.val

        with open(file_name, 'w', encoding="utf-8") as outfile:
            outfile.write(str(json.dumps(_model, ensure_ascii=False)))







    @classmethod
    def load_model(cls,file_name):
        jet = cls()
        with open(file_name, 'r') as infile:
            _model = json.load(infile)

        print ('_model',_model)

        jet.model_type = 'jet'

        jet.init_BlazarSED()
        jet.parameters = ModelParameterArray()

        jet.set_electron_distribution(str(_model['electron_distribution']))
        _l=jet.get_spectral_component_names_list()
        #print ('->',_model['EC_components_list'],_model['model_spectral_components'],_l,self._allowed_EC_components_list)
        for c in _model['model_spectral_components']:
            print (c)
            if str(c) not  in _l:

                #print ('test',c.replace('EC_',''), _model['EC_components_list'],c.replace('EC_','') in _model['EC_components_list'])
                if c.replace('EC_','') in _model['EC_components_list']:
                    print ('add EC',c.replace('EC_',''))
                    jet.add_EC_component(str(c.replace('EC_','')))
                else:
                    jet.add_spectral_component(str(c))

        jet.SED = jet.get_spectral_component_by_name('Sum').SED

        jet.set_emitting_region(str(_model['beaming_expr']))
        jet.set_electron_distribution(str(_model['electron_distribution']))
        _par_dict=_model['pars']
        jet.show_pars()
        for k in _par_dict.keys():
            print ('set', k,_par_dict[k])
            jet.set_par(par_name=str(k),val=_par_dict[str(k)])

        return jet

    def set_electron_distribution(self,name):
        self.electron_distribution = ElectronDistribution(name, self)

        self._electron_distribution_name=name

        elec_models_list = ['lp', 'pl', 'lppl', 'lpep', 'plc', 'bkn','spitkov','lppl_pile_up','bkn_pile_up']

        if name not in elec_models_list:
            print ("electron distribution model %s not allowed" % name)
            print ("please choose among: ", elec_models_list)
            return

        if self._electron_distribution_dic is not None:
            self.del_par_from_dic(self._electron_distribution_dic)


        self._electron_distribution_dic = self.electron_distribution._build_electron_distribution_dic(self._electron_distribution_name)

        self.add_par_from_dic(self._electron_distribution_dic)


    def show_spectral_components(self):

        print ("Spectral components for Jet model:%s"%(self.name))

        for comp in self.spectral_components:

            print ("comp: %s "%(comp.name))

        print




    def get_spectral_component_by_name(self,name,verbose=True):
        for i in range(len(self.spectral_components)):
            if self.spectral_components[i].name==name:
                return self.spectral_components[i]
        else:
            if verbose==True:
                print ("no spectral components with name %s found"%name)
                print ("names in array are:")
                self.show_spectral_components()

            return None



    def list_spectral_components(self):
        for i in range(len(self.spectral_components)):
            print (self.spectral_components[i].name)

    def get_spectral_component_names_list(self):
        _l=[]
        for i in range(len(self.spectral_components)):
            _l.append(self.spectral_components[i].name)
        return _l




    def del_spectral_component(self,name,verbose=True):

        comp=self.get_spectral_component_by_name(name,verbose=verbose)

        if comp is not None:
            self.spectral_components.remove(comp)






    def add_spectral_component(self,name):

        self.spectral_components.append(JetSpecComponent(name,self._blob))







    def set_emitting_region(self,beaming_expr):

        if  self._emitting_region_dic is not None:
            self.del_par_from_dic(self._emitting_region_dic)

        self._beaming_expr=beaming_expr

        set_str_attr(self._blob,'BEAMING_EXPR',beaming_expr)

        #setattr(self._blob,'BEAMING_EXPR',beaming_expr)

        self._emitting_region_dic=build_emitting_region_dic(beaming_expr=beaming_expr)

        self.add_par_from_dic(self._emitting_region_dic)




    def set_N_from_nuLnu(self,L_0, nu_0):
        """
        sets the normalization of N to match the rest frame luminosity L_0, at a given frequency nu_0
        """
        N = 1.0
        setattr(self._blob, 'N', N)
        gamma_grid_size = self._blob.gamma_grid_size

        self._blob.gamma_grid_size = 100

        z = self._blob.z_cosm

        self.init_BlazarSED()

        delta = self._blob.beam_obj
        nu_blob = nu_0 / delta
        # print 'delta for set N form L' ,delta

        L_out = BlazarSED.Lum_Sync_at_nu(self._blob, nu_blob) * delta ** 4
        # if L_out <= 0.0:
        #    L_out = BlazarSED.Lum_Sync_at_nu(blob, nu_blob) * delta ** 4

        ratio = (L_out / L_0)

        self._blob.gamma_grid_size = gamma_grid_size

        #print 'N',N/ratio
        self.set_par('N', val=N/ratio)


    def set_N_from_nuFnu(self, nuFnu_obs, nu_obs):
        """
        sets the normalization of N to match the observed flux nu0F_nu0 at a given frequency nu_0
        """

        self.init_BlazarSED()
        DL = self._blob.dist

        L = nuFnu_obs * DL * DL * 4.0 * np.pi

        nu_rest = nu_obs * (1 + self._blob.z_cosm)

        self.set_N_from_nuLnu( L, nu_rest)



    def set_B_eq(self, nuFnu_obs, nu_obs, B_min=1E-9,B_max=1.0,N_pts=20,plot=False):
        """
        returns equipartiont B
        """
        # print B_min


        b_grid = np.logspace(np.log10(B_min), np.log10(B_max), N_pts)
        print ('B grid min ',B_min)
        print ('B grid max ',B_max)
        print ('grid points',N_pts)
        U_e = np.zeros(N_pts)
        U_B = np.zeros(N_pts)
        N = np.zeros(N_pts)
        self.set_par('B', b_grid[0])
        self.init_BlazarSED()

        for ID, b in enumerate(b_grid):
            self.set_par('B', b)
            self.set_par('N', 1.0)
            # print 'B_eq',ID
            self.set_N_from_nuFnu(nuFnu_obs, nu_obs)
            N[ID]=self.get_par_by_name('N').val
            self.init_BlazarSED()
            #
            U_e[ID] = self._blob.U_e
            U_B[ID] = self._blob.UB
            # delta=Jet.get_beaming()
            # print "check L_in=%4.4e L_out=%4.4e"%(L_0,(L_0/delta**4)/BlazarSED.Power_Sync_Electron(Jet._Jet__blob))

        ID_min = np.argmin(U_B + U_e)

        if plot==True:
            import  pylab as plt
            fig, ax = plt.subplots()
            ax.plot(b_grid,U_e)
            ax.plot(b_grid,U_B)
            ax.plot(b_grid,U_B+U_e)
            ax.scatter(b_grid[ID_min],U_e[ID_min]+U_B[ID_min])

            ax.semilogy()
            ax.semilogx()
            plt.show()


        self.set_par('B', val=b_grid[ID_min])
        self.set_par('N', val=N[ID_min])
        print ('setting B to ',b_grid[ID_min])
        print ('setting N to ',N[ID_min])
        return b_grid[ID_min],b_grid,U_B,U_e







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
                    self._blob.do_EC_Disk=0
                    self.del_spectral_component('EC_Disk',verbose=False)
                    self.del_spectral_component('Disk',verbose=False)
                    self.EC_components_list.remove(EC_component)

            if EC_component=='BLR':
                if self.get_spectral_component_by_name('BLR', verbose=False) is None:
                    self._blob.do_EC_BLR=0

                    self.del_spectral_component('EC_BLR',verbose=False)
                    self.EC_components_list.remove(EC_component)

            if EC_component=='DT':
                if self.get_spectral_component_by_name('DT', verbose=False) is None:
                    self._blob.do_EC_DT=0
                    self.del_spectral_component('EC_DT',verbose=False)
                    self.del_spectral_component('DT',verbose=False)
                    self.EC_components_list.remove(EC_component)

            if EC_component=='CMB':
                if self.get_spectral_component_by_name('CMB', verbose=False) is None:
                    self._blob.do_EC_CMB=0
                    self.del_spectral_component('EC_CMB',verbose=False)
                    self.EC_components_list.remove(EC_component)

            if EC_component=='CMB_stat':
                if self.get_spectral_component_by_name('CMB_stat', verbose=False) is None:
                    self._blob.do_EC_CMB_stat=0
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
                self._blob.do_EC_Disk=1
                if self.get_spectral_component_by_name('EC_Disk',verbose=False) is None:
                    self.add_spectral_component('EC_Disk')
                    self.EC_components_list.append(EC_component)

                if self.get_spectral_component_by_name('Disk',verbose=False) is None:
                    self.add_spectral_component('Disk')


            if EC_component=='BLR':
                self._blob.do_EC_BLR=1
                if self.get_spectral_component_by_name('EC_BLR',verbose=False) is None:
                    self.add_spectral_component('EC_BLR')
                    self.EC_components_list.append(EC_component)

                if self.get_spectral_component_by_name('Disk',verbose=False) is None:
                    self.add_spectral_component('Disk')
                    self.EC_components_list.append('Disk')

            if EC_component=='DT':
                self._blob.do_EC_DT=1
                if self.get_spectral_component_by_name('EC_DT',verbose=False) is None:
                    self.add_spectral_component('EC_DT')
                    self.EC_components_list.append(EC_component)
                if self.get_spectral_component_by_name('DT',verbose=False) is None:
                    self.add_spectral_component('DT')
                    self.EC_components_list.append('DT')

            if EC_component=='CMB':
                self._blob.do_EC_CMB=1
                if self.get_spectral_component_by_name('EC_CMB',verbose=False) is None:
                    self.add_spectral_component('EC_CMB')
                    self.EC_components_list.append(EC_component)

            if EC_component == 'CMB_stat':
                self._blob.do_EC_CMB_stat = 1
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
                pval=getattr(self._blob,pname)
                ptype=model_dic[key][0]
                vmin=model_dic[key][1]
                vmax=model_dic[key][2]
                punit=model_dic[key][3]

                froz=False

                if len(model_dic[key])>4:
                    froz=model_dic[key][4]


                self.parameters.add_par(JetParameter(self._blob,name=pname,par_type=ptype,val=pval,val_min=vmin,val_max=vmax,units=punit,frozen=froz))



    #PROTECTED MEMBER ACCESS
    #!! controlla i paramteri ....
    def get_electron_distribution_name(self):
        return self.electron_distribution._name



    def get_DL_cm(self):
        self.init_BlazarSED()
        return self._blob.dist



    def get_beaming(self,theta=None,bulk_factor=None,beaming=None):

        if theta is None:
            theta=self._blob.theta

        if bulk_factor is None:
            bulk_factor=self._blob.BulkFactor

        if beaming is None:
            beaming=self._blob.beam_obj

        BlazarSED.SetBeaming(self._blob)

        return self._blob.beam_obj

    def set_flag(self,flag):
        self._blob.STEM=flag

    def get_flag(self):
        return self._blob.STEM

    def get_path(self):
        return self._blob.path

    def set_path(self,path,clean_work_dir=True):
        set_str_attr(self._blob,'path',path)
        #set_str(self._blob.path,path)
        makedir(path,clean_work_dir=clean_work_dir)

    def set_SSC_mode(self,val):
        self._blob.do_SSC=val

    def get_SSC_mode(self):
        return self._blob.do_SSC

    def set_IC_mode(self,val):
        self._blob.do_IC=val

    def get_IC_mode(self):
        return self._blob.do_IC

    def set_sync_mode(self,val):
        self._blob.do_Sync=val

    def get_sync_mode(self):
        return self._blob.do_Sync

    def set_IC_nu_size(self,val):
        self._blob.nu_IC_size=val

    def get_IC_nu_size(self):
        return self._blob.nu_IC_size

    def set_seed_nu_size(self,val):
        self._blob.nu_seed_size=val

    def get_seed_nu_size(self):
        return self._blob.nu_seed_size

    def set_gamma_grid_size(self,val):
        self._blob.gamma_grid_size=val

    def get_gamma_grid_size(self):
        return self._blob.gamma_grid_size


    def set_verbosity(self,val):
        self._blob.verbose=val

    def get_verbosity(self):
        return  self._blob.verbose

    def debug_synch(self):
        print ("nu stop synch", self._blob.nu_stop_Sync)
        print ("nu stop synch ssc", self._blob.nu_stop_Sync_ssc)
        print ("ID MAX SYNCH", self._blob.NU_INT_STOP_Sync_SSC)

    def debug_SSC(self):
        print ("nu start SSC", self._blob.nu_start_SSC)
        print ("nu stop SSC", self._blob.nu_stop_SSC)
        print ("ID MAX SSC", self._blob.NU_INT_STOP_COMPTON_SSC)

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

        print ("-----------------------------------------------------------------------------------------")
        print ("model parameters for jet model:")
        print


        print ("electron distribution type = %s  "%(self._electron_distribution_name))

        self.parameters.show_pars()

        print ("-----------------------------------------------------------------------------------------")


    def plot_model(self,plot_obj=None,clean=False,autoscale=True,label=None,comp=None,sed_data=None):
        if plot_obj is None:
            plot_obj=Plot(sed_data=sed_data)


        if clean==True:
            plot_obj.clean_model_lines()

        line_style='-'

        if comp is not None:
            c = self.get_spectral_component_by_name(comp)
            if label is not None:
                comp_label = label
            else:
                comp_label = c.name
            plot_obj.add_model_plot(c.SED, autoscale=autoscale, line_style=line_style, label=comp_label)

        else:
            for c in self.spectral_components:
                comp_label = c.name
                if c.name=='Sum':
                    if label is not None:
                        comp_label=label

                plot_obj.add_model_plot(c.SED, autoscale=autoscale, line_style=line_style, label=comp_label)

        return plot_obj



    def init_BlazarSED(self):
        BlazarSED.Init(self._blob)




    def eval(self,init=True,fill_SED=True,nu=None,get_model=False,loglog=False,plot=None,label=None,phys_output=False):
        """
        Runs the BlazarSED  code for the current `JetModel` instance.

        :param init: (boolean), "defualt=True" initializes the BlazarSED code
            for the current `Jet` instance parameters values.

        """




        if init==True:

            BlazarSED.Init(self._blob)



        BlazarSED.Run_SED(self._blob)


        if phys_output==True:
            BlazarSED.EnergeticOutput(self._blob)

        nu_sed_sum,nuFnu_sed_sum= self.get_SED_points()
        #print('nu_sed_sum,nuFnu_sed_sum',nu_sed_sum,nuFnu_sed_sum)


        if fill_SED==True:

            self.SED.fill(nu=nu_sed_sum,nuFnu=nuFnu_sed_sum)

            for i in range(len(self.spectral_components)):

                #print ('fill name',self.spectral_components[i].name)
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

            size=self._blob.spec_array_size
            x=zeros(size)
            y=zeros(size)

            for i in range(size):
                x[i]=BlazarSED.get_freq_array(nu_ptr,self._blob,i)
                y[i]=BlazarSED.get_freq_array(nuFnu_ptr,self._blob,i)

                #print"%s %e %e"%(name,x[i],y[i])

            msk=y<self.flux_plot_lim


            y[msk]=self.flux_plot_lim



            if log_log==True:
                msk = y < 0
                y[msk] = self.flux_plot_lim

                #x=x[msk]
                #    y=y[msk]

                x=log10(x)
                y=log10(y)



            return x,y

        except:
            raise RuntimeError ('model evaluation failed in get_SED_points')


    def energetic_report(self,write_file=False,getstring=True,wd=None,name=None):

        _energetic = BlazarSED.EnergeticOutput(self._blob,0)
        _par_array=ModelParameterArray()

        _name = [i for i in _energetic.__class__.__dict__.keys() if i[:1] != '_']
        for _n in _name:
            if _n[0]=='L':
                par_type='Lum. rest. frme.'
                units='erg/s'
            if _n[0] == 'U':
                par_type = 'Energy dens.'
                units = 'erg/cm^3'
            if _n[0] == 'j':
                par_type = 'jet Lum.'
                units = 'erg/s'
            _par_array.add_par(ModelParameter(name=_n, val=getattr(_energetic, _n), units=units,par_type=par_type))

        print("-----------------------------------------------------------------------------------------")
        print("jet eneregetic report:")
        self._energetic_report = _par_array.show_pars(getstring=False)
        self._energetic_report=_par_array.show_pars(getstring=getstring)
        print("-----------------------------------------------------------------------------------------")

        if write_file==True:

            if wd is None:
                wd = self.wd

            if name is None:
                name = 'energetic_report_%s' % self.name + '.txt'

            outname = '%s/%s' % (wd, name)

            outfile = open(outname, 'w')

            for text in self._energetic_report:
                print(text, file=outfile)

            outfile.close()



    def get_SED_peak(self,peak_name=None,freq_range=None,log_log=False):

        if peak_name is not None and freq_range is not None:
            print ("either you provide peak_name or freq_range")
            raise ValueError
        elif   peak_name is not None:
            try:
                if log_log==False:
                    return getattr(self._blob,peak_name)
                else:
                     return np.log10(getattr(self._blob,peak_name) )
            except:
                print ("peak name %s, not found, check name"%peak_name)
                raise ValueError

        else:
            x,y=self.get_SED_points(log_log=log_log)
            msk1=x>freq_range[0]
            msk2=x<freq_range[1]
            #print x,freq_range,x[msk1*msk2]
            y_m= y[msk1*msk2].max()
            x_id= np.argmax(y[msk1*msk2])
            return x[msk1*msk2][x_id],y_m



def set_str(blob_attr,val):
    print ('set',blob_attr,'to',val)
    try:
        try:
            blob_attr = val
        except:
            blob_attr = val.encode('ascii')

    except Exception as e:
        raise RuntimeError(e)

def set_str_attr(obj,name,val):
    #print('set obj', obj,'name',name ,'to', val)
    try:

        try:
            setattr(obj, name,val)
        except:
            setattr(obj, name, val.encode('ascii'))
    except Exception as e:
        raise RuntimeError('error setting attr',name,'execption:',e)







        
