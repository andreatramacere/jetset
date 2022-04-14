__author__ = "Andrea Tramacere"


import json
import dill as pickle
import six
import numpy as np
import copy
import warnings

from .jet_spectral_components import JetSpecComponent, SpecCompList

from .model_parameters import ModelParameterArray, ModelParameter, _show_table
from .base_model import Model
from .output import makedir,WorkPlace
from  .plot_sedfit import PlotSED,plt
from .cosmo_tools import Cosmo
from .utils import safe_run,set_str_attr, old_model_warning, get_info, clean_var_name
from .jet_paramters import *
from .jet_emitters import *
from .jet_emitters_factory import EmittersFactory
from .jet_tools import *
from .mathkernel_helper import bessel_table_file_path


# on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

# if on_rtd is True:
#     try:
#         from .jetkernel import jetkernel as BlazarSED
#     except ImportError:
#         from .mock import jetkernel as BlazarSED
# else:

from .jetkernel import jetkernel as BlazarSED



__all__=['Jet','JetBase']


class JetBase(Model):
    """ JetBase class.
    This class allows to build a ``Jet`` model providing the interface to the
    C code, giving  full access to the physical parameters and
    providing the methods to run the code.
    A :class:`Jet` object  will store the
    the physical parameters in  the ::py:attr:`Jet.parameters`  that is :class:`.ModelParameterArray` class,
    i.e. a collection of :class:`JetParameter` objects.
    All the physical parameters are  also accessible as attributes of
    the  ::py:attr:`Jet.parameters`
    """
    def __repr__(self):
        return str(self.show_model())

    

    def __init__(self,
                 cosmo=None,
                 name='test',
                 emitters_type='electrons',
                 emitters_distribution='pl',
                 emitters_distribution_log_values=False,
                 beaming_expr='delta',
                 jet_workplace=None,
                 verbose=None,
                 nu_size=500,
                 clean_work_dir=True,
                 **keywords):

        """

        Parameters
        ----------
        cosmo
        name
        emitters_type
        emitters_distribution
        emitters_distribution_log_values
        beaming_expr
        jet_workplace
        verbose
        nu_size
        clean_work_dir
        """
        super(JetBase,self).__init__(  **keywords)




        if cosmo is not None:
            self.cosmo=cosmo
        else:
            self.cosmo= Cosmo()
        #print('cosmo', self.cosmo)
        self.name = clean_var_name(name)

        self.model_type='jet'
        self._emitters_type=emitters_type
        self._scale='lin-lin'

        self._blob = self.build_blob(verbose=verbose)
        self._static_spec_arr_grid_size = BlazarSED.static_spec_arr_grid_size
        self._nu_static_size = BlazarSED.static_spec_arr_size
        self.nu_size = nu_size
        self.nu_grid_size=self._get_nu_grid_size_blob()
        if jet_workplace is None:
            jet_workplace=WorkPlace()
            out_dir= jet_workplace.out_dir + '/' + self.name + '_jet_prod/'
        else:
            out_dir=jet_workplace.out_dir+'/'+self.name+'_jet_prod/'

        self.set_path(out_dir,clean_work_dir=clean_work_dir)

        self.set_flag(self.name)


        self._allowed_EC_components_list=['EC_BLR',
                                          'DT',
                                          'EC_DT',
                                          'Star',
                                          'EC_Star',
                                          #CMB',
                                          'EC_CMB',
                                          'Disk',
                                          'Disk_MultiBB',
                                          'Disk_Mono',
                                          'EC_Disk',
                                          'All']

        self._allwed_disk_type =['BB', 'Mono', 'MultiBB']
        self.EC_components_list =[]
        self.spectral_components_list=[]


        self.spectral_components= SpecCompList(self.spectral_components_list)

        self.add_basic_components()
        self.SED=self.get_spectral_component_by_name('Sum').SED

        self.parameters = JetModelParameterArray(model=self)

        self._emitting_region_dic = None
        self._electron_distribution_dic= None
        self._external_photon_fields_dic= None
        self._original_emitters_distr = None

        self._setup(emitters_distribution,emitters_distribution_log_values,beaming_expr,emitters_type)


    def _setup(self, emitters_distribution, emitters_distribution_log_values, beaming_expr, emitters_type):
        self.EC_components_list = []
        self.spectral_components_list = []

        self.spectral_components = SpecCompList(self.spectral_components_list)
        self.add_basic_components()

        self.SED = self.get_spectral_component_by_name('Sum').SED

        self.parameters = JetModelParameterArray(model=self)

        self._emitting_region_dic = None
        self._emitters_distribution_dic = None
        self._external_photon_fields_dic = None



        self._blob.IC_adaptive_e_binning = 0
        self._blob.do_IC_down_scattering = 0

        self.set_emitting_region(beaming_expr)

        self.flux_plot_lim=1E-30
        self.set_emiss_lim(1E-120)

        self._IC_states = {}
        self._IC_states['on'] = 1
        self._IC_states['off'] = 0

        self._external_field_transf = {}
        self._external_field_transf['blob'] = 0
        self._external_field_transf['disk'] = 1
        self._jetkernel_interp='linear'
        self.set_emitters_distribution(emitters_distribution, emitters_distribution_log_values, emitters_type,
                                       init=False)
        self.set_blob()


    def __getstate__(self):
        return  self._serialize_model()

    def __setstate__(self,state):
        self.__init__()
        self._decode_model(state)

    def _serialize_model(self):

        _model = {}
        _model['version']=get_info()['version']
        _model['name'] = self.name
        if isinstance(self.emitters_distribution,JetkernelEmittersDistribution):
            _model['emitters_distribution'] = self._emitters_distribution_name
            _model['emitters_distribution_log_values'] = self._emitters_distribution_log_values
            _model['emitters_type'] = self._emitters_type
            _model['emitters_distribution_class']='JetkernelEmittersDistribution'
        elif isinstance(self.emitters_distribution,EmittersDistribution):
            self._original_emitters_distr._copy_from_jet(self)
            _model['custom_emitters_distribution'] =self._original_emitters_distr
            _model['emitters_distribution_class'] = 'EmittersDistribution'
        else:
            raise  RuntimeError('emitters distribuion type not valid',type(self._emitters_distribution))

        if hasattr(self,'T_esc_e_second'):
            _model['T_esc_e_second']=self.T_esc_e_second

        _model['beaming_expr'] = self._beaming_expr
        _model['spectral_components_name'] = self.get_spectral_component_names_list()
        _model['spectral_components_state'] = [c.state for c in self.spectral_components_list]
        _model['EC_components_name'] = self.EC_components_list
        _model['basic_components_name'] = self.basic_components_list
        _model['cosmo'] = self.cosmo

        _model['pars'] = {}
        _model['pars']=self.parameters._serialize_pars()

        _model['internal_pars'] = {}
        _model['internal_pars']['nu_size'] = self.nu_size
        _model['internal_pars']['nu_seed_size'] = self.nu_seed_size
        _model['internal_pars']['gamma_grid_size'] = self.gamma_grid_size
        _model['internal_pars']['IC_nu_size']=self.IC_nu_size
        _model['internal_pars']['nu_min']=self.nu_min
        _model['internal_pars']['nu_max']=self.nu_max
        _model['internal_pars']['Norm_distr'] = self.Norm_distr
        return _model

    def save_model(self,file_name):
        pickle.dump(self._serialize_model(), open(file_name, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

    @classmethod
    @safe_run
    def load_model(cls, file_name):
        try:
            _model = pickle.load(open(file_name, "rb"))
        except Exception as e:
            raise RuntimeError('The model you loaded is not valid please check the file name', e)

        try:
            jet = cls(name='no_name')
            jet._decode_model(_model)
            jet._fix_par_dep_on_load()
            jet.show_pars()
            jet.eval()
            return jet
        except Exception as e:
            if 'version' in _model.keys():
                v=_model['version']
            else:
                v='unknown(<1.1.2)'
            msg = "trying  to load a model saved with jetset version %s \n "%v
            msg += "currenlty using jetset version %s \n "%get_info()['version']
            msg += "caused the following problem:\n %s\n"%(repr(e))

            raise RuntimeError (msg)




    def _decode_model(self,_model):

        if 'version' in _model.keys():
            self._set_version(_model['version'])
        else:
            self._set_version(v='unknown(<1.1.2)')

        self.cosmo = _model['cosmo']
        self.model_type = 'jet'
        self.name = _model['name']
        self.set_blob()
        self.parameters = JetModelParameterArray(model=self)


        if 'emitters_type' in _model.keys():
            emitters_type=str(_model['emitters_type'])
        else:
            emitters_type = 'electrons'

        if 'electron_distribution' in _model.keys():
            _v=_model['electron_distribution']
            del(_model['electron_distribution'])
            _model['emitters_distribution']=_v

        if 'electron_distribution_log_values' in _model.keys():
            _v=_model['electron_distribution_log_values']
            del(_model['electron_distribution_log_values'])
            _model['emitters_distribution_log_values']=_v

        if _model['emitters_distribution_class'] == 'JetkernelEmittersDistribution':
            self.set_emitters_distribution(distr=_model['emitters_distribution'],
                                           log_values=_model['emitters_distribution_log_values'],
                                           emitters_type=emitters_type,
                                           init=False)
        elif _model['emitters_distribution_class'] == 'EmittersDistribution':
            self.set_emitters_distribution(distr=_model['custom_emitters_distribution'], init=False)
        else:
            raise RuntimeError('emitters distribuion type not valid', type(self._emitters_distribution))

        for c in self.basic_components_list:
            if c not in _model['basic_components_name']:
                self.del_spectral_component(c)

        if self.emitters_distribution.emitters_type == 'protons':

            self.add_pp_gamma_component()
            self.add_pp_neutrino_component()
            self.add_bremss_ep_component()
            self.T_esc_e_second=_model['T_esc_e_second']

        #for k in _model['pars'].keys():
            #print ('-->',k,_model['pars'][k])
        if 'disk_type' in _model['pars'].keys():
            disk_type=_model['pars']['disk_type']['val']
        else:
            disk_type=None
        self.add_EC_component(_model['EC_components_name'],disk_type=disk_type)

        for ID, c in enumerate(_model['spectral_components_name']):
            comp = getattr(self.spectral_components, c)
            if comp._state_dict != {}:
                comp.state = _model['spectral_components_state'][ID]

        self.SED = self.get_spectral_component_by_name('Sum').SED

        self.set_emitting_region(str(_model['beaming_expr']))
        #self.set_electron_distribution(str(_model['electron_distribution']))

        _par_dict = _model['pars']
        _non_user_dict={}
        _user_dict={}
        for k, v in _model['pars'].items():
            if v['par_type'] == 'user_defined':
                _user_dict[k]=v
            else:
                _non_user_dict[k]=v
        self.parameters._decode_pars(_non_user_dict)

        for k, v in _user_dict.items():
            #print('==>',v)
            v['name']=k
            self.parameters.add_par(ModelParameter(**v))

        _par_dict = _model['internal_pars']
        for k in _par_dict.keys():
            #print ('set', k,_par_dict[k])
            setattr(self,k,_par_dict[str(k)])
            #self.set_par(par_name=str(k), val=_par_dict[str(k)])

    @classmethod
    def load_old_model(cls, file_name):

        old_model_warning()
        jet = cls()
        with open(file_name, 'r') as infile:
            _model = json.load(infile)

        # print ('_model',_model)

        jet.model_type = 'jet'

        jet.init_BlazarSED()

        jet.parameters = JetModelParameterArray(model=jet)

        if 'electron_distribution' in _model.keys():
            _v = _model['electron_distribution']
            del (_model['electron_distribution'])
            _model['emitters_distribution'] = _v

        if 'electron_distribution_log_values' in _model.keys():
            _v = _model['electron_distribution_log_values']
            del (_model['electron_distribution_log_values'])
            _model['emitters_distribution_log_values'] = _v

        jet.set_emitters_distribution(name=str(_model['emitters_distribution']),
                                       log_values=_model['emitters_distribution_log_values'],
                                       emitters_type='electrons',
                                       init=False)

        for c in jet.basic_components_list:
            if c not in _model['basic_components_name']:
                jet.del_spectral_component(c)


        jet.add_EC_component(_model['EC_components_name'])

        for ID, c in enumerate(_model['spectral_components_name']):
            comp = getattr(jet.spectral_components, c)
            if comp._state_dict != {}:
                comp.state = _model['spectral_components_state'][ID]

        jet.SED = jet.get_spectral_component_by_name('Sum').SED

        jet.set_emitting_region(str(_model['beaming_expr']))

        _par_dict = _model['pars']
        jet.show_pars()
        for k in _par_dict.keys():
            # print ('set', k,_par_dict[k])
            jet.set_par(par_name=str(k), val=_par_dict[str(k)])

        jet.eval()
        return jet



    def build_blob(self, verbose=None):

        blob = BlazarSED.MakeBlob()

        blob.x_Bessel_min = 1E-17
        blob.x_Bessel_max = 7.2E2

        blob.x_ave_Bessel_min = 1E-16
        blob.x_ave_Bessel_max = 3.5E2

        blob.log_x_Bessel_min = np.log10( blob.x_Bessel_min)
        blob.log_x_Bessel_max = np.log10( blob.x_Bessel_max)


        blob.log_x_ave_Bessel_min = np.log10( blob.x_ave_Bessel_min)
        blob.log_x_ave_Bessel_max = np.log10( blob.x_ave_Bessel_max)




        F_Sync_x_ptr = getattr(blob, 'F_Sync_x')

        F_Sync_y_ptr = getattr(blob,  'F_Sync_y')

        F_ave_Sync_x_ptr = getattr(blob,  'F_ave_Sync_x')

        F_ave_Sync_y_ptr = getattr(blob, 'F_ave_Sync_y')

        log_F_Sync_x_ptr = getattr(blob, 'log_F_Sync_x')

        log_F_Sync_y_ptr = getattr(blob, 'log_F_Sync_y')

        log_F_ave_Sync_x_ptr = getattr(blob, 'log_F_ave_Sync_x')

        log_F_ave_Sync_y_ptr = getattr(blob, 'log_F_ave_Sync_y')


        d = np.genfromtxt(bessel_table_file_path)
        log_F_Sync_x=np.log10(d[:,0])
        log_F_Sync_y=np.log10(d[:,1])
        log_F_ave_Sync_x=np.log10(d[:,2])
        log_F_ave_Sync_y=np.log10(d[:,3])

        for ID, l in enumerate(d):
            BlazarSED.set_bessel_table(F_Sync_x_ptr,blob, l[0], ID)
            BlazarSED.set_bessel_table(F_Sync_y_ptr, blob,  l[1], ID)
            BlazarSED.set_bessel_table(F_ave_Sync_x_ptr, blob,  l[2], ID)
            BlazarSED.set_bessel_table(F_ave_Sync_y_ptr, blob,  l[3], ID)

            BlazarSED.set_bessel_table(log_F_Sync_x_ptr, blob, log_F_Sync_x[ID], ID)
            BlazarSED.set_bessel_table(log_F_Sync_y_ptr, blob, log_F_Sync_y[ID], ID)
            BlazarSED.set_bessel_table(log_F_ave_Sync_x_ptr, blob, log_F_ave_Sync_x[ID], ID)
            BlazarSED.set_bessel_table(log_F_ave_Sync_y_ptr, blob, log_F_ave_Sync_y[ID], ID)

        blob.BESSEL_TABLE_DONE=1

        if verbose is None:
            blob.verbose = 0
        else:
            blob.verbose = verbose

        set_str_attr(blob, 'path', './')

        set_str_attr(blob, 'MODE', 'custom')

        blob.gamma_grid_size = 200

        blob.nu_IC_size = 100
        blob.nu_seed_size = 100

        blob.nu_grid_size= 1000

        blob.do_Sync = 2

        blob.do_SSC = 1

        blob.R = 5.0e15

        blob.B = 0.1

        blob.z_cosm = 0.1

        blob.BulkFactor = 10
        blob.theta = 0.1

        blob.N = 100


        blob.NH_pp = 1

        blob.L_Disk = 1E45

        blob.L_DT = 1E45

        blob.gmin = 2

        blob.gmax = 1e6

        blob.nu_start_Sync = 1e6
        blob.nu_stop_Sync = 1e20

        blob.nu_start_SSC = 1e14
        blob.nu_stop_SSC = 1e30

        blob.nu_start_grid = 1e6
        blob.nu_stop_grid = 1e30

        return blob

    def set_emitting_region(self,beaming_expr):

        if  self._emitting_region_dic is not None:
            self.del_par_from_dic(self._emitting_region_dic)

        self._beaming_expr=beaming_expr

        set_str_attr(self._blob,'BEAMING_EXPR',beaming_expr)

        self._emitting_region_dic=build_emitting_region_dic(self.cosmo,beaming_expr=beaming_expr)

        self.parameters.add_par_from_dict(self._emitting_region_dic,self,'_blob',JetParameter)

    @property
    def IC_adaptive_e_binning(self,):
        return np.int(self._blob.IC_adaptive_e_binning)

    @IC_adaptive_e_binning.setter
    def IC_adaptive_e_binning(self,state):
        if type(state) == bool:
            pass
        else:
            raise RuntimeError('state has to be boolean')
        self._blob.IC_adaptive_e_binning=np.int(state)

    @staticmethod
    def available_emitters_distributions():
        EmittersFactory.available_distributions()

    def set_emitters_distribution(self, distr=None, log_values=False, emitters_type='electrons', init=True):

        if init is True:
            self.set_blob()
        self._emitters_distribution_log_values = log_values

        if self._emitters_distribution_dic is not None:
            self.del_par_from_dic(self._emitters_distribution_dic)

        if isinstance(distr, ArrayDistribution):
            self._emitters_distribution_name = 'from_array'
            self.emitters_distribution = JetkernelEmittersDistribution.from_array(self, distr, emitters_type=emitters_type)

        elif isinstance(distr, JetkernelEmittersDistribution):
            self.emitters_distribution = distr

            self._emitters_distribution_name = self.emitters_distribution.name
            self._emitters_distribution_dic = self.emitters_distribution._parameters_dict

            self.parameters.add_par_from_dict(self._emitters_distribution_dic,self,'_blob',JetParameter)


        elif isinstance(distr, EmittersDistribution):
            self._original_emitters_distr = distr
            self.emitters_distribution = copy.deepcopy(distr)
            self._update_emitters_pars_dependence()
            self.emitters_distribution.set_jet(self)
            self.emitters_distribution._update_parameters_dict()
            self._emitters_distribution_name = self.emitters_distribution.name
            self._emitters_distribution_dic = self.emitters_distribution._parameters_dict
            self.parameters.add_par_from_dict(self._emitters_distribution_dic, self, '_blob', JetParameter)
            self._attach_pars_to_jet(preserve_value_emitters=True)

            self.emitters_distribution.update()
        elif isinstance(distr, str):
            nf=EmittersFactory()
            self.emitters_distribution = nf.create_emitters(distr, log_values=log_values, emitters_type=emitters_type)
            self._original_emitters_distr = copy.deepcopy(self.emitters_distribution)
            self.emitters_distribution.set_jet(self)
            self.emitters_distribution._update_parameters_dict()
            self._emitters_distribution_name = self.emitters_distribution.name
            self._emitters_distribution_dic = self.emitters_distribution._parameters_dict
            self.parameters.add_par_from_dict(self._emitters_distribution_dic, self, '_blob', JetParameter)
            self._attach_pars_to_jet(preserve_value_emitters=True)
            self._update_emitters_pars_dependence()
            self.emitters_distribution.update()
        else:
            raise RuntimeError('distr',type(distr),'not valid should be a string or an',type(EmittersDistribution),'instance')

    def _attach_pars_to_jet(self, preserve_value_emitters=False):

        v_dict={}
        for par in self.emitters_distribution.parameters.par_array[::]:
            self.emitters_distribution.parameters.del_par(par)
            v_dict[par.name]=par.val
        for k in self._emitters_distribution_dic.keys():
            par = self.parameters.get_par_by_name(k)
            if preserve_value_emitters is True:
                par.set(val=v_dict[k],skip_dep_par_warning=True)
            self.emitters_distribution.parameters.add_par(par)


    def _update_emitters_pars_dependence(self):
        for par in self.emitters_distribution.parameters.par_array:
            if par.immutable is True:
                warnings.warn('parameter dependence has to be reassigned after emitters distribution is assigned to a jet, now will be reset no dep')
                self.emitters_distribution.parameters.reset_dependencies()
                break

    def get_emitters_distribution_name(self):
        return self.emitters_distribution.name


    def show_spectral_components(self):

        print ("Spectral components for Jet model:%s"%(self.name))

        for comp in self.spectral_components_list:

            print ("comp: %s "%(comp.name))

        print()


    def get_spectral_component_by_name(self,name,verbose=True):
        for i in range(len(self.spectral_components_list)):
            if self.spectral_components_list[i].name==name:
                return self.spectral_components_list[i]
        else:
            if verbose==True:
                print ("no spectral components with name %s found"%name)
                print ("names in array are:")
                self.show_spectral_components()

            return None

    def list_spectral_components(self):
        for i in range(len(self.spectral_components_list)):
            print (self.spectral_components_list[i].name)

    def get_spectral_component_names_list(self):
        _l=[]
        for i in range(len(self.spectral_components_list)):
            _l.append(self.spectral_components_list[i].name)
        return _l

    def del_spectral_component(self,name):
        print('deleting spectral component', name)
        if name in self.EC_components_list:
            self.del_EC_component(name)
        else:
            self._del_spectral_component(name)
        if name in self.basic_components_list:
            self.basic_components_list.remove(name)

    def _del_spectral_component(self, name, verbose=True):

        comp=self.get_spectral_component_by_name(name,verbose=verbose)

        if comp._state_dict!={}:

            comp.state='off'

        if comp is not None:

            self.spectral_components_list.remove(comp)

            self.spectral_components.__delattr__(name)



    def _add_spectral_component(self, name, var_name=None, state_dict=None,state=None):

        self.spectral_components_list.append(
            JetSpecComponent(self, name, self._blob, var_name=var_name, state_dict=state_dict, state=state))
        setattr(self.spectral_components,name,self.spectral_components_list[-1])

    def _update_spectral_components(self):
        _l=[]
        for ID,s in enumerate(self.spectral_components_list):
            self.spectral_components_list[ID]= JetSpecComponent(self, s.name, self._blob, var_name=s._var_name, state_dict=s._state_dict, state=s.state)
            setattr(self.spectral_components, self.spectral_components_list[ID].name, self.spectral_components_list[ID])


    def add_basic_components(self):
        self.basic_components_list=['Sum','Sync','SSC']

        self._add_spectral_component('Sum')
        self._add_spectral_component('Sync', var_name='do_Sync', state_dict=dict((('on', 1), ('off', 0), ('self-abs', 2))),state='self-abs')

        self._add_spectral_component('SSC', var_name='do_SSC', state_dict=dict((('on', 1), ('off', 0))))


    def add_sync_component(self,state='self-abs'):
        self._add_spectral_component('Sync', var_name='do_Sync',
                                     state_dict=dict((('on', 1), ('off', 0), ('self-abs', 2))), state=state)

    def add_SSC_component(self,state='on'):
        self._add_spectral_component('SSC', var_name='do_SSC', state_dict=dict((('on', 1), ('off', 0))),state=state)

    def del_EC_component(self,EC_components_list, disk_type='BB'):
        if isinstance(EC_components_list, six.string_types):
            EC_components_list = [EC_components_list]

        if 'All' in EC_components_list:
            EC_components_list=self._allowed_EC_components_list[::]


        for EC_component in EC_components_list:


            if EC_component not in self._allowed_EC_components_list:
                raise RuntimeError("EC_component %s not allowed" % EC_component, "please choose among ",
                                   self._allowed_EC_components_list)

            if EC_component=='Disk':
                if self.get_spectral_component_by_name('Disk', verbose=False) is not None:
                    self._del_spectral_component('Disk', verbose=False)
                    self.EC_components_list.remove('Disk')

                if self.get_spectral_component_by_name('EC_Disk',verbose=False) is not None:
                    self._del_spectral_component('EC_Disk')
                    self._blob.do_EC_Disk = 0
                    self.EC_components_list.remove('EC_Disk')

                if self.get_spectral_component_by_name('EC_BLR', verbose=False) is not None:
                    self._blob.do_EC_BLR=0
                    self._del_spectral_component('EC_BLR', verbose=False)
                    self.EC_components_list.remove('EC_BLR')


            if EC_component=='EC_Disk':
                if self.get_spectral_component_by_name('EC_Disk', verbose=False) is not None:
                    self._blob.do_EC_Disk=0
                    self._del_spectral_component('EC_Disk', verbose=False)
                    self.EC_components_list.remove('EC_Disk')


            if EC_component=='EC_BLR':
                if self.get_spectral_component_by_name('EC_BLR', verbose=False) is not None:
                    self._blob.do_EC_BLR=0
                    self._del_spectral_component('EC_BLR', verbose=False)
                    self.EC_components_list.remove('EC_BLR')

            if EC_component=='DT':
                if self.get_spectral_component_by_name('DT', verbose=False) is not None:
                    self._del_spectral_component('DT', verbose=False)
                    self.EC_components_list.remove('DT')
                if self.get_spectral_component_by_name('EC_DT', verbose=False) is not None:
                    self._blob.do_EC_DT = 0
                    self._del_spectral_component('EC_DT', verbose=False)
                    self.EC_components_list.remove('EC_DT')

            if EC_component=='EC_DT':
                if self.get_spectral_component_by_name('EC_DT', verbose=False) is not None:
                    self._blob.do_EC_DT=0
                    self._del_spectral_component('EC_DT', verbose=False)
                    self.EC_components_list.remove('EC_DT')

            if EC_component=='EC_CMB':
                if self.get_spectral_component_by_name('EC_CMB', verbose=False) is not None:
                    self._blob.do_EC_CMB=0
                    self._del_spectral_component('EC_CMB', verbose=False)
                    self.EC_components_list.remove('EC_CMB')

            if EC_component=='Star':
                if self.get_spectral_component_by_name('Star', verbose=False) is not None:
                    self._blob.do_star=0
                    self._del_spectral_component('Star', verbose=False)
                    self.EC_components_list.remove('Star')

        self.del_par_from_dic(build_ExtFields_dic(EC_components_list,disk_type))



    def add_EC_component(self,EC_components_list=[],disk_type='BB'):

        if disk_type is not None:
            if disk_type not in self._allwed_disk_type:
                raise RuntimeError('disk type',disk_type,'not in allwowed', self._allwed_disk_type)
        else:
            if self.get_par_by_name('disk_type') is not None:
                disk_type=self.get_par_by_name('disk_type').val


        if isinstance(EC_components_list, six.string_types):
            EC_components_list=[EC_components_list]

        if 'All' in EC_components_list:
            EC_components_list=self._allowed_EC_components_list[::]

        for EC_component in EC_components_list:
            if EC_component not in self._allowed_EC_components_list:
                raise RuntimeError("EC_component %s not allowed" % EC_component, "please choose among ",
                                   self._allowed_EC_components_list)

            if EC_component == 'Disk':
                if self.get_spectral_component_by_name('Disk',verbose=False) is None:
                    self._add_spectral_component('Disk',var_name='do_Disk', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('Disk')

            if EC_component=='EC_Disk':
                #self._blob.do_EC_Disk=1
                if self.get_spectral_component_by_name('EC_Disk',verbose=False) is None:
                    self._add_spectral_component('EC_Disk', var_name='do_EC_Disk', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('EC_Disk')

                if self.get_spectral_component_by_name('Disk',verbose=False) is None:
                    self._add_spectral_component('Disk',var_name='do_Disk', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('Disk')

            if EC_component=='EC_BLR':
                #self._blob.do_EC_BLR=1
                if self.get_spectral_component_by_name('EC_BLR',verbose=False) is None:
                    self._add_spectral_component('EC_BLR', var_name='do_EC_BLR', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('EC_BLR')

                if self.get_spectral_component_by_name('Disk',verbose=False) is None:
                    self._add_spectral_component('Disk',var_name='do_Disk', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('Disk')


            if EC_component == 'DT':
                if self.get_spectral_component_by_name('DT',verbose=False) is None:
                    self._add_spectral_component('DT',var_name='do_DT', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('DT')

                if self.get_spectral_component_by_name('Disk',verbose=False) is None:
                    self._add_spectral_component('Disk',var_name='do_Disk', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('Disk')



            if EC_component == 'Star':
                if self.get_spectral_component_by_name('Star',verbose=False) is None:
                    self._add_spectral_component('Star',var_name='do_Star', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('Star')




            if EC_component=='EC_DT':
                #self._blob.do_EC_DT=1
                if self.get_spectral_component_by_name('EC_DT',verbose=False) is None:
                    self._add_spectral_component('EC_DT', var_name='do_EC_DT', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('EC_DT')

                if self.get_spectral_component_by_name('DT',verbose=False) is None:
                    self._add_spectral_component('DT',var_name='do_DT', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('DT')

                if self.get_spectral_component_by_name('Disk',verbose=False) is None:
                    self._add_spectral_component('Disk',var_name='do_Disk', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('Disk')

            if EC_component=='EC_CMB':
                #self._blob.do_EC_CMB=1
                if self.get_spectral_component_by_name('EC_CMB',verbose=False) is None:
                    self._add_spectral_component('EC_CMB', var_name='do_EC_CMB', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('EC_CMB')

        #IF disk_type is already a parameter it has to be updated here to make it effective
        self.parameters.add_par_from_dict(build_ExtFields_dic(self.EC_components_list, disk_type),self,'_blob',JetParameter)
        #print('--> disk_type',disk_type)
        if disk_type is not None:

            #remove the old disk type parameters
            if  self.parameters.get_par_by_name('disk_type') is not None:
                self.del_par_from_dic(build_ExtFields_dic(['Disk'],self.parameters.get_par_by_name('disk_type').val))
            else:
                self.EC_components_list.append('Disk')
            self.parameters.add_par_from_dict(build_ExtFields_dic(self.EC_components_list, disk_type),self,'_blob',JetParameter)
            self.parameters.disk_type.val = disk_type



    def del_par_from_dic(self,model_dic):
        """
        """
        for key in model_dic.keys():

            par=self.parameters.get_par_by_name(key)

            if par is not None:
                self.parameters.del_par(par)



    def get_DL_cm(self,eval=False):

        if eval is True:
            self.set_blob()

        if self.cosmo._c is not None:
            return self.cosmo.get_DL_cm(self.parameters.z_cosm.val)
        else:
            return self.cosmo.get_DL_cm()


    @safe_run
    def get_beaming(self,):

        BlazarSED.SetBeaming(self._blob)
        return self._blob.beam_obj


    def set_flag(self,flag):
        self._blob.STEM=flag

    def get_flag(self):
        return self._blob.STEM

    def get_path(self):
        return self._blob.path

    def set_path(self,path,clean_work_dir=True):
        if path.endswith('/'):
            pass
        else:
            path+='/'

        set_str_attr(self._blob,'path',path)
        makedir(path,clean_work_dir=clean_work_dir)

    def get_IC_mode(self):
        return dict(map(reversed, self._IC_states.items()))[self._blob.do_IC]



    def set_external_field_transf(self,val):
        if val not in self._external_field_transf.keys():
            raise RuntimeError('val',val,'not in allowed values',self._external_field_transf.keys())

        self._blob.EC_stat=self._external_field_transf[val]

    def get_external_field_transf(self):
        return dict(map(reversed, self._external_field_transf.items()))[self._blob.EC_stat]

    def set_emiss_lim(self,val):
        self._blob.emiss_lim=val

    def get_emiss_lim(self):
        return self._blob.emiss_lim


    @property
    def IC_nu_size(self):
        return self._blob.nu_IC_size

    @IC_nu_size.setter
    def IC_nu_size(self, val):
        self.set_IC_nu_size(val)

    def get_IC_nu_size(self):
        return self._blob.nu_IC_size

    def set_IC_nu_size(self, val):
        if val > self._nu_static_size:
            raise RuntimeError('value can not exceed',self._nu_static_size)
        self._blob.nu_IC_size = val

    @property
    def nu_seed_size(self):
        return self._blob.nu_seed_size

    def get_seed_nu_size(self):
        return self._blob.nu_seed_size

    @nu_seed_size.setter
    def nu_seed_size(self,val):
        self.set_seed_nu_size(val)

    def set_seed_nu_size(self,val):
        if val>self._nu_static_size:
            raise RuntimeError('value can not exceed',self._nu_static_size)
        self._blob.nu_seed_size=val



    def set_gamma_grid_size(self,val):
        self.emitters_distribution.set_grid_size(gamma_grid_size=val)

    @property
    def gamma_grid_size(self):
        return self._blob.gamma_grid_size

    @gamma_grid_size.setter
    def gamma_grid_size(self,val):
        self.emitters_distribution.set_grid_size(gamma_grid_size=val)

    @property
    def nu_min(self):
        return self._get_nu_min_grid()

    @nu_min.setter
    def nu_min(self, val):
        if hasattr(self, '_blob'):
            self._set_nu_min_grid(val)

    def _get_nu_min_grid(self):
        return  self._blob.nu_start_grid


    def _set_nu_min_grid(self, val):
        self._blob.nu_start_grid=val

    @property
    def Norm_distr(self):
        if hasattr(self,'emitters_distribution'):
            if self.emitters_distribution._user_defined is False:
                return self._blob.Norm_distr
            else:
                return self.emitters_distribution.normalize
        else:
            return None


    @Norm_distr.setter
    def Norm_distr(self, val):
        if hasattr(self, 'emitters_distribution'):

                if val == 1 or val is True:
                    if self.emitters_distribution._user_defined is False:
                        self._blob.Norm_distr = 1
                    else:
                        self.emitters_distribution.normalize = val
                    self.parameters.N.par_type='emitters_density'
                    self.parameters.N.units ='1/cm3'
                elif val == 0 or val is False:
                    if self.emitters_distribution._user_defined is False:
                        self._blob.Norm_distr = 0
                    else:
                        self.emitters_distribution.normalize = val
                        self.parameters.N.par_type = 'scaling_factor'
                else:
                    raise RuntimeError('value', val, 'not allowed, allowed 0/1 or False/True')




    def switch_Norm_distr_ON(self):
        self.Norm_distr=True

    def switch_Norm_distr_OFF(self):
        self.Norm_distr = False

    @property
    def nu_max(self):
        return self._get_nu_max_grid()

    @nu_max.setter
    def nu_max(self, val):
        if hasattr(self, '_blob'):
            self._set_nu_max_grid(val)

    def _set_nu_max_grid(self, val):
        self._blob.nu_stop_grid=val

    def _get_nu_max_grid(self):
        return  self._blob.nu_stop_grid

    @property
    def nu_size(self):
        return self._nu_size

    @nu_size.setter
    def nu_size(self, size):
        self._nu_size = size

    def set_nu_grid_size(self, val):
        self._set_nu_grid_size_blob(val)

    @property
    def nu_grid_size(self):
        return self._get_nu_grid_size_blob()

    @nu_grid_size.setter
    def nu_grid_size(self, val):
        self._set_nu_grid_size_blob(val)


    def _set_nu_grid_size_blob(self, val):
        if val < 100:
            val = 100
        if val > self._static_spec_arr_grid_size:
            raise RuntimeError('value can not exceed', self._nu_static_size)
        self._blob.nu_grid_size=val

    def _get_nu_grid_size_blob(self):
        return  self._blob.nu_grid_size



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


    def show_emitters_distribution(self):
        print('-'*80)
        print('%s distribution:'%self._emitters_type)
        print(" type: %s  " % (self._emitters_distribution_name))
        print(" gamma energy grid size: ", self.gamma_grid_size)
        print(" gmin grid : %e" % self._blob.
              gmin_griglia)
        print(" gmax grid : %e" % self._blob.gmax_griglia)
        print(" normalization ", self.Norm_distr)
        print(" log-values ", self._emitters_distribution_log_values)
        print('')
        self.parameters.par_array.sort(key=lambda x: x.name, reverse=False)

        self.parameters.show_pars(names_list=self._emitters_distribution_dic.keys())

    def show_model(self):
        """
        shortcut to :class:`ModelParametersArray.show_pars` method
        shows all the paramters in the model

        """
        print("")
        print('-'*80)

        print ("jet model description")
        print('-'*80)
        print("name: %s  " % (self.name))
        print('')
        print('%s distribution:'%self._emitters_type)
        print(" type: %s  " % (self._emitters_distribution_name))
        print (" gamma energy grid size: ",self.gamma_grid_size)
        print (" gmin grid : %e"%self._blob.gmin_griglia)
        print (" gmax grid : %e"%self._blob.gmax_griglia)
        print(" normalization ", self.Norm_distr)
        print(" log-values ", self._emitters_distribution_log_values)
        print('')
        if 'Disk' in self.EC_components_list:
            print('accretion disk:')
            BlazarSED.set_Disk(self._blob)
            print(' disk Type: %s'%self.parameters.get_par_by_name('disk_type').val)
            print(' L disk: %e (erg/s)'%self.parameters.get_par_by_name('L_Disk').val)
            print(' T disk: %e (K)'%self._blob.T_Disk)
            print(' nu peak disk: %e (Hz)'%BlazarSED.eval_nu_peak_Disk(self._blob.T_Disk))
            if self.parameters.get_par_by_name('disk_type').val == 'MultiBB':
                print(' Sw radius %e (cm)'%self._blob.R_Sw)
                print(' L Edd. %e (erg/s)'%self._blob.L_Edd)
                yr=86400*365
                print(' accr_rate: %e (M_sun/yr)'%(yr*self._blob.accr_rate/BlazarSED.m_sun))
                print(' accr_rate Edd.: %e (M_sun/yr)'%(yr*self._blob.accr_Edd/BlazarSED.m_sun))


            print('')
        print('radiative fields:')
        print (" seed photons grid size: ", self.nu_seed_size)
        print (" IC emission grid size: ", self.get_IC_nu_size())
        print (' source emissivity lower bound :  %e' % self._blob.emiss_lim)
        print (' spectral components:')
        for _s in self.spectral_components_list:
            print("   name:%s,"%_s.name, 'state:', _s.state)
        print('external fields transformation method:', self.get_external_field_transf())
        print ('')
        print ('SED info:')
        print (' nu grid size jetkernel: %d' % self.nu_grid_size)
        print (' nu size: %d' % self.nu_size)
        print (' nu mix (Hz): %e' % self._get_nu_min_grid())
        print (' nu max (Hz): %e' % self._get_nu_max_grid())
        print('')
        print('flux plot lower bound   :  %e' % self.flux_plot_lim)
        print('')
        print('-'*80)
        self.show_pars()
        print('-' * 80)

    def plot_model(self,plot_obj=None,clean=False,label=None,comp=None,sed_data=None,color=None,auto_label=True,line_style='-',frame='obs', density=False):
        plot_obj=self._set_up_plot(plot_obj,sed_data,frame,density)

        if clean is True:
            plot_obj.clean_model_lines()

        if comp is not None:
            c = self.get_spectral_component_by_name(comp)

            if label is not None:
                comp_label = label
            elif label is None and auto_label is True:
                comp_label = c.name
            else:
                comp_label=None
            if c.state!='off':
                plot_obj.add_model_plot(c.SED, line_style=line_style, label=comp_label,flim=self.flux_plot_lim,color=color,auto_label=auto_label, density=density,update=False, frame=frame)

        else:
            for c in self.spectral_components_list:
                comp_label = c.name
                if auto_label is not True:
                    comp_label=label
                #print('comp label',comp_label)
                if c.state != 'off' and c.name!='Sum':
                    plot_obj.add_model_plot(c.SED, line_style=line_style, label=comp_label,flim=self.flux_plot_lim,auto_label=auto_label,color=color, density=density,update=False, frame=frame)

            c=self.get_spectral_component_by_name('Sum')
            if label is not None:
                comp_label = label
            else:
                comp_label='Sum'

            plot_obj.add_model_plot(c.SED, line_style='--', label=comp_label, flim=self.flux_plot_lim,color=color, density=density,update=False)

        plot_obj.update_plot()

        return plot_obj

    @safe_run
    def set_blob(self):
        if hasattr(self, 'T_esc_e_second'):
            if self.T_esc_e_second is None:
                self._blob.T_esc_e_second = self.parameters.R.val / BlazarSED.vluce_cm
            else:
                self._blob.T_esc_e_second = self.T_esc_e_second
        if self.emitters_distribution._user_defined is True:
            self.emitters_distribution._fill()
        BlazarSED.Init(self._blob, self.get_DL_cm())
        if self.emitters_distribution._user_defined is True:
            self.emitters_distribution._set_blob()


    @safe_run
    def set_external_fields(self):
        self.set_blob()
        BlazarSED.spectra_External_Fields(1,self._blob)

    def lin_func(self, lin_nu, init, phys_output=False, update_emitters=True):
        if self.emitters_distribution is None:
            raise RuntimeError('emitters distribution not defined')

        if init is True:
            self.set_blob()
            self._update_spectral_components()
        BlazarSED.Run_SED(self._blob)

        if phys_output==True:
            BlazarSED.EnergeticOutput(self._blob,0)

        if self.emitters_distribution._user_defined is False:
            if update_emitters is True:
                self.emitters_distribution.update()
            else:
                self.emitters_distribution._fill()

        nu_sed_sum, nuFnu_sed_sum = self.spectral_components.Sum.get_SED_points(lin_nu=lin_nu,log_log=False,interp=self._jetkernel_interp)
        return nu_sed_sum, nuFnu_sed_sum

    def _eval_model(self, lin_nu, log_nu, init, loglog, phys_output=False, update_emitters=True):
        log_model = None
        lin_nu, lin_model = self.lin_func(lin_nu, init, phys_output, update_emitters)
        if loglog is True:
            log_model = np.log10(lin_model)

        return lin_model, log_model


    def _prepare_nu_model(self, nu, loglog):
        if nu is None:
            lin_nu = np.logspace(np.log10(self.nu_min), np.log10(self.nu_max), self.nu_size)
            log_nu = np.log10(lin_nu)
        else:
            if np.shape(nu) == ():
                nu = np.array([nu])

            if loglog is True:
                lin_nu = np.power(10., nu)
                log_nu = nu
            else:
                log_nu = np.log10(nu)
                lin_nu = nu

        return lin_nu, log_nu

    @safe_run
    def eval(self,
             init=True,
             fill_SED=True,
             nu=None,
             get_model=False,
             loglog=False,
             plot=None,
             label=None,
             phys_output=False,
             update_emitters=True):
        out_model = None
        lin_nu, log_nu = self._prepare_nu_model(nu, loglog)

        lin_model, log_model= self._eval_model(lin_nu, log_nu ,init, loglog, phys_output=phys_output,
                                                update_emitters=update_emitters)
        #print('-->',lin_nu.min(),lin_nu.max())
        if fill_SED is True:
            self._fill(lin_nu,lin_model)

        if get_model is True:

            if loglog is True:
                out_model = log_model
            else:
                out_model = lin_model

        return out_model


    def energetic_report(self,verbose=True):
        self._build_energetic_report()            
        if verbose is True:
            _show_table(self.energetic_report_table)

    def _build_energetic_report(self,):
        self.energetic_dict={}
        BlazarSED.SetBeaming(self._blob)
        _energetic = BlazarSED.EnergeticOutput(self._blob,0)
        _par_array=ModelParameterArray()

        _name = [i for i in _energetic.__class__.__dict__.keys() if i[:1] != '_']
        _par_array.add_par(ModelParameter(name='BulkLorentzFactor', val=self._blob.BulkFactor, units='',par_type=''))
        self.energetic_dict['BulkLorentzFactor']= self._blob.BulkFactor
        try:
            for _n in _name:
                units = 'skip_this'
                if _n[0]=='L':
                    par_type='Lum. blob rest. frme.'
                    units='erg/s'
                elif _n[0] == 'U' and 'DRF' not in _n:
                    par_type = 'Energy dens. blob rest. frame'
                    units = 'erg/cm^3'
                elif _n[0] == 'U' and 'DRF'  in _n:
                    par_type = 'Energy dens. disk rest. frame'
                    units = 'erg/cm^3'
                elif _n[0] == 'j':
                    par_type = 'jet Lum.'
                    units = 'erg/s'
                else:
                    if _n !='thisown':
                        warnings.warn('energetic name %s not understood'%_n)

                if units == 'skip_this':
                    pass
                else:
                    self.energetic_dict[_n]=getattr(_energetic, _n)

                    _par_array.add_par(ModelParameter(name=_n, val=getattr(_energetic, _n), units=units,par_type=par_type))


            if self.emitters_distribution.emitters_type=='electrons':
                _i=_par_array.par_table['name']=='U_p_target'
                _par_array.par_table.remove_rows(_i)
                _=self.energetic_dict.pop('U_p_target')

            if self.emitters_distribution.emitters_type=='protons':
                _i=_par_array.par_table['name']=='U_p_cold'
                _par_array.par_table.remove_rows(_i)
                _=self.energetic_dict.pop('U_p_cold')

            self.energetic_report_table = _par_array.par_table
            self.energetic_report_table.remove_columns(['log','frozen','phys. bound. min','phys. bound. max'])
            self.energetic_report_table.rename_column('par type','type')

        except Exception as e:
            print('_energetic',_energetic)
            raise RuntimeError('energetic_report failed',e)







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
            x,y=self.spectral_components.Sum.get_SED_points(log_log=log_log)
            msk1=x>freq_range[0]
            msk2=x<freq_range[1]
            y_m= y[msk1*msk2].max()
            x_id= np.argmax(y[msk1*msk2])
            return x[msk1*msk2][x_id],y_m

    def get_component_peak(self,comp_name=None,log_log=False):
        comp = self.get_spectral_component_by_name(comp_name)

        ID = np.argmax(comp.SED.nuFnu.value)
        x_p, y_p = comp.SED.nu[ID].value, comp.SED.nuFnu[ID].value

        if log_log is True:
            x_p, y_p=np.log10([x_p,y_p])

        return x_p, y_p




class Jet(JetBase):
    """ Jet class

    """

    def __init__(self,
                 cosmo=None,
                 name=None,
                 emitters_type='electrons',
                 emitters_distribution='plc',
                 emitters_distribution_log_values=False,
                 beaming_expr='delta',
                 T_esc_e_second=None,
                 jet_workplace=None,
                 verbose=None,
                 clean_work_dir=True,
                 electron_distribution=None,
                 proton_distribution=None,
                 electron_distribution_log_values=None,
                 proton_distribution_log_values=None):
        """

        Parameters
        ----------
        cosmo
        name
        emitters_type
        emitters_distribution
        emitters_distribution_log_values
        beaming_expr
        T_esc_e_second
        jet_workplace
        verbose
        clean_work_dir
        electron_distribution
        proton_distribution
        electron_distribution_log_values
        proton_distribution_log_values
        """
        if electron_distribution is not None:
            emitters_type = 'electrons'
            emitters_distribution= electron_distribution
        if electron_distribution_log_values is not None:
            emitters_distribution_log_values = electron_distribution_log_values

        if proton_distribution is not None:
            emitters_type = 'protons'
            emitters_distribution= proton_distribution
        if proton_distribution_log_values is not None:
            emitters_distribution_log_values = proton_distribution_log_values


        if name is None:
            _name = 'jet'
        else:
            _name = name

        super(Jet,self).__init__(cosmo=cosmo,
                                 name=_name,
                                 emitters_type=emitters_type,
                                 emitters_distribution=emitters_distribution,
                                 emitters_distribution_log_values=emitters_distribution_log_values,
                                 beaming_expr=beaming_expr,
                                 jet_workplace=jet_workplace,
                                 verbose=verbose,
                                 clean_work_dir=clean_work_dir)

        if name is None or name == '':
            if self.emitters_distribution.emitters_type == 'electrons':
                name = 'jet_leptonic'
            elif self.emitters_distribution.emitters_type == 'protons':
                name = 'jet_hadronic_pp'
            else:
                name = 'jet'

        self.name = clean_var_name(name)

        if self.emitters_distribution.emitters_type == 'protons':
            self.T_esc_e_second=T_esc_e_second
            self.add_pp_gamma_component()
            self.add_pp_neutrino_component()
            self.add_bremss_ep_component()

        if self.emitters_distribution.emitters_type == 'electrons':
            self.electron_distribution=self.emitters_distribution

        if self.emitters_distribution.emitters_type == 'protons':
            self.protons_distribution=self.emitters_distribution

    @staticmethod
    def available_electron_distributions():
        JetBase.available_emitters_distributions()

    @staticmethod
    def available_proton_distributions():
        JetBase.available_emitters_distributions()

    def get_proton_distribution_name(self):
        return self.get_emitters_distribution_name()

    def get_electron_distribution_name(self):
        return self.get_emitters_distribution_name()

    def show_proton_distribution(self):
        self.show_emitters_distribution()

    def show_electron_distribution(self):
        self.show_emitters_distribution()

    def add_bremss_ep_component(self):
        self._add_spectral_component('Bremss_ep', var_name='do_bremss_ep', state_dict=dict((('on', 1), ('off', 0))))

    def add_pp_gamma_component(self):
        self._add_spectral_component('PP_gamma', var_name='do_pp_gamma', state_dict=dict((('on', 1), ('off', 0))))

    def add_pp_neutrino_component(self):
        self._add_spectral_component('PP_neutrino_tot', var_name='do_pp_neutrino',
                                     state_dict=dict((('on', 1), ('off', 0))))
        self._add_spectral_component('PP_neutrino_mu', var_name='do_pp_neutrino',
                                     state_dict=dict((('on', 1), ('off', 0))))
        self._add_spectral_component('PP_neutrino_e', var_name='do_pp_neutrino',
                                     state_dict=dict((('on', 1), ('off', 0))))





    def set_N_from_U_emitters(self,U, gmin=None, gmax=None):
        """ Sets the normalization of N to match the energy density of the primary emitters
        Parameters
        ----------
        U: float, (erg/cm3)
        gmin: float, optional,
            minimum value to evaluate the integral

        gmax: float, optional,
            maximum value to evaluate the integral
        Returns
        -------

        """
        N = self.parameters.N.val
        #gamma_grid_size = self._blob.gamma_grid_size
        #self.emitters_distribution.set_grid_size(1000)
        #self.emitters_distribution._fill()
        #self.set_blob()
        #BlazarSED.EvalU_e(self._blob)
        ratio = U/self.emitters_distribution.eval_U(gmin=gmin, gmax=gmax)
        #self.emitters_distribution.set_grid_size(gamma_grid_size)
        self.emitters_distribution._fill()
        self.set_par('N', val=N *ratio)
        self.set_blob()

    def set_N_from_U_vol_emitters(self, U_vol, gmin=None, gmax=None):
        """Sets the normalization of N to match the volume integrated energy of the primary emitters

        Parameters
        ----------
        U_vol: float (erg)

        gmin: float, optional,
            minimum value to evaluate the integral

        gmax: float, optional,
            maximum value to evaluate the integral

        Returns
        -------

        """
        #gamma_grid_size = self._blob.gamma_grid_size
        #self.emitters_distribution.set_grid_size(1000)
        #self.set_blob()
        U=U_vol/ self._blob.Vol_sphere
        #self.emitters_distribution.set_grid_size(gamma_grid_size)
        self.set_N_from_U_emitters(U, gmin=gmin, gmax=gmax)


    def set_N_from_L_sync(self,L_sync):
        """Sets the normalization of N to match the src integrated Luminosity of the   synchrotron emission

        Parameters
        ----------
        L_sync : float (erg/s)

        Returns
        -------

        """
        self.set_par('N', val=1.0)
        #gamma_grid_size = self._blob.gamma_grid_size
        #self.emitters_distribution.set_grid_size(100)
        #self.set_blob()
        delta = self._blob.beam_obj
        ratio = L_sync/(BlazarSED.Power_Sync_Electron(self._blob)* delta ** 4)
        #self.emitters_distribution.set_grid_size(gamma_grid_size)
        self.set_par('N', val=ratio)


    def set_N_from_F_sync(self, F_sync):
        """Sets the normalization of N to match the observed integrated synchrotron flux

        Parameters
        ----------
        F_sync : float, (erg cm-1 s-1 Hz-1)
            observed integrated synchrotron flux

        Returns
        -------

        """
        DL = self.get_DL_cm()
        L = F_sync * DL * DL * 4.0 * np.pi
        self.set_N_from_L_sync(L)

    def set_N_from_nuLnu(self,nuLnu_src, nu_src):
        """Sets the normalization of N to match the src Luminosity of the   synchrotron emission at src frequency nu

        Parameters
        ----------
        nuLnu_src : float, (erg/s)
            Luminosity of the   synchrotron emission at src frequency nu

        nu_src: float (Hz)
            synchrotron emission at src frequency nu (Hz)

        Returns
        -------

        """
        self.set_par('N',val=1.0)
        #gamma_grid_size = self._blob.gamma_grid_size
        #self.emitters_distribution.set_grid_size(100)
        self.set_blob()
        delta = self._blob.beam_obj
        nu_blob = nu_src / delta
        L_out = BlazarSED.Lum_Sync_at_nu(self._blob, nu_blob) * delta ** 4
        N_out = nuLnu_src / L_out
        #self.emitters_distribution.set_grid_size(gamma_grid_size)
        self.set_par('N', val=N_out)


    def set_N_from_nuFnu(self, nuFnu_obs, nu_obs):
        """Sets the normalization of N to match the observed flux nuFnu_obs at a given frequency nu_obs

        Parameters
        ----------
        nuFnu_obs: float, (erg cm-1 s-1 Hz-1)
            observed differential synchrotron flux

        nu_obs: float, (Hz)
            synchrotron emission at src frequency nu (Hz)

        Returns
        -------

        """

        self.set_blob()
        DL =  self.get_DL_cm()
        L = nuFnu_obs * DL * DL * 4.0 * np.pi
        nu_rest = nu_obs * (1 + self.parameters.z_cosm.val)
        self.set_N_from_nuLnu( L, nu_rest)



    def set_B_eq(self, nuFnu_obs, nu_obs, B_min=1E-9,B_max=1.0,N_pts=20,plot=False):
        """Sets the magnetic field (B) equipartition from numerical minimization over a logarithmic grid  of B values,
        for a given observed flux of the  synchrotron emission (nuFnu_obs) at a given observed frequency (nu_obs)

        Parameters
        ----------
        nuFnu_obs: float, (erg cm-1 s-1 Hz-1)
            observed differential synchrotron flux

        nu_obs: float, (Hz)
            synchrotron emission at src frequency nu (Hz)

        B_min: float, (Gauss), optional
            lower bound for B-grid

        B_max:  float, (Gauss), optional
            upper bound for B-grid

        N_pts: int, optional
            Other number of points to build the B-grid

        plot: book, optional, default=False
            if True plots the numerical grid

        Returns
        -------

        """

        b_grid = np.logspace(np.log10(B_min), np.log10(B_max), N_pts)
        print ('B grid min ',B_min)
        print ('B grid max ',B_max)
        print ('grid points',N_pts)
        U_e = np.zeros(N_pts)
        U_B = np.zeros(N_pts)
        N = np.zeros(N_pts)
        self.set_par('B', b_grid[0])
        self.set_blob()

        for ID, b in enumerate(b_grid):
            self.set_par('B', b)
            self.set_par('N', 1.0)
            # print 'B_eq',ID
            self.set_N_from_nuFnu(nuFnu_obs, nu_obs)
            N[ID]=self.get_par_by_name('N').val
            self.set_blob()
            #
            U_e[ID] = self._blob.U_e
            U_B[ID] = self._blob.UB
            # delta=Jet.get_beaming()
            # print "check L_in=%4.4e L_out=%4.4e"%(L_0,(L_0/delta**4)/BlazarSED.Power_Sync_Electron(Jet._Jet__blob))

        ID_min = np.argmin(U_B + U_e)

        if plot==True:
            #import  pylab as plt
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

