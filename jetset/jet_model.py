"""Jet model classes and utilities for spectral energy distribution calculations."""

__author__ = "Andrea Tramacere"


#import json
#import dill as pickle
import six
import numpy as np
import copy
import warnings
import os
from astropy import units as u
from astropy import constants
#from contextlib import redirect_stdout
#import io
import multiprocessing

from .jet_spectral_components import JetSpecComponent, SpecCompList

from .model_parameters import ModelParameterArray, ModelParameter, _show_table
from .base_model import Model
from .output import makedir,WorkPlace
from  .plot_sedfit import plt
from .cosmo_tools import Cosmo
from .utils import set_str_attr, set_particle_attr, old_model_warning, get_info, clean_var_name, get_nested_attr
from .jet_paramters import *
from .jet_emitters import *
from .jet_emitters_factory import EmittersFactory, InjEmittersFactory
from .jet_kernel_tools import set_emitters_c_array1d_fast as set_emitters
from .jet_kernel_tools import get_emitters_c_array1d_fast as get_emitters
from .jet_tools import *
from .mathkernel_helper import bessel_table_file_path
from .internal_absorption import InternalAbsorption

on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if on_rtd is True:
    from .mock import jetkernel as BlazarSED
else:
    from .jetkernel import jetkernel as BlazarSED



__all__=['Jet','JetBase','GalacticBeamed','GalacticUnbeamed']


class _JetEmittersDistributionProxy(object):
    """Proxy used by Jet equilibrium mode to hide carrier parameters."""

    def __init__(self, target):
        object.__setattr__(self, '_target', target)

    def _parameters_access_error(self):
        return AttributeError(
            'jet.emitters_distribution.parameters is not available when using '
            'leptonic equilibrium on Jet; use jet.parameters '

        )

    def unwrap(self):
        """Return the wrapped distribution object."""
        return self._target

    def __getattr__(self, name):
        if name == 'parameters':
            raise self._parameters_access_error()
        return getattr(self._target, name)

    def __setattr__(self, name, value):
        if name == 'parameters':
            raise self._parameters_access_error()
        setattr(self._target, name, value)

    def __dir__(self):
        names = set(dir(self._target))
        names.discard('parameters')
        return sorted(names)

    def __repr__(self):
        return repr(self._target)

    def __str__(self):
        return str(self._target)


class JetBase(Model):
    """Base jet model wrapping the jetkernel backend.

    The class exposes a Python interface to configure physical parameters,
    emitter distributions, spectral components, and evaluation options.
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
                 verbose=False,
                 nu_size=500,
                 clean_work_dir=True,
                 geometry='spherical',
                 **keywords):
        """Build a jet model and initialize the backend state.

        Parameters
        ----------
        cosmo : Cosmo, optional
            Cosmology helper. If omitted, a default :class:`Cosmo` instance is used.
        name : str, optional
            Model name.
        emitters_type : str, optional
            Primary emitter type, typically ``'electrons'`` or ``'protons'``.
        emitters_distribution : str or BaseEmittersDistribution, optional
            Emitter distribution identifier or distribution object.
        emitters_distribution_log_values : bool, optional
            If ``True``, distribution parameters are interpreted in log space.
        beaming_expr : str, optional
            Beaming parametrization (for example ``'delta'`` or ``'bulk_theta'``).
        jet_workplace : WorkPlace, optional
            Output workspace used to store model products.
        verbose : bool, optional
            Verbosity flag forwarded to the backend.
        nu_size : int, optional
            Number of frequency points used by model evaluation.
        clean_work_dir : bool, optional
            If ``True``, clean existing output directory before use.
        geometry : str, optional
            Emission-region geometry. Allowed values are
            ``'spherical'`` and ``'spherical_shell'``.
        """

        super(JetBase,self).__init__( cosmo=cosmo, **keywords)

        self.name = clean_var_name(name)

        self.model_type='jet'
        #self._emitters_type=emitters_type
        self._scale='lin-lin'
        self.verbose=verbose
        self._blob = self._build_blob(verbose=verbose)
        self._static_spec_arr_grid_size = BlazarSED.static_spec_arr_grid_size
        self._nu_static_size = BlazarSED.static_spec_arr_size
        self.nu_size = nu_size
        self.nu_grid_size=self._get_nu_grid_size_blob()
        N=multiprocessing.cpu_count()
        if N>20:
            print("*** limiting number of C threads to 20, use set_num_c_threads to override this configuration ***")
            N=20
            
        self.set_num_c_threads(N=N)
        if jet_workplace is None:
            jet_workplace=WorkPlace()
            out_dir= jet_workplace.out_dir + '/' + self.name + '_jet_prod/'
        else:
            out_dir=jet_workplace.out_dir+'/'+self.name+'_jet_prod/'

        self._allowed_geometry=['spherical','spherical_shell']

        self.geometry=geometry
        self.set_path(out_dir,clean_work_dir=clean_work_dir)

        self.set_flag(self.name)

        

        self._allowed_EC_components_list=['EC_BLR',
                                          'DT',
                                          'EC_DT',
                                          'Corona',
                                          'EC_Corona',
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
        self._spectral_components_list=[]
        self._hidden_spectral_components_list = []
        self._internal_absorption_comp={}
        self.spectral_components= SpecCompList(self._spectral_components_list)

        self.add_basic_components()
        self.SED=self.get_spectral_component_by_name('Sum').SED

        self.parameters = JetModelParameterArray(model=self)

        self._emitting_region_dict = None
        self._electron_distribution_dic= None
        self._external_photon_fields_dic= None
        self._original_inj_emitters_distr = None
        self._inj_emitters_distribution = None
        self._leptonic_equilibrium = False
        self._energetic = None
        self._setup(emitters_distribution,emitters_distribution_log_values,beaming_expr,emitters_type)


    def _setup(self, emitters_distribution, emitters_distribution_log_values, beaming_expr, emitters_type):
        """Initialize internal model state and bind the selected emitter distribution.

        Parameters
        ----------
        emitters_distribution : str or BaseEmittersDistribution
            Emitter distribution definition used to configure the model.
        emitters_distribution_log_values : bool
            If ``True``, interpret distribution parameters provided in logarithmic space.
        beaming_expr : str
            Beaming parameterization used by the emitting region (for example ``'delta'``).
        emitters_type : str
            Primary particle type handled by the model (for example ``'electrons'``).
        """
        self.EC_components_list = []
        self._spectral_components_list = []

        self.spectral_components = SpecCompList(self._spectral_components_list)
        self.add_basic_components()

        self.SED = self.get_spectral_component_by_name('Sum').SED

        self.parameters = JetModelParameterArray(model=self)

        self._emitting_region_dict = None
        self._emitters_distribution_dic = None
        self._external_photon_fields_dic = None



        self._blob.core.IC_adaptive_e_binning = 0
        self._blob.core.do_IC_down_scattering = 0
        if hasattr(emitters_distribution,'emitters_type'):
            emitters_type=emitters_distribution.emitters_type

        self.set_emitting_region(beaming_expr,emitters_type)

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
        j= self._serialize_model()
        return j 

    def __setstate__(self,state):
        self.__init__()
        self._decode_model(state)
        self._fix_par_dep_on_load(verbose=False)
        self._internal_absorption_comp = {}
        if '_internal_absorption_comp' in state:
            for _, p in self._extract_internal_abs_pars(state['_internal_absorption_comp']).items():
                self.enable_internal_absorption(**p)

    @staticmethod
    def _extract_internal_abs_pars(serialized_internal_abs):
        out = {}
        if not isinstance(serialized_internal_abs, dict):
            return out

        for comp, item in serialized_internal_abs.items():
            if isinstance(item, dict) and isinstance(item.get('pars'), dict):
                out[comp] = dict(item['pars'])
                continue
            if isinstance(item, dict) and ('comp' in item):
                out[comp] = dict(item)
        return out
        
    @staticmethod
    def _copy_emitters_parameter_state(source_distr, target_distr):
        for source_par in source_distr.parameters.par_array:
            target_par = target_distr.parameters.get_par_by_name(source_par.name)
            if target_par is None:
                target_distr.add_par(
                    source_par.name,
                    par_type=source_par.par_type,
                    val=source_par.val,
                    vmin=source_par.val_min,
                    vmax=source_par.val_max,
                    unit=source_par.units,
                    log=source_par.islog,
                    frozen=source_par.frozen,
                )
                target_par = target_distr.parameters.get_par_by_name(source_par.name)

            target_par.set(
                val=source_par.val,
                val_min=source_par.val_min,
                val_max=source_par.val_max,
                units=source_par.units,
                frozen=source_par.frozen,
                log=source_par.islog,
                skip_dep_par_warning=True,
            )

    def _build_serializable_emitters_distribution(self):
        src = self.emitters_distribution

        if isinstance(src, EmittersArrayDistribution):
            snapshot = EmittersArrayDistribution(
                name=src.name,
                emitters_type=src.emitters_type,
                normalize=src.normalize,
                gamma_array=np.asarray(src._array_gamma, dtype=np.float64).copy(),
                n_gamma_array=np.asarray(src._array_n_gamma, dtype=np.float64).copy(),
                gamma_grid_size=int(src._gamma_grid_size),
            )
            self._copy_emitters_parameter_state(src, snapshot)
            snapshot._update_parameters_dict()
            return snapshot

        if isinstance(src, EmittersDistribution):
            if hasattr(self,'_emitters_from_factory'):
                is_factory_distribution =self._emitters_from_factory
            else:
                available = set(EmittersFactory.available_distributions_list())
                is_factory_distribution = src.name in available
            if is_factory_distribution:
                snapshot = EmittersFactory().create_emitters(
                    src.name,
                    gamma_grid_size=int(src._gamma_grid_size),
                    log_values=src._log_values,
                    emitters_type=src.emitters_type,
                    normalize=src.normalize,
                )
            else:
                snapshot = EmittersDistribution(
                    name=src.name,
                    spectral_type=src.spectral_type,
                    gamma_grid_size=int(src._gamma_grid_size),
                    log_values=src._log_values,
                    emitters_type=src.emitters_type,
                    normalize=src.normalize,
                )

            self._copy_emitters_parameter_state(src, snapshot)

            if is_factory_distribution is False or src.spectral_type == 'user_defined':
                distr_func = getattr(src, '_py_distr_func', None)
                if distr_func is None:
                    distr_func = getattr(src, 'distr_func', None)
                    if hasattr(distr_func, 'py_func'):
                        distr_func = distr_func.py_func
                if distr_func is not None:
                    snapshot.set_distr_func(distr_func)

            snapshot._update_parameters_dict()
            return snapshot

        raise RuntimeError('emitters distribution type not valid', type(src))

    def _serialize_model(self):
        _model = {}
        _model['version']=get_info()['version']
        _model['name'] = self.name
        _model['emitters_type'] = self.emitters_distribution.emitters_type
       
        if self._inj_emitters_distribution is None:
            if isinstance(self.emitters_distribution,EmittersDistribution):
                _model['custom_emitters_distribution'] = self._build_serializable_emitters_distribution()
                clean_numba(_model['custom_emitters_distribution'])
                _model['emitters_distribution_class'] = 'EmittersDistribution'
            else:
                raise RuntimeError('emitters distribution type not valid', type(self.emitters_distribution))

        else:
            if isinstance(self._inj_emitters_distribution ,InjEmittersDistribution):
                _model['custom_emitters_distribution']=self._inj_emitters_distribution
                clean_numba(_model['custom_emitters_distribution'])
                _model['emitters_distribution_class'] = 'InjEmittersDistribution'
            else:
                raise  RuntimeError('inj_emitters_distribution distribution type not valid',type(self._inj_emitters_distribution))
        
        if hasattr(self,'geometry'):
            _model['geometry']=self.geometry

        _model['beaming_expr'] = self._beaming_expr
        _model['spectral_components_name'] = self.get_spectral_component_names_list()
        _model['spectral_components_state'] = [c.state for c in self._spectral_components_list]
        _model['EC_components_name'] = self.EC_components_list
        _model['basic_components_name'] = self.basic_components_list
        _model['cosmo'] = self.cosmo._serialize_model()
        _model['pars'] = {}
        _model['pars']=self.parameters._serialize_pars()
        _model['external_field_transf']=self.get_external_field_transf()
        _model['_internal_absorption_comp'] = {}
        for comp, item in self._internal_absorption_comp.items():
            if isinstance(item, dict) and isinstance(item.get('pars'), dict):
                _model['_internal_absorption_comp'][comp] = {'pars': dict(item['pars'])}
        _model['internal_pars'] = {}
        _model['internal_pars']['nu_size'] = self.nu_size
        _model['internal_pars']['nu_seed_size'] = self.nu_seed_size
        _model['internal_pars']['gamma_grid_size'] = self.gamma_grid_size
        _model['internal_pars']['IC_nu_size']=self.IC_nu_size
        _model['internal_pars']['nu_min']=self.nu_min
        _model['internal_pars']['nu_max']=self.nu_max
        _model['internal_pars']['Norm_distr'] = self.Norm_distr
        return _model

    #def save_model(self,file_name,to_string=False):
    #    pickle.dump(self, open(file_name, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
      

    @classmethod
    #@safe_run
    def load_model(cls, file_name_or_obj, from_string=False):
        """Load a serialized jet model.

        Parameters
        ----------
        file_name_or_obj : str or bytes or file-like
            Serialized model source accepted by the internal pickle loader.
        from_string : bool, optional
            If ``True``, interpret ``file_name_or_obj`` as in-memory content
            instead of a filesystem path.

        Returns
        -------
        JetBase
            Reconstructed model with backend state refreshed via ``set_blob``.
        """
        try:
            jet=cls._load_pickle(file_name_or_obj,from_string=from_string)
            jet.set_blob()
            jet._update_spectral_components()
            return jet
        except Exception as e:
            raise RuntimeError('The model you loaded is not valid please check the file name', e)



    def _decode_model(self,_model):
        
        if 'version' in _model.keys():
            self._set_version(_model['version'])
        else:
            self._set_version(v='unknown(<1.1.2)')
        self.cosmo = Cosmo.from_model(_model['cosmo'])
        self.model_type = 'jet'
        self.name = _model['name']
        if 'geometry' in _model.keys():
            self.geometry=_model['geometry']
        else:
             warnings.warn('the model you are loading has not geometry specification. Assuming it is spherical!')
             self.geometry='spherical'
             
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


        if _model['emitters_distribution_class'] == 'EmittersDistribution' or _model['emitters_distribution_class'] == 'InjEmittersDistribution':
            self.set_emitters_distribution(distr=_model['custom_emitters_distribution'], init=False)
        else:
            raise RuntimeError(
                'emitters distribution type not valid',
                type(_model.get('custom_emitters_distribution'))
            )

    
        for c in self.basic_components_list:
            if c not in _model['basic_components_name']:
                self.del_spectral_component(c)

        if self.emitters_distribution.emitters_type == 'protons':

            self.add_pp_gamma_component()
            self.add_pp_neutrino_component()
            self.add_bremss_ep_component()

       
        if 'disk_type' in _model['pars'].keys():
            disk_type=_model['pars']['disk_type']['val']
        else:
            disk_type=None

        self.add_EC_component(_model['EC_components_name'],disk_type=disk_type)
        self.set_external_field_transf(_model['external_field_transf'])
        
        for ID, c in enumerate(_model['spectral_components_name']):
            comp = getattr(self.spectral_components, c)
            if comp._state_dict != {}:
                comp.state = _model['spectral_components_state'][ID]

        self.SED = self.get_spectral_component_by_name('Sum').SED
       
        self.set_emitting_region(str(_model['beaming_expr']),self.emitters_distribution.emitters_type)
        if isinstance(self,GalacticBeamed):
            self._handle_z(d=self.cosmo._DL_cm)
        
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
            v['name']=k
            self.parameters.add_par(ModelParameter(**v))

        _par_dict = _model['internal_pars']
        for k in _par_dict.keys():
            setattr(self,k,_par_dict[str(k)])


    def _build_blob(self, verbose=False):

        blob = BlazarSED.MakeBlob()

        blob.Sync.x_Bessel_min = 1E-17
        blob.Sync.x_Bessel_max = 7.2E2

        blob.Sync.x_ave_Bessel_min = 1E-16
        blob.Sync.x_ave_Bessel_max = 3.5E2

        blob.Sync.log_x_Bessel_min = np.log10( blob.Sync.x_Bessel_min)
        blob.Sync.log_x_Bessel_max = np.log10( blob.Sync.x_Bessel_max)


        blob.Sync.log_x_ave_Bessel_min = np.log10( blob.Sync.x_ave_Bessel_min)
        blob.Sync.log_x_ave_Bessel_max = np.log10( blob.Sync.x_ave_Bessel_max)




        F_Sync_x_ptr = get_nested_attr(blob, 'Sync.F_Sync_x')

        F_Sync_y_ptr = get_nested_attr(blob, 'Sync.F_Sync_y')

        G_Sync_x_ptr = get_nested_attr(blob, 'Sync.G_Sync_x')

        G_Sync_y_ptr = get_nested_attr(blob, 'Sync.G_Sync_y')


        F_ave_Sync_x_ptr = get_nested_attr(blob, 'Sync.F_ave_Sync_x')

        F_ave_Sync_y_ptr = get_nested_attr(blob, 'Sync.F_ave_Sync_y')

        log_F_Sync_x_ptr = get_nested_attr(blob, 'Sync.log_F_Sync_x')

        log_F_Sync_y_ptr = get_nested_attr(blob, 'Sync.log_F_Sync_y')

        log_F_ave_Sync_x_ptr = get_nested_attr(blob, 'Sync.log_F_ave_Sync_x')

        log_F_ave_Sync_y_ptr = get_nested_attr(blob, 'Sync.log_F_ave_Sync_y')

        log_G_Sync_x_ptr = get_nested_attr(blob, 'Sync.log_G_Sync_x')

        log_G_Sync_y_ptr = get_nested_attr(blob, 'Sync.log_G_Sync_y')
        d = np.genfromtxt(bessel_table_file_path,comments='#')
        log_F_Sync_x=np.log10(d[:,0])
        log_F_Sync_y=np.log10(d[:,1])
        log_F_ave_Sync_x=np.log10(d[:,2])
        log_F_ave_Sync_y=np.log10(d[:,3])
        log_G_Sync_x=np.log10(d[:,4])
        log_G_Sync_y=np.log10(d[:,5])

        for ID, l in enumerate(d):
            BlazarSED.set_bessel_table(F_Sync_x_ptr,blob, l[0], ID)
            BlazarSED.set_bessel_table(F_Sync_y_ptr, blob,  l[1], ID)
            BlazarSED.set_bessel_table(F_ave_Sync_x_ptr, blob,  l[2], ID)
            BlazarSED.set_bessel_table(F_ave_Sync_y_ptr, blob,  l[3], ID)
            BlazarSED.set_bessel_table(G_Sync_x_ptr,blob, l[4], ID)
            BlazarSED.set_bessel_table(G_Sync_y_ptr, blob,  l[5], ID)


            BlazarSED.set_bessel_table(log_F_Sync_x_ptr, blob, log_F_Sync_x[ID], ID)
            BlazarSED.set_bessel_table(log_F_Sync_y_ptr, blob, log_F_Sync_y[ID], ID)
            BlazarSED.set_bessel_table(log_F_ave_Sync_x_ptr, blob, log_F_ave_Sync_x[ID], ID)
            BlazarSED.set_bessel_table(log_F_ave_Sync_y_ptr, blob, log_F_ave_Sync_y[ID], ID)
            BlazarSED.set_bessel_table(log_G_Sync_x_ptr, blob, log_G_Sync_x[ID], ID)
            BlazarSED.set_bessel_table(log_G_Sync_y_ptr, blob, log_G_Sync_y[ID], ID)

        blob.core.BESSEL_TABLE_DONE=1

        if verbose is False:
            blob.core.verbose = 0
        else:
            blob.core.verbose = 1

        set_str_attr(blob, 'core.path', './')

        set_str_attr(blob, 'core.MODE', 'custom')

        blob.emitters.gamma_grid_size = 200

        blob.core.nu_IC_size = 100
        blob.core.nu_seed_size = 100

        blob.core.nu_grid_size= 1000

        blob.core.do_Sync = 2

        blob.core.do_SSC = 1

        blob.core.R = 5.0e15
        blob.core.h_sh=0.1

        blob.core.B = 0.1

        blob.core.z_cosm = 0.1

        blob.core.BulkFactor = 10
        blob.core.theta = 0.1

        blob.emitters.N = 100

        blob.PP_gamma.NH_pp = 1

        blob.emitters.NH_cold_to_rel_e = 1.0       

        blob.Disk.L_Disk = 1E45

        blob.DT.L_DT = 1E45
        blob.Corona.L_Corona = 1E44
        blob.Corona.R_Corona = 1E13
        blob.Corona.R_H_Corona = 0.0
        blob.Corona.alpha_Corona = 1.0
        blob.Corona.nu_cut_low_Corona = 0.0
        blob.Corona.nu_cut_Corona = 1E20

        blob.emitters.gmin = 2

        blob.emitters.gmax = 1e6

        blob.Sync.spec.nu_min = 1e6
        blob.Sync.spec.nu_max = 1e20

        blob.SSC.spec.nu_min = 1e14
        blob.SSC.spec.nu_max = 1e30

        blob.core.nu_start_grid = 1e6
        blob.core.nu_stop_grid = 1e30
        
        return blob

    def set_emitting_region(self,beaming_expr,emitters_type):
        """Configure emitting-region parameters on the backend.

        Parameters
        ----------
        beaming_expr : str
            Beaming parametrization used to build region parameters.
        emitters_type : str
            Emitter family used to select the corresponding parameter set.
        """
        if  self._emitting_region_dict is not None:
            self.del_par_from_dic(self._emitting_region_dict)

        self._beaming_expr=beaming_expr

        set_str_attr(self._blob,'core.BEAMING_EXPR',beaming_expr)

        self._emitting_region_dict=build_emitting_region_dict(self.cosmo,beaming_expr=beaming_expr,emitters_type=emitters_type)
        self._set_geometry()
        self.parameters.add_par_from_dict(self._emitting_region_dict,self,'_blob',JetParameter)
        

    def _set_geometry(self):
        self._blob.core.GEOMETRY=self.geometry
        if self.geometry == 'spherical':
            pass
        elif self.geometry == 'spherical_shell':
            self._emitting_region_dict.pop('R')
            self._emitting_region_dict['R_sh'] = JetModelDictionaryPar(ptype='region_size', vmin=1E3, vmax=1E30,
                                                                       punit='cm', froz=False, log=False,
                                                                       jetkernel_par_name='core.R_sh')
            self._emitting_region_dict['h_sh'] = JetModelDictionaryPar(ptype='scaling_factor', val=0.1, vmin=0,
                                                                       punit='', vmax=1, froz=False, log=False,
                                                                       jetkernel_par_name='core.h_sh')

    
    @property
    def spectral_components_list(self):
        """Return visible spectral components.

        Returns
        -------
        list
            Spectral components with ``hidden == False``.
        """
        return [s for s in self._spectral_components_list if s.hidden is False]

    @property
    def geometry(self,):
        """Return the emitting-region geometry."""
        return self._geometry

    @geometry.setter
    def geometry(self,geometry):
        """Set the emitting-region geometry."""
        if geometry in self._allowed_geometry:
            self._geometry=geometry
        else:
            raise RuntimeError("geometry %s not allowed" % geometry, "please choose among ",
                                   self._allowed_geometry)


    def make_conical_jet(self, R=None,theta_open=5.,theta_open_min=1,theta_open_max=10):
        """Convert the model to a conical geometry parametrization.

        The method introduces a user parameter ``theta_open`` and sets
        ``R = tan(theta_open) * R_H`` as a dependent parameter relation.

        Parameters
        ----------
        R : float, optional
            Target region size used to initialize ``R_H`` after linking.
            If omitted, the current ``R`` value is used.
        theta_open : float, optional
            Semi-opening angle in degrees.
        theta_open_min : int, optional
            Lower bound for ``theta_open``.
        theta_open_max : int, optional
            Upper bound for ``theta_open``.

        Raises
        ------
        RuntimeError
            If dependency setup fails (for example because parameters already
            exist with incompatible dependency state).
        """
        if R is None:
            R=self.parameters.R.val
        try:
            self.add_user_par(name='theta_open',val=theta_open,units='deg',val_min=theta_open_min,val_max=theta_open_max)
            self.make_dependent_par(par='R', depends_on=['R_H','theta_open'],
                              par_expr='np.tan(np.radians(theta_open))*R_H')
            self.parameters.R_H.free()
           
            R_H=R/(np.tan(np.radians(self.parameters.theta_open.val)))
            print('setting R_H to',R_H)
            self.parameters.R_H.val=R_H
        except RuntimeError as e:
            msg='\n'
            msg+='something wrong, probably theta_open has already been added, or R is already dependant\n'
            msg+='please, try to reset the dependency of R: e.g my_jet.parameters.R.reset_dependencies()\n'
            msg+='or, try to remove theta_open: e.g my_jet.parameters.del_par(my_jet.parameters.theta_open)\n'
            msg+='the error is in:\n'
            msg+=str(e)

            raise RuntimeError(msg)
        
        except Exception as e:
            raise RuntimeError('unexpected exception',str(e))
    

    def set_EC_dependencies(self):
        """Attach standard EC size-luminosity relations for BLR and DT.

        When disk parameters are available, this sets dependent-parameter
        expressions for ``R_BLR_in``, ``R_BLR_out``, and ``R_DT`` as
        functions of ``L_Disk``.
        """
        try:
            if  self.parameters.get_par_by_name('L_Disk') is not None and self.parameters.get_par_by_name('R_BLR_in') is not None:
                try:            
                    self.make_dependent_par(par='R_BLR_in', depends_on=['L_Disk'], par_expr='3E17*(L_Disk/1E46)**0.5')
                    self.make_dependent_par(par='R_BLR_out', depends_on=['R_BLR_in'], par_expr='R_BLR_in*1.1')
                except RuntimeError as e:
                    print('functional dependencies already invoked for BLR',str(e))
            if  self.parameters.get_par_by_name('L_Disk') is not None and self.parameters.get_par_by_name('R_DT') is not None:
                try:           
                    self.make_dependent_par(par='R_DT', depends_on=['L_Disk'], par_expr='2E19*(L_Disk/1E46)**0.5')
                except RuntimeError as e:
                    print('functional dependencies already invoked for DT',str(e)) 

        except Exception as e:
            print('unexpected exception',str(e))

    @property
    def IC_adaptive_e_binning(self,):
        """Return whether adaptive electron binning is enabled for IC."""
        return np.intc(self._blob.core.IC_adaptive_e_binning)

    @IC_adaptive_e_binning.setter
    def IC_adaptive_e_binning(self,state):
        """Enable or disable adaptive electron binning for IC."""
        if type(state) == bool:
            pass
        else:
            raise RuntimeError('state has to be boolean')
        self._blob.core.IC_adaptive_e_binning=np.intc(state)

    @staticmethod
    def available_emitters_distributions():
        """Print available emitter distributions from the factory."""
        EmittersFactory.available_distributions()

    @staticmethod
    def _coerce_inj_emitters_distribution(q_inj):
        if isinstance(q_inj, InjEmittersDistribution) is False:
            raise RuntimeError('q_inj has to be an instance of InjEmittersDistribution')
        if q_inj.emitters_type != 'electrons':
            raise RuntimeError('q_inj emitters_type has to be electrons')

        if hasattr(q_inj, 'distr_func') is True and q_inj.distr_func is not None:
            return q_inj

        try:
            q_inj_from_factory = InjEmittersFactory().create_inj_emitters(name=q_inj.spectral_type,
                                                                           gamma_grid_size=q_inj._gamma_grid_size,
                                                                           log_values=q_inj._log_values,
                                                                           emitters_type=q_inj.emitters_type,
                                                                           normalize=q_inj.normalize)
            for p in q_inj.parameters.par_array:
                _p = q_inj_from_factory.parameters.get_par_by_name(p.name)
                if _p is not None:
                    _p.set(val=p.val, skip_dep_par_warning=True)
            return q_inj_from_factory
        except Exception as e:
            raise RuntimeError('q_inj has no distribution function; use InjEmittersFactory().create_inj_emitters(...) or set_distr_func(...)') from e

    def _disable_leptonic_equilibrium(self, remove_parameters=False):
        self._leptonic_equilibrium = False
        self._inj_emitters_distribution = None
        self._original_inj_emitters_distr = None
        self._blob.emitters.do_equilibrium = 0
        if remove_parameters is True:
            for par_name in ('T_esc_e_primaries','L_inj'):
                p = self.parameters.get_par_by_name(par_name)
                if p is not None:
                    self.parameters.del_par(p)

    def _ensure_leptonic_equilibrium_parameters(self):
        if self.parameters.get_par_by_name('T_esc_e_primaries') is None:
            _d={}
            _d['T_esc_e_primaries']=JetModelDictionaryPar(ptype='escape_time', 
                                                          vmin=1, 
                                                          vmax=None, 
                                                          punit='R/c',
                                                          froz=False, log=False,
                                                          val=1.0,
                                                          jetkernel_par_name='emitters.T_esc_e_primaries')

            self.parameters.add_par_from_dict(_d,self,'_blob',JetParameter)
       

    def _sync_inj_emitters_distribution_from_jet_parameters(self):
        if self._inj_emitters_distribution is None:
            return
        for inj_par in self._inj_emitters_distribution.parameters.par_array:
            jet_par = self.parameters.get_par_by_name(inj_par.name)
            if jet_par is not None:
                inj_par.set(val=jet_par.val, skip_dep_par_warning=True)

    def _sync_jet_parameters_from_inj_emitters_distribution(self):
        if self._inj_emitters_distribution is None:
            return
        for inj_par in self._inj_emitters_distribution.parameters.par_array:
            jet_par = self.parameters.get_par_by_name(inj_par.name)
            if jet_par is not None:
                jet_par.set(val=inj_par.val, skip_dep_par_warning=True)

    def _build_equilibrium_carrier_distribution(self):
        size = int(self._inj_emitters_distribution._gamma_grid_size)
        gmax = self._inj_emitters_distribution.parameters.get_par_by_name('gmax').val_lin
        gmax = max(gmax, 1.0)
        gamma_array = np.logspace(0.0, np.log10(gmax), size)
        n_gamma_array = np.zeros(gamma_array.size, dtype=np.float64)
        return EmittersArrayDistribution(name=f'{self._inj_emitters_distribution.name}_eq_carrier',
                                         emitters_type='electrons',
                                         gamma_array=gamma_array,
                                         n_gamma_array=n_gamma_array,
                                         normalize=False)

    def _set_equilibrium_injection_on_blob(self):
        if self._inj_emitters_distribution is None:
            raise RuntimeError('leptonic equilibrium mode requires an InjEmittersDistribution')

        self._sync_inj_emitters_distribution_from_jet_parameters()

        p_gmin = self.parameters.get_par_by_name('gmin')
        p_gmax = self.parameters.get_par_by_name('gmax')
        if p_gmin is None or p_gmax is None:
            raise RuntimeError('gmin/gmax parameters are required to initialize leptonic equilibrium')

        self._blob.emitters.gmin = 1
        #NOTE: safe boundary for eq. evolution
        self._blob.emitters.gmax = p_gmax.val_lin*2
        self._blob.emitters.gamma_grid_size = int(self._inj_emitters_distribution._gamma_grid_size)

        set_str_attr(self._blob, 'core.DISTR', 'jetset')
        set_particle_attr(self._blob, 'electrons')

        BlazarSED.setNgrid(self._blob)
        BlazarSED.build_Ne_jetset(self._blob)
        size = int(self._blob.emitters.gamma_grid_size)
        gamma_ptr = get_nested_attr(self._blob, 'emitters.griglia_gamma_jetset_Ne_log')
        ne_ptr = get_nested_attr(self._blob, 'emitters.Ne_jetset')
        gamma_blob, _ = get_emitters(gamma_ptr, ne_ptr, self._blob, size)

        self._inj_emitters_distribution._gamma_grid_size = size
        self._inj_emitters_distribution._gamma_grid = np.asarray(gamma_blob, dtype=np.float64)
        f = self._inj_emitters_distribution._eval_func(gamma=self._inj_emitters_distribution._gamma_grid)
        gmin_inj = self._inj_emitters_distribution.parameters.get_par_by_name('gmin').val_lin
        gmax_inj = self._inj_emitters_distribution.parameters.get_par_by_name('gmax').val_lin
        f = np.asarray(f, dtype=np.float64)
        f[self._inj_emitters_distribution._gamma_grid < gmin_inj] = 0.0
        f[self._inj_emitters_distribution._gamma_grid > gmax_inj] = 0.0
        norm = 1.0
        if self._inj_emitters_distribution.normalize is True:
            integ = np.trapezoid(f, self._inj_emitters_distribution._gamma_grid)
            if integ > 0:
                norm = 1.0 / integ
        BlazarSED.InitRadiative(self._blob,0)
        volume=self._blob.core.Vol_region
        self._inj_emitters_distribution._set_L_inj(self.parameters.get_par_by_name('L_inj').val,volume )
        q_inj = f * norm * self._inj_emitters_distribution.parameters.get_par_by_name('Q').val
        q_inj = np.asarray(q_inj, dtype=np.float64)
        q_inj[np.isnan(q_inj)] = 0
        q_inj[np.isinf(q_inj)] = 0

        self._inj_emitters_distribution.f = q_inj
        self._inj_emitters_distribution.gamma_e = self._inj_emitters_distribution._gamma_grid.copy()
        self._inj_emitters_distribution.n_gamma_e = q_inj.copy()

        # Safety/clarity: clear transport buffer before writing q_inj.
        set_emitters(ne_ptr, self._blob, size, np.zeros(size, dtype=np.float64))
        set_emitters(ne_ptr, self._blob, size, q_inj)

    def set_emitters_distribution(self, distr=None, log_values=False, emitters_type='electrons', init=True):
        """Set or replace the emitter distribution used by the model.

        This method supports analytic distributions (by name), array-based
        distributions, and injection distributions for leptonic-equilibrium
        mode. It updates Jet parameters and backend bindings accordingly.

        Parameters
        ----------
        distr : str or EmittersDistribution or ArrayDistribution or InjEmittersDistribution
            Distribution specification.
        log_values : bool, optional
            If ``True``, interpret analytic distribution parameters in log space.
        emitters_type : str, optional
            Emitter type used when creating distributions from names/arrays.
        init : bool, optional
            If ``True``, initialize backend state before replacing the
            distribution.
        """

        self._emitters_from_factory=False
        if init is True:
            self.set_blob()
        self._emitters_distribution_log_values = log_values

        if self._emitters_distribution_dic is not None:
            self.del_par_from_dic(self._emitters_distribution_dic)

        if isinstance(distr, InjEmittersDistribution):
            q_inj = self._coerce_inj_emitters_distribution(distr)
            if hasattr(q_inj, '_activate_numba'):
                q_inj._activate_numba()

            self._leptonic_equilibrium = True
            self._blob.emitters.do_equilibrium = 1
            self._inj_emitters_distribution = copy.deepcopy(q_inj)
            self._original_inj_emitters_distr = copy.deepcopy(q_inj)

            self._emitters_distribution_dic = copy.deepcopy(self._inj_emitters_distribution._parameters_dict)
            if 'gmin' in self._emitters_distribution_dic:
                # gmin drives only q_inj in equilibrium mode, not the solved Ne grid.
                self._emitters_distribution_dic['gmin'].is_in_jetkernel = False
                self._emitters_distribution_dic['gmin'].jetkernel_par_name = None
            if 'gmax' in self._emitters_distribution_dic:
                self._emitters_distribution_dic['gmax'].is_in_jetkernel = True
                self._emitters_distribution_dic['gmax'].jetkernel_par_name = 'emitters.gmax'

            self.parameters.add_par_from_dict(self._emitters_distribution_dic, self, '_blob', JetParameter)
            self._sync_jet_parameters_from_inj_emitters_distribution()
            self._ensure_leptonic_equilibrium_parameters()

            self.emitters_distribution = self._build_equilibrium_carrier_distribution()
            self.emitters_distribution.set_jet(self)
            self._sync_jet_parameters_from_inj_emitters_distribution()
            self.emitters_distribution._update_parameters_dict()
            self._emitters_distribution_name = self.emitters_distribution.name
            if self.parameters.get_par_by_name('L_inj') is None:
                self.parameters.add_par( ModelParameter(name='L_inj', par_type='L_inj', val=1E-3, val_min=0, val_max=None, units='erg/s'))
            if self.parameters.get_par_by_name('Q') is not None:
                self.parameters.del_par(self.parameters.get_par_by_name('Q'))
            
        elif isinstance(distr, ArrayDistribution):
            self._disable_leptonic_equilibrium(remove_parameters=True)
            self._emitters_distribution_name = 'from_array'
            self.emitters_distribution = EmittersDistribution.from_array(self, distr, emitters_type=emitters_type)
            self._emitters_distribution_dic = self.emitters_distribution._parameters_dict
            self.parameters.add_par_from_dict(self._emitters_distribution_dic, self, '_blob', JetParameter)
            self._attach_emitters_pars_to_jet(preserve_value_emitters=True)
            self.emitters_distribution.update()

        elif isinstance(distr, EmittersDistribution):
            self._disable_leptonic_equilibrium(remove_parameters=True)
            if hasattr(distr,'_activate_numba'):
                distr._activate_numba()
            self.emitters_distribution = copy.deepcopy(distr)
            self._update_emitters_pars_dependence()
            self.emitters_distribution.set_jet(self)
            self.emitters_distribution._update_parameters_dict()
            self._emitters_distribution_name = self.emitters_distribution.name
            self._emitters_distribution_dic = self.emitters_distribution._parameters_dict
            self.parameters.add_par_from_dict(self._emitters_distribution_dic, self, '_blob', JetParameter)
            self._attach_emitters_pars_to_jet(preserve_value_emitters=True)
            
            self.emitters_distribution.update()
            
        elif isinstance(distr, str):
            self._disable_leptonic_equilibrium(remove_parameters=True)
            nf=EmittersFactory()
            self._emitters_from_factory=True
            self.emitters_distribution = nf.create_emitters(distr, log_values=log_values, emitters_type=emitters_type)
            if hasattr( self.emitters_distribution,'_activate_numba'):
                self.emitters_distribution._activate_numba()
            self.emitters_distribution.set_jet(self)
            self.emitters_distribution._update_parameters_dict()
            self._emitters_distribution_name = self.emitters_distribution.name
            self._emitters_distribution_dic = self.emitters_distribution._parameters_dict
            self.parameters.add_par_from_dict(self._emitters_distribution_dic, self, '_blob', JetParameter)
            self._attach_emitters_pars_to_jet(preserve_value_emitters=True)
            self._update_emitters_pars_dependence()
            self.emitters_distribution.update()
        else:
            raise RuntimeError('distr',type(distr),'not valid should be a string or an',type(EmittersDistribution),'instance')

    def _attach_emitters_pars_to_jet(self, preserve_value_emitters=False, skip=None):

        if skip is None:
            skip = []

        v_dict={}
        for par in self.emitters_distribution.parameters.par_array[::]:
            if par.name not in skip:
                self.emitters_distribution.parameters.del_par(par)
                v_dict[par.name]=par.val
        for k in self._emitters_distribution_dic.keys():
            par = self.parameters.get_par_by_name(k)
            if par is None:
                continue
            if preserve_value_emitters is True and k in v_dict:
                par.set(val=v_dict[k],skip_dep_par_warning=True)
            self.emitters_distribution.parameters.add_par(par)

  


    def _update_emitters_pars_dependence(self):
        for par in self.emitters_distribution.parameters.par_array:
            if par.immutable is True:
                warnings.warn('parameter dependence has to be reassigned after emitters distribution is assigned to a jet, now will be reset no dep')
                self.emitters_distribution.parameters.reset_dependencies()
                break

    
    def get_emitters_distribution_name(self):
        """Return the active emitter-distribution name."""
        return self.emitters_distribution.name


    def show_spectral_components(self):
        """Print currently registered spectral components."""
        print ("Spectral components for Jet model:%s"%(self.name))

        for comp in self._spectral_components_list:

            print ("comp: %s "%(comp.name))

        print()


   
    def sed_table(self, restframe='obs'):
        """Provide the astropy table with SED spectral components

        Parameters
        ----------
        restframe : str, optional
            , by default 'obs'

        Returns
        -------
        astropy table
        """
        try:
            self.spectral_components.build_table(restframe=restframe)
        except:
            self.eval()
            
        self._SED_table = self.spectral_components.table
        return  self._SED_table

    def get_spectral_component_by_name(self,name,verbose=True):
        """Return a spectral component by name.

        Parameters
        ----------
        name : str
            Component name.
        verbose : bool, optional
            If ``True``, print available names when not found.

        Returns
        -------
        JetSpecComponent or None
            Matching component, or ``None`` if missing.
        """
        for i in range(len(self._spectral_components_list)):
            if self._spectral_components_list[i].name==name:
                return self._spectral_components_list[i]
        else:
            if verbose==True:
                print ("no spectral components with name %s found"%name)
                print ("names in array are:")
                self.show_spectral_components()

            return None

    def list_spectral_components(self):
        """Print the list of spectral-component names."""
        for i in range(len(self._spectral_components_list)):
            print (self._spectral_components_list[i].name)

    def get_spectral_component_names_list(self):
        """Return spectral-component names as a list."""
        _l=[]
        for i in range(len(self._spectral_components_list)):
            _l.append(self._spectral_components_list[i].name)
        return _l

    def del_spectral_component(self,name):
        """Delete a spectral component from the model."""
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

            self._spectral_components_list.remove(comp)

            self.spectral_components.__delattr__(name)



    def _add_spectral_component(self, name, var_name=None, state_dict=None,state=None):
        #print("==> adding",  name, var_name)
        self._spectral_components_list.append(
            JetSpecComponent(self, name, self._blob, var_name=var_name, state_dict=state_dict, state=state))
        setattr(self.spectral_components,name,self._spectral_components_list[-1])

    def _update_spectral_components(self,tau=None):
        for ID,s in enumerate(self._spectral_components_list):
            self._spectral_components_list[ID]= JetSpecComponent(self, s.name, self._blob, var_name=s._var_name, state_dict=s._state_dict, state=s.state, tau=tau)
            self._spectral_components_list[ID].hidden = s.hidden
            setattr(self.spectral_components, self._spectral_components_list[ID].name, self._spectral_components_list[ID])
            

    
    def add_basic_components(self):
        """Add default ``Sum``, ``Sync``, and ``SSC`` components."""
        self.basic_components_list=['Sum','Sync','SSC']

        self._add_spectral_component('Sum')
        self._add_spectral_component('Sync', var_name='core.do_Sync', state_dict=dict((('on', 1), ('off', 0), ('self-abs', 2))),state='self-abs')

        self._add_spectral_component('SSC', var_name='core.do_SSC', state_dict=dict((('on', 1), ('off', 0))))


    def add_sync_component(self,state='self-abs'):
        """Add a synchrotron component with the requested state."""
        self._add_spectral_component('Sync', var_name='core,do_Sync',
                                     state_dict=dict((('on', 1), ('off', 0), ('self-abs', 2))), state=state)

    def add_SSC_component(self,state='on'):
        """Add a synchrotron self-Compton component."""
        self._add_spectral_component('SSC', var_name='core.do_SSC', state_dict=dict((('on', 1), ('off', 0))),state=state)

    def del_EC_component(self,EC_components_list, disk_type='BB'):
        """Remove one or more external-Compton related components.
        
        Parameters
        ----------
        EC_components_list : str or list of str, optional

        disk_type : str, optional
            Disk template type for disk-dependent EC components.
        """
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
                    self._blob.core.do_EC_Disk = 0
                    self.EC_components_list.remove('EC_Disk')

                if self.get_spectral_component_by_name('EC_BLR', verbose=False) is not None:
                    self._blob.core.do_EC_BLR=0
                    self._del_spectral_component('EC_BLR', verbose=False)
                    self.EC_components_list.remove('EC_BLR')


            if EC_component=='EC_Disk':
                if self.get_spectral_component_by_name('EC_Disk', verbose=False) is not None:
                    self._blob.core.do_EC_Disk=0
                    self._del_spectral_component('EC_Disk', verbose=False)
                    self.EC_components_list.remove('EC_Disk')


            if EC_component=='EC_BLR':
                if self.get_spectral_component_by_name('EC_BLR', verbose=False) is not None:
                    self._blob.core.do_EC_BLR=0
                    self._del_spectral_component('EC_BLR', verbose=False)
                    self.EC_components_list.remove('EC_BLR')

            if EC_component=='DT':
                if self.get_spectral_component_by_name('DT', verbose=False) is not None:
                    self._del_spectral_component('DT', verbose=False)
                    self.EC_components_list.remove('DT')
                if self.get_spectral_component_by_name('EC_DT', verbose=False) is not None:
                    self._blob.core.do_EC_DT = 0
                    self._del_spectral_component('EC_DT', verbose=False)
                    self.EC_components_list.remove('EC_DT')

            if EC_component=='EC_DT':
                if self.get_spectral_component_by_name('EC_DT', verbose=False) is not None:
                    self._blob.core.do_EC_DT=0
                    self._del_spectral_component('EC_DT', verbose=False)
                    self.EC_components_list.remove('EC_DT')

            if EC_component=='Corona':
                if self.get_spectral_component_by_name('Corona', verbose=False) is not None:
                    self._blob.core.do_Corona=0
                    self._del_spectral_component('Corona', verbose=False)
                    self.EC_components_list.remove('Corona')
                if self.get_spectral_component_by_name('EC_Corona', verbose=False) is not None:
                    self._blob.core.do_EC_Corona=0
                    self._del_spectral_component('EC_Corona', verbose=False)
                    self.EC_components_list.remove('EC_Corona')

            if EC_component=='EC_Corona':
                if self.get_spectral_component_by_name('EC_Corona', verbose=False) is not None:
                    self._blob.core.do_EC_Corona=0
                    self._del_spectral_component('EC_Corona', verbose=False)
                    self.EC_components_list.remove('EC_Corona')

            if EC_component=='EC_CMB':
                if self.get_spectral_component_by_name('EC_CMB', verbose=False) is not None:
                    self._blob.core.do_EC_CMB=0
                    self._del_spectral_component('EC_CMB', verbose=False)
                    self.EC_components_list.remove('EC_CMB')

            if EC_component=='Star':
                if self.get_spectral_component_by_name('Star', verbose=False) is not None:
                    self._blob.core.do_Star=0
                    self._del_spectral_component('Star', verbose=False)
                    self.EC_components_list.remove('Star')

            if EC_component=='EC_Star':
                if self.get_spectral_component_by_name('EC_Star', verbose=False) is not None:
                    self._blob.core.do_EC_Star=0
                    self._del_spectral_component('EC_Star', verbose=False)
                    self.EC_components_list.remove('EC_Star')

        self.del_par_from_dic(build_ExtFields_dic(EC_components_list,disk_type))



    def add_EC_component(self,EC_components_list=[],disk_type=None):
        """Add one or more external-Compton related components.

        Parameters
        ----------
        EC_components_list : str or list of str, optional
            Component names to add. ``'All'`` enables all supported EC entries.
        disk_type : str, optional
            Disk template type for disk-dependent EC components.
        """

        if disk_type is not None:
            if disk_type not in self._allwed_disk_type:
                raise RuntimeError('disk type',disk_type,'not in allwowed', self._allwed_disk_type)
        else:
            if self.get_par_by_name('disk_type') is not None:
                disk_type=self.get_par_by_name('disk_type').val
            else:
                disk_type='BB'

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
                    self._add_spectral_component('Disk',var_name='core.do_Disk', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('Disk')

            if EC_component=='EC_Disk':
                #self._blob.core.do_EC_Disk=1
                if self.get_spectral_component_by_name('EC_Disk',verbose=False) is None:
                    self._add_spectral_component('EC_Disk', var_name='core.do_EC_Disk', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('EC_Disk')

                if self.get_spectral_component_by_name('Disk',verbose=False) is None:
                    self._add_spectral_component('Disk',var_name='core.do_Disk', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('Disk')

            if EC_component=='EC_BLR':
                #self._blob.core.do_EC_BLR=1
                if self.get_spectral_component_by_name('EC_BLR',verbose=False) is None:
                    self._add_spectral_component('EC_BLR', var_name='core.do_EC_BLR', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('EC_BLR')

                if self.get_spectral_component_by_name('Disk',verbose=False) is None:
                    self._add_spectral_component('Disk',var_name='core.do_Disk', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('Disk')


            if EC_component == 'DT':
                if self.get_spectral_component_by_name('DT',verbose=False) is None:
                    self._add_spectral_component('DT',var_name='core.do_DT', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('DT')

                if self.get_spectral_component_by_name('Disk',verbose=False) is None:
                    self._add_spectral_component('Disk',var_name='core.do_Disk', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('Disk')



            if EC_component == 'Star':
                if self.get_spectral_component_by_name('Star',verbose=False) is None:
                    self._add_spectral_component('Star',var_name='core.do_Star', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('Star')


            if EC_component == 'EC_Star':
                if self.get_spectral_component_by_name('Star',verbose=False) is None:
                    self._add_spectral_component('Star',var_name='core.do_Star', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('Star')

                if self.get_spectral_component_by_name('EC_Star',verbose=False) is None:
                    self._add_spectral_component('EC_Star', var_name='core.do_EC_Star', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('EC_Star')


            if EC_component=='EC_DT':
                #self._blob.core.do_EC_DT=1
                if self.get_spectral_component_by_name('EC_DT',verbose=False) is None:
                    self._add_spectral_component('EC_DT', var_name='core.do_EC_DT', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('EC_DT')

                if self.get_spectral_component_by_name('DT',verbose=False) is None:
                    self._add_spectral_component('DT',var_name='core.do_DT', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('DT')

                if self.get_spectral_component_by_name('Disk',verbose=False) is None:
                    self._add_spectral_component('Disk',var_name='core.do_Disk', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('Disk')

            if EC_component == 'Corona':
                if self.get_spectral_component_by_name('Corona',verbose=False) is None:
                    self._add_spectral_component('Corona',var_name='core.do_Corona', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('Corona')

            if EC_component == 'EC_Corona':
                if self.get_spectral_component_by_name('EC_Corona',verbose=False) is None:
                    self._add_spectral_component('EC_Corona', var_name='core.do_EC_Corona', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('EC_Corona')

                if self.get_spectral_component_by_name('Corona',verbose=False) is None:
                    self._add_spectral_component('Corona',var_name='core.do_Corona', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('Corona')

            if EC_component=='EC_CMB':
                #self._blob.core.do_EC_CMB=1
                if self.get_spectral_component_by_name('EC_CMB',verbose=False) is None:
                    self._add_spectral_component('EC_CMB', var_name='core.do_EC_CMB', state_dict=dict((('on', 1), ('off', 0))))
                    self.EC_components_list.append('EC_CMB')

        #IF disk_type is already a parameter it has to be updated here to make it effective
        self.parameters.add_par_from_dict(build_ExtFields_dic(self.EC_components_list, disk_type),self,'_blob',JetParameter)
        #print('--> disk_type',disk_type)
        
        if disk_type is not None and 'Disk' in self.EC_components_list:
            #remove the old disk type parameters
            if  self.parameters.get_par_by_name('disk_type') is not None:
                self.del_par_from_dic(build_ExtFields_dic(['Disk'],self.parameters.get_par_by_name('disk_type').val))
            else:
                self.EC_components_list.append('Disk')
            self.parameters.add_par_from_dict(build_ExtFields_dic(self.EC_components_list, disk_type),self,'_blob',JetParameter)
            self.parameters.disk_type.val = disk_type


    def enable_internal_absorption(self,
                                comp,
                                nu_min=None,
                                N_soft=50,
                                N_hard=20,
                                N_R_H=20,
                                N_theta=30,
                                use_R_H_profile_extrapolation=False,
                                use_sigma_gamma_gamma_fast=False):
        
        """Enable internal gamma-gamma absorption for a seed component.

        Parameters
        ----------
        comp : str
            Seed-photon component name used to compute optical depth.
        nu_min : float, optional
            Minimum frequency used for the internal absorption solver (Hz).
        N_soft, N_hard, N_R_H, N_theta : int, optional
            Numerical grid sizes used by the absorption solver.
        use_R_H_profile_extrapolation : bool, optional
            If ``True``, extrapolate the radial profile when required.
        use_sigma_gamma_gamma_fast : bool, optional
            If ``True``, enable tabulated/interpolated sigma_gamma_gamma in C backend.
        """
        self._internal_absorption_comp[comp]={}
        self._internal_absorption_comp[comp]['pars']=dict(N_hard=N_hard,
                                                          N_soft=N_soft,
                                                          N_theta=N_theta,
                                                          N_R_H=N_R_H,
                                                          nu_min=nu_min,
                                                          comp=comp,
                                                          use_R_H_profile_extrapolation=use_R_H_profile_extrapolation,
                                                          use_sigma_gamma_gamma_fast=use_sigma_gamma_gamma_fast)
        
        self._internal_absorption_comp[comp]['obj']=InternalAbsorption(jet=self,
                                                                nu_min=nu_min,
                                                                seed_photons_name=comp,
                                                                N_soft=N_soft,
                                                                N_hard=N_hard,
                                                                N_R_H=N_R_H,
                                                                N_theta=N_theta,
                                                                use_R_H_profile_extrapolation=use_R_H_profile_extrapolation,
                                                                use_sigma_gamma_gamma_fast=use_sigma_gamma_gamma_fast)
        self._configure_internal_absorption_on_blob(comp, self._internal_absorption_comp[comp]['pars'])

    def remove_internal_absorption(self,comp):
        """Disable internal absorption for a component."""
        if comp in self._internal_absorption_comp.keys():
            self._disable_internal_absorption_on_blob(comp)
            del self._internal_absorption_comp[comp]
            if len(self._internal_absorption_comp.keys()) > 0:
                BlazarSED.recompute_internal_absorption_tau(self._blob)
    
    def show_internal_absorption_components(self):
        """Print enabled internal-absorption components and settings."""
        if len(self._internal_absorption_comp.keys())>0:
            for comp in  self._internal_absorption_comp.keys():
                print('internal absorption  for component:', comp)
                for p in self._internal_absorption_comp[comp]['pars'].items():
                    print(p)
                print()
        else:
            print('internal absorption not enabled in this jet model')
            

    def eval_internal_absorption(self, comp, peak=False, nu=None):
        """Evaluate internal absorption for a component.

        Parameters
        ----------
        comp : str
            Enabled internal-absorption component name.
        peak : bool, optional
            If ``True``, use peak-optimized seed-photon sampling.
        nu : array-like or None, optional
            Source-frame frequency grid in Hz. If omitted, a default model
            grid is used even when the model has not been evaluated yet.

        Notes
        -----
        Internal absorption is evaluated through the isolated C path only
        (worker-blob evaluation, then merge of IA output only).

        Returns
        -------
        tuple
            ``(tau, nu_src)`` if available, otherwise ``(None, None)``.
        """
        if comp in self._internal_absorption_comp.keys():
            return self._internal_absorption_comp[comp]['obj'].eval(get_tau=True,
                                                                     lin_nu=nu,
                                                                     peak=peak)
        return None,None
    
    def del_par_from_dic(self,model_dic):
        """Delete parameters listed in a model-parameter dictionary."""
        for key in model_dic.keys():

            par=self.parameters.get_par_by_name(key)

            if par is not None:
                self.parameters.del_par(par)



    def get_DL_cm(self,eval_model=False):
        """Return luminosity distance in centimeters."""
        if eval_model is True:
            self.set_blob()

        if self.cosmo._c is not None:
            return self.cosmo.get_DL_cm(self.parameters.z_cosm.val)
        else:
            return self.cosmo.get_DL_cm()


    #@safe_run
    def get_beaming(self,):
        """Return the current beaming factor from the backend."""
        BlazarSED.SetBeaming(self._blob)
        return self._blob.core.beam_obj


    def set_flag(self,flag):
        """Set the backend STEM flag used by jetkernel output helpers."""
        self._blob.core.STEM=flag

    def get_flag(self):
        """Return the backend STEM flag."""
        return self._blob.core.STEM

    def get_path(self):
        """Return the backend working path."""
        return self._blob.core.path

    def set_path(self,path,clean_work_dir=True):
        """Set and create the backend working path."""
        if path.endswith('/'):
            pass
        else:
            path+='/'

        set_str_attr(self._blob,'core.path',path)
        makedir(path,clean_work_dir=clean_work_dir)

    def get_IC_mode(self):
        """Return inverse-Compton mode label (``'on'`` or ``'off'``)."""
        return dict(map(reversed, self._IC_states.items()))[self._blob.core.do_IC]



    def set_external_field_transf(self,val):
        """Set external-field transformation frame.

        Parameters
        ----------
        val : str
            Allowed values are ``'blob'`` and ``'disk'``.
        """
        if val not in self._external_field_transf.keys():
            raise RuntimeError('val',val,'not in allowed values',self._external_field_transf.keys())
        self._blob.core.EC_stat=self._external_field_transf[val]
        self._blob.core.EC_stat_orig=self._external_field_transf[val]

    def get_external_field_transf(self):
        """Return external-field transformation frame label."""
        return dict(map(reversed, self._external_field_transf.items()))[self._blob.core.EC_stat]

    def set_emiss_lim(self,val):
        """Set lower emissivity bound used by backend calculations."""
        self._blob.core.emiss_lim=val

    def get_emiss_lim(self):
        """Return lower emissivity bound used by backend calculations."""
        return self._blob.core.emiss_lim


    @property
    def IC_nu_size(self):
        """Return inverse-Compton spectral-grid size."""
        return self._blob.core.nu_IC_size

    @IC_nu_size.setter
    def IC_nu_size(self, val):
        """Set inverse-Compton spectral-grid size."""
        self.set_IC_nu_size(val)

    def get_IC_nu_size(self):
        """Return inverse-Compton spectral-grid size."""
        return self._blob.core.nu_IC_size

    def set_IC_nu_size(self, val):
        """Set inverse-Compton spectral-grid size."""
        if val > self._nu_static_size:
            raise RuntimeError('value can not exceed',self._nu_static_size)
        self._blob.core.nu_IC_size = val

    @property
    def nu_seed_size(self):
        """Return seed-photon spectral-grid size."""
        return self._blob.core.nu_seed_size

    def get_seed_nu_size(self):
        """Return seed-photon spectral-grid size."""
        return self._blob.core.nu_seed_size

    @nu_seed_size.setter
    def nu_seed_size(self,val):
        """Set seed-photon spectral-grid size."""
        self.set_seed_nu_size(val)

    def set_seed_nu_size(self,val):
        """Set seed-photon spectral-grid size."""
        if val>self._nu_static_size:
            raise RuntimeError('value can not exceed',self._nu_static_size)
        self._blob.core.nu_seed_size=val



    def set_gamma_grid_size(self,val):
        """Set particle Lorentz-factor grid size for emitters."""
        self.emitters_distribution.set_grid_size(gamma_grid_size=val)

    @property
    def gamma_grid_size(self):
        """Return particle Lorentz-factor grid size."""
        return self._blob.emitters.gamma_grid_size

    @gamma_grid_size.setter
    def gamma_grid_size(self,val):
        """Set particle Lorentz-factor grid size."""
        self.emitters_distribution.set_grid_size(gamma_grid_size=val)

    @property
    def nu_min(self):
        """Return minimum frequency of the model grid (Hz)."""
        return self._get_nu_min_grid()

    @nu_min.setter
    def nu_min(self, val):
        """Set minimum frequency of the model grid (Hz)."""
        if hasattr(self, '_blob'):
            self._set_nu_min_grid(val)

    def _get_nu_min_grid(self):
        return  self._blob.core.nu_start_grid


    def _set_nu_min_grid(self, val):
        self._blob.core.nu_start_grid=val

    @property
    def Norm_distr(self):
        """Return normalization mode of the emitter distribution.

        Returns
        -------
        bool or int
            ``True``/``1`` means physical-density normalization,
            ``False``/``0`` means scaling-factor normalization.
        """
        if hasattr(self,'emitters_distribution'):
            if self.emitters_distribution._user_defined is False:
                return self._blob.emitters.Norm_distr
            else:
                return self.emitters_distribution.normalize
        else:
            return None


    @Norm_distr.setter
    def Norm_distr(self, val):
        """Set normalization mode of the emitter distribution."""
        if hasattr(self, 'emitters_distribution') and self.get_par_by_name('N') is not None:

            if val == 1 or val is True:
                if self.emitters_distribution._user_defined is False:
                    self._blob.emitters.Norm_distr = 1
                else:
                    self.emitters_distribution.normalize = val
                self.parameters.N.par_type='emitters_density'
                self.parameters.N.units ='1/cm3'
            elif val == 0 or val is False:
                if self.emitters_distribution._user_defined is False:
                    self._blob.emitters.Norm_distr = 0
                else:
                    self.emitters_distribution.normalize = val
                    self.parameters.N.par_type = 'scaling_factor'
            else:
                raise RuntimeError('value', val, 'not allowed, allowed 0/1 or False/True')



    def switch_Norm_distr_ON(self):
        """Enable physical-density normalization mode."""
        self.Norm_distr=True

    def switch_Norm_distr_OFF(self):
        """Enable scaling-factor normalization mode."""
        self.Norm_distr = False

    @property
    def nu_max(self):
        """Return maximum frequency of the model grid (Hz)."""
        return self._get_nu_max_grid()

    @nu_max.setter
    def nu_max(self, val):
        """Set maximum frequency of the model grid (Hz)."""
        if hasattr(self, '_blob'):
            self._set_nu_max_grid(val)

    def _set_nu_max_grid(self, val):
        self._blob.core.nu_stop_grid=val

    def _get_nu_max_grid(self):
        return  self._blob.core.nu_stop_grid

    @property
    def nu_size(self):
        """Return number of sampled frequencies in Python-side evaluation."""
        return self._nu_size

    @nu_size.setter
    def nu_size(self, size):
        """Set number of sampled frequencies in Python-side evaluation."""
        self._nu_size = size

    def set_nu_grid_size(self, val):
        """Set jetkernel frequency-grid size."""
        self._set_nu_grid_size_blob(val)

    @property
    def nu_grid_size(self):
        """Return jetkernel frequency-grid size."""
        return self._get_nu_grid_size_blob()

    @nu_grid_size.setter
    def nu_grid_size(self, val):
        """Set jetkernel frequency-grid size."""
        self._set_nu_grid_size_blob(val)


    def _set_nu_grid_size_blob(self, val):
        if val < 100:
            val = 100
        if val > self._static_spec_arr_grid_size:
            raise RuntimeError('value can not exceed', self._nu_static_size)
        self._blob.core.nu_grid_size=val

    def _get_nu_grid_size_blob(self):
        return  self._blob.core.nu_grid_size



    def set_verbosity(self,val):
        """Set backend verbosity level."""
        self._blob.core.verbose=val

    def get_verbosity(self):
        """Return backend verbosity level."""
        return  self._blob.core.verbose

    def debug_synch(self):
        """Print synchrotron integration debug information."""
        print ("nu stop synch", self._blob.Sync.spec.nu_max)
        print ("nu stop synch ssc", self._blob.Sync.nu_stop_Sync_ssc)
        print ("ID MAX SYNCH", self._blob.Sync.NU_INT_STOP_Sync_SSC)

    def debug_SSC(self):
        """Print SSC integration debug information."""
        print ("nu start SSC", self._blob.SSC.spec.nu_min)
        print ("nu stop SSC", self._blob.SSC.spec.nu_max)
        print ("ID MAX SSC", self._blob.SSC.NU_INT_STOP_COMPTON_SSC)


    def show_emitters_distribution(self):
        """Print emitter-distribution settings and parameters."""
        print('-'*80)
        print('%s distribution:'%self.emitters_distribution.emitters_type)
        print(" type: %s  " % (self._emitters_distribution_name))
        print(" gamma energy grid size: ", self.gamma_grid_size)
        print(" gmin grid : %e" % self._blob.emitters.gmin_griglia)
        print(" gmax grid : %e" % self._blob.emitters.gmax_griglia)
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
        print("model description: ")
        print('-'*80)
        print("type:",self.__class__.__name__)
        print("name: %s  " % (self.name))
        print("geometry: %s  " % (self.geometry))
        print('')
        print('%s distribution:'%self.emitters_distribution.emitters_type)
        print(" type: %s  " % (self._emitters_distribution_name))
        print (" gamma energy grid size: ",self.gamma_grid_size)
        print (" gmin grid : %e"%self._blob.emitters.gmin_griglia)
        print (" gmax grid : %e"%self._blob.emitters.gmax_griglia)
        print(" normalization: ", self.Norm_distr)
        print(" log-values: ", self._emitters_distribution_log_values)
        #print(" leptonic equilibrium: ", getattr(self, '_leptonic_equilibrium', False))
        _p = self.parameters.get_par_by_name('NH_cold_to_rel_e') 
        if _p is not None:
            print(' ratio of cold protons to relativistic electrons: %e'%_p.val)
        print('')
        if 'Disk' in self.EC_components_list:
            print('accretion disk:')
            BlazarSED.set_Disk(self._blob)
            print(' disk Type: %s'%self.parameters.get_par_by_name('disk_type').val)
            print(' L disk: %e (erg/s)'%self.parameters.get_par_by_name('L_Disk').val)
            print(' T disk: %e (K)'%self._blob.Disk.T_Disk)
            print(' nu peak disk: %e (Hz)'%BlazarSED.eval_nu_peak_Disk(self._blob.Disk.T_Disk))
            if self.parameters.get_par_by_name('disk_type').val == 'MultiBB':
                print(' Sw radius %e (cm)'%self._blob.Disk.R_Sw)
                print(' L Edd. %e (erg/s)'%self._blob.Disk.L_Edd)
                yr=86400*365
                print(' accr_rate: %e (M_sun/yr)'%(yr*self._blob.Disk.accr_rate/BlazarSED.m_sun))
                print(' accr_rate Edd.: %e (M_sun/yr)'%(yr*self._blob.Disk.accr_Edd/BlazarSED.m_sun))

        if 'EC_Star' in self.EC_components_list:
            print('Star:')
            print(' Star radius %e (cm)'%self._blob.Star.R_Star, ',%e (R_Sun)'%(self._blob.Star.R_Star/constants.R_sun.to('cm').value))
            print('')
        print('radiative fields:')
        print (" seed photons grid size: ", self.nu_seed_size)
        print (" IC emission grid size: ", self.get_IC_nu_size())
        print (' source emissivity lower bound :  %e' % self._blob.core.emiss_lim)
        print (' spectral components:')
        for _s in self._spectral_components_list:
            print("   name:%s,"%_s.name, 'state:', _s.state)
            print("   name:%s,"%_s.name, 'hidden:', _s.hidden)
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
        """Plot model spectral components on a :class:`~jetset.plot_sedfit.PlotSED`.

        Parameters
        ----------
        plot_obj : jetset.plot_sedfit.PlotSED, optional
            Existing plot object. If omitted, a new one is created.
        clean : bool, optional
            If ``True``, remove previously plotted model lines before adding new
            ones.
        label : str, optional
            Label to use for the plotted line(s). If not provided, component
            names are used when ``auto_label=True``.
        comp : str, optional
            Name of a single spectral component to plot. If omitted, all active
            components plus ``Sum`` are plotted.
        sed_data : jetset.data_loader.ObsData, optional
            Observational data to add to the same plot.
        color : str, optional
            Matplotlib color for model line(s).
        auto_label : bool, optional
            If ``True``, automatically derive labels from component names.
        line_style : str, optional
            Line style used for non-sum components.
        frame : {"obs", "src", "blob"}, optional
            Output frame used for plotting.
        density : bool, optional
            If ``True``, plot density representation.

        Returns
        -------
        jetset.plot_sedfit.PlotSED
            Plot object with added model traces.
        """
        plot_obj=self._set_up_plot(plot_obj,sed_data,frame,density)

        if clean is True:
            plot_obj.clean_model_lines()

        if comp is not None:
            c = self.get_spectral_component_by_name(comp)
            if c is None:
                raise RuntimeError('spectral component', comp, 'not found')

            if label is not None:
                comp_label = label
            elif label is None and auto_label is True:
                comp_label = c.name
            else:
                comp_label=None
            if c.state!='off':
                plot_obj.add_model_plot(c.SED, line_style=line_style, label=comp_label,flim=self.flux_plot_lim,color=color,auto_label=auto_label,update=False, frame=frame)

        else:
            for c in self.spectral_components_list:
                comp_label = c.name
                if auto_label is not True:
                    comp_label=label
                #print('comp label',comp_label)
                if c.state != 'off' and c.name!='Sum':
                    plot_obj.add_model_plot(c.SED, line_style=line_style, label=comp_label,flim=self.flux_plot_lim,auto_label=auto_label,color=color,update=False, frame=frame)

            c=self.get_spectral_component_by_name('Sum')
            if label is not None:
                comp_label = label
            else:
                comp_label='Sum'

            plot_obj.add_model_plot(c.SED, line_style='--', label=comp_label, flim=self.flux_plot_lim,color=color,update=False, frame=frame)

        plot_obj.update_plot()

        return plot_obj

    #@safe_run
    def set_blob(self):
        """Synchronize Python parameters/distributions to the backend blob."""
        if self._leptonic_equilibrium is True:
            self._blob.emitters.do_equilibrium = 1
            self._set_equilibrium_injection_on_blob()
        else:
            self._blob.emitters.do_equilibrium = 0
            if self.emitters_distribution._user_defined is True:
                self.emitters_distribution._fill()


        BlazarSED.Init(self._blob, self.get_DL_cm())
        if self.emitters_distribution._user_defined is True:
            self.emitters_distribution._set_blob()
        for comp, item in self._internal_absorption_comp.items():
            if isinstance(item, dict) and isinstance(item.get('pars'), dict):
                self._configure_internal_absorption_on_blob(comp, item['pars'])

    #@safe_run
    def set_external_fields(self):
        """Update external radiation fields on the backend."""
        self.set_blob()
        BlazarSED.spectra_External_Fields(1,self._blob,1)

    def _get_internal_abs_component_on_blob(self, comp):
        if comp == 'BLR':
            return self._blob.core.internal_abs.BLR
        if comp == 'DT':
            return self._blob.core.internal_abs.DT
        if comp == 'Corona':
            return self._blob.core.internal_abs.Corona
        raise RuntimeError('internal absorption component %s not valid' % comp)

    def _configure_internal_absorption_on_blob(self, comp, pars):
        c_comp = self._get_internal_abs_component_on_blob(comp)
        new_use_rh = int(bool(pars.get('use_R_H_profile_extrapolation', False)))
        new_use_sigma_fast = int(bool(pars.get('use_sigma_gamma_gamma_fast', False)))
        new_n_soft = int(pars.get('N_soft', 50))
        new_n_hard = int(pars.get('N_hard', 50))
        new_n_r_h = int(pars.get('N_R_H', 50))
        new_n_theta = int(pars.get('N_theta', 50))
        nu_min = pars.get('nu_min', None)
        new_nu_min = -1.0 if nu_min is None else float(nu_min)

        config_changed = (
            int(c_comp.is_enabled) != 1 or
            int(c_comp.use_R_H_profile_extrapolation) != new_use_rh or
            int(c_comp.use_sigma_gamma_gamma_fast) != new_use_sigma_fast or
            int(c_comp.peak_mode) != 0 or
            int(c_comp.N_soft) != new_n_soft or
            int(c_comp.N_hard) != new_n_hard or
            int(c_comp.N_R_H) != new_n_r_h or
            int(c_comp.N_theta) != new_n_theta or
            float(c_comp.nu_min) != float(new_nu_min)
        )

        c_comp.is_enabled = 1
        if config_changed:
            c_comp.is_valid = 0
        c_comp.use_R_H_profile_extrapolation = new_use_rh
        c_comp.use_sigma_gamma_gamma_fast = new_use_sigma_fast
        c_comp.peak_mode = 0
        c_comp.N_soft = new_n_soft
        c_comp.N_hard = new_n_hard
        c_comp.N_R_H = new_n_r_h
        c_comp.N_theta = new_n_theta
        c_comp.nu_min = new_nu_min
        

    def _disable_internal_absorption_on_blob(self, comp):
        c_comp = self._get_internal_abs_component_on_blob(comp)
        c_comp.is_enabled = 0
        c_comp.is_valid = 0
   

    def lin_func(self, lin_nu, init, phys_output=False, update_emitters=True):
        """Evaluate model spectrum in linear units on ``lin_nu``.

        Internal absorption is applied at C level during ``Run_SED``.
        """
        if self.emitters_distribution is None:
            raise RuntimeError('emitters distribution not defined')

        if init is True:
            self.set_blob()
            self._update_spectral_components()
        BlazarSED.Run_SED(self._blob)

        if phys_output==True:
            BlazarSED.EnergeticOutput(self._blob)

        if self.emitters_distribution._user_defined is False:
            if update_emitters is True:
                self.emitters_distribution.update()
            else:
                self.emitters_distribution._fill()

        nu_sed_sum, nuFnu_sed_sum = self.spectral_components.Sum.get_SED_points(lin_nu=lin_nu,log_log=False,interp=self._jetkernel_interp)
        
        nuFnu_sed_hidden=0.
        for _sc in self._spectral_components_list:
            if _sc.hidden is True:
                _, _nuFnu_sed= _sc.get_SED_points(lin_nu=lin_nu,log_log=False,interp=self._jetkernel_interp)
                nuFnu_sed_hidden += _nuFnu_sed
        
        nuFnu_sed_sum = nuFnu_sed_sum - nuFnu_sed_hidden
       
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

    #@safe_run
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
        """Evaluate the jet model.

        Parameters are compatible with :meth:`Model.eval`, with additional
        control over backend initialization and emitter updates.
        """
        out_model = None
        lin_nu, log_nu = self._prepare_nu_model(nu, loglog)
        
        lin_model, log_model= self._eval_model(lin_nu, log_nu ,init, loglog, phys_output=phys_output,
                                                update_emitters=update_emitters)

        if fill_SED is True:
            self._fill(lin_nu,lin_model)
            self.spectral_components.Sum.SED.nuFnu=lin_model
        
        if get_model is True:

            if loglog is True:
                out_model = log_model
            else:
                out_model = lin_model

        return out_model


    def energetic_report(self,verbose=True,):
        """Build and optionally print energetic quantities table."""
        self._build_energetic_report()            
        if verbose is True:
            _show_table(self.energetic_report_table)

    def _build_energetic_dict(self):
        self.energetic_dict={}
        BlazarSED.SetBeaming(self._blob)

        #NOTE: workaround for the users be sure they  evaluate correctly the jet energetics for 'disk' frame
        _change_frame=False
        if self.get_external_field_transf()=='disk':
            self.set_external_field_transf('blob')
            self.eval()
            _change_frame=True
        self._energetic = BlazarSED.EnergeticOutput(self._blob)
        
        if _change_frame is True:
            self.set_external_field_transf('disk')
            self.eval()

        _par_array=ModelParameterArray()

        _name = [i for i in self._energetic.__class__.__dict__.keys() if i[:1] != '_']
        _par_array.add_par(ModelParameter(name='BulkLorentzFactor', val=self._blob.core.BulkFactor, units='',par_type='jet-bulk-factor'))
        self.energetic_dict['BulkLorentzFactor']= self._blob.core.BulkFactor
        try:
            for _n in _name:
                units = 'skip_this'
                if _n[0]=='L':
                    par_type='Lum. blob rest. frame.'
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
                    self.energetic_dict[_n]=getattr(self._energetic, _n)

                    _par_array.add_par(ModelParameter(name=_n, val=getattr(self._energetic, _n), units=units,par_type=par_type))
        except Exception as e:
            print('_energetic',self._energetic)
            raise RuntimeError('energetic_report failed',e)

        return  _par_array

    def _set_energetic_report(self,_par_array):
        if self.emitters_distribution.emitters_type=='electrons':
                _par_array.add_par(ModelParameter(name='NH_cold_to_rel_e', val=self.parameters.NH_cold_to_rel_e.val, units='',par_type='cold_p_to_rel_e_ratio'))
                self.energetic_dict['NH_cold_to_rel_e'] = self.parameters.NH_cold_to_rel_e.val

        self.energetic_report_table = _par_array.par_table
        self.energetic_report_table.remove_columns(['log','frozen','phys. bound. min','phys. bound. max'])
        self.energetic_report_table.rename_column('par type','type')

        if self.emitters_distribution.emitters_type=='electrons':
            _i=self.energetic_report_table['name']=='U_p_target'
            _i=np.argwhere(_i==True).flatten()[0]
            self.energetic_report_table.remove_row(_i)
            _=self.energetic_dict.pop('U_p_target')
            
            _i=self.energetic_report_table['name']=='U_p'
            _i=np.argwhere(_i==True).flatten()[0]
            self.energetic_report_table.remove_row(_i)
            _=self.energetic_dict.pop('U_p')

            _i=self.energetic_report_table['name']=='L_pp_gamma_rf'	
            _i=np.argwhere(_i==True).flatten()[0]
            self.energetic_report_table.remove_row(_i)
            _=self.energetic_dict.pop('L_pp_gamma_rf')

            _i=self.energetic_report_table['name']=='jet_L_p'	
            _i=np.argwhere(_i==True).flatten()[0]
            self.energetic_report_table.remove_row(_i)
            _=self.energetic_dict.pop('jet_L_p')

        if self.emitters_distribution.emitters_type=='protons':
            _i=self.energetic_report_table['name']=='U_p_cold'
            _i=np.argwhere(_i==True).flatten()[0]
            self.energetic_report_table.remove_row(_i)
            _=self.energetic_dict.pop('U_p_cold')

            _i=self.energetic_report_table['name']=='jet_L_p_cold'
            _i=np.argwhere(_i==True).flatten()[0]
            self.energetic_report_table.remove_row(_i)
            _=self.energetic_dict.pop('jet_L_p_cold')

        if isinstance(self,GalacticUnbeamed):
            i=[]
            for n in self.energetic_report_table['name']:
                if 'jet' in n:
                    i.append(True)
                    _=self.energetic_dict.pop(n)
                else:
                    i.append(False)
            i=np.array(i)
            i=np.argwhere(i==True)
            self.energetic_report_table.remove_rows(i)

            _d_l=['U_Disk','U_BLR','U_DT','U_Corona','U_CMB','U_Synch_DRF','U_Star']
            i=[]
            for n in self.energetic_report_table['name']:
                if n in _d_l:
                    i.append(True)
                    _=self.energetic_dict.pop(n)
                else:
                    i.append(False)
            i=np.array(i)
            i=np.argwhere(i==True)
            self.energetic_report_table.remove_rows(i)

            for n in self.energetic_report_table['name']:
                if n.endswith('_DRF'):
                    _i=self.energetic_report_table['name']==n
                    _i=np.argwhere(_i==True).flatten()[0]
                    nn=self.energetic_report_table[_i]['name'].replace('_DRF','')
                    self.energetic_report_table[_i]['name']=nn
                    self.energetic_dict[nn]=self.energetic_dict.pop(n)
                if n.endswith('_rf'):
                    _i=self.energetic_report_table['name']==n
                    _i=np.argwhere(_i==True).flatten()[0]
                    nn=self.energetic_report_table[_i]['name'].replace('_rf','')
                    self.energetic_report_table[_i]['name']=nn
                    self.energetic_dict[nn]=self.energetic_dict.pop(n)

            for r in self.energetic_report_table:
                if 'blob' in r['type']:
                    nn=r['type'].replace('blob','')
                    r['type']=nn
                    
                if 'disk' in r['type']:
                    nn=r['type'].replace('disk','')
                    r['type']=nn

    def _build_energetic_report(self,):
        try:
            _par_array=self._build_energetic_dict()
            self._set_energetic_report(_par_array)
        except Exception as e:
            raise RuntimeError('energetic_report failed',e)







    def get_SED_peak(self,peak_name=None,freq_range=None,log_log=False):
        """Return SED peak from backend attributes or a frequency interval.

        Parameters
        ----------
        peak_name : str, optional
            Backend peak attribute name (for example ``'nu_p_Sync_obs'``).
        freq_range : tuple, optional
            ``(nu_min, nu_max)`` range used to search the peak numerically.
        log_log : bool, optional
            If ``True``, return values in log10 scale.
        """
        if peak_name is not None and freq_range is not None:
            print ("either you provide peak_name or freq_range")
            raise ValueError
        elif   peak_name is not None:
            try:
                if log_log==False:
                    return get_nested_attr(self._blob, peak_name)
                else:
                     return np.log10(get_nested_attr(self._blob, peak_name))
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
        """Return peak frequency and flux for a spectral component."""
        comp = self.get_spectral_component_by_name(comp_name)

        ID = np.argmax(comp.SED.nuFnu.value)
        x_p, y_p = comp.SED.nu[ID].value, comp.SED.nuFnu[ID].value

        if log_log is True:
            x_p, y_p=np.log10([x_p,y_p])

        return x_p, y_p

    def set_num_c_threads(self,N=None,verbose=False):
        """Configure the number of C backend threads used by jetkernel.

        Parameters
        ----------
        N : int, optional
            Number of threads to assign to ``self._blob.core.N_THREADS``.
            If ``None``, use ``multiprocessing.cpu_count()``.
        verbose : bool, optional
            If ``True``, print the selected thread count.

        Raises
        ------
        RuntimeError
            If ``N`` is not an integer value.
        """

        if N is None:
            N=multiprocessing.cpu_count()
        if verbose:
            print("===> setting C threads to",N)
        if isinstance(N,int):
            self._blob.core.N_THREADS=N
        else:
            raise RuntimeError('N must be integer')



class Jet(JetBase):
    """Concrete jet model for leptonic or hadronic emitter populations."""

    def __init__(self,
                 cosmo=None,
                 name=None,
                 emitters_type='electrons',
                 emitters_distribution='plc',
                 emitters_distribution_log_values=False,
                 beaming_expr='delta',
                 jet_workplace=None,
                 verbose=False,
                 clean_work_dir=True,
                 electron_distribution=None,
                 proton_distribution=None,
                 electron_distribution_log_values=None,
                 proton_distribution_log_values=None,
                 geometry='spherical'):
        """Initialize a :class:`Jet` model.

        Parameters
        ----------
        cosmo : Cosmo, optional
            Cosmology helper object.
        name : str, optional
            Model name.
        emitters_type : str, optional
            Emitter type if no explicit electron/proton distribution alias is used.
        emitters_distribution : str or distribution object, optional
            Default distribution configuration.
        emitters_distribution_log_values : bool, optional
            If ``True``, interpret analytic distribution parameters in log space.
        beaming_expr : str, optional
            Beaming parametrization.
        jet_workplace : WorkPlace, optional
            Output workspace.
        verbose : bool, optional
            Verbosity flag.
        clean_work_dir : bool, optional
            If ``True``, clean model output directory before use.
        electron_distribution : str or distribution object, optional
            Convenience alias to select an electron distribution.
        proton_distribution : str or distribution object, optional
            Convenience alias to select a proton distribution.
        electron_distribution_log_values : bool, optional
            Log-value mode for the electron distribution alias.
        proton_distribution_log_values : bool, optional
            Log-value mode for the proton distribution alias.
        geometry : str, optional
            Emission-region geometry.
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


        if electron_distribution is not None and emitters_type == 'protons':
            raise RuntimeError("if you provide an electron_distribution, you can't set emitters_type=",emitters_type, 'use emitters_distribution instead')
        
        if proton_distribution is not None and emitters_type == 'electrons':
            raise RuntimeError("if you provide a proton_distribution, you can't set emitters_type=",emitters_type, 'use emitters_distribution instead')

        


        if name is None:
            _name = 'jet'
        else:
            _name = name = clean_var_name(name)

        super(Jet,self).__init__(cosmo=cosmo,
                                 name=_name,
                                 emitters_type=emitters_type,
                                 emitters_distribution=emitters_distribution,
                                 emitters_distribution_log_values=emitters_distribution_log_values,
                                 beaming_expr=beaming_expr,
                                 jet_workplace=jet_workplace,
                                 verbose=verbose,
                                 clean_work_dir=clean_work_dir,
                                 geometry=geometry)

        if name is None or name == '':
            if self.emitters_distribution.emitters_type == 'electrons':
                name = 'jet_leptonic'
            elif self.emitters_distribution.emitters_type == 'protons':
                name = 'jet_hadronic_pp'
            else:
                name = 'jet'

        self.name = clean_var_name(name)

        if self.emitters_distribution.emitters_type == 'protons':
            self.add_pp_gamma_component()
            self.add_pp_neutrino_component()
            self.add_bremss_ep_component()

        if self.emitters_distribution.emitters_type == 'electrons':
            self.electron_distribution=self.emitters_distribution
            self.parameters.NH_cold_to_rel_e.hidden=False

        if self.emitters_distribution.emitters_type == 'protons':
            self.protons_distribution=self.emitters_distribution

    def set_emitters_distribution(self, distr=None, log_values=False, emitters_type='electrons', init=True):
        """Set emitters distribution and hide carrier parameters in Jet equilibrium mode."""
        super(Jet, self).set_emitters_distribution(distr=distr,
                                                   log_values=log_values,
                                                   emitters_type=emitters_type,
                                                   init=init)
        if isinstance(distr, InjEmittersDistribution):
            self.emitters_distribution = _JetEmittersDistributionProxy(self.emitters_distribution)
        elif isinstance(self.emitters_distribution, _JetEmittersDistributionProxy):
            self.emitters_distribution = self.emitters_distribution.unwrap()

    @staticmethod
    def available_electron_distributions():
        """Print available electron-distribution names."""
        JetBase.available_emitters_distributions()

    @staticmethod
    def available_proton_distributions():
        """Print available proton-distribution names."""
        JetBase.available_emitters_distributions()

    def get_proton_distribution_name(self):
        """Return the active proton-distribution name."""
        return self.get_emitters_distribution_name()

    def get_electron_distribution_name(self):
        """Return the active electron-distribution name."""
        return self.get_emitters_distribution_name()

    def show_proton_distribution(self):
        """Print proton-distribution configuration."""
        self.show_emitters_distribution()

    def show_electron_distribution(self):
        """Print electron-distribution configuration."""
        self.show_emitters_distribution()

    def add_bremss_ep_component(self):
        """Add electron-proton bremsstrahlung spectral component."""
        self._add_spectral_component('Bremss_ep', var_name='Bremss_ep.do_bremss_ep', state_dict=dict((('on', 1), ('off', 0))))

    def add_pp_gamma_component(self):
        """Add proton-proton gamma-ray spectral component."""
        self._add_spectral_component('PP_gamma', var_name='PP_gamma.do_pp_gamma', state_dict=dict((('on', 1), ('off', 0))))

    def add_pp_neutrino_component(self):
        """Add proton-proton neutrino spectral components."""
        self._add_spectral_component('PP_neutrino_tot', var_name='PP_neutrino.do_pp_neutrino',
                                     state_dict=dict((('on', 1), ('off', 0))))
        self._add_spectral_component('PP_neutrino_mu', var_name='PP_neutrino.do_pp_neutrino',
                                     state_dict=dict((('on', 1), ('off', 0))))
        self._add_spectral_component('PP_neutrino_e', var_name='PP_neutrino.do_pp_neutrino',
                                     state_dict=dict((('on', 1), ('off', 0))))


    def set_N_from_U_emitters(self,U, gmin=None, gmax=None):
        """Set the normalization  of N or L_inj (leptonic eq.), to match target emitter energy density.

        Parameters
        ----------
        U : float
            Target energy density in ``erg cm^-3``.
        gmin, gmax : float, optional
            Integration bounds in Lorentz factor used for ``eval_U``.
        """
        
        ratio = U/self.emitters_distribution.eval_U(gmin=gmin, gmax=gmax)
        
        if self._leptonic_equilibrium:
            L_inj=self.parameters.L_inj.val
            self.set_par('L_inj', val=L_inj *ratio)
            self.set_blob()
        else:
            self.emitters_distribution._fill()
            N = self.parameters.N.val
            self.set_par('N', val=N *ratio)
            self.set_blob()

    def set_N_from_U_vol_emitters(self, U_vol, gmin=None, gmax=None):
        """Set the normalization  of N or L_inj (leptonic eq.), to match target integrated emitter energy.

        Parameters
        ----------
        U_vol : float
            Target integrated energy in ``erg`` over emitting volume.
        gmin, gmax : float, optional
            Integration bounds in Lorentz factor used for ``eval_U``.
        """
        U=U_vol/ self._blob.core.Vol_region
        self.set_N_from_U_emitters(U, gmin=gmin, gmax=gmax)


    def set_N_from_L_sync(self,L_sync):
        """Set the normalization  of N or L_inj (leptonic eq.), to match target integrated synchrotron luminosity.

        Parameters
        ----------
        L_sync : float
            Target source-frame synchrotron luminosity in ``erg s^-1``.
        """
        if self._leptonic_equilibrium:
            self.set_blob()
            self.set_par('L_inj', val=1E40)
            self.set_blob()
            delta = self.get_beaming()
            ratio = L_sync/(BlazarSED.Power_Sync_Electron(self._blob)* delta ** 4)*1E40
            self.set_par('L_inj', val=ratio)
        
        else:
            self.set_par('N', val=1.0)
            
            self.set_blob()
            delta = self.get_beaming()
            ratio = L_sync/(BlazarSED.Power_Sync_Electron(self._blob)* delta ** 4)
            self.set_par('N', val=ratio)
       
    def set_N_from_F_sync(self, F_sync):
        """Set the normalization  of N or L_inj (leptonic eq.), to match target observed integrated synchrotron flux.

        Parameters
        ----------
        F_sync : float
            Target observed integrated synchrotron flux in ``erg cm^-2 s^-1``.
        """
        DL = self.get_DL_cm()
        L = F_sync * DL * DL * 4.0 * np.pi
        self.set_N_from_L_sync(L)

    def set_N_from_nuLnu(self,nuLnu_src, nu_src):
        """Set the normalization  of N or L_inj (leptonic eq.), to match target source-frame ``nuLnu`` at ``nu_src``.

        Parameters
        ----------
        nuLnu_src : float
            Target source-frame luminosity at ``nu_src`` in ``erg s^-1``.
        nu_src : float
            Source-frame frequency in ``Hz``.
        """
        if self._leptonic_equilibrium:
            self.set_par('L_inj',val=1E40)
            
            self.set_blob()
            delta = self._blob.core.beam_obj
            nu_blob = nu_src / delta
            L_out = BlazarSED.Lum_Sync_at_nu(self._blob, nu_blob) * delta ** 4
            L_out = nuLnu_src / L_out*1E40
            self.set_par('L_inj', val=L_out)
        
        else:
            self.set_par('N',val=1.0)
            self.set_blob()
            delta = self._blob.core.beam_obj
            nu_blob = nu_src / delta
            L_out = BlazarSED.Lum_Sync_at_nu(self._blob, nu_blob) * delta ** 4
            N_out = nuLnu_src / L_out
            self.set_par('N', val=N_out)


    def set_N_from_nuFnu(self, nuFnu_obs, nu_obs):
        """Set the normalization  of N or L_inj (leptonic eq.), to match target observed ``nuFnu`` at ``nu_obs``.

        Parameters
        ----------
        nuFnu_obs : float
            Target observed differential flux in ``erg cm^-2 s^-1``.
        nu_obs : float
            Observed frequency in ``Hz``.
        """

        self.set_blob()
        DL =  self.get_DL_cm()
        L = nuFnu_obs * DL * DL * 4.0 * np.pi
        nu_rest = nu_obs * (1 + self.parameters.z_cosm.val)
        self.set_N_from_nuLnu( L, nu_rest)



    def set_B_eq(self, nuFnu_obs, nu_obs, B_min=1E-9,B_max=1.0,N_pts=20,plot=False):
        """Estimate equipartition magnetic field by scanning a logarithmic B grid, 
        for a given observed flux of the  synchrotron emission (nuFnu_obs) at a given observed frequency (nu_obs)


        Parameters
        ----------
        nuFnu_obs : float
            Target observed differential synchrotron flux in ``erg cm^-2 s^-1``.
        nu_obs : float
            Observed frequency in ``Hz``.
        B_min, B_max : float, optional
            Scan bounds for magnetic field in Gauss.
        N_pts : int, optional
            Number of logarithmically spaced trial values.
        plot : bool, optional
            If ``True``, show diagnostic curves for ``U_e`` and ``U_B``.

        Returns
        -------
        tuple
            ``(B_best, B_grid, U_B, U_e)``.
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
            if self._leptonic_equilibrium:
                self.set_par('L_inj', 1E40)
            else:
                self.set_par('N', 1.0)
            # print 'B_eq',ID
            self.set_N_from_nuFnu(nuFnu_obs, nu_obs)
            if self._leptonic_equilibrium:
                pass
            else:
                N[ID]=self.get_par_by_name('N').val
            self.set_blob()
            #
            U_e[ID] = self._blob.emitters.U_e
            U_B[ID] = self._blob.Sync.UB
            # delta=Jet.get_beaming()
            # print "check L_in=%4.4e L_out=%4.4e"%(L_0,(L_0/delta**4)/BlazarSED.Power_Sync_Electron(Jet._Jet__blob))

        ID_min = np.argmin(U_B + U_e)

        if plot==True:
            #import  pylab as plt
            fig, ax = plt.subplots()
            ax.plot(b_grid,U_e,label=r'$u_e$')
            ax.plot(b_grid,U_B,label=r'$u_B$')
            ax.plot(b_grid,U_B+U_e)
            ax.scatter(b_grid[ID_min],U_e[ID_min]+U_B[ID_min])

            ax.semilogy()
            ax.semilogx()
            ax.set_xlabel(r'$B$ (G)')
            ax.set_ylabel(r'$u$ erg cm$-^3$')
            ax.legend
            plt.show()
           


        self.set_par('B', val=b_grid[ID_min])
        self.set_par('N', val=N[ID_min])
        print ('setting B to ',b_grid[ID_min])
        print ('setting N to ',N[ID_min])
        return b_grid[ID_min],b_grid,U_B,U_e

    def eval_synch_pol(self,nu_range_obs,sin_theta_B_prime=None):
        """Evaluate synchrotron polarization in the observer frame.

        Returns
        -------
        tuple of ndarray
            Polarization degree and corresponding ``nuFnu`` values.
        """
        Sync_kernel_initial=self._blob.core.Sync_kernel
        self._blob.core.Sync_kernel=0
        if sin_theta_B_prime is not None:
            self._blob.Sync.sin_psi=sin_theta_B_prime
            #BlazarSED.InitRadiative(self._blob,0)
        nuF_nu=self.eval(get_model=True,nu=nu_range_obs)
    
        pol_nu=np.zeros(nu_range_obs.size)
        #TODO: this will be removed when eval_Sync_polarization will follow the same pattern of synch flux
        nu_range_pol_blob=nu_range_obs/self.get_beaming()*(1+self.parameters.z_cosm.val)
        for ID,nu in enumerate(nu_range_pol_blob):
            pol_nu[ID]=BlazarSED.eval_Sync_polarization(self._blob,nu)
        
        m=np.logical_or(nuF_nu<=0,np.isnan(pol_nu))
        pol_nu[m]=0
        nuF_nu[m]=0

        Sync_kernel_initial=Sync_kernel_initial
        return pol_nu,nuF_nu
    
    def eval_synch_pol_blob(self,nu_range_blob,sin_theta_B_prime=None):
        """Evaluate synchrotron polarization in the blob frame.

        Parameters
        ----------
        nu_range_blob : array-like or float
            Blob-frame frequency values in Hz. A scalar is accepted and is
            internally promoted to a one-dimensional array.

        Returns
        -------
        pol_nu_blob : ndarray
            Synchrotron polarization degree evaluated at ``nu_range_blob``.
        nuLnu_blob : ndarray
            Synchrotron ``nuLnu`` values in the blob frame at the same
            frequencies.
        """
        Sync_kernel_initial=self._blob.core.Sync_kernel
        if sin_theta_B_prime is not None:
            self._blob.Sync.sin_psi=sin_theta_B_prime
        #TODO: this will be removed when eval_Sync_polarization will follow the same pattern of synch flux
        nu_range_blob = np.atleast_1d(np.asarray(nu_range_blob, dtype=np.float64))
        nu_range_obs=nu_range_blob*self.get_beaming()/(1+self.parameters.z_cosm.val)
        self.eval(nu=nu_range_obs)
        # Work on a local copy to avoid mutating cached SED arrays.
        nuLnu_blob=np.asarray(self.spectral_components.Sync.SED.nuLnu_blob, dtype=np.float64).copy()
        pol_nu_blob=np.zeros(nu_range_blob.size)
        #TODO: this will be removed when eval_Sync_polarization will follow the same pattern of synch flux
        for ID,nu in enumerate(nu_range_blob):
            pol_nu_blob[ID]=BlazarSED.eval_Sync_polarization(self._blob,nu)
        
        m=np.logical_or(nuLnu_blob<=0,~np.isfinite(pol_nu_blob))
        pol_nu_blob[m]=0
        nuLnu_blob[m]=0

        Sync_kernel_initial=Sync_kernel_initial
        return pol_nu_blob,nuLnu_blob


class GalacticBeamed(Jet):

    """Beamed galactic-source specialization of :class:`Jet`."""
    def __init__(self,
                 distance=3*u.kpc,
                 name=None,
                 emitters_type='electrons',
                 emitters_distribution='plc',
                 emitters_distribution_log_values=False,
                 beaming_expr='delta',
                 jet_workplace=None,
                 verbose=False,
                 clean_work_dir=True,
                 electron_distribution=None,
                 proton_distribution=None,
                 electron_distribution_log_values=None,
                 proton_distribution_log_values=None,
                 geometry='spherical'):
        """Initialize a galactic beamed model with fixed luminosity distance."""
        if name is None:
            _name = 'unbeamed'
        else:
            _name = name = clean_var_name(name)

        if np.isscalar(distance):
            _d=distance
        else:
            try:
                _d=distance.to('cm')
            except:
                raise RuntimeError('distance must be either an astropy unit object convertible to cm, or a scalar in cm')

        
        cosmo=Cosmo(DL_cm=_d,verbose=False)
        super(GalacticBeamed,self).__init__(cosmo=cosmo,
                                    name=_name,
                                    emitters_type=emitters_type,
                                    emitters_distribution=emitters_distribution,
                                    emitters_distribution_log_values=emitters_distribution_log_values,
                                    electron_distribution=electron_distribution,
                                    proton_distribution=proton_distribution,
                                    electron_distribution_log_values=electron_distribution_log_values,
                                    proton_distribution_log_values=proton_distribution_log_values,
                                    beaming_expr=beaming_expr,
                                    jet_workplace=jet_workplace,
                                    verbose=verbose,
                                    clean_work_dir=clean_work_dir,
                                    geometry=geometry)

        if name is None or name == '':
            if self.emitters_distribution.emitters_type == 'electrons':
                name = 'galactic_beamed_leptonic'
            elif self.emitters_distribution.emitters_type == 'protons':
                name = 'galactic_beamed_hadronic_pp'
            else:
                name = 'galactic_beamed'

        self.name = clean_var_name(name)
        self._handle_z(d=_d)

    def _dummy_z_par_func(self,DL_cm):
            self.cosmo._DL_cm=DL_cm*u.cm
            return 0

    def _handle_z(self,d=u.kpc*1):
        self._blob.core.z_cosm=0
        self.parameters.z_cosm.val=0
        self.parameters.z_cosm.hidden=True
        self.parameters.z_cosm.frozen=True
        _dmax=1000*u.kpc.to('cm')
        self.add_user_par(name='DL_cm',units='cm',val=d.value,val_min=0,val_max=_dmax)
        p=self.parameters.get_par_by_name('DL_cm')
        p.par_type='distance'


class GalacticUnbeamed(GalacticBeamed):

    """Unbeamed galactic-source specialization of :class:`GalacticBeamed`."""
    def __init__(self,
                 distance=3*u.kpc,
                 name=None,
                 emitters_type='electrons',
                 emitters_distribution='plc',
                 emitters_distribution_log_values=False,
                 jet_workplace=None,
                 verbose=False,
                 clean_work_dir=True,
                 electron_distribution=None,
                 proton_distribution=None,
                 electron_distribution_log_values=None,
                 proton_distribution_log_values=None,
                 geometry='spherical'):
        """Initialize a galactic unbeamed model (beam factor fixed to 1)."""
        super(GalacticUnbeamed,self).__init__(distance=distance,
                                            name=name,
                                            emitters_type=emitters_type,
                                            emitters_distribution=emitters_distribution,
                                            emitters_distribution_log_values=emitters_distribution_log_values,
                                            beaming_expr='delta',
                                            jet_workplace=jet_workplace,
                                            verbose=verbose,
                                            clean_work_dir=clean_work_dir,
                                            electron_distribution=electron_distribution,
                                            proton_distribution=proton_distribution,
                                            electron_distribution_log_values=electron_distribution_log_values,
                                            proton_distribution_log_values=proton_distribution_log_values,
                                            geometry=geometry)

        if name is None or name == '':
            if self.emitters_distribution.emitters_type == 'electrons':
                name = 'galactic_unbeamed_leptonic'
            elif self.emitters_distribution.emitters_type == 'protons':
                name = 'galactic_unbeamed_hadronic_pp'
            else:
                name = 'galactic_unbeamed'
        
        self.name=name


        self.parameters.beam_obj.val=1
        self.parameters.beam_obj.hidden=True
        self.parameters.R_H.hidden=True
    
        self.set_external_field_transf('disk')
  
