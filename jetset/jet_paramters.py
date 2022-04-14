
__author__ = "Andrea Tramacere"

from .model_parameters import ModelParameterArray, ModelParameter
from .utils import safe_run
import  numpy as np


__all__=['JetParameter','JetModelDictionaryPar','JetModelParameterArray']


class JetModelDictionaryPar(object):

    def __init__(self,
                 ptype=None,
                 vmin=None,
                 vmax=None,
                 punit=None,
                 froz=False,
                 log=False,
                 val=None,
                 jetkernel_par_name=None,
                 is_in_jetkernel=True,
                 allowed_values=None ):

        self.ptype =ptype
        self.val=val
        self.vmin =vmin
        self.vmax =vmax
        self.punit =punit
        self.froz =froz
        self.log =log
        self.jetkernel_par_name = jetkernel_par_name
        self.is_in_jetkernel=is_in_jetkernel
        self.allowed_values =allowed_values


class JetParameter(ModelParameter):
    """
    This class is a subclass of the :class:`.ModelParameter` class,
    extending the base class to  handles SSC/EC parameters,
    overriding the :meth:`.ModelParameter.set` in order to propagate the
    parameter value to the BlazarSED object instance
    """
    def __init__(self,
                 model,
                 jetkernel_struct_name,
                 jetkernel_parameter_name,
                 is_in_jetkernel=True,
                 **keywords):

        self._model = model
        #self._jetkernel_attr_name = jetkernel_attr_name
        self._jetkernel_parameter_name=  jetkernel_parameter_name
        self._jetkernel_struct_name =   jetkernel_struct_name

        _allowed_par_types = ['scaling_factor',
                              'acceleration_time',
                              'cooling_time',
                              'escape_time',
                              'gamma_grid',
                              'time_grid',
                              'accelerator_width',
                              'fp_coeff_index',
                              'turbulence_scale',
                              'inj_luminosity',
                              'acc_energy',
                              'time_ev_output',
                              'exp_start_time',
                              'emitters_density',
                              'target_density',
                              'Disk',
                              'BLR',
                              'DT',
                              'Star',
                              'region_size',
                              'radius_expansion_index',
                              'region_position',
                              'beta_expansion',
                              'jet_base_height',
                              'electron_energy',
                              'LE_spectral_slope',
                              'HE_spectral_slope',
                              'high-energy-cut-off',
                              'low-energy-cut-off',
                              'spectral_curvature',
                              'turn-over-energy',
                              'magnetic_field',
                              'magnetic_field_index',
                              'jet_opening_angle',
                              'beaming',
                              'jet-viewing-angle',
                              'jet-bulk-factor',
                              'redshift']

        self._is_in_jetkernel=is_in_jetkernel
        super(JetParameter,self).__init__(allowed_par_types=_allowed_par_types,  **keywords)

        if 'val' in keywords.keys():
            val=keywords['val']
            self.assign_val_to_jetkernel(self.name,val)

    @safe_run
    def set(self,**keywords):
        """
        overrides the  :meth:`.ModelParameter.set` method in order to propagate the
        parameter value to the BlazarSED object instance
        """
        super(JetParameter,self).set(**keywords )

        if 'val' in keywords.keys():
            #print('val',keywords['val'],'to',(self.name))
            self.assign_val_to_jetkernel(self.name,keywords['val'])

            #This instruction is not requested since it is performed in
            #parent class
            #if self._depending_par is not None:
            #    self._depending_par.set(val = self._depending_par._func(self._val.val))


    def assign_val_to_jetkernel(self, name, val):
        """
        """
        if self._jetkernel_parameter_name is not None:
            name=self._jetkernel_parameter_name

        #_jetkernel = getattr(self._model, self._jetkernel_attr)
        _jetkernel_struct = getattr(self._model, self._jetkernel_struct_name)

        if self._is_in_jetkernel is True:
            if hasattr(_jetkernel_struct,name):
                b=getattr(_jetkernel_struct,name)



                if type(b)==int:
                    if self._val.islog is True:
                        val=10**val
                    val=int(val)

                elif type(b)==float:
                    if self._val._islog is True:
                        val=10**val
                    val=float(val)

                elif type(b)==str:
                    val=val
                setattr(_jetkernel_struct,name,val)



class JetModelParameterArray(ModelParameterArray):

    def __init__(self,model=None):
        super(JetModelParameterArray, self).__init__(model=model)



    def add_par_from_dict(self,
                          model_dic,
                          model ,
                          struct_name,
                          parameter_class):
        """

        """
        #print('--->',model_dic)
        for key in model_dic.keys():

            pname = key

            jetkernel_par_name = model_dic[key].jetkernel_par_name
            if jetkernel_par_name is  None:
                jetkernel_par_name=pname

            p_test = self.get_par_by_name(pname)


            _jetkernel_struct=getattr(model,struct_name)
            if p_test is None:
                if hasattr(_jetkernel_struct, jetkernel_par_name) and model_dic[key].is_in_jetkernel is True:
                    pval = getattr(_jetkernel_struct, jetkernel_par_name)
                elif model_dic[key].val is not None:
                    pval = model_dic[key].val
                else:
                    raise RuntimeError('par', pname, 'not found in jetkernel and model dict')

                log = model_dic[key].log

                allowed_values = model_dic[key].allowed_values

                # IMPORTANT
                # This has to stay here, because the parameter, even if log, is stored as linear
                if log is True:
                    pval = np.log10(pval)

                ptype = model_dic[key].ptype
                vmin = model_dic[key].vmin
                vmax = model_dic[key].vmax
                punit = model_dic[key].punit

                froz = model_dic[key].froz
                is_in_jetkernel = model_dic[key].is_in_jetkernel

                if (hasattr(model_dic[key],'_is_dependent')):
                    _is_dependent=model_dic[key]._is_dependent
                else:
                    _is_dependent=False

                if (hasattr(model_dic[key],'_depending_pars')):
                    _depending_pars=model_dic[key]._depending_pars
                else:
                    _depending_pars=[]

                if (hasattr(model_dic[key],'_func')):
                    _func=model_dic[key]._func
                else:
                    _func=None

                if (hasattr(model_dic[key],'_master_pars')):
                    _master_pars=model_dic[key]._master_pars
                else:
                    _master_pars=[]

                #print('==>, p',pname,'log_val',log,pval)
                self.add_par(parameter_class(model,
                                             struct_name,
                                             jetkernel_par_name,
                                             name=pname,
                                             par_type=ptype,
                                             val=pval,
                                             val_min=vmin,
                                             val_max=vmax,
                                             units=punit,
                                             frozen=froz,
                                             log=log,
                                             is_in_jetkernel=is_in_jetkernel,
                                             allowed_values=model_dic[key].allowed_values,
                                             _is_dependent=_is_dependent,
                                             _depending_pars= _depending_pars,
                                             _func=_func,
                                             _master_pars=_master_pars))