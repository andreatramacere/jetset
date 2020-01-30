from __future__ import absolute_import, division, print_function

from builtins import (str, open, super, range,
                      object, map)




__author__ = "Andrea Tramacere"



import os


on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if on_rtd == True:
    try:
        from jetkernel import jetkernel as BlazarSED
    except ImportError:
        from .mock import jetkernel as BlazarSED
else:
    from jetkernel import jetkernel as BlazarSED


from .jetkernel_models_dic import available_N_distr,N_distr_descr,n_seed_dic,allowed_disk_type

from .jet_paramters import  *

from .jet_emitters import *

__all__=[ 'build_emitting_region_dic','build_ExtFields_dic','BLR_constraints','DT_constraints']








def build_emitting_region_dic(beaming_expr='delta'):
    """

    Builds a dictionary to init the :class:`.JetParameter`
    objects   for the emitting region:

        - **R**, the radius of the emitting region in cm

        - **B**, the magnetic field in G

        - **beaming**, the beaming factor

        - **z**, the redshift


    """

    model_dic = {}
    model_dic['R'] = JetModelDictionaryPar(ptype='region_size', vmin=1E3, vmax=1E30, punit='cm', froz=False, log=False)
    #    ['region_size',0,30,'cm',False,True]
    model_dic['R_H'] = JetModelDictionaryPar(ptype='region_position', vmin=0, vmax=None, punit='cm', froz=True)
    # ['region_position', 0, None, 'cm']
    model_dic['B'] = JetModelDictionaryPar(ptype='magnetic_field', vmin=0, vmax=None, punit='G')
    # ['magnetic_field',0,None,'G']

    if beaming_expr == 'bulk_theta':
        model_dic['theta'] = JetModelDictionaryPar(ptype='jet-viewing-angle', vmin=0, vmax=None, punit='deg')
        # ['jet-viewing-angle',0.0,None,'deg']
        model_dic['BulkFactor'] = JetModelDictionaryPar(ptype='jet-bulk-factor', vmin=1.0, vmax=None,
                                                        punit='Lorentz-factor')
        # ['jet-bulk-factor',1.0,None,'Lorentz-factor']
    elif beaming_expr == 'delta' or beaming_expr == '':
        model_dic['beam_obj'] = JetModelDictionaryPar(ptype='beaming', vmin=1E-4, vmax=None, punit='Lorentz-factor')
        # ['beaming', 1, None, '']
    else:
        raise RuntimeError('''wrong beaming_expr="%s" value, allowed 'delta' or 'bulk_theta' ''' % (beaming_expr))

    model_dic['z_cosm'] = JetModelDictionaryPar(ptype='redshift', vmin=0, vmax=None, punit='')
    # ['redshift',0,None,'']

    return model_dic


def BLR_constraints(L_Disk):
    r1_min = 1E17 * (L_Disk / 1E45) ** 0.5
    r2_min = r1_min * 1.01
    r2_max = r1_min * 2
    print('-> BLR constraints', L_Disk, r1_min, r2_min, r2_max)
    return r1_min, r2_min, r2_max


def DT_constraints(L_Disk):
    return 1E18 * (L_Disk / 1E45) ** 0.5


def build_ExtFields_dic(EC_model_list, allowed_EC_components_list, ):
    """

        """

    model_dic = {}

    for EC_model in EC_model_list:

        # print('EC_model',EC_model)
        # if EC_model not in allowed_EC_components_list:
        #   raise RuntimeError("EC model %s not allowed"%EC_model,"please choose among ", allowed_EC_components_list)

        if 'Disk' in EC_model:
            model_dic['L_Disk'] = JetModelDictionaryPar(ptype='Disk', vmin=0, vmax=None, punit='erg/s')
            # ['Disk',0,None,'erg/s']
            model_dic['R_inner_Sw'] = JetModelDictionaryPar(ptype='Disk', vmin=0, vmax=None, punit='Sw. radii')
            # ['Disk',0,None,'Sw. radii']
            model_dic['R_ext_Sw'] = JetModelDictionaryPar(ptype='Disk', vmin=0, vmax=None, punit='Sw. radii')
            # ['Disk',0,None,'Sw. radii']
            model_dic['T_Disk'] = JetModelDictionaryPar(ptype='Disk', vmin=0, vmax=None, punit='K')
            # ['Disk',0,None,'K']
            model_dic['accr_eff'] = JetModelDictionaryPar(ptype='Disk', vmin=0, vmax=None, punit='')
            # ['Disk',0,None,'']
            model_dic['disk_type'] = JetModelDictionaryPar(ptype='Disk', vmin=None, vmax=None, punit='', froz=True,
                                                           allowed_values=allowed_disk_type)
            # ['Disk',None,None,'',True,False,allowed_disk_type]
            model_dic['M_BH'] = JetModelDictionaryPar(ptype='Disk', vmin=0, vmax=None, punit='M_sun')
            # ['Disk', 0, None, 'M_sun']

        if 'BLR' in EC_model:
            # r1_BLR_min, r2_BLR_min, r2_BLR_max = BLR_constraints(L_Disk)
            model_dic['tau_BLR'] = JetModelDictionaryPar(ptype='BLR', vmin=0, vmax=1.0, punit='')
            # ['BLR',0.0,1.0,'']
            model_dic['R_BLR_in'] = JetModelDictionaryPar(ptype='BLR', vmin=0, vmax=None, punit='cm', froz=True)
            # ['BLR',0,None,'cm',True]
            model_dic['R_BLR_out'] = JetModelDictionaryPar(ptype='BLR', vmin=0, vmax=None, punit='cm', froz=True)
            # ['BLR',0,None,'cm',True]

        if 'DT' in EC_model:
            model_dic['T_DT'] = JetModelDictionaryPar(ptype='DT', vmin=0, vmax=None, punit='K')
            # ['DT',0.0,None,'K']
            model_dic['R_DT'] = JetModelDictionaryPar(ptype='DT', vmin=0, vmax=None, punit='cm')
            # ['DT',0,None,'cm',True]
            model_dic['tau_DT'] = JetModelDictionaryPar(ptype='DT', vmin=0, vmax=1.0, punit='')
            # ['DT',0.0,1.0,'']

    return model_dic
