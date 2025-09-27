__author__ = "Andrea Tramacere"

import os
import numpy as np
#import copy
#import ctypes

on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if on_rtd == True:
    try:
        from .jetkernel import jetkernel as BlazarSED
    except ImportError:
        from .mock import jetkernel as BlazarSED
else:
    from .jetkernel import jetkernel as BlazarSED

from .jetkernel_models_dic import allowed_disk_type

from .jet_paramters import *


__all__=[ 'get_spectral_c_array','get_emitters_c_array1d']



def get_spectral_c_array(x_ptr, y_ptr, blob_object, size):
    x = np.zeros(size)
    y = np.zeros(size)

    for i in range(size):
        x[i] = BlazarSED.get_spectral_array(x_ptr, blob_object, i)
        y[i] = BlazarSED.get_spectral_array(y_ptr, blob_object, i)


    #if deep_copy is True:
    #    x = (ctypes.c_double * size).from_address(int(x_ptr))
    #    y = (ctypes.c_double * size).from_address(int(y_ptr))

    #x=copy.deepcopy(np.asarray(x))
    #y=copy.deepcopy(np.asarray(y)) 
    return x,y
    
    #x = BlazarSED.get_spectral_array_np(x_ptr, blob_object)
    #y = BlazarSED.get_spectral_array_np(y_ptr, blob_object)

    #return np.array(x_ptr),np.array(y_ptr)

def get_emitters_c_array1d(gamma_prt,n_ptr,blob_object,size):
    x = np.zeros(size)
    y = np.zeros(size)
    for ID in range(size):
        x[ID] = BlazarSED.get_elec_array(gamma_prt, blob_object, ID)
        y[ID] = BlazarSED.get_elec_array(n_ptr, blob_object, ID)

    return x,y

    #return np.array(gamma_prt),np.array(n_ptr)