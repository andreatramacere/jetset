__author__ = "Andrea Tramacere"

import os
import numpy as np
import ctypes

on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if on_rtd:
    from .mock import jetkernel as BlazarSED
else:
    from .jetkernel import jetkernel as BlazarSED


from .jet_paramters import *


__all__=[ 'get_spectral_c_array_read_only','get_emitters_c_array1d','get_emitters_c_array1d_fast','set_emitters_c_array1d','set_emitters_c_array1d_fast']



def get_spectral_c_array_read_only(x_ptr, y_ptr, size):
    

    """Return spectral c array read only.
    
    Parameters
    ----------
    x_ptr : object
        Parameter controlling x ptr.
    y_ptr : object
        Parameter controlling y ptr.
    size : object
        Parameter controlling size.
    
    Returns
    -------
    object
        Requested value.
    """
    x_ctype = (ctypes.c_double * size).from_address(int(x_ptr))
    y_ctype = (ctypes.c_double * size).from_address(int(y_ptr))

    x = np.ctypeslib.as_array(x_ctype)
    y = np.ctypeslib.as_array(y_ctype)
     
    return x.copy(),y.copy()
    

    #x = np.zeros(size)
    #y = np.zeros(size)

    #for i in range(size):
    #    x[i] = BlazarSED.get_spectral_array(x_ptr, blob_object, i)
    #    y[i] = BlazarSED.get_spectral_array(y_ptr, blob_object, i)


    #if deep_copy is True:
    #    x = (ctypes.c_double * size).from_address(int(x_ptr))
    #    y = (ctypes.c_double * size).from_address(int(y_ptr))

    #x=copy.deepcopy(np.asarray(x))
    #y=copy.deepcopy(np.asarray(y)) 
   
    #x = BlazarSED.get_spectral_array_np(x_ptr, blob_object)
    #y = BlazarSED.get_spectral_array_np(y_ptr, blob_object)

    #return np.array(x_ptr),np.array(y_ptr)



def get_emitters_c_array1d(gamma_prt, n_ptr, blob_object, size):
    """Return emitters c array1d.
    
    Parameters
    ----------
    gamma_prt : object
        Frequency/energy control value for gamma prt.
    n_ptr : object
        Parameter controlling n ptr.
    blob_object : object
        Parameter controlling blob object.
    size : object
        Parameter controlling size.
    
    Returns
    -------
    object
        Requested value.
    """
    x = np.zeros(size)
    y = np.zeros(size)
    if size != int(blob_object.emitters.gamma_grid_size):
        raise RuntimeError("mismatch between expected and actual c-array size")
    
    for ID in range(size):
        x[ID] = BlazarSED.get_elec_array(gamma_prt, blob_object, ID)
        y[ID] = BlazarSED.get_elec_array(n_ptr, blob_object, ID)

    return x,y

def get_emitters_c_array1d_fast(gamma_ptr, n_ptr, blob_object, size):
    
    """Return emitters c array1d fast.
    
    Parameters
    ----------
    gamma_ptr : object
        Frequency/energy control value for gamma ptr.
    n_ptr : object
        Parameter controlling n ptr.
    blob_object : object
        Parameter controlling blob object.
    size : object
        Parameter controlling size.
    
    Returns
    -------
    object
        Requested value.
    """
    if size != int(blob_object.emitters.gamma_grid_size):
        raise RuntimeError("mismatch between expected and actual c-array size")

    # hard guard: NULL pointers
    if int(gamma_ptr) == 0 or int(n_ptr) == 0:
        raise RuntimeError("emitters arrays not allocated yet")

    gamma_ct = (ctypes.c_double * size).from_address(int(gamma_ptr))
    n_ct     = (ctypes.c_double * size).from_address(int(n_ptr))

    gamma = np.ctypeslib.as_array(gamma_ct).copy()
    n     = np.ctypeslib.as_array(n_ct).copy()
    return gamma, n



def set_emitters_c_array1d(n_ptr, blob_object, size, values):
    """Set emitters c array1d.
    
    Parameters
    ----------
    n_ptr : object
        Parameter controlling n ptr.
    blob_object : object
        Parameter controlling blob object.
    size : object
        Parameter controlling size.
    values : object
        Parameter controlling values.
    """
    if size != int(blob_object.emitters.gamma_grid_size):
        raise RuntimeError("mismatch between expected and actual c-array size")
    for idx in range(size):
        BlazarSED.set_elec_array(n_ptr,blob_object,values[idx], idx)


def set_emitters_c_array1d_fast(n_ptr, blob_object, size, values):
    """Set emitters c array1d fast.
    
    Parameters
    ----------
    n_ptr : object
        Parameter controlling n ptr.
    blob_object : object
        Parameter controlling blob object.
    size : object
        Parameter controlling size.
    values : object
        Parameter controlling values.
    """
    if size != int(blob_object.emitters.gamma_grid_size):
        raise RuntimeError("mismatch between expected and actual c-array size")
    if int(n_ptr) == 0:
        raise RuntimeError("emitters array not allocated yet")

    arr = np.ascontiguousarray(values, dtype=np.float64)
    if arr.size != size:
        raise ValueError("size mismatch")

    dest = (ctypes.c_double * size).from_address(int(n_ptr))
    ctypes.memmove(dest, arr.ctypes.data, arr.nbytes)


