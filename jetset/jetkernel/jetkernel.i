/* jetkernel.i */
%module jetkernel
%{
    #define SWIG_FILE_WITH_INIT
#include <stddef.h>
#include <numpy/arrayobject.h>
    /* Put header files here or function declarations like below */
    #include "../../jetkernel_src/include/Blazar_SED.h"
%}

%init %{
import_array();
%}

%inline %{
static PyObject *_jetkernel_numpy_from_double(void *owner,
                                              swig_type_info *owner_type,
                                              double *data,
                                              size_t size) {
    npy_intp dims[1] = { (npy_intp)size };
    return PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, (void *)data);
}

PyObject *get_spectral_array_np(double *arr, struct blob *pt) {
    return _jetkernel_numpy_from_double(pt, SWIGTYPE_p_blob, arr, pt->nu_grid_size);
}

PyObject *get_elec_array_np(double *arr, struct blob *pt) {
    return _jetkernel_numpy_from_double(pt, SWIGTYPE_p_blob, arr, pt->gamma_grid_size);
}

PyObject *get_temp_ev_gamma_array_np(double *arr, struct temp_ev *pt_ev) {
    return _jetkernel_numpy_from_double(pt_ev, SWIGTYPE_p_temp_ev, arr, pt_ev->gamma_grid_size);
}

PyObject *get_temp_ev_time_array_np(double *arr, struct temp_ev *pt_ev) {
    return _jetkernel_numpy_from_double(pt_ev, SWIGTYPE_p_temp_ev, arr, pt_ev->NUM_SET);
}

PyObject *get_temp_ev_N_gamma_array_np(double *arr, struct temp_ev *pt_ev) {
    npy_intp dims[2] = { (npy_intp)pt_ev->NUM_SET, (npy_intp)pt_ev->gamma_grid_size };
    return PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, (void *)arr);
}
%}

/* Parse the header file to generate wrappers */
%include "../../jetkernel_src/include/Blazar_SED.h"
