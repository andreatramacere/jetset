%module jetkernel

%{
    #define SWIG_FILE_WITH_INIT
    // Include your C headers
    #include "../../jetkernel_src/include/Blazar_SED.h"

    // Include NumPy C API header so npy_intp, NPY_DOUBLE, etc. are declared
    #include <numpy/arrayobject.h>
%}

// Bring in NumPy SWIG interface helpers
%include "numpy.i"

// Initialize NumPy C-API when the module loads
%init %{
    import_array();  // REQUIRED: initializes NumPy C API once
%}

// ---------- INPUT TYPEMAPS (Python → C) ----------

// Accept numpy.ndarray (float64) wherever a double* is expected
%typemap(in) double * {
    void *argp = NULL;

    // First try: see if it's already a SWIG pointer to double
    if (SWIG_ConvertPtr($input, &argp, SWIGTYPE_p_double, 0) == SWIG_OK) {
        $1 = (double *) argp;
    }
    else if (PyArray_Check($input)) {
        // Otherwise, try to treat as numpy array
        PyArrayObject *arr = (PyArrayObject *) $input;
        if (PyArray_TYPE(arr) != NPY_DOUBLE) {
            PyErr_SetString(PyExc_TypeError,
                            "Expected numpy array of dtype float64");
            SWIG_fail;
        }
        if (!PyArray_ISCARRAY(arr)) {
            PyErr_SetString(PyExc_TypeError,
                            "Expected C-contiguous numpy array");
            SWIG_fail;
        }
        $1 = (double *) PyArray_DATA(arr);
    }
    else {
        PyErr_SetString(PyExc_TypeError,
                        "Expected a numpy.ndarray or double * pointer");
        SWIG_fail;
    }
}
// ---------- OUTPUT TYPEMAPS (C → Python) ----------

// Apply to ALL fixed-size double arrays in structs
%typemap(out) double [ANY] {
    npy_intp dims[1] = { $1_dim0 };  // number of elements
    $result = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, $1);
}

// Apply to ALL fixed-size int arrays in structs
%typemap(out) int [ANY] {
    npy_intp dims[1] = { $1_dim0 };
    $result = PyArray_SimpleNewFromData(1, dims, NPY_INT, $1);
}

// Include the header with the actual API declarations
%include "../../jetkernel_src/include/Blazar_SED.h"
