#!/usr/bin/env python

"""
===================================================================
Moudule: frame_converter
===================================================================

This module contains helper fucntions necessary to perform restframes
conversions




Classes and Inheritance Structure
-------------------------------------------------------------------

.. inheritance-diagram:: BlazarSEDFit.frame_converter
   


Classes relations
----------------------------------------------

.. figure::  classes_frame_converter.png    
   :align:   center     


  
.. autosummary::
   
   
    
Module API
-------------------------------------------------------------------
"""

import math as m
import numpy as np

__all__=['convert_nu_to_blob','convert_nu_to_src']

def convert_nu_to_blob(nu,in_frame,delta,z):
    """converts Energies/frequencies from in_frame  to blob rest frame


        Args:
            delta: beaming factor
            z: redshift
            in_frame: startign reference frame, obs=observer (earth), src=corrected for redshift, blob=corrected for redshift and beaming

        Returns:
            E/nu in blob rest frame

    """
    if in_frame=="obs":
        z_c=z
        delta_c=delta
    elif in_frame=="src":
        z_c=0
        delta_c=delta
    elif in_frame=="blob":
        delta_c=1
        z_c=0
    else:
        raise RuntimeError, "reference frame keyword not valid %s"%in_frame

    return nu*(1+z_c)/delta_c




def convert_nu_to_src(nu,z,in_frame):
    """converts Energies/frequencies from in_frame  to blob rest frame


        Args:
            delta: beaming factor
            z: redshift
            in_frame: startign reference frame, obs=observer (earth), src=corrected for redshift, blob=corrected for redshift and beaming
    

        Returns:
            E/nu in blob rest frame

    """
    if in_frame=="obs":
        z_c=z
    elif in_frame=="src":
        z_c=0
    else:
        raise RuntimeError, "reference frame keyword not valid %s"%in_frame

    return nu*(1+z_c)


def convert_nuFnu_to_nuLnu_src(nuFnu,z,in_frame,dl):
    if in_frame == "obs":
        z_c = z
    elif in_frame == "src":
        z_c = 0
    else:
        raise RuntimeError, "reference frame keyword not valid %s" % in_frame

    return 4.0*np.pi*dl**2



    
    
    
    