
__author__ = "Andrea Tramacere"

import numpy as np
from .utils import *

__all__=['convert_nu_to_blob','convert_nu_to_src','convert_nuFnu_to_nuLnu_src']

def convert_nu_to_blob(nu,in_frame,delta,z):
    """converts Energies/frequencies from in_frame  to blob rest frame


        Args:
            delta: beaming factor
            z: redshift
            in_frame: startign reference frame, obs=observer (earth), src=corrected for redshift, blob=corrected for redshift and beaming

        Returns:
            E/nu in blob rest frame

    """
    check_frame(frame=in_frame)

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
        unexpected_behaviour()

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
    check_frame(in_frame)
    if in_frame=="obs":
        z_c=z
    elif in_frame=="src":
        z_c=0
    else:
         unexpected_behaviour()

    return nu*(1+z_c)


def convert_nuFnu_to_nuLnu_src(nuFnu,z,in_frame,dl):
    check_frame(in_frame)
    if in_frame == "obs":
       c=4.0*np.pi*dl**2
    elif in_frame == "src":
        c=1.0
    else:
        unexpected_behaviour()

    return c*nuFnu


def convert_nuLnu_src_to_nuLnu_blob(nuLnu,beaming,in_frame):
    check_frame(in_frame)
    c=1/(beaming**4)
    return c*nuLnu


def convert_nu_src_to_nu_blob(nu_src,beaming,in_frame):
    check_frame(in_frame)
    return nu_src/beaming

    
    
    