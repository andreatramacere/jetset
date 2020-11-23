#from __future__ import absolute_import, division, print_function

#from builtins import (bytes, str, open, super, range,
#                      zip, round, input, int, pow, object, map, zip)

__author__ = "Andrea Tramacere"

import os
import glob

__all__=[]

_kernel_dir = os.path.join(os.path.dirname(__file__), 'mathkernel')
bessel_table_file_path = os.path.join(_kernel_dir, 'F_Sync.dat')
