
__author__ = "Andrea Tramacere"

import os
import glob

__all__=[]

_kernel_dir = os.path.join(os.path.dirname(__file__), 'mathkernel')
bessel_table_file_path = os.path.join(_kernel_dir, 'F_Sync.dat')
