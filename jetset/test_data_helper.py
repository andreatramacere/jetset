from __future__ import absolute_import, division, print_function

from builtins import (bytes, str, open, super, range,
                      zip, round, input, int, pow, object, map, zip)

__author__ = "Andrea Tramacere"

import os
import glob

__all__=[]

test_data_path=os.path.dirname(__file__)+'/test_data/'
test_SED_data_path=os.path.dirname(__file__)+'/test_data/SEDs_data'
test_SEDs=[]

for test_SED in glob.glob(test_SED_data_path+'/*.ecsv'):
    test_SEDs.append(test_SED)

test_SEDs.sort()