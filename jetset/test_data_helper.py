
import os
import glob

__all__=[]

test_data_path=os.path.dirname(__file__)+'/test_data/'
test_SED_data_path=os.path.dirname(__file__)+'/test_data/SEDs_data'
test_SEDs=[]

for test_SED in glob.glob(test_SED_data_path+'/*.dat'):
    test_SEDs.append(test_SED)