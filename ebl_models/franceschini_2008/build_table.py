from astropy.table import Table
from astropy.units import Unit as u
from astropy.units import spectral
import glob
import numpy as np
from astropy.io import fits


_c_tau=[]
_n_tau=[]
_tau=[]
z_array=[]
data = np.loadtxt('tau_fran08.dat',usecols=(0,2))
E =  data[0:50,0]* u('TeV')
_n_tau.append('energies')
_c_tau.append(E)
for i in range(int(len(data[:,1])/len(E))):
    tau=data[i*50:i*50+50,1]
    _c_tau.append(tau.flatten())
    z= 1e-3*(i+1.)
    z_array.append( z)
    _n_tau.append('z_%2.8f' % z)
    #print(z,E,tau)

#MUST use ascii, fits forma can not handle table with more than 1000 cols
t_tau=Table(data=_c_tau,names=_n_tau,meta={'redshift':z_array})
t_tau.write('tau_franceschini_2008.dat',format='ascii.ecsv',overwrite=True,)


