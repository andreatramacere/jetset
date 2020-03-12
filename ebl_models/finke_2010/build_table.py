from astropy.table import Table
from astropy.units import Unit as u
from astropy.units import spectral
import glob
import numpy as np
from astropy.io import fits

f_list=glob.glob('tau_modelC_total*.dat')

_c_tau=[]
_n_tau=[]
_tau=[]
z_array=[]
for ID,f_name in enumerate(f_list):
    z=f_name.split('_z')[1]
    z=z.replace('z','')
    z=z.replace('.dat','')
    d = np.loadtxt(f_name)

    if ID==0:
        E = d[:, 0] * u('TeV')
        _c_tau.append(E)
        _n_tau.append('energies')

    _n_tau.append('z_%s'%z)
    z_array.append(z)
    _c_tau.append(d[:,1])


t_tau=Table(data=_c_tau,names=_n_tau,meta={'redshift':z_array})
t_tau.write('tau_finke_2010.fits',format='fits',overwrite=True)


