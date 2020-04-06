from astropy.table import Table
from astropy.units import Unit as u
from astropy.units import spectral
import glob
import numpy as np
from astropy.io import fits


_c_tau=[]
_n_tau=[]
_tau=[]

d = np.loadtxt('tau_dominguez10.dat')
z_array= [ 0.01, 0.02526316, 0.04052632, 0.05578947, 0.07105263, 0.08631579, 0.10157895, 0.11684211, 0.13210526,  0.14736842,    0.16263158, 0.17789474,  0.19315789,  0.20842105,  0.22368421,        0.23894737, 0.25421053,  0.26947368,  0.28473684, 0.3 ,  0.35,  0.4 , 0.45,  0.5 ,  0.55,  0.6 ,  0.65,  0.7 ,  0.75,  0.8 ,  0.85,  0.9 ,   0.95,  1.,1.2,1.4,1.6,1.8,2.]
E = d[:, 0] * u('TeV')
print(E,z_array)

_n_tau.append('energies')
_c_tau.append(E)
for ID,z in enumerate(z_array):
    _n_tau.append('z_%1.8f'%z)
    _c_tau.append(d[:,ID+1])

print(len(_n_tau))
t_tau=Table(data=_c_tau,names=_n_tau,meta={'redshift':z_array})
t_tau.write('tau_dominguez_2010.fits',format='fits',overwrite=True)


