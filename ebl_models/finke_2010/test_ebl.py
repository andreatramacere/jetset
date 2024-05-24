from astropy.table import Table
from astropy.units import Unit as u
from astropy.units import spectral
import glob
import numpy as np
from astropy.io import fits
import pylab as plt
from scipy import interpolate

data=Table.read('tau_finke_2010.fits')
nu=np.logspace(-1.5,3,30)

def get_tau_funtcion(data_table):
    y = data_table['energies']
    x = np.array(data_table.meta['REDSHIFT'], dtype=np.float64)
    cn = [name for name in data_table.colnames if name != 'energies']
    z = np.array([data_table[n].data for n in cn])
    return interpolate.RectBivariateSpline(x, y, z)



def test_scipy(data,z,plot=False):
    f_tau = get_tau_funtcion(data)
    tau = f_tau(z,nu).T
    print(tau)
    if plot is True:
        plt.ion()
        E = nu * u('Hz').to('TeV', equivalencies=spectral())
        plt.loglog(E,np.exp(-tau),'-')
        plt.grid()
        plt.show()