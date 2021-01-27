#! /usr/bin/env python27

import numpy as np
from scipy.optimize import leastsq
from .leastsqbound import leastsqbound

def func(p,x):
    """model data as y = m*x+b """
    m,b = p
    return m*np.array(x)+b

def err(p,y,x):
    return y-func(p,x)

# extract data
temp = np.genfromtxt("sample_data.dat")
x = temp[:,0]
y = temp[:,1]

# perform unbounded least squares fitting
p0 = [1.0,0.0]
p, cov_x, infodic, mesg, ier = leastsq(err, p0, args=(y, x), full_output=True)

# print out results 
# print "Standard Least Squares fitting results:"
# print "p:", p
# print "cov_x:", cov_x
# print "infodic['nfev']:", infodic['nfev']
# print "infodic['fvec']:", infodic['fvec']
# print "infodic['fjac']:", infodic['fjac']
# print "infodic['ipvt']:", infodic['ipvt']
# print "infodic['qtf']:", infodic['qtf']
# print "mesg:", mesg
# print "ier:", ier
# print ""

# same as above using no bounds
p0 = [1.0,0.0]
p, cov_x, infodic, mesg, ier = leastsqbound(err, p0, args=(y, x), 
                                                full_output=True)

# print out results 
# print "Bounded Least Squares fitting with no bounds results:"
# print "p:", p
# print "cov_x:", cov_x
# print "infodic['nfev']:", infodic['nfev']
# print "infodic['fvec']:", infodic['fvec']
# print "infodic['fjac']:", infodic['fjac']
# print "infodic['ipvt']:", infodic['ipvt']
# print "infodic['qtf']:", infodic['qtf']
# print "mesg:", mesg
# print "ier:", ier
# print ""


# perform bounded least squares fitting
p0 = [1.0,0.0]
bounds = [(0.0, 2.0), (-10.0, 10.0)] 
p, cov_x, infodic, mesg, ier = leastsqbound(err, p0, args=(y, x),
                                            bounds = bounds, full_output=True)

# print out results 
# print "Bounded Least Squares fitting results:"
# print "p:", p
# print "cov_x:", cov_x
# print "infodic['nfev']:", infodic['nfev']
# print "infodic['fvec']:", infodic['fvec']
# print "infodic['fjac']:", infodic['fjac']
# print "infodic['ipvt']:", infodic['ipvt']
# print "infodic['qtf']:", infodic['qtf']
# print "mesg:", mesg
# print "ier:", ier
