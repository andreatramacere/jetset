from __future__ import absolute_import, division, print_function

from builtins import (bytes, str, open, super, range,
                      zip, round, input, int, pow, object, map, zip)

__author__ = "Andrea Tramacere"

from .minimizer import  log_like

import emcee
from itertools import cycle

import numpy as np
import scipy as sp
from scipy import stats
from .plot_sedfit import  plt
import corner
import pickle
from multiprocessing import Pool,cpu_count
import warnings
import  time

__all__=['McmcSampler']



class Counter(object):


    def __init__(self,count_tot):
        self.count = 0
        self.count_OK = 0
        self.count_tot = count_tot
        self._progress_iter = cycle(['|', '/', '-', '\\'])

class McmcSampler(object):

    def __init__(self,model_minimizer):

        self.model_minimizer=model_minimizer
        self.ndim = self.model_minimizer.free_pars
        self._progress_iter = cycle(['|', '/', '-', '\\'])


    def run_sampler(self,nwalkers=500,steps=100,pos=None,burnin=50,use_UL=False, threads=1):


        if threads>cpu_count():
            threads=cpu_count()
            warnings.warn('number of threads has been reduced to cpu_count',cpu_count())

        counter=Counter(nwalkers*steps)


        if pos is None:
            pos = emcee.utils.sample_ball(np.array([p.best_fit_val for p in self.model_minimizer.fit_par_free]),
                                          np.array([p.best_fit_val * 0.005 for p in self.model_minimizer.fit_par_free]),
                                          nwalkers)

        self.calls = 0
        self.calls_OK = 0
        self.use_UL = use_UL

        self.pos = pos
        self.nwalkers = nwalkers
        self.steps = steps
        self.calls_tot = nwalkers * steps
        self.labels = [par.name for par in self.model_minimizer.fit_par_free]
        self.labels_units = [par.units for par in self.model_minimizer.fit_par_free]

        print('mcmc run starting')
        start = time.time()
        if emcee.__version__ >="3":
            with Pool(processes=threads) as pool:

                self.sampler = emcee.EnsembleSampler(nwalkers, self.ndim, log_prob,threads=threads,args=(self.model_minimizer,use_UL,counter),pool=pool)

                self.sampler.run_mcmc(pos,steps,progress=True)
        else:
            self.sampler = emcee.EnsembleSampler(nwalkers, self.ndim, log_prob, threads=threads,
                                                 args=(self.model_minimizer, use_UL, counter), n=pool)

            self.sampler.run_mcmc(pos, steps, progress=True)

        end = time.time()
        comp_time = end - start
        print("mcmc run done, with %d threads took %2.2f seconds"%(threads,comp_time))

        self.samples = self.sampler.chain[:, burnin:, :].reshape((-1, self.ndim))

        self.sampler_out=SamplerOutput(self.samples,
                                       self.sampler.chain,
                                       nwalkers,
                                       steps,
                                       burnin,
                                       labels=self.labels,
                                       labels_units=self.labels_units)

        self.model_minimizer.reset_to_best_fit()

    def plot_par(self,p=None,nbins=20,log_plot=False):
        return self.sampler_out.plot_par(p=p,nbins=nbins,log_plot=log_plot)

    def corner_plot(self, ):
        return self.sampler_out.corner_plot()

    def get_par(self, p=None):
        return self.sampler_out.get_par(p)

    def seve_run(self,name):
        self.sampler_out.save(name)



    def _progess_bar(self,):
        if np.mod(self.calls, 10) == 0 and self.calls != 0:
            print("\r%s progress=%3.3f%% calls=%d accepted=%d" % (next(self._progress_iter),float(100*self.calls)/(self.calls_tot),self.calls,self.calls_OK), end="")



def emcee_log_like(theta,model_minimizer,counter,use_UL):
    _warn = False
    for pi in range(len(theta)):
        model_minimizer.fit_par_free[pi].set(val=theta[pi])
        if np.isnan(theta[pi]):
            _warn=True





    _m = model_minimizer.fit_Model.eval(nu=model_minimizer.nu_fit, fill_SED=False, get_model=True, loglog=model_minimizer.loglog)

    _res_sum, _res, _res_UL= log_like(model_minimizer.nuFnu_fit,
                    _m,
                    model_minimizer.err_nuFnu_fit,
                    model_minimizer.UL,
                    use_UL=use_UL)

    #_progess_bar(counter)
    return  _res_sum


def _progess_bar(counter):

    if np.mod(counter.count, 10) == 0 and counter.count != 0:
        print("\r%s progress=%3.3f%% calls=%d accepted=%d" % (next( counter._progress_iter),float(100* counter.count)/( counter.count_tot),counter.count,counter.count_OK), end="")



def log_prob(theta,model_minimizer,use_UL,counter):
    lp = log_prior(theta,model_minimizer)
    counter.count += 1
    if not np.isfinite(lp):
        return -np.inf
    counter.count_OK += 1
    return lp + emcee_log_like(theta,model_minimizer,counter,use_UL)

def log_prior(theta,model_minimizer):
    _r=0.
    bounds = [(par.fit_range_min, par.fit_range_max) for par in model_minimizer.fit_par_free]
    for pi in range(len(theta)):
        if bounds[pi][1] is not None:
            if theta[pi]<bounds[pi][1]:
                pass
            else:
                _r = -np.inf
        if bounds[pi][0] is not None:
            if theta[pi]>bounds[pi][0]:
                pass
            else:
                _r=-np.inf

    return _r


class SamplerOutput(object):

    def __init__(self,
                 samples,
                 chain,
                 nwalkers,
                 steps,
                 burnin,
                 labels=None,
                 labels_units=None):

        self.samples=samples
        self.chain=chain
        self.labels=labels
        self.nwalkers=nwalkers
        self.steps=steps
        self.burnin=burnin
        self.labels_units=labels_units

    def save(self, name):
        with open(name, 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)

    @classmethod
    def from_file(self,name):
        with open(name, 'rb') as input:

            return pickle.load(input)

    def corner_plot(self,labels=None):
        if labels is not None:
            _id=[]
            if type(labels)==list:
                pass
            else:
                labels=[labels]

            for l in labels:
                _id.append(self.labels.index(l))
            f = corner.corner(self.samples[:,_id], labels=self.labels[_id])
        else:
            f = corner.corner(self.samples, labels=self.labels)
        return f

    def get_par(self,p):
        if type(p) == int:
            pass

        else:
            try:
                p = self.labels.index(p)
            except:
                raise RuntimeError('paramter p', p, 'not found')

        if p > len(self.labels):
            raise RuntimeError('label id larger then labels size')

        return self.samples[:, p].flatten(),p

    def plot_par(self,p,nbins=20,log_plot=False):


        _d,idx = self.get_par(p)

        n = self.labels[idx]

        if self.labels_units is not None:
            if self.labels_units[idx] is not None:
                n += ' (%s)' % self.labels_units[idx]

        f = plt.figure()
        ax=f.add_subplot(111)

        if log_plot==True:
            if  np.any(_d<=0):
                raise RuntimeWarning('negative values in p')
            else:
                _d=np.log10(_d)

        ax.hist(_d,
                bins=nbins,
                density=True)

        ax.set_xlabel(n)

        return f