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
from .plot_sedfit import  plt, PlotSED
import corner
import pickle
from multiprocessing import cpu_count, Pool
import warnings
import  time
import copy

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

        self._progress_iter = cycle(['|', '/', '-', '\\'])

    def run_sampler(self,nwalkers=500,steps=100,pos=None,burnin=50,use_UL=False,bound=0.2,bound_rel=False,threads=None,walker_start_bound=0.002,labels=None):




        counter=Counter(nwalkers*steps)




        self.calls = 0
        self.calls_OK = 0
        self.use_UL = use_UL

        self.pos = pos
        self.nwalkers = nwalkers
        self.steps = steps
        self.calls_tot = nwalkers * steps
        if labels is None:
            self.par_array = copy.deepcopy(self.model_minimizer.fit_par_free)

            self.labels = [par.name for par in self.par_array]
            self.labels_units = [par.units for par in self.par_array]
            self.labels_start_val = [p.best_fit_val for p in self.par_array]
        else:
            self.labels=[]
            self.par_array=[]
            self.labels_units=[]
            self.labels_start_val=[]
            for l in labels:
                p=copy.deepcopy(self.model_minimizer.fit_Model.parameters.get_par_by_name(l))
                self.par_array.append(p)
                self.labels.append( p.name)
                self.labels_units.append(p.units)
                self.labels_start_val.append(p.best_fit_val)

        #print ('->',self.labels)
        self.ndim = len(self.labels)
        if pos is None:
            pos = emcee.utils.sample_ball(np.array([p.best_fit_val for p in self.par_array]),
                                          np.array([p.best_fit_val * walker_start_bound for p in self.par_array]),
                                          nwalkers)

        self._build_bounds(bound=bound,bound_rel=bound_rel)
        print('mcmc run starting')
        start = time.time()
        if threads is not None and threads>1:
            if threads > cpu_count():
                threads = cpu_count()
                warnings.warn('number of threads has been reduced to cpu_count=%d'%cpu_count())

            if emcee.__version__ >="3":
                with Pool(processes=threads) as pool:

                    self.sampler = emcee.EnsembleSampler(nwalkers, self.ndim, log_prob,threads=threads,args=(self.model_minimizer,use_UL,counter,self._bounds, self.par_array),pool=pool)
                    self.sampler.run_mcmc(pos, steps, progress=True)

            else:
                self.sampler = emcee.EnsembleSampler(nwalkers, self.ndim, log_prob, threads=threads,
                                                     args=(self.model_minimizer, use_UL, counter,self._bounds, self.par_array))

                self.sampler.run_mcmc(pos, steps, progress=True)
        else:
            threads=1
            self.sampler = emcee.EnsembleSampler(nwalkers, self.ndim, log_prob,args=(self.model_minimizer, use_UL, counter,self._bounds, self.par_array))
            self.sampler.run_mcmc(pos, steps, progress=True)



        end = time.time()
        comp_time = end - start
        print("mcmc run done, with %d threads took %2.2f seconds"%(threads,comp_time))

        self.samples = self.sampler.chain[:, burnin:, :].reshape((-1, self.ndim))
        self.acceptance_fraction=np.mean(self.sampler.acceptance_fraction)

        self._reset_to_best_fit()

    def _reset_to_best_fit(self):
        for par in self.par_array:
            self.model_minimizer.fit_Model.set(par_name=par.name, val=par.best_fit_val)

    def _build_bounds(self, bound=0.2,bound_rel=True):

        self._bounds=[]

        if np.shape(bound) == ():
            bound=[bound,bound]
        elif np.shape(bound) == (2,):
            pass
        else:
            raise RuntimeError('bound shape', np.shape(bound), 'it is wrong, has to be scalare or (2,)')




        for par in self.par_array:
            if  bound_rel is False and par.best_fit_err is not None:

                _min =  par.best_fit_val - par.best_fit_err * bound[0]
                _max =  par.best_fit_val + par.best_fit_err * bound[1]
            else:
                _min = par.best_fit_val * (1.0  - bound[0])
                _max = par.best_fit_val * (1.0  + bound[1])


            if par.fit_range_min is not None:
                _min= max(_min, par.fit_range_min)

            if par.fit_range_max is not None:
                _max= min(_max, par.fit_range_max)

            self._bounds.append([_min, _max])



    def corner_plot(self, labels=None, quantiles=(0.16, 0.5, 0.84)):
        if labels is not None:
            _id = []
            if type(labels) == list:
                pass
            else:
                labels = [labels]

            for l in labels:
                _id.append(self.labels.index(l))
            f = corner.corner(self.samples[:, _id], quantiles=quantiles, labels=self.labels[_id],
                              truths=self.labels_start_val[_id])
        else:
            f = corner.corner(self.samples, quantiles=quantiles, labels=self.labels, truths=self.labels_start_val)
        return f

    def get_par(self, p):
        if type(p) == int:
            pass

        else:
            try:
                p = self.labels.index(p)
            except:
                raise RuntimeError('paramter p', p, 'not found')

        if p > len(self.labels):
            raise RuntimeError('label id larger then labels size')

        return self.samples[:, p].flatten(), p

    def plot_chain(self, p, log_plot=False):
        _d, idx = self.get_par(p)

        n = self.labels[idx]

        if self.labels_units is not None:
            if self.labels_units[idx] is not None:
                n += ' (%s)' % self.labels_units[idx]

        f = plt.figure()
        ax = f.add_subplot(111)

        if log_plot == True:
            if np.any(_d <= 0):
                raise RuntimeWarning('negative values in p')
            else:
                _d = np.log10(_d)

        ax.plot(_d, '-', color='k', alpha=0.3)
        alpha_true = np.median(_d)
        ax.axhline(alpha_true, color='blue')

        ax.set_ylabel(n)

        return f

    def plot_par(self, p, nbins=20, log_plot=False):

        _d, idx = self.get_par(p)

        n = self.labels[idx]

        if self.labels_units is not None:
            if self.labels_units[idx] is not None:
                n += ' (%s)' % self.labels_units[idx]

        f = plt.figure()
        ax = f.add_subplot(111)

        if log_plot == True:
            if np.any(_d <= 0):
                raise RuntimeWarning('negative values in p')
            else:
                _d = np.log10(_d)

        ax.hist(_d,
                bins=nbins,
                density=True)

        ax.set_xlabel(n)

        return f

    def plot_model(self, sed_data, fit_range=None, size=10, labels=None,frame='obs'):
        if labels is None:
            labels = self.labels

        p = PlotSED()
        p.add_data_plot(sed_data, fit_range=fit_range)

        self._reset_to_best_fit()
        self.model_minimizer.fit_Model.eval(fill_SED=True)

        x, y = self.model_minimizer.fit_Model.SED.get_model_points(log_log=True, frame=frame)
        y = np.zeros((size,x.size))

        for ID,ID_rand in enumerate(np.random.randint(len(self.samples), size=size)):
            for par_name in labels:
                self.model_minimizer.fit_Model.set(par_name=par_name, val=self.get_par(par_name)[0][ID_rand])

            self.model_minimizer.fit_Model.eval()
            x, y[ID] = self.model_minimizer.fit_Model.SED.get_model_points(log_log=True, frame=frame)


        y_min=np.amin(y, axis=0)
        y_max=np.amax(y, axis=0)
        p.sedplot.fill_between(x,y_max,y_min,color='gray',alpha=0.3,label='mcmc model range')

        self._reset_to_best_fit()
        self.model_minimizer.fit_Model.eval(fill_SED=True)
        p.add_model_plot(self.model_minimizer.fit_Model, color='red')
        p.add_residual_plot(self.model_minimizer.fit_Model, sed_data, fit_range=fit_range, color='red')

        return p


    def _progess_bar(self,):
        if np.mod(self.calls, 10) == 0 and self.calls != 0:
            print("\r%s progress=%3.3f%% calls=%d accepted=%d" % (next(self._progress_iter),float(100*self.calls)/(self.calls_tot),self.calls,self.calls_OK), end="")

    def save(self, name):
        with open(name, 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)

    @classmethod
    def load(self, name):
        with open(name, 'rb') as input:
            return pickle.load(input)



def emcee_log_like(theta,model_minimizer,counter,use_UL,par_array):
    _warn = False
    for pi in range(len(theta)):
        par_array[pi].set(val=theta[pi])
        if np.isnan(theta[pi]):
            _warn=True





    _m = model_minimizer.fit_Model.eval(nu=model_minimizer.nu_fit, fill_SED=False, get_model=True, loglog=model_minimizer.loglog)

    _res_sum, _res, _res_UL = log_like(model_minimizer.nuFnu_fit,
                                       _m,
                                       model_minimizer.err_nuFnu_fit,
                                       model_minimizer.UL,
                                       use_UL=use_UL)

    #_progess_bar(counter)
    return  _res_sum


def _progess_bar(counter):

    if np.mod(counter.count, 10) == 0 and counter.count != 0:
        print("\r%s progress=%3.3f%% calls=%d accepted=%d" % (next( counter._progress_iter),float(100* counter.count)/( counter.count_tot),counter.count,counter.count_OK), end="")



def log_prob(theta,model_minimizer,use_UL,counter,bounds,par_array):
    lp = log_prior(theta,bounds)
    counter.count += 1
    if not np.isfinite(lp):
        pass
    else:
        lp += emcee_log_like(theta,model_minimizer,counter,use_UL,par_array)
    counter.count_OK += 1
    return lp




def log_prior(theta,bounds):
    _r=0.
    #bounds = [(par.fit_range_min, par.fit_range_max) for par in model_minimizer.fit_par_free]
    skip=False
    for pi in range(len(theta)):

        if bounds[pi][1] is not None:
            if theta[pi]>bounds[pi][1]:
                skip=True
            else:
                pass
                #_r = -np.inf
        if bounds[pi][0] is not None:
            if theta[pi]<bounds[pi][0]:
                skip=True
            else:
                pass
                #_r=-np.inf

    if skip is True:
        _r=-np.inf

    return _r


# class SamplerOutput(object):
#
#     def __init__(self,
#                  samples,
#                  chain,
#                  nwalkers,
#                  steps,
#                  burnin,
#                  labels_start_val,
#                  acceptance_fraction,
#                  labels,
#                  labels_units=None,):
#
#         self.samples=samples
#         self.chain=chain
#         self.labels=labels
#         self.labels_start_val=labels_start_val
#         self.nwalkers=nwalkers
#         self.steps=steps
#         self.burnin=burnin
#         self.labels_units=labels_units
#         self.acceptance_fraction=acceptance_fraction
#
#     def save(self, name):
#         with open(name, 'wb') as output:
#             pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)
#
#     @classmethod
#     def from_file(self,name):
#         with open(name, 'rb') as input:
#
#             return pickle.load(input)
#
#     def corner_plot(self,labels=None,quantiles=(0.16, 0.5, 0.84)):
#         if labels is not None:
#             _id=[]
#             if type(labels)==list:
#                 pass
#             else:
#                 labels=[labels]
#
#             for l in labels:
#                 _id.append(self.labels.index(l))
#             f = corner.corner(self.samples[:,_id], quantiles=quantiles, labels=self.labels[_id],truths=self.labels_start_val[_id])
#         else:
#             f = corner.corner(self.samples, quantiles=quantiles,labels=self.labels, truths=self.labels_start_val)
#         return f
#
#     def get_par(self,p):
#         if type(p) == int:
#             pass
#
#         else:
#             try:
#                 p = self.labels.index(p)
#             except:
#                 raise RuntimeError('paramter p', p, 'not found')
#
#         if p > len(self.labels):
#             raise RuntimeError('label id larger then labels size')
#
#         return self.samples[:, p].flatten(),p
#
#     def plot_chain(self,p,log_plot=False):
#         _d, idx = self.get_par(p)
#
#         n = self.labels[idx]
#
#         if self.labels_units is not None:
#             if self.labels_units[idx] is not None:
#                 n += ' (%s)' % self.labels_units[idx]
#
#         f = plt.figure()
#         ax = f.add_subplot(111)
#
#         if log_plot == True:
#             if np.any(_d <= 0):
#                 raise RuntimeWarning('negative values in p')
#             else:
#                 _d = np.log10(_d)
#
#         ax.plot(_d, '-', color='k', alpha=0.3)
#         alpha_true=np.median(_d)
#         ax.axhline(alpha_true, color='blue')
#
#         ax.set_ylabel(n)
#
#         return f
#
#
#     def plot_par(self,p,nbins=20,log_plot=False):
#
#
#         _d,idx = self.get_par(p)
#
#         n = self.labels[idx]
#
#         if self.labels_units is not None:
#             if self.labels_units[idx] is not None:
#                 n += ' (%s)' % self.labels_units[idx]
#
#         f = plt.figure()
#         ax=f.add_subplot(111)
#
#         if log_plot==True:
#             if  np.any(_d<=0):
#                 raise RuntimeWarning('negative values in p')
#             else:
#                 _d=np.log10(_d)
#
#         ax.hist(_d,
#                 bins=nbins,
#                 density=True)
#
#         ax.set_xlabel(n)
#
#         return f
#
#     def plot_model(self,sed_data,model_minimizer,fit_range=None,size=100):
#
#         p = PlotSED()
#         p.add_data_plot(sed_data, fit_range=fit_range)
#         p.add_model_plot(model_minimizer.fit_Model, color='red')
#         p.add_residual_plot(model_minimizer.fit_Model, sed_data, fit_range=fit_range, color='red')
#         for ID in np.random.randint(len(self.samples), size=size):
#             for par_name in range(len(self.labels)):
#                 model_minimizer.fit_Model.set(par_name=par_name, value=self.get_par(par_name)[ID])
#
#             model_minimizer.fit_Model.eval()
#             p.add_model_plot(model_minimizer.fit_Model, color='gray')
#
#         p.rescale(y_min=-13, x_min=6, x_max=28.5)
#
#         model_minimizer.reset_to_best_fit()
#
#         return p