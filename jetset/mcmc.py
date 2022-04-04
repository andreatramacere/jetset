__author__ = "Andrea Tramacere"

from .minimizer import  _eval_res

import emcee
from itertools import cycle


import numpy as np
import scipy as sp
from scipy import stats
import corner
import dill as pickle
from multiprocessing import cpu_count, Pool
import multiprocessing as mp
import warnings
import  time
import copy
import threading
from .plot_sedfit import  plt, PlotSED, set_mpl

__all__=['McmcSampler']



class Counter(object):


    def __init__(self,count_tot):
        self.count = 0
        self.count_OK = 0
        self.count_tot = count_tot
        self._progress_iter = cycle(['|', '/', '-', '\\'])


class RunThread(object):
    def __init__(self, target_class):
        self.target_class = target_class

    def run(self):
        self.target_class.model=self.target_class.model.clone()
        self.target_class.sampler.run_mcmc(self.target_class._pos, self.target_class._npernode, rstate0=np.random.get_state(), progress=True,store = True)


class McmcSampler(object):

    def __init__(self,model_minimizer):

        self.model = model_minimizer.fit_model
        self.data = model_minimizer.data
        self.fit_par_free = model_minimizer.fit_par_free

        self._progress_iter = cycle(['|', '/', '-', '\\'])



    def run_sampler(self,
                    nwalkers=500,
                    steps=100,
                    pos=None,
                    burnin=50,
                    use_UL=False,
                    bound=0.2,
                    bound_rel=False,
                    threads=None,
                    walker_start_bound=0.002,
                    use_labels_dict=None,
                    loglog = False,
                    progress='notebook'):

        counter=Counter(nwalkers*steps)
        self.calls = 0
        self.calls_OK = 0
        self.use_UL = use_UL

        self.pos = pos
        self.nwalkers = nwalkers
        self.steps = steps
        self.calls_tot = nwalkers * steps

        if use_labels_dict is None:
            self.par_array =  self.fit_par_free

            self.labels = [par.name for par in self.par_array]
            self.labels_units = [par.units for par in self.par_array]
            self.labels_start_val = [p.best_fit_val for p in self.par_array]
        else:
            self.labels=[]
            self.par_array=[]
            self.labels_units=[]
            self.labels_start_val=[]
            for model_name in use_labels_dict.keys():
                for par_name in use_labels_dict[model_name]:
                    p= self.model.parameters.get_par_by_name(model_name,par_name)
                    if p is not None:
                        self.par_array.append(p)
                        self.labels.append( p.name)
                        self.labels_units.append(p.units)
                        self.labels_start_val.append(p.best_fit_val)
                    else:
                        warnings.warn('par %s'%par_name+' not present in model, will be sckipped')

        self.par_array_best_fit=copy.deepcopy(self.par_array)
        self.ndim = len(self.labels)
        if pos is None:
            pos = emcee.utils.sample_ball(np.array([p.best_fit_val for p in self.par_array]),
                                          np.array([p.best_fit_val * walker_start_bound for p in self.par_array]),
                                          nwalkers)

        self._build_bounds(bound=bound,bound_rel=bound_rel)
        if emcee.__version__ < "3":
            raise RuntimeError('Please update to emcee v>=3.0.0')

        print('mcmc run starting')
        print('')
        start = time.time()
        if progress == 'notebook':
            pass

        if progress is True:
            tqdm=None


        if threads is not None and threads>1:
            # if threads > cpu_count():
            #     threads = cpu_count()
            #     warnings.warn('number of threads has been reduced to cpu_count=%d'%cpu_count())
            #


            # with Pool(processes=threads) as pool:
            #
            #     self.sampler = emcee.EnsembleSampler(nwalkers, self.ndim, log_prob, threads=threads,
            #                                          args=(self.model, self.data, use_UL, counter,self._bounds, self.par_array, loglog),pool=pool)
            #     self.sampler.run_mcmc(pos, steps, progress=True)



            warnings.warn('multithreadign not implemented yet')
            threads=1
            self.sampler = emcee.EnsembleSampler(nwalkers, self.ndim, log_prob, args=(self.model, self.data, use_UL, counter, self._bounds, self.par_array, loglog))
            self.sampler.run_mcmc(pos, steps, progress=progress)
        else:
            threads=1
            self.sampler = emcee.EnsembleSampler(nwalkers, self.ndim, log_prob,args=(self.model, self.data, use_UL, counter,self._bounds, self.par_array, loglog))
            self.sampler.run_mcmc(pos, steps, progress=progress)



        end = time.time()
        comp_time = end - start
        print("mcmc run done, with %d threads took %2.2f seconds"%(threads,comp_time))

        #self.samples = self.sampler.chain[:, burnin:, :].reshape((-1, self.ndim))
        self.samples = self.sampler.get_chain(flat=True,discard=burnin)
        self.samples_log_prob  = self.sampler.get_log_prob(flat=True,discard=burnin)
        self.acceptance_fraction=np.mean(self.sampler.acceptance_fraction)

        self.reset_to_best_fit()



    def get_par_quantiles(self,p,quantiles=(0.16,0.5,0.84)):
        _d, idx=self.get_par(p)

        return np.array(np.quantile(_d,quantiles))

    def reset_to_best_fit(self):
        for ID,par in enumerate(self.par_array):
            par.val = self.par_array_best_fit[ID].val

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

                #_min =  par.best_fit_val - par.best_fit_err * bound[0]
                #_max =  par.best_fit_val + par.best_fit_err * bound[1]
                delta_p = par.best_fit_err * bound[1]
                delta_m = par.best_fit_err * bound[0]

            else:
                delta_p = np.fabs(par.best_fit_val)*bound[1]
                delta_m = np.fabs(par.best_fit_val)*bound[0]
                #_min = par.best_fit_val * (1.0  - bound[0])
                #_max = par.best_fit_val * (1.0  + bound[1])

            _min = par.best_fit_val - delta_m
            _max = par.best_fit_val + delta_p


            if par.fit_range_min is not None:
                _min= max(_min, par.fit_range_min)

            if par.fit_range_max is not None:
                _max= min(_max, par.fit_range_max)

            self._bounds.append([_min, _max])


    def show_pars(self):
        pass

    def corner_plot(self, labels = None, quantiles = (0.16, 0.5, 0.84), levels = None, title_kwargs = {}):
        _id = []

        if labels is None:
            labels=self.labels

        if type(labels) == list:
            pass
        else:
            labels = [labels]

        for l in labels:
            _id.append(self.labels.index(l))

        f = corner.corner(self.samples[:, _id],
                          quantiles=quantiles, labels=[self.labels[i] for i in _id],
                          truths=[self.labels_start_val[i] for i in _id],
                          title_kwargs=title_kwargs,show_titles = True,
                          levels = levels)

        title = 'quantiles ='+str(quantiles)
        f.suptitle(title,y=1.0)
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

        traces=self.sampler.chain[:, :, idx]

        if self.labels_units is not None:
            if self.labels_units[idx] is not None:
                n += ' (%s)' % self.labels_units[idx]

        alpha_true = np.median(_d)

        f = plt.figure()
        ax = f.add_subplot(111)

        if log_plot == True:
            n = 'log10(%s)'%n
            if np.any(_d <= 0):
                raise RuntimeWarning('negative values in p')
            else:
                traces = np.log10(traces)
                alpha_true = np.log10(alpha_true)
        for t in traces:
            ax.plot(t, '-', color='k', alpha=0.5)



        ax.axhline(alpha_true, color='blue')

        ax.set_ylabel(n)
        ax.set_xlabel('steps')
        return f

    def plot_par(self, p, nbins=20, log_plot=False,quantiles=(0.16,0.5,0.84),figsize=None):
        set_mpl()

        _d, idx = self.get_par(p)

        par_name = self.labels[idx]

        x_name = par_name
        if self.labels_units is not None:
            if self.labels_units[idx] is not None and str(self.labels_units[idx]).strip() != '' :
                x_name += ' (%s)' % self.labels_units[idx]

        f = plt.figure(figsize=figsize)
        ax = f.add_subplot(111)

        if log_plot == True:
            x_name = 'log10(%s)' % x_name
            if np.any(_d <= 0):
                raise RuntimeWarning('negative values in p')
            else:
                _d = np.log10(_d)

        q_vals=self.get_par_quantiles(p,quantiles=quantiles)

        q_diff = np.diff(q_vals)

        ax.hist(_d,
                bins=nbins,
                density=True,
                alpha=0.5,
                label=r'%s = $%.3e^{+%.3e}_{-%.3e} $'%(par_name,q_vals[0],q_diff[1],q_diff[0]))

        for q in q_vals:
            if log_plot is True:
                q=np.log10(q)
            ax.axvline(q,c='black',ls='--',lw=0.5)

        ax.set_xlabel(x_name)
        f.suptitle('quantiles = %s'%str(quantiles))
        ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), ncol=1)
        return f

    def plot_model(self, sed_data=None, fit_range=None, size=100, frame='obs', density=False):

        if sed_data is None:
            sed_data=self.sed_data

        if fit_range is None:
            fit_range = [self.model.nu_min_fit, self.model.nu_max_fit]

        p = self.model._set_up_plot(None, sed_data, frame, density)

        self.reset_to_best_fit()
        self.model.eval(fill_SED=True)

        x, y = self.model.SED.get_model_points(log_log=False, frame=frame)
        #if density is True:
        #    y=y-x
        if size is None:
            size = len(self.samples)
            ID_mcmc = np.arange(size)
        else:
            size = min(len(self.samples), int(size))
            ID_mcmc = np.random.randint(len(self.samples), size=size)

        y = np.zeros((size,x.size))

        for ID,ID_rand in enumerate(ID_mcmc):

            for ID_par, pi in enumerate(self.par_array):
                pi.set(val=self.get_par(ID_par)[0][ID_rand])

            self.model.eval(fill_SED=True)
            x, y[ID] = self.model.SED.get_model_points(log_log=False, frame=frame)

        if density is True:
            y=y/x
        y_min=np.amin(y, axis=0)
        y_max=np.amax(y, axis=0)
        msk = y_min > np.log10(self.model.flux_plot_lim)
        p.sedplot.fill_between(x[msk],y_max[msk],y_min[msk],color='gray',alpha=0.3,label='mcmc model range')

        self.reset_to_best_fit()
        self.model.eval(fill_SED=True)

        p.add_model_plot(self.model, color='red',fit_range = fit_range,density=density,flim=self.model.flux_plot_lim)
        p.add_model_residual_plot(model = self.model, data = sed_data, fit_range =  fit_range, color='red')


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



def emcee_log_like(theta,fit_model,data,use_UL,par_array,loglog):
    _warn = False
    for pi in range(len(theta)):
        par_array[pi].set(val=theta[pi])
        if np.isnan(theta[pi]):
            _warn=True

    _m = fit_model.eval(nu=data['x'], fill_SED=False, get_model=True, loglog=loglog)

    _res_sum, _res, _res_UL = _eval_res(data['y'],
                                        _m,
                                        data['dy'],
                                        data['UL'],
                                        use_UL=use_UL)
    #_progess_bar(counter)
    return  _res_sum *-0.5


def _progess_bar(counter):

    if np.mod(counter.count, 10) == 0 and counter.count != 0:
        print("\r%s progress=%3.3f%% calls=%d accepted=%d" % (next( counter._progress_iter),float(100* counter.count)/( counter.count_tot),counter.count,counter.count_OK), end="")



def log_prob(theta,fit_model,data,use_UL,counter,bounds,par_array,loglog):
    lp = log_prior(theta,bounds)
    counter.count += 1
    if not np.isfinite(lp):
        lp = -np.inf
    else:
        lp += emcee_log_like(theta,fit_model,data,use_UL,par_array,loglog)
    counter.count_OK += 1
    return lp




def log_prior(theta,bounds):
    _r=0.
    #bounds = [(par.fit_range_min, par.fit_range_max) for par in model_minimizer.fit_par_free]
    #skip=False
    for pi in range(len(theta)):

        if bounds[pi][1] is not None:
            if theta[pi]>bounds[pi][1]:
                _r=-np.inf

        if bounds[pi][0] is not None:
            if theta[pi]<bounds[pi][0]:
                _r=-np.inf


    #if skip is True:
    #_r=-np.inf
    #print(theta[pi],bounds)
    return _r


