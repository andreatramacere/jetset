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


#class RunThread(object):
#    def __init__(self, target_class):
#        self.target_class = target_class

#    def run(self):
#        self.target_class.model=self.target_class.model.clone()
#        self.target_class.sampler.run_mcmc(self.target_class._pos, self.target_class._npernode, rstate0=np.random.get_state(), progress=True,store = True)

#to prevent from deprecation to error in emcee 
def sample_ball(p0, std, size=1):
    """
    Produce a ball of walkers around an initial parameter value.

    :param p0: The initial parameter value.
    :param std: The axis-aligned standard deviation.
    :param size: The number of samples to produce.

    """
    assert len(p0) == len(std)
    return np.vstack(
        [p0 + std * np.random.normal(size=len(p0)) for i in range(size)]
    )

class McmcSampler(object):

    def __init__(self,model_minimizer):
        if emcee.__version__ < "3":
            raise RuntimeError('Please update to emcee v>=3.0.0')
        
        self.model = model_minimizer.fit_model
        self.data = model_minimizer.data
        self.fit_par_free = model_minimizer.fit_par_free

        self._progress_iter = cycle(['|', '/', '-', '\\'])
        self.labels=None
        self.par_array=None
        self._bounds=None


    def set_labels(self,use_labels_dict=None):
        """_summary_

        Parameters
        ----------
        use_labels_dict : _type_, optional
            _description_, by default None
        """
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
                        warnings.warn('par %s'%par_name+' not present in model, will be skipped')


    def set_bounds(self,bound=0.2,bound_rel=False,preserve_fit_range=True):

        self._build_bounds(bound=bound,bound_rel=bound_rel,preserve_fit_range=preserve_fit_range)
   

        
    def run_sampler(self,
                    nwalkers=None,
                    steps=100,
                    pos=None,
                    burnin=50,
                    use_UL=False,
                    threads=None,
                    walker_start_bound=0.005,
                    loglog = False,
                    progress='notebook',
                    use_labels_dict=None,
                    bound=None,
                    bound_rel=False):

        
        self.calls = 0
        self.calls_OK = 0
        self.use_UL = use_UL

    

        if use_labels_dict is not None or bound is not None:
            warnings.warn('use_labels_dict, and bounds in run_sampler are deprecated and will result in an error in the next version, please use the .set_labels and .set_bounds methods as explained in the documentation')
            self.set_labels(use_labels_dict=use_labels_dict)
            self.set_bounds(bound=bound,bound_rel=bound_rel)
        
        
        self.plot_labels=copy.deepcopy(self.labels)
        self.par_array_best_fit=copy.deepcopy(self.par_array)
        self.ndim = len(self.labels)
        self.pos = pos
        self.burnin=burnin
       
        if nwalkers is None:
            self.nwalkers = 4*len(self.labels)
            print('setting nwalkers to:', self.nwalkers)
        else:
            self.nwalkers = nwalkers
        
        if self.nwalkers < 2*len(self.labels):
            raise RuntimeError("numbers of walkers has to be at least two times the number of sampling pars")

        counter=Counter( self.nwalkers*steps)
        self.steps = steps
        self.calls_tot = self.nwalkers * steps

        if pos is None:
            pos = sample_ball(np.array([p.best_fit_val for p in self.par_array]),
                              np.array([p.best_fit_val * walker_start_bound for p in self.par_array]),
                              self.nwalkers)

        
        print('mcmc run starting')
        print('')
        start = time.time()
        if threads is not None and threads>1:
            warnings.warn('python multithreading is not effective, JetSeT uses C threads to speedup computation')
            threads=1
            self.sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim, log_prob, args=(self.model, self.data, use_UL, counter, self._bounds, self.par_array, loglog))
            self.sampler.run_mcmc(pos, steps, progress=progress)
        else:
            threads=1
            self.sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim, log_prob, args=(self.model, self.data, use_UL, counter, self._bounds, self.par_array, loglog))
            self.sampler.run_mcmc(pos, steps, progress=progress)



        end = time.time()
        comp_time = end - start
        print("mcmc run done, with %d threads took %2.2f seconds"%(threads,comp_time))

        
        self.samples = self.sampler.get_chain(flat=True,discard=burnin)
        self.samples_log_prob  = self.sampler.get_log_prob(flat=True,discard=burnin)
        self.acceptance_fraction=np.mean(self.sampler.acceptance_fraction)
        self.reset_to_best_fit()


    def set_plot_label(self,label,plot_label):
        self.plot_labels[self.labels.index(label)]=plot_label


    def get_par_quantiles(self,p,quantiles=(0.16,0.5,0.84)):
        _d, idx=self.get_par(p)

        return np.array(np.quantile(_d,quantiles))

    def reset_to_best_fit(self):
        for ID,par in enumerate(self.par_array):
            par.val = self.par_array_best_fit[ID].val

    #def reset_to_mcmc(self,quantile=0.5):
    #    for ID,par in enumerate(self.par_array):
    #        q_vals=self.get_par_quantiles(ID,quantiles=quantile)
    #        par.val = q_vals

    def reset_to_mcmc_best_fit(self):
        _prob_max = np.argmax(self.samples_log_prob)
        _id_prob_max = np.unravel_index(_prob_max, self.samples_log_prob.shape)
        
        print("----------------------------")
        print("MCMC best fit solution")
        for par in self.par_array:
            _sample, _idx = self.get_par(par.name)
            par.val = _sample[_id_prob_max]
            print(f"{par.name}: {_sample[_id_prob_max]:}")
        print("----------------------------")

    def _build_bounds(self, bound=0.2,bound_rel=True,preserve_fit_range=True):

        self._bounds=[]

        if np.shape(bound) == ():
            bound=[bound,bound]
        elif np.shape(bound) == (2,):
            pass
        else:
            raise RuntimeError('bound shape', np.shape(bound), 'it is wrong, has to be a scalar or (2,)')

        for par in self.par_array:
            if  not bound_rel  and par.best_fit_err is not None:
                #print('1',par.best_fit_err, bound[1])
                #print('1',par.best_fit_err, bound[0])
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


            if par.fit_range_min is not None and preserve_fit_range is True:
                _min= max(_min, par.fit_range_min)
            elif par.val_min is not None:
                _min= max(_min, par.val_min)

            if par.fit_range_max is not None:
                _max= min(_max, par.fit_range_max)
            elif par.val_max is not None:
                _max= min(_max, par.val_max)
            
            print('par:',par.name,' best fit value: ',par.best_fit_val,' mcmc bounds:',[_min, _max])
            self._bounds.append([_min, _max])



    def corner_plot(self, labels = None, quantiles = (0.16, 0.5, 0.84), levels = None, title_kwargs = {}, **kwargs):
        """_summary_

        Parameters
        ----------
        labels : _type_, optional
            _description_, by default None
        quantiles : tuple, optional
            _description_, by default (0.16, 0.5, 0.84)
        levels :levels=(1 - np.exp(-0.5),)
            check for proper levels definition for 2d: https://corner.readthedocs.io/en/latest/pages/sigmas/
        title_kwargs : dict, optional
            _description_, by default {}

        Returns
        -------
        _type_
            _description_
        """
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
                          quantiles=quantiles, labels=self.plot_labels,
                          truths=[self.labels_start_val[i] for i in _id],
                          title_kwargs=title_kwargs,show_titles = True,
                          levels = levels,**kwargs)

        title = 'quantiles ='+str(quantiles)
        f.suptitle(title,y=1.0)
        return f

    def get_par(self, p,):
        if type(p) == int:
            pass

        else:
            try:
                p = self.labels.index(p)
            except:
                raise RuntimeError('parameter p', p, 'not found')

        if p > len(self.labels):
            raise RuntimeError('label id larger then labels size')

        return self.samples[:, p].flatten(), p


    def plot_chain(self,p=None,log_plot=False):
        if p is None:
            p=self.labels
        else:
            p=np.atleast_1d(p)

        f, axes = plt.subplots(len(p), sharex=True)
        axes=np.atleast_1d(axes)
        for ID,_p in enumerate(p):
            self._plot_chain(_p,axes[ID],log_plot=log_plot)
        
        axes[-1].set_xlabel('steps')
        return f

    def _plot_chain(self, p,ax, log_plot=False):
        _d, idx = self.get_par(p)

        n = self.plot_labels[idx]

        traces=self.sampler.chain[:, :, idx]

        if self.labels_units is not None:
            if self.labels_units[idx] is not None:
                n += ' (%s)' % self.labels_units[idx]

        alpha_true = np.median(_d)

       

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
        if hasattr(self,'burnin'):
            ax.axvline(self.burnin, ls='--',color='orange',alpha=0.5)

        ax.set_ylabel(n)
        
        

    def plot_par(self, p, nbins=20, log_plot=False,quantiles=(0.16,0.5,0.84),figsize=None):
        set_mpl()

        _d, idx = self.get_par(p)

        par_name = self.plot_labels[idx]

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

    def plot_model(self, sed_data=None, fit_range=None, size=100, frame='obs', density=False,quantiles=None, get_model=False, plot_mcmc_best_fit_model=False,rnd_seed=0):
        """_summary_

        Parameters
        ----------
        sed_data : _type_, optional
            _description_, by default None
        fit_range : _type_, optional
            _description_, by default None
        size : int, optional
            _description_, by default 100
        frame : str, optional
            _description_, by default 'obs'
        density : bool, optional
            _description_, by default False
        quantiles : _type_, optional
            _description_, by default None
        get_model : bool, optional
            _description_, by default False
        plot_mcmc_best_fit_model : bool, optional
            _description_, by default False
        rnd_seed : int, optional
            _description_, by default 0

        Returns
        -------
        _type_
            _description_
        """
        if sed_data is None:
            sed_data=self.sed_data

        if fit_range is None:
            fit_range = [self.model.nu_min_fit, self.model.nu_max_fit]

        p = self.model._set_up_plot(None, sed_data, frame, density)
        x,y=self._get_model_samples(size=size,rnd_seed=rnd_seed,frame=frame)
       
        if density is True:
            y=y/x
        
        if quantiles is None:
            y_min=np.amin(y, axis=0)
            y_max=np.amax(y, axis=0)
            _l='mcmc model range'
            #msk = y_min > self.model.flux_plot_lim
            #l=p.sedplot.fill_between(x[msk],y_max[msk],y_min[msk],color='gray',alpha=0.3,label='mcmc model range')
        else:
            _l='mcmc model conf. %s'%quantiles
            y_min,y_max=np.quantile(y, quantiles, axis=0)
        
        msk = y_min > self.model.flux_plot_lim
        l=p.sedplot.fill_between(x[msk],y_max[msk],y_min[msk],color='gray',alpha=0.3,label=_l)

        p.lines_model_list.append(l)
        msk = y_min > self.model.flux_plot_lim
        if plot_mcmc_best_fit_model is False:
            self.reset_to_best_fit()
            label=None
        else:
            label='mcmc best fit'
            self.reset_to_mcmc_best_fit()
        
        self.model.eval(fill_SED=True)
        p.add_model_plot(self.model, color='red',fit_range = fit_range,density=density,flim=self.model.flux_plot_lim,label=label)
        p.add_model_residual_plot(model = self.model, data = sed_data, fit_range =  fit_range, color='red')
        
        if get_model is True:
            return p, [x[msk],y_min[msk],y_max[msk]]
        else:
            return p


    def _get_model_samples(self,size,rnd_seed,frame):
        x, _y = self.model.SED.get_model_points(log_log=False, frame=frame)
        
        if size is None:
            size = len(self.samples)
        else:
            size = min(len(self.samples), int(size))
            
        rng = np.random.default_rng(rnd_seed)
        ID_mcmc = rng.integers(0,len(self.samples),size=size)

        y = np.zeros((size,x.size))

        for ID,ID_rand in enumerate(ID_mcmc):
            for ID_par, pi in enumerate(self.par_array):
                pi.set(val=self.get_par(ID_par)[0][ID_rand])

            self.model.eval(fill_SED=True)
            x, y[ID] = self.model.SED.get_model_points(log_log=False, frame=frame)
        
        return x,y

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

    _model = fit_model.eval(nu=data['x'], fill_SED=False, get_model=True, loglog=loglog)

    _res_sum, _res, _res_UL = _eval_res(data['y'],
                                        _model,
                                        data['dy'],
                                        data['UL'],
                                        use_UL=use_UL)
    return  _res_sum *-0.5


def _progess_bar(counter):

    if np.mod(counter.count, 10) == 0 and counter.count != 0:
        print("\r%s progress=%3.3f%% calls=%d accepted=%d" % (next( counter._progress_iter),float(100* counter.count)/( counter.count_tot),counter.count,counter.count_OK), end="")



def log_prob(theta,fit_model,data,use_UL,counter,bounds,par_array,loglog):
    lp = log_prior(theta,bounds)
    counter.count += 1
    if not np.isfinite(lp):
        res = -np.inf
    else:
        ll= emcee_log_like(theta,fit_model,data,use_UL,par_array,loglog)
        res= lp+ll
    counter.count_OK += 1
    return res




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


