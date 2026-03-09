"""MCMC sampling utilities built on emcee for JetSeT fit models."""

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
from collections import OrderedDict
from astropy.table import Table
from .plot_sedfit import  plt, PlotSED, set_mpl
from .model_parameters import _show_table

__all__=['McmcSampler']



class Counter(object):


    """Progress counter used during MCMC sampling loops.

    Notes
    -----
    Stores total calls, accepted calls, and a spinner iterator used by
    text-based progress feedback.
    """
    def __init__(self,count_tot):
        """Create a new `Counter` instance.
        
        Parameters
        ----------
        count_tot : object
            Total number of expected iterations/samples.
        """
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
    """Sample ball.
    
    Parameters
    ----------
    p0 : object
        Initial parameter vector.
    std : object
        Per-parameter standard deviations for sampling.
    size : int, optional
        Number of samples or sample size.
    
    Returns
    -------
    object
        Computed value.
    """
    assert len(p0) == len(std)
    return np.vstack(
        [p0 + std * np.random.normal(size=len(p0)) for i in range(size)]
    )

class McmcSampler(object):

    """Run and manage ``emcee`` sampling for a fitted JetSeT model.

    Notes
    -----
    Wraps an input :class:`~jetset.minimizer.ModelMinimizer`, clones its model
    state, stores sampler outputs, and provides serialization-safe state
    handling for chain analysis and plotting.
    """
    def __init__(self,model_minimizer):
        """Create a new `McmcSampler` instance.
        
        Parameters
        ----------
        model_minimizer : object
            Initialized model-minimizer object.
        """
        if emcee.__version__ < "3":
            raise RuntimeError('Please update to emcee v>=3.0.0')
        #self.model_minimizer
        self.model = model_minimizer.fit_model.clone()
        self.data = model_minimizer.data
        #self._fit_par_free = self.model_minimizer.fit_par_free
        self._par_array=None
        self._bounds=None
        
        
        self._progress_iter = cycle(['|', '/', '-', '\\'])
       
    @staticmethod
    def _new_progress_iter():
        return cycle(['|', '/', '-', '\\'])

    def __getstate__(self):
        # Keep serialized state free of runtime-only objects not stable across Python versions.
        state = self.__dict__.copy()
        state['_progress_iter'] = None
        state['sampler'] = None
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._progress_iter = self._new_progress_iter()
        if hasattr(self, 'chain') and self.chain is not None:
            self.chain = self._as_walker_first_chain(self.chain)
        if hasattr(self, 'log_prob_chain') and self.log_prob_chain is not None:
            self.log_prob_chain = self._as_walker_first_log_prob(self.log_prob_chain)

    def _as_walker_first_chain(self, chain):
        _chain = np.asarray(chain)
        if _chain.ndim != 3:
            return _chain

        if hasattr(self, 'nwalkers'):
            if _chain.shape[0] == self.nwalkers:
                return _chain
            if _chain.shape[1] == self.nwalkers:
                return np.swapaxes(_chain, 0, 1)
        return _chain

    def _as_walker_first_log_prob(self, log_prob_chain):
        _logp = np.asarray(log_prob_chain)
        if _logp.ndim != 2:
            return _logp

        if hasattr(self, 'nwalkers'):
            if _logp.shape[0] == self.nwalkers:
                return _logp
            if _logp.shape[1] == self.nwalkers:
                return np.swapaxes(_logp, 0, 1)
        return _logp

    def _cache_sampler_chains(self):
        # Persist full (pre-burnin) traces in walker-first layout for plotting after reload.
        _chain = None
        try:
            _chain = self.sampler.get_chain()
        except Exception:
            _chain = self.sampler.get_chain(flat=False)
        self.chain = self._as_walker_first_chain(_chain)

        try:
            _logp = self.sampler.get_log_prob()
        except Exception:
            _logp = self.sampler.get_log_prob(flat=False)
        self.log_prob_chain = self._as_walker_first_log_prob(_logp)


    def set_labels(self,use_labels_dict=None):
        """Set labels.
        
        Parameters
        ----------
        use_labels_dict : object, optional
            If ``True``, enable labels dict.
        """

        self._par_array=[]
        self._bounds=None
        if use_labels_dict is None:
            for comp in  self.model.components._components_list:
                mod = getattr(self.model,comp.name)
                for p in mod.parameters.par_array:
                    if p.frozen == False:
                        self._append_label_dict(comp.name,p)
        else:
            for model_name in use_labels_dict.keys():
                for par_name in use_labels_dict[model_name]:
                    p= self.model.parameters.get_par_by_name(model_name,par_name)
                    if p is not None:
                        if p.frozen == False:
                            self._append_label_dict(model_name,p)

    def _append_label_dict(self,comp_name,p):
        self._par_array.append(OrderedDict(comp_name=comp_name,
                                            name=p.name,
                                            full_name=f"{comp_name}.{p.name}",
                                            val=p.val,
                                            val_mix=p.val_min,
                                            val_max=p.val_max,
                                            labels_start_val=p.val,
                                            minimizer_best_fit_val=p.best_fit_val,
                                            minimizer_best_fit_err=p.best_fit_err,
                                            minimizer_fit_range_min=p.fit_range_min,
                                            minimizer_fit_range_max=p.fit_range_max,
                                            units=p.units,
                                            plot_label=p.name,
                                            bounds=[None,None])) 

    @property
    def par_table(self):
        """Par table.
        
        Returns
        -------
        object
            Requested value.
        """
        self._build_par_table()
        return self._par_table

    def _build_par_table(self, names_list=None):
        if hasattr(self,'samples_log_prob') and self.samples_log_prob is not None:
            self.reset_to_mcmc_best_fit(verbose=False)
            _prob_max = np.argmax(self.samples_log_prob)
            _id_prob_max = np.unravel_index(_prob_max, self.samples_log_prob.shape)
                
        if self._par_array is None:
            raise RuntimeError('please set labels, using .set_labels, before showing parameters')

        _rows = []
        for idx, par in enumerate(self._par_array):
            if names_list is not None:
                if par['name'] not in names_list and par['full_name'] not in names_list:
                    continue

            _bounds = par.get('bounds', [None, None])
            if isinstance(_bounds, (list, tuple, np.ndarray)) and len(_bounds) == 2:
                _bound_min, _bound_max = _bounds[0], _bounds[1]
            else:
                _bound_min, _bound_max = None, None

            if hasattr(self,'samples_log_prob') and self.samples_log_prob is not None:
                mcmc_best_fit_val = self.get_sample(idx)[_id_prob_max]
                q_016,q_05,q_084=self.get_par_quantiles( par['name'] )
            else:
                mcmc_best_fit_val=None
                q_016,q_05,q_084=[None,None,None]
            _rows.append((
                idx,
                par.get('comp_name'),
                par.get('name'),
                #par.get('full_name'),
                par.get('val'),
                mcmc_best_fit_val,
                q_016,
                q_05,
                q_084,
                #par.get('minimizer_best_fit_val'),
                #par.get('minimizer_best_fit_err'),
                #par.get('minimizer_fit_range_min'),
                #par.get('minimizer_fit_range_max'),
                par.get('val_min', par.get('val_mix')),
                par.get('val_max'),
                _bound_min,
                _bound_max,
                str(par.get('units')),
                par.get('plot_label'),
            ))

        _names = [
            'idx',
            'model name',
            'name',
            #'full name',
            'current val',
            'mcmc best fit val',
            'quantile 0.16',
            'quantile 0.50',
            'quantile 0.84',
            'val min',
            'val max',
            'mcmc bound min',
            'mcmc bound max',
            'units',
            'plot label',
        ]
        self._par_table = Table(rows=_rows, names=_names, masked=False)

    def show_pars(self, getstring=False, names_list=None, sort_key=None):
        """Display pars.
        
        Parameters
        ----------
        getstring : bool, optional
            If ``True``, return text output instead of printing.
        names_list : object, optional
            Ordered list of parameter/component names.
        sort_key : object, optional
            Key used to sort table-like outputs.
        
        Returns
        -------
        object
            Computed value.
        """
        self._build_par_table(names_list=names_list)
        if sort_key is not None:
            self.par_table.sort(sort_key)

        if getstring is True:
            return self.par_table.pformat_all()
        else:
            _show_table(self.par_table)

    @property
    def labels(self):
        """Labels.
        
        Returns
        -------
        object
            Requested value.
        """
        return self.par_table
    
    def set_bounds(self,bound=0.2,bound_rel=False,preserve_fit_range=True):
        """Set bounds.
        
        Parameters
        ----------
        bound : float, optional
            Absolute parameter-bound span.
        bound_rel : bool, optional
            Relative parameter-bound span.
        preserve_fit_range : bool, optional
            Range for preserve fit.
        """
        self._set_bounds(bound=bound,bound_rel=bound_rel,preserve_fit_range=preserve_fit_range)
    
    def _set_bounds(self, bound=0.2,bound_rel=True,preserve_fit_range=True):

        self._bounds=[]

        if np.shape(bound) == ():
            bound=[bound,bound]
        elif np.shape(bound) == (2,):
            pass
        else:
            raise RuntimeError('bound shape', np.shape(bound), 'it is wrong, has to be a scalar or (2,)')

        for par in self._par_array:
            if  not bound_rel  and par['best_fit_err'] is not None:
                delta_p = par['minimizer_best_fit_err'] * bound[1]
                delta_m = par['minimizer_best_fit_err'] * bound[0]

            else:
                delta_p = np.fabs(par['minimizer_best_fit_val'])*bound[1]
                delta_m = np.fabs(par['minimizer_best_fit_val'])*bound[0]
        
            _min = par['minimizer_best_fit_val'] - delta_m
            _max = par['minimizer_best_fit_val'] + delta_p


            if par['minimizer_fit_range_min'] is not None and preserve_fit_range is True:
                _min= max(_min, par['minimizer_fit_range_min'] )
            elif par['val_min'] is not None:
                _min= max(_min, par['val_min'])

            if par['minimizer_fit_range_max'] is not None:
                _max= min(_max, par['minimizer_fit_range_max'])
            elif par['val_max'] is not None:
                _max= min(_max, par['val_max'])
            
            print('par:',par['name'],' best fit value: ',par['minimizer_best_fit_val'],' mcmc bounds:',[_min, _max])
            par['bounds']=[_min, _max]

    def _build_sampler_bounds(self):
        self._bounds=[]
        for par in self._par_array:
            self._bounds.append(par['bounds'])

    def get_par(self, par_name_or_idx, comp_name=None, get_index=False):
        """Return par.
        
        Parameters
        ----------
        par_name_or_idx : object
            Parameter identifier by name or index.
        comp_name : object, optional
            Model-component name.
        get_index : bool, optional
            If ``True``, also return the index of the selected item.
        
        Returns
        -------
        object
            Requested value.
        """
        if type(par_name_or_idx) == int:
            par_idx=par_name_or_idx

        else:
            par_name=par_name_or_idx
            try:
                if comp_name is None:
                    par_idx = [par['name'] for par in self._par_array].index(par_name)
                else:
                    par_idx = [par['name']+par['comp_name'] for par in self._par_array].index(par_name+comp_name)
            except:
                raise RuntimeError('parameter p', par_name, 'not found')

        if par_idx > len(self._par_array):
            raise RuntimeError('label id larger then labels size')

        if get_index is True:
            return self._par_array[par_idx], par_idx
        else:
            return self._par_array[par_idx]
        


    def set_plot_label(self,par_name,plot_label,comp_name=None):
        """Set plot label.
        
        Parameters
        ----------
        par_name : object
            Parameter name.
        plot_label : object
            Custom label used in plots.
        comp_name : object, optional
            Model-component name.
        """
        p=self.get_par(par_name,comp_name=comp_name)
        p['plot_label']=plot_label

    def reset_to_minimizer_best_fit(self):
        """Reset sampled parameter values to minimizer best-fit values.

        Notes
        -----
        This updates only the internal parameter dictionary used by the
        sampler helper; it does not run a new minimization.
        """
        for par in  self._par_array:
            par['val'] = par['minimizer_best_fit_val']



    def reset_to_mcmc_best_fit(self,verbose=True):
        """Reset to mcmc best fit.
        
        Parameters
        ----------
        verbose : bool, optional
            If ``True``, print additional information.
        """
        _prob_max = np.argmax(self.samples_log_prob)
        _id_prob_max = np.unravel_index(_prob_max, self.samples_log_prob.shape)
        
        if verbose:
            print("----------------------------")
            print("MCMC best fit solution")
            for ID,par in enumerate(self._par_array):
                par['val'] = self.get_sample(ID)[_id_prob_max]
                print(f"{par['name']}: {par['val']}")
            print("----------------------------")

        
    def run_sampler(self,
                    nwalkers=None,
                    steps=100,
                    pos=None,
                    burnin=50,
                    use_UL=False,
                    threads=None,
                    walker_start_bound=0.005,
                    loglog = False,
                    progress='notebook'):


        """Run sampler.
        
        Parameters
        ----------
        nwalkers : int, optional
            Number of MCMC walkers.
        steps : int, optional
            Number of MCMC steps.
        pos : object, optional
            Initial walker positions.
        burnin : int, optional
            Number of burn-in steps to discard.
        use_UL : bool, optional
            If ``True``, enable ul.
        threads : object, optional
            Number of worker threads/processes.
        walker_start_bound : float, optional
            Initial spread factor for walker starting points.
        loglog : bool, optional
            If ``True``, operate in log10 space.
        progress : str, optional
            If ``True``, display sampling progress.
        """
        if self._par_array is None:
            raise RuntimeError('please set the labels, using .set_labels, before running the sampler')

        self._build_sampler_bounds()

        self.calls = 0
        self.calls_OK = 0
        self.use_UL = use_UL
    
       
        self.ndim = len(self._par_array)
        self.pos = pos
        self.burnin=burnin
       
        if nwalkers is None:
            self.nwalkers = 4*len(self._par_array)
            print('setting nwalkers to:', self.nwalkers)
        else:
            self.nwalkers = nwalkers
        
        if self.nwalkers < 2*len(self._par_array):
            raise RuntimeError("numbers of walkers has to be at least two times the number of sampling pars")

        counter=Counter( self.nwalkers*steps)
        self.steps = steps
        self.calls_tot = self.nwalkers * steps

        if pos is None:
            pos = sample_ball(np.array([p['minimizer_best_fit_val'] for p in self._par_array]),
                              np.array([p['minimizer_best_fit_val'] * walker_start_bound for p in self._par_array]),
                              self.nwalkers)

        
        print('mcmc run starting')
        print('')
        start = time.time()
        if threads is not None and threads>1:
            warnings.warn('python multithreading is not effective, JetSeT uses C threads to speedup computation')
            threads=1
            self.sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim, log_prob, args=(self.model, self.data, use_UL, counter, self._bounds, self._par_array, loglog))
            self.sampler.run_mcmc(pos, steps, progress=progress)
        else:
            threads=1
            self.sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim, log_prob, args=(self.model, self.data, use_UL, counter, self._bounds, self._par_array, loglog))
            self.sampler.run_mcmc(pos, steps, progress=progress)



        end = time.time()
        comp_time = end - start
        print("mcmc run done, with %d threads took %2.2f seconds"%(threads,comp_time))

        
        self.samples = self.sampler.get_chain(flat=True,discard=burnin)
        self.samples_log_prob  = self.sampler.get_log_prob(flat=True,discard=burnin)
        self._cache_sampler_chains()
        self.acceptance_fraction=np.mean(self.sampler.acceptance_fraction)
        self.reset_to_mcmc_best_fit()


    def get_par_quantiles(self,par_name,comp_name=None,quantiles=(0.16,0.5,0.84)):
        """Return par quantiles.
        
        Parameters
        ----------
        par_name : str
            Parameter name.
        comp_name : object, optional
            Model-component name.
        quantiles : tuple, optional
            Quantiles to evaluate/report.
        
        Returns
        -------
        object
            Requested value.
        """
        return np.array(np.quantile(self.get_sample(par_name,comp_name=comp_name),quantiles))

    

    def corner_plot(self, comp_name=None,quantiles = (0.16, 0.5, 0.84), levels = None, title_kwargs = {}, **kwargs):
        """Corner plot.
        
        Parameters
        ----------
        comp_name : str, optional
            Model-component name.
        quantiles : tuple, optional
            Quantiles to evaluate/report.
        levels : object, optional
            Contour levels for corner plots.
        title_kwargs : dict, optional
            Keyword arguments forwarded to plot-title rendering.
        **kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        object
            Computed value.
        """
       
        if comp_name is None:
            components=np.unique([p['comp_name'] for  p in self._par_array])
        else:
            components=[comp_name]
        f_list=[]
        for c in components:
            _idxs = []
            truths = []
            
            
            msk=np.array([p['comp_name']==c for  p in self._par_array])

            #if labels is None:
            #    plot_labels= [p['plot_label'] for p in np.array(self._par_array)[msk]]
            
            #elif type(labels) == list:
            #    plot_labels=labels
            #else:
            #    plot_labels = [labels]
            names=[p['name'] for p in np.array(self._par_array)[msk]]
            plot_labels=[]
            if msk.sum()>0:
                for name in names:
                    
                    _idxs.append(self.get_par(name,comp_name=c,get_index=True)[1])

    
                for _idx in _idxs:
                    truths.append(self.get_par(_idx)['minimizer_best_fit_val'])
                    plot_labels.append(self.get_par(_idx)['plot_label'])

                f = corner.corner(self.samples[:, _idxs],
                                quantiles=quantiles, 
                                labels=plot_labels,
                                truths=truths,
                                title_kwargs=title_kwargs,
                                show_titles = True,
                                levels = levels,**kwargs)

                
                #print(c,str(quantiles))
                title = c + ' quantiles ='+str(quantiles)

                f.suptitle(title,y=1.0)
                f_list.append(f)
        return f_list

    
    def plot_chain(self,par_name=None, comp_name=None,log_plot=False):
        """Plot chain.
        
        Parameters
        ----------
        par_name : str, optional
            Parameter name.
        comp_name : object, optional
            Model-component name.
        log_plot : bool, optional
            If ``True``, use logarithmic plot scaling where applicable.
        
        Returns
        -------
        object
            Plot object or generated visualization.
        """
        f_list=[]
        if comp_name is None:
            components=np.unique([p['comp_name'] for  p in self._par_array])
        else:
            components=[comp_name]
        for c in components:
            msk=np.array([p['comp_name']==c for  p in self._par_array])
            if par_name is None:
                
                par_names=[par['name'] for par in  np.array(self._par_array)[msk]]
            else:
                par_names=np.atleast_1d(par_name)

            f, axes = plt.subplots(len(par_names), sharex=True)
            axes=np.atleast_1d(axes)
            for ID,_p_name in enumerate(par_names):
                self._plot_chain(_p_name,axes[ID],comp_name=comp_name,log_plot=log_plot)
            
            axes[-1].set_xlabel('steps')
            f_list.append(f)
            
        return f_list

    def _plot_chain(self, par_name,ax,comp_name=None, log_plot=False):
        par = self.get_par(par_name,comp_name=comp_name)

        n = par['plot_label']

        traces=self.get_trace(par_name,comp_name=comp_name)

        if par['units'] is not None:
           n += ' (%s)' % par['units']

        _s=self.get_sample(par_name,comp_name=comp_name)
        alpha_true = np.median(_s)       

        if log_plot == True:
            n = 'log10(%s)'%n
            if np.any(_s <= 0):
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
        
      
    
    def get_trace(self, par_name,comp_name=None):
        """Return trace.
        
        Parameters
        ----------
        par_name : str
            Parameter name.
        comp_name : object, optional
            Model-component name.
        
        Returns
        -------
        object
            Requested value.
        """
        _p,p_idx=self.get_par(par_name,comp_name=comp_name,get_index=True)
        if hasattr(self, 'chain') and self.chain is not None:
            return self.chain[:, :, p_idx]
        if hasattr(self, 'sampler') and self.sampler is not None:
            try:
                _chain = self._as_walker_first_chain(self.sampler.get_chain(flat=False))
                return _chain[:, :, p_idx]
            except Exception:
                _chain = self._as_walker_first_chain(self.sampler.get_chain())
                return _chain[:, :, p_idx]
        raise RuntimeError('MCMC traces are not available in this sampler')


    def plot_par(self, par_name, comp_name=None,nbins=20, log_plot=False,quantiles=(0.16,0.5,0.84),figsize=None):
        """Plot par.
        
        Parameters
        ----------
        par_name : str
            Parameter name.
        comp_name : object, optional
            Model-component name.
        nbins : int, optional
            Number of bins for histogram estimates.
        log_plot : bool, optional
            If ``True``, use logarithmic plot scaling where applicable.
        quantiles : tuple, optional
            Quantiles to evaluate/report.
        figsize : object, optional
            Matplotlib figure size.
        
        Returns
        -------
        object
            Plot object or generated visualization.
        """
        set_mpl()

        par = self.get_par(par_name,comp_name=comp_name)

        par_name = par['name']

        x_name = par_name
        if par['units'] is not None:
            x_name += ' (%s)' % par['units']

        _d=self.get_sample(par_name,comp_name=comp_name)
        

        f = plt.figure(figsize=figsize)
        ax = f.add_subplot(111)
        
        if log_plot == True:
            x_name = 'log10(%s)' % x_name
            if np.any(_d <= 0):
                raise RuntimeWarning('negative values in p')
            else:
                _d = np.log10(_d)

        q_vals=self.get_par_quantiles(par_name,comp_name=comp_name,quantiles=quantiles)

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
    
    def get_sample(self, par_name,comp_name=None):
        """Return sample.
        
        Parameters
        ----------
        par_name : str
            Parameter name.
        comp_name : object, optional
            Model-component name.
        
        Returns
        -------
        object
            Requested value.
        """
        _p,p_idx=self.get_par(par_name,comp_name=comp_name,get_index=True)
        return self.samples[:,p_idx]

    def plot_model(self, sed_data=None, fit_range=None, size=100, frame='obs', density=False,quantiles=None, get_model=False, plot_mcmc_best_fit_model=False,rnd_seed=0):
        """Plot model.
        
        Parameters
        ----------
        sed_data : object, optional
            Observational SED data container.
        fit_range : [float,float], optional
            Range for fit.
        size : int, optional
            Number of samples or sample size.
        frame : str, optional
            Reference frame for data/model values.
        density : bool, optional
            If ``True``, use density representation instead of integrated quantity.
        quantiles : object, optional
            Quantiles to evaluate/report.
        get_model : bool, optional
            If ``True``, return model values.
        plot_mcmc_best_fit_model : bool, optional
            If ``True``, overlay MCMC best-fit model in plots.
        rnd_seed : int, optional
            Random seed used for reproducible sampling.
        
        Returns
        -------
        object
            Plot object or generated visualization.
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
            self.reset_to_minimizer_best_fit()
            label=None
        else:
            label='mcmc best fit'
            self.reset_to_mcmc_best_fit()
        
        self.model.eval(fill_SED=True)
        p.add_model_plot(self.model, color='red',fit_range = fit_range,flim=self.model.flux_plot_lim,label=label)
        p.add_model_residual_plot(model = self.model, data = sed_data, fit_range =  fit_range, color='red')
        self.reset_to_mcmc_best_fit(verbose=False)
        if get_model is True:
            return p, [x[msk],y_min[msk],y_max[msk]]
        else:
            return p


    def _get_model_samples(self,size,rnd_seed,frame):
        self.model.eval()
        x, _y = self.model.SED.get_model_points(log_log=False, frame=frame)
        
        if size is None:
            size = len(self.samples)
        else:
            size = min(len(self.samples), int(size))
            
        rng = np.random.default_rng(rnd_seed)
        ID_mcmc = rng.integers(0,len(self.samples),size=size)

        y = np.zeros((size,x.size))

        for ID,ID_rand in enumerate(ID_mcmc):
            for id_p,par in enumerate(self._par_array):
                par['val']=self.get_sample(id_p)[ID_rand]
                self.model.parameters.set_par(model_name= par['comp_name'],par_name=par['name'],val=par['val'])
                
            self.model.eval(fill_SED=True)
            x, y[ID] = self.model.SED.get_model_points(log_log=False, frame=frame)
        
        return x,y

    def _progess_bar(self,):
        if np.mod(self.calls, 10) == 0 and self.calls != 0:
            print("\r%s progress=%3.3f%% calls=%d accepted=%d" % (next(self._progress_iter),float(100*self.calls)/(self.calls_tot),self.calls,self.calls_OK), end="")

    def save(self, name):
        """Save object state to disk.
        
        Parameters
        ----------
        name : object
            Name identifier.
        """
        with open(name, 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)

    @classmethod
    #def load(self, name):
    #    with open(name, 'rb') as input:
    #        return pickle.load(input)
    def load(cls, file_name):

        """Load object state from disk.
        
        Parameters
        ----------
        file_name : object
            Input/output file path.
        
        Returns
        -------
        object
            Loaded object.
        """
        try:
            c = pickle.load(open(file_name, "rb"))
            if isinstance(c, McmcSampler):
                #c.__init__(c.minimizer)
                if hasattr(c,'model'):
                    c.model=c.model._build_model(c.model)
                return  c
            else:
                raise RuntimeError('The model you loaded is not valid please check the file name')

        except Exception as e:
            raise RuntimeError(e)



def emcee_log_like(theta,fit_model,data,use_UL,par_array,loglog):
    """Emcee log like.
    
    Parameters
    ----------
    theta : object
        Parameter vector sampled by MCMC.
    fit_model : object
        Model instance used for fitting.
    data : object
        Input data table or array.
    use_UL : object
        If ``True``, enable ul.
    par_array : object
        Array/list of model parameters.
    loglog : object
        If ``True``, operate in log10 space.
    
    Returns
    -------
    object
        Computed value.
    """
    _warn = False
    for pi in range(len(theta)):
        
        
        par_array[pi]['val']=theta[pi]
        fit_model.parameters.set_par(model_name= par_array[pi]['comp_name'],par_name=par_array[pi]['name'],val=par_array[pi]['val'])
        if np.isnan(theta[pi]):
            _warn=True

    _model = fit_model.eval(nu=data['x'], fill_SED=False, get_model=True, loglog=loglog)

    _res_sum, _res, _res_UL = _eval_res(data['y'],
                                        _model,
                                        data['dy'],
                                        data['UL'],
                                        use_UL=use_UL)
    return  _res_sum *-0.5


#def _progess_bar(counter):
#
#    if np.mod(counter.count, 10) == 0 and counter.count != 0:
#        print("\r%s progress=%3.3f%% calls=%d accepted=%d" % (next( counter._progress_iter),float(100* counter.count)/( counter.count_tot),counter.count,counter.count_OK), end="")



def log_prob(theta,fit_model,data,use_UL,counter,bounds,par_array,loglog):
    """Log prob.
    
    Parameters
    ----------
    theta : object
        Parameter vector sampled by MCMC.
    fit_model : object
        Model instance used for fitting.
    data : object
        Input data table or array.
    use_UL : object
        If ``True``, enable ul.
    counter : object
        Sampling progress counter object.
    bounds : object
        Bounds for sampled/fitted parameters.
    par_array : object
        Array/list of model parameters.
    loglog : object
        If ``True``, operate in log10 space.
    
    Returns
    -------
    object
        Computed value.
    """
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
    """Log prior.
    
    Parameters
    ----------
    theta : object
        Parameter vector sampled by MCMC.
    bounds : object
        Bounds for sampled/fitted parameters.
    
    Returns
    -------
    object
        Computed value.
    """
    _r=0.

    for pi in range(len(theta)):

        if bounds[pi][1] is not None:
            if theta[pi]>bounds[pi][1]:
                _r=-np.inf

        if bounds[pi][0] is not None:
            if theta[pi]<bounds[pi][0]:
                _r=-np.inf

    return _r
