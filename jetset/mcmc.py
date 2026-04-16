"""MCMC sampling utilities built on emcee for JetSeT fit models."""

__author__ = "Andrea Tramacere"

from .minimizer import  _eval_res

import emcee
from itertools import cycle
import types


import numpy as np
import corner
import dill as pickle

import warnings
import  time

from .plot_sedfit import  plt, set_mpl
from .mcmc_parameters import(
    McmcCompositeModelParameterArray,
    get_mcmc_bound_max,
    get_mcmc_bound_min,
    set_mcmc_bound_max,
    set_mcmc_bound_min,
    _check_par_mcmc_bounds
)

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
    def __init__(self,model_minimizer,build_mcmc_parameters=True):
        """Create a new `McmcSampler` instance.
        
        Parameters
        ----------
        model_minimizer : object
            Initialized model-minimizer object.
        """
        if emcee.__version__ < "3":
            raise RuntimeError('Please update to emcee v>=3.0.0')
        self.model = model_minimizer.fit_model.clone()
        self.data = model_minimizer.data
        self._bounds_sampler=[]
        
        self.model.parameters.__class__=McmcCompositeModelParameterArray
        self._progress_iter = cycle(['|', '/', '-', '\\'])
        if build_mcmc_parameters:
            self._bulild_mcmc_paramters()
    
    @property
    def _par_array_sampler(self):
        return [p for p in self.model.parameters.par_array if p.frozen is False]

    @staticmethod
    def _new_progress_iter():
        return cycle(['|', '/', '-', '\\'])

    def __getstate__(self):
        # Keep serialized state free of runtime-only objects not stable across Python versions.
        state = self.__dict__.copy()
        state['_progress_iter'] = None
        state['sampler'] = None
        state['sate_of_mcmc_parameters']={}
        for model in self.model.components.components_list:
            state['sate_of_mcmc_parameters'][model.name]={}
            for p in self.model.parameters.par_array:
                state['sate_of_mcmc_parameters'][model.name][p.name]={}
                state['sate_of_mcmc_parameters'][model.name][p.name]['best_fit_mcmc_val']=p.best_fit_mcmc_val
                state['sate_of_mcmc_parameters'][model.name][p.name]['q_16']=p.q_16
                state['sate_of_mcmc_parameters'][model.name][p.name]['q_50']=p.q_50
                state['sate_of_mcmc_parameters'][model.name][p.name]['q_84']=p.q_84
                state['sate_of_mcmc_parameters'][model.name][p.name]['plot_label']=p.plot_label
                state['sate_of_mcmc_parameters'][model.name][p.name]['mcmc_bound_min']=p.mcmc_bound_min
                state['sate_of_mcmc_parameters'][model.name][p.name]['mcmc_bound_max']=p.mcmc_bound_max
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


    def _bulild_mcmc_paramters(self):
        """dynamically create McmcParameter with proper inheritance"""

        for p in self.model.parameters.par_array:
            Base = p.__class__
            if not Base.__name__.startswith("McmcParameter_"):
                McmcParameter = types.new_class(f"McmcParameter_{Base.__name__}", (Base,))
                McmcParameter.x = property(set_mcmc_bound_max,get_mcmc_bound_max)
                McmcParameter.x = property(set_mcmc_bound_min,get_mcmc_bound_min)
                p.__class__ = McmcParameter
                p._check_par_mcmc_bounds = types.MethodType(_check_par_mcmc_bounds, p)
               

                p.best_fit_mcmc_val=p.val
                p.q_16=None
                p.q_50=None
                p.q_84=None
                if not p.frozen:
                    if p.fit_range_min is not None:
                        p.mcmc_bound_min=p.fit_range_min
                    else:
                        p.mcmc_bound_min=None
                    if  p.fit_range_max is not None:
                        p.mcmc_bound_max= p.fit_range_max 
                    else:
                        p.mcmc_bound_max=None
                else:
                    p.mcmc_bound_min=p.val_min
                    p.mcmc_bound_max=p.val_max
                p.plot_label=p.name
    
   
    @property
    def best_fit_par_table(self):
        return self.model.parameters.best_fit_par_table
    
    
    @property
    def sampler_parameters(self):
        """sampler table.
        
        Returns
        -------
        object
            Requested value.
        """
    
 
        return self.model.parameters._build_sampler_par_table()
        


    def set_bounds(self,par_name=None,comp_name=None,bound=0.2,bound_rel=False,preserve_fit_range=True,par_bounds=None):
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
        #all pars
        if comp_name is None and par_name is None:
            if par_bounds is not None:
                raise RuntimeError('if you pass par_bounds, please provide both model comp_name and par_name')
        
            for par in self._par_array_sampler:
                self._set_bounds(par=par,bound=bound,bound_rel=bound_rel,preserve_fit_range=preserve_fit_range,par_bounds=par_bounds)
            return
        
        elif comp_name is not None and par_name is not None:
            if np.shape(par_bounds)!=(2,):
                raise RuntimeError('please provide par_bounds as [min_bound, max_bound], with min_bound<max_bound')
            if par_bounds[0]>=par_bounds[1]:
                raise RuntimeError('please provide par_bounds as [min_bound, max_bound], with min_bound<max_bound')

            par=self.get_par(par_name,comp_name=comp_name)
            self._set_bounds(par=par,bound=bound,bound_rel=bound_rel,preserve_fit_range=preserve_fit_range,par_bounds=par_bounds)
            return

        else:
            raise RuntimeError('please, provide both par_name and comp_name')
        

    
       
    def _set_bounds(self, par, bound=0.2,bound_rel=True,preserve_fit_range=True,par_bounds=None):

        if par.best_fit_val is None:
            ref_val=par.val
        else:
            ref_val=par.best_fit_val 

        if par_bounds == [] or par_bounds is None:
            if np.shape(bound) == ():
                bound=[bound,bound]
            elif np.shape(bound) == (2,):
                pass
            else:
                raise RuntimeError('bound shape', np.shape(bound), 'it is wrong, has to be a scalar or (2,)')
            
            if  not bound_rel:
                delta_p = ref_val * bound[1]
                delta_m = ref_val  * bound[0]

            else:
                delta_p = np.fabs(ref_val)*bound[1]
                delta_m = np.fabs(ref_val)*bound[0]
        
            _min = ref_val - delta_m
            _max = ref_val + delta_p
        else:
            _min,_max=par_bounds

        if par.fit_range_min is not None and preserve_fit_range is True:
            _min= max(_min, par.fit_range_min )
        elif par.val_min is not None:
            _min= max(_min, par.val_min)

        if par.fit_range_max is not None:
            _max= min(_max, par.fit_range_max )
        elif par.val_max is not None:
            _max= min(_max, par.val_max)
        
        if ref_val>=_max or ref_val<=_min:
            raise RuntimeError(f'please set bounds for par: {par.name} of model comp: {par.model.name} such that  bound_min<{ref_val}<bound_max')
        print('par:',par.name,' ref value: ',ref_val,' mcmc bounds:',[_min, _max])
        par.mcmc_bound_max=_max
        par.mcmc_bound_min=_min

    def _build_sampler_bounds(self):
        self._bounds_sampler=[]
        for par in self._par_array_sampler:
            self._bounds_sampler.append([par.mcmc_bound_min,par.mcmc_bound_max])

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
                    par_idx = [par.name for par in self._par_array_sampler].index(par_name)
                else:
                    par_idx = [par.name+par.model.name for par in self._par_array_sampler].index(par_name+comp_name)
            except:
                raise RuntimeError('parameter p', par_name, 'not found')

        if par_idx > len(self._par_array_sampler):
            raise RuntimeError('label id larger then labels size')

        if get_index is True:
            return self._par_array_sampler[par_idx], par_idx
        else:
            return self._par_array_sampler[par_idx]
        


    def set_plot_label(self,par_name,plot_label,comp_name):
        """Set the display label used for a sampled parameter in plots.

        Parameters
        ----------
        par_name : str
            Name of the parameter whose label must be updated.
        plot_label : str
            New text label used in corner and chain plots.
        comp_name : str
            Name of the model component that owns ``par_name``.

        Raises
        ------
        RuntimeError
            If the requested parameter/component pair is not found.
        """
        p=self.get_par(par_name,comp_name=comp_name)
        p.plot_label=plot_label

    def reset_to_minimizer_best_fit(self):
        """Reset sampled parameter values to minimizer best-fit values.

        Notes
        -----
        This updates only the internal parameter dictionary used by the
        sampler helper; it does not run a new minimization.
        """
       
        for par in self._par_array_sampler:
            if par.best_fit_val is not None:
                par.val = par.best_fit_val


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
        for ID,par in enumerate(self._par_array_sampler):
            par.val= self.get_sample(ID)[_id_prob_max]
            par.best_fit_mcmc_val= par.val
            

            quantiles=self.get_par_quantiles(par_name=par.name,comp_name=par.model.name,quantiles=(0.16,0.5,0.84))
            par.q_16=quantiles[0]
            par.q_50=quantiles[1]
            par.q_84=quantiles[2]
            if verbose:
                print(f"comp: {par.model.name} par: {par.name}  mcmc best fit val: {par.val} quantiles(0.16,0.5,0.84): {quantiles} ")
        if verbose:
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
        to_fix=False
        err_str='\n'
        for par in self._par_array_sampler:
            if par.mcmc_bound_min is None or par.mcmc_bound_max is None:
                err_str+=f"please set mcmc_bound_min and mcmc_bound_max for par: {par.name} in model component: {par.model.name}\n"
                to_fix=True
        if to_fix:
            raise RuntimeError(f"can not run sampler if you do not set bounds for these parameters: {err_str}")

        self._build_sampler_bounds()

        self.calls = 0
        self.calls_OK = 0
        self.use_UL = use_UL
    
       
        self.ndim = len(self._par_array_sampler)
        self.pos = pos
        self.burnin=burnin
       
        if nwalkers is None:
            self.nwalkers = 4*len(self._par_array_sampler)
            print('setting nwalkers to:', self.nwalkers)
        else:
            self.nwalkers = nwalkers
        
        if self.nwalkers < 2*len(self._par_array_sampler):
            raise RuntimeError("numbers of walkers has to be at least two times the number of sampling pars")

        counter=Counter( self.nwalkers*steps)
        self.steps = steps
        self.calls_tot = self.nwalkers * steps

        if pos is None:
            pos = sample_ball(np.array([p.val for p in self._par_array_sampler]),
                              np.array([p.val * walker_start_bound for p in self._par_array_sampler]),
                              self.nwalkers)

        
        print('mcmc run starting')
        print('')
        start = time.time()
        if threads is not None and threads>1:
            warnings.warn('python multithreading is not effective, JetSeT uses C threads to speedup computation')
            threads=1
            self.sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim, log_prob, args=(self.model, self.data, use_UL, counter, self._bounds_sampler, self._par_array_sampler, loglog))
            self.sampler.run_mcmc(pos, steps, progress=progress)
        else:
            threads=1
            self.sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim, log_prob, args=(self.model, self.data, use_UL, counter, self._bounds_sampler, self._par_array_sampler, loglog))
            self.sampler.run_mcmc(pos, steps, progress=progress)



        end = time.time()
        comp_time = end - start
        print("mcmc run done, with %d threads took %2.2f seconds"%(threads,comp_time))

        
        self.samples = self.sampler.get_chain(flat=True,discard=burnin)
        self.samples_log_prob  = self.sampler.get_log_prob(flat=True,discard=burnin)
        self.posterior_weights=None
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

    

    def corner_plot(self, comp_name=None,quantiles = (0.16, 0.5, 0.84), levels = None, title_kwargs = {}, per_component=True, **kwargs):
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
        if per_component:
            if comp_name is None:
                components=np.unique([p.model.name for  p in self._par_array_sampler])
            else:
                components=[comp_name]
            f_list=[]
            for c in components:
                _idxs = []
                truths = []                
                msk=np.array([p.model.name==c for  p in self._par_array_sampler])
                names=[p.name for p in np.array(self._par_array_sampler)[msk]]
                plot_labels=[]
                if msk.sum()>0:
                    for name in names:
                        
                        _idx=self.get_par(name,comp_name=c,get_index=True)[1]        
                        _idxs.append(_idx)
        
                        p=self.get_par(_idx)
                        truths.append(p.best_fit_mcmc_val)
                        if hasattr(p,'plot_label'):
                            plot_labels.append(p.plot_label)
                        else:
                            plot_labels.append(p.name)

                    f=self._do_corner_plot(c,
                                        self.samples,
                                        _idxs,
                                        quantiles,
                                        plot_labels,
                                        truths,
                                        title_kwargs,
                                        levels,
                                        weights=self.posterior_weights,
                                        **kwargs)
                    f_list.append(f)
            return f_list
        
        else:
            f_list=[]
            _idxs = []
            truths = []
            plot_labels=[]
            components=np.unique([p.model.name for  p in self._par_array_sampler])
            for c in components:
                msk=np.array([p.model.name==c for  p in self._par_array_sampler])
                if msk.sum()>0:
                    _names=([p.name for p in np.array(self._par_array_sampler)[msk]])
                    for name in _names:
                        _idx=self.get_par(name,comp_name=c,get_index=True)[1]        
                        _idxs.append(_idx)

                        p=self.get_par(_idx)
                        truths.append(p.best_fit_mcmc_val)
                        if hasattr(p,'plot_label'):
                            plot_labels.append(f'{p.model.name}\n{p.plot_label}')
                        else:
                            plot_labels.append(f'{p.model.name}\n{p.name}')

            f=self._do_corner_plot(self.model.name,
                                self.samples,
                                _idxs,
                                quantiles,
                                plot_labels,
                                truths,
                                title_kwargs,
                                levels,
                                weights=self.posterior_weights,
                                **kwargs)
            f_list.append(f)
            return f_list



    def _do_corner_plot(self, c,samples,_idxs,quantiles,plot_labels,truths,title_kwargs,levels,**kwargs):
        f = corner.corner(samples[:, _idxs],
                        quantiles=quantiles, 
                        labels=plot_labels,
                        truths=truths,
                        title_kwargs=title_kwargs,
                        show_titles = True,
                        levels = levels,
                        **kwargs)

        title = c + ' quantiles ='+str(quantiles)

        f.suptitle(title,y=1.0)

        return f

    
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
            components=np.unique([p.model.name for  p in self._par_array_sampler])
        else:
            components=[comp_name]
        for c in components:
            msk=np.array([p.model.name==c for  p in self._par_array_sampler])
            if par_name is None:
                
                par_names=[par.name for par in  np.array(self._par_array_sampler)[msk]]
            else:
                par_names=np.atleast_1d(par_name)

            f, axes = plt.subplots(len(par_names),figsize=(10, 3*len(par_names)), sharex=True)
            axes=np.atleast_1d(axes)
            for ID,_p_name in enumerate(par_names):
                self._plot_chain(_p_name,axes[ID],comp_name=comp_name,log_plot=log_plot)
            
            axes[-1].set_xlabel('steps')
            f_list.append(f)
            
        return f_list

    def _plot_chain(self, par_name,ax,comp_name=None, log_plot=False):
        par = self.get_par(par_name,comp_name=comp_name)
      
        n = par.plot_label

        traces=self.get_trace(par_name,comp_name=comp_name)

        if par.units is not None:
           n += ' (%s)' % par.units

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

        par_name = par.name

        x_name = par_name
        if par.units is not None:
            x_name += ' (%s)' % par.units

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

    def plot_model(self, sed_data=None, fit_range=None, size=100, frame='obs', density=False,quantiles=None, get_model=False, plot_mcmc_best_fit_model=True,rnd_seed=0,plot_components=False):
        """Plot posterior model envelope and a best-fit reference curve.

        Parameters
        ----------
        sed_data : ObsData, optional
            Observational SED data used for plotting and residuals. If ``None``,
            ``self.sed_data`` is used.
        fit_range : sequence of float, optional
            Two-element fit interval ``[nu_min, nu_max]`` used for model and
            residual overlays. If ``None``, ``[self.model.nu_min_fit,
            self.model.nu_max_fit]`` is used.
        size : int, optional
            Number of posterior samples used to build the shaded model envelope.
        frame : {'obs', 'src'}, optional
            Frame used to evaluate and display model SED values.
        density : bool, optional
            If ``True``, convert sampled ``nuFnu`` curves to ``Fnu`` by dividing
            by frequency before computing the envelope.
        quantiles : tuple of float, optional
            Lower/upper quantiles for the shaded envelope (for example
            ``(0.16, 0.84)``). If ``None``, use the full min/max range.
        get_model : bool, optional
            If ``True``, also return the sampled envelope arrays
            ``[x, y_min, y_max]`` after flux-limit masking.
        plot_mcmc_best_fit_model : bool, optional
            If ``True``, overlay the MCMC best-fit model. If ``False``, overlay
            the minimizer best-fit model.
        rnd_seed : int, optional
            Seed used to draw posterior samples reproducibly.
        plot_components : bool, optional
            If ``True``, plot component/sub-component curves before the final
            reference curve.

        Returns
        -------
        PlotSED or tuple
            Plot object. If ``get_model`` is ``True``, returns
            ``(plot_obj, [x, y_min, y_max])``.
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
        
        if not plot_mcmc_best_fit_model :
            self.reset_to_minimizer_best_fit()
            label=None
        else:
            label='mcmc best fit'
            self.reset_to_mcmc_best_fit(verbose=False)
       
        self.model.eval()

        if plot_components:            
            self.model.plot_model(sed_data=sed_data,plot_obj=p,only_components=True)
       
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
            for id_p,par in enumerate(self._par_array_sampler):
                par.val=self.get_sample(id_p)[ID_rand]
                self.model.parameters.set_par(model_name= par.model.name,par_name=par.name,val=par.val)
                
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
                c.model.parameters.__class__=McmcCompositeModelParameterArray
                if hasattr(c,'model'):
                    c.model=c.model._build_model(c.model)
                for p in c.model.parameters.par_array:
                    Base = p.__class__
                    if not Base.__name__.startswith("McmcParameter_"):
                        McmcParameter = types.new_class(f"McmcParameter_{Base.__name__}", (Base,))
                        McmcParameter.x = property(set_mcmc_bound_max,get_mcmc_bound_max)
                        McmcParameter.x = property(set_mcmc_bound_min,get_mcmc_bound_min)
                        p.__class__ = McmcParameter
                        p._check_par_mcmc_bounds = types.MethodType(_check_par_mcmc_bounds, p)
                        p.best_fit_mcmc_val=c.sate_of_mcmc_parameters[p.model.name][p.name]['best_fit_mcmc_val']
                        p.q_16=c.sate_of_mcmc_parameters[p.model.name][p.name]['q_16']
                        p.q_50=c.sate_of_mcmc_parameters[p.model.name][p.name]['q_50']
                        p.q_84=c.sate_of_mcmc_parameters[p.model.name][p.name]['q_84']
                        p.plot_label=c.sate_of_mcmc_parameters[p.model.name][p.name]['plot_label']
                        p.mcmc_bound_min=c.sate_of_mcmc_parameters[p.model.name][p.name]['mcmc_bound_min']
                        p.mcmc_bound_max=c.sate_of_mcmc_parameters[p.model.name][p.name]['mcmc_bound_max']

                delattr(c,'sate_of_mcmc_parameters')
               
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

    theta = np.asarray(theta, dtype=float)

    # Fast reject before mutating model
    if not np.all(np.isfinite(theta)):
        return -1e100

    try:
        for pi in range(len(theta)):
            
            
            par_array[pi].val=theta[pi]
            fit_model.parameters.set_par(model_name= par_array[pi].model.name,par_name=par_array[pi].name,val=par_array[pi].val)
            if np.isnan(theta[pi]):
                _warn=True
    except Exception:
        # Parameter-setting failure
        return -1e100

    try:
        _model = fit_model.eval(nu=data['x'], fill_SED=False, get_model=True, loglog=loglog)
    except Exception:
        # Parameter-setting failure
        return -1e100

    _model = np.asarray(_model, dtype=float)
    if not np.all(np.isfinite(_model)):
        return -1e100
    try:
        _res_sum, _res, _res_UL = _eval_res(data['y'],
                                            _model,
                                            data['dy'],
                                            data['UL'],
                                            use_UL=use_UL)
    except Exception:
        return -1e100

    if not np.isfinite(_res_sum):
        return -1e100
    
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
