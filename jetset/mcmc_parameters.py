"""MCMC sampling paramters ."""

__author__ = "Andrea Tramacere"


import ast
from astropy.table import Table
from .model_parameters import _show_table,CompositeModelParameterArray


class McmcCompositeModelParameterArray(CompositeModelParameterArray):

    def _build_sampler_par_table(self, names_list=None):

        _name=[]
        _model_name = []
        _best_fit_mcmc_val=[]
        _q16=[]
        _q50=[]
        _q84=[]
        _val_start=[]
        _mcmc_bound_min=[]
        _mcmc_bound_max=[]
        _frozen=[]
        _val=[]

        _fields=[_model_name,_name,_val,_best_fit_mcmc_val,_q16,_q50,_q84,_mcmc_bound_min,_mcmc_bound_max,_frozen]
        _names=['model name','name','val','bestfit mcmc','q16','q50','q84','mcmc bound min','mcmc bound max','frozen']
      
        #if self.model is None:
        #   _fields.pop(0)
        #   _names.pop(0)

        for par in self.par_array:

            append=False

            if names_list is not None:
                if par.name in names_list:
                    append=True
            else:
                append = True

            if par.hidden is True:
                append=False

            if type(par.val_start) is str:
                try:
                    ast.literal_eval(par.val_start)
                except:
                    append= False

            if append:
                #if self.model is not None:
                _model_name.append(par.model.name)
                #else:
                #    _model_name.append('no_name')

                if par._linked is True and par._linked_root_model != par.model:
                    _p_name = par.name + '(L,%s)' % par._linked_root_model.name

                    _val_start.append(None)
                    _best_fit_mcmc_val.append(None)
                    _q16.append(None)
                    _q50.append(None)
                    _q84.append(None)

                    _mcmc_bound_min.append(par.mcmc_bound_min)
                    _mcmc_bound_max.append(par.mcmc_bound_max)

                else:
                    if par._linked is False and par._depending_pars!=[]:
                        _p_name = par.name + '(M)'
                    else:
                        _p_name = par.name
                    if par.frozen == True:
                        best_fit_mcmc_val = None
                        q16 = None
                        q50 = None
                        q84 = None
                    else:
                        best_fit_mcmc_val=par.best_fit_mcmc_val
                        q16 = par.q_16
                        q50 = par.q_50
                        q84 = par.q_84

                    _val_start.append(par.val_start)
                    _best_fit_mcmc_val.append(best_fit_mcmc_val)
                    _q16.append(q16)
                    _q50.append(q50)
                    _q84.append(q84)

                    _mcmc_bound_min.append(par.mcmc_bound_min)
                    _mcmc_bound_max.append(par.mcmc_bound_max)

                if par._is_dependent is True and par._linked is False:
                    for p in par._master_pars:
                        _p_name = '*' + par.name + '(D,%s)' % p.name
                    if par._linked_root_model is not None:
                        if par._linked_root_model != par.model:
                            _p_name += par._linked_root_model.name

                _name.append(_p_name)
                _frozen.append(par.frozen)
                _val.append(par.val)

        t = Table(_fields, names=_names, masked=False)

       # self._fromat_column_entry(t)

        return t

def get_mcmc_bound_max(self):
    return self._mcmc_bound_max

def set_mcmc_bound_max(self,v):
    self._check_par_mcmc_bounds(v)
    self._mcmc_bound_max=v

def get_mcmc_bound_min(self):
    return self._mcmc_bound_min

def set_mcmc_bound_min(self,v):
    self._check_par_mcmc_bounds(v)
    self._mcmc_bound_min=v

def _check_par_mcmc_bounds(self,v):
    if v is not None:
        if self.val_max is not None:
            if v>self.val_max :
                raise RuntimeError(f'par {self.name} mcmc_bound_max {v} > val_max {self.val_max}')
        if self.val_min is not None:
            if v<self.val_min :
                raise RuntimeError(f'par {self.name} mcmc_bound_max {v} < val_min {self.val_min}')