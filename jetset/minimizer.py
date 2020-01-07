from __future__ import absolute_import, division, print_function

from builtins import (bytes, str, open, super, range,
                      zip, round, input, int, pow, object, map, zip)

__author__ = "Andrea Tramacere"

from itertools import cycle


import scipy as sp

import numpy as np

import os

import sys

from scipy.stats import chi2

try:
    import iminuit
    minuit_installed=True

except:
    minuit_installed=False
#from iminuit.frontends import ConsoleFrontend
#from iminuit.frontends import console
import  pylab as plt
from scipy.optimize import leastsq

from scipy.optimize import least_squares,curve_fit

from leastsqbound.leastsqbound import  leastsqbound

from .output import section_separator,WorkPlace,makedir





__all__=['FitResults','fit_SED','Minimizer','LSMinimizer','LSBMinimizer','MinutiMinimizer','ModelMinimizer']





class FitResults(object):
    """
    Class to store the fit results 
    
    Parameters
     
    
    Members
    :ivar fit_par: fit_par
    :ivar info: info
    :ivar mesg: mesg
    :ivar success: success
    :ivar chisq: chisq
    :ivar dof: dof
    :ivar chisq_red: chisq_red
    :ivar null_hyp_sig: null_hyp_sig
    :ivar fit_report: ivar get_report()
    -------
    
    """
    
    def __init__(self,
                 name,
                 mm,
                 parameters,
                 calls,
                 mesg,
                 success,
                 chisq,
                 dof,
                 chisq_red,
                 null_hyp_sig,
                 wd,
                 chisq_no_UL=None,
                 dof_no_UL=None,
                 chisq_red_no_UL=None,
                 null_hyp_sig_no_UL=None):

        self.name=name
        self.parameters=parameters
        self.mm=mm
        self.calls=calls
        self.mesg=mesg
        self.success=success
        self.chisq=chisq
        self.dof=dof
        self.chisq_red=chisq_red
        self.null_hyp_sig=null_hyp_sig

        self.chisq_no_UL = chisq_no_UL
        self.dof_no_UL = dof_no_UL
        self.chisq_red_no_UL = chisq_red_no_UL
        self.null_hyp_sig_no_UL = null_hyp_sig_no_UL

        self.fit_report=self.get_report()
        self.wd=wd
         

    def get_report(self):
        out=[]
        out.append("")        
        out.append("**************************************************************************************************")
        out.append("Fit report")
        out.append("")
        out.append("Model: %s"%self.name)
        #pars_rep=self.parameters.show_pars(getstring=True)
        #for string in pars_rep:
        #    out.append(string)
        self.parameters._build_par_table()
        self.model_table=self.parameters._par_table
        self.parameters._build_best_fit_par_table()
        self.bestfit_table = self.parameters._best_fit_par_table
        for _l in self.model_table.pformat_all():
            out.append(_l)
        out.append("")
        out.append("converged=%s"%self.success)
        out.append("calls=%d"%self.calls)
        try:
            out.append("mesg=",self.mesg)
        except:
            out.append(self.mesg)
        out.append("dof=%d"%self.dof)
        out.append("chisq=%f, chisq/red=%f null hypothesis sig=%f"%(self.chisq,self.chisq_red,self.null_hyp_sig))
        if self.dof_no_UL is not None:
            out.append("")
            out.append("stats without the UL")
            out.append("dof  UL=%d" % self.dof_no_UL)
            out.append(
                "chisq=%f, chisq/red=%f null hypothesis sig=%f" % (self.chisq_no_UL, self.chisq_red_no_UL, self.null_hyp_sig_no_UL))
            out.append("")
        out.append("")
        out.append("best fit pars")
        for _l in self.bestfit_table.pformat_all():
            out.append(_l)
        #pars_rep=self.parameters.show_best_fit_pars(getstring=True)
        #for string in pars_rep:
        #    out.append(string)
            

        out.append("**************************************************************************************************")    
        out.append("")
        return out
        
    def _update_asymm_errors(self):
        for pi in range(len(self.mm.fit_par_free)):
            if  hasattr(self.mm.minimizer,'asymm_errors'):
                self.mm.fit_par_free[pi].err_p=self.mm.minimizer.asymm_errors[pi][0]
                self.mm.fit_par_free[pi].err_m=self.mm.minimizer.asymm_errors[pi][1]

    def show_report(self):
        self.fit_report=self.get_report()
        for text in self.fit_report:
            try:
                print (text)
            except Exception as e:
                raise(RuntimeWarning,'problem in formatting text for report',e)

    
    def save_report(self,wd=None,name=None):
        if wd is None:
            wd=self.wd
        
        if name is None:
            name='best_fit_report_%s'%self.name+'.txt'
        
            
        outname='%s/%s'%(wd,name)
         
        outfile=open(outname,'w')
    
        
        for text in self.fit_report:
        
            print(text,file=outfile)
            
        outfile.close()










class ModelMinimizer(object):


    def __init__(self,minimizer_type):
        __accepted__ = ['lsb', 'minuit', 'ls']


        if minimizer_type == 'lsb':
            self.minimizer=LSBMinimizer(self)

        elif minimizer_type=='ls':
            self.minimizer=LSMinimizer(self)

        elif minimizer_type=='minuit':
            self.minimizer=MinutiMinimizer(self)

        elif minimizer_type not in __accepted__:
            raise RuntimeError('minimizer ', minimizer_type, 'not accepted, please choose among', __accepted__)

        else:
            raise RuntimeError('minimizer factory failed')

        #print('minimizer',minimizer_type)


    def _prepare_fit(self,fit_Model,sed_data,nu_fit_start,nu_fit_stop,fitname=None,fit_workplace=None,loglog=False,silent=False,get_conf_int=False,use_facke_err=False,use_UL=False):

        if fitname is None:
            fitname = fit_Model.name + '_' + WorkPlace.flag


        if fit_workplace is None:
            fit_workplace = WorkPlace()
            out_dir = fit_workplace.out_dir + '/' + fitname + '/'
        else:
            out_dir = fit_workplace.out_dir + '/' + fitname + '/'

        makedir(out_dir)

        for model in fit_Model.components:
            if model.model_type == 'jet':
                model.set_path(out_dir)

        if sed_data.data['dnuFnu_data'] is None:
            sed_data.data['dnuFnu_data'] = np.ones(sed_data.data['nu_data'].size)

        # filter data points
        msk1 = sed_data.data['nu_data'] > nu_fit_start
        msk2 = sed_data.data['nu_data'] < nu_fit_stop
        msk_zero_error = sed_data.data['dnuFnu_data'] > 0.0
        # msk = s.array([(el>nu_fit_start) and (el<nu_fit_stop) for el in SEDdata.data['nu_data']])
        # print msk1.size,msk2.size,SEDdata.data['UL'].size
        if use_UL == True:
            msk = msk1 * msk2 * msk_zero_error

        else:
            msk = msk1 * msk2 * msk_zero_error*~sed_data.data['UL']

        UL = sed_data.data['UL'][msk]

        if loglog == False:
            nu_fit = sed_data.data['nu_data'][msk]
            nuFnu_fit = sed_data.data['nuFnu_data'][msk]
            if use_facke_err == False:
                err_nuFnu_fit = sed_data.data['dnuFnu_data'][msk]
            else:
                err_nuFnu_fit = sed_data.data['dnuFnu_facke'][msk]
        else:
            nu_fit = sed_data.data['nu_data_log'][msk]
            nuFnu_fit = sed_data.data['nuFnu_data_log'][msk]
            err_nuFnu_fit = sed_data.data['dnuFnu_data_log'][msk]
            if use_facke_err == False:
                err_nuFnu_fit = sed_data.data['dnuFnu_data_log'][msk]
            else:
                err_nuFnu_fit = sed_data.data['dnuFnu_facke_log'][msk]


        if silent == False:
            print("filtering data in fit range = [%e,%e]" % (nu_fit_start, nu_fit_stop))
            print("data length", nu_fit.size)


        # print nu_fit,len(nu_fit)

        # set starting value of parameters
        for par in fit_Model.parameters.par_array:
            if get_conf_int == True:
                par.set_fit_initial_value(par.best_fit_val)
            else:
                par.set_fit_initial_value(par.val)

        fit_par_free = [par for par in fit_Model.parameters.par_array if par.frozen == False]

        pinit = [par.get_fit_initial_value() for par in fit_par_free]

        # bounds
        free_pars = 0

        for pi in range(len(fit_Model.parameters.par_array)):

            if fit_Model.parameters.par_array[pi].frozen == False:
                free_pars += 1

        if silent == False:
            print(section_separator)
            print("*** start fit process ***")
            print("initial pars: ")

            fit_Model.parameters.show_pars()
            print("----- ")

        self.out_dir=out_dir
        self.pinit=pinit
        self.free_pars=free_pars
        self.fit_par_free=fit_par_free
        self.nu_fit=nu_fit
        self.nuFnu_fit=nuFnu_fit
        self.err_nuFnu_fit=err_nuFnu_fit
        self.fit_Model=fit_Model
        self.loglog=loglog
        self.UL=UL
        self.fit_Model.nu_min_fit = nu_fit_start
        self.fit_Model.nu_max_fit = nu_fit_stop
        #print('-->nu_fit_start A %e' % nu_fit_start)
        #print('-->nu_fit_start B %e' % self.fit_Model.nu_min_fit)

    def fit(self,fit_Model,
            sed_data,
            nu_fit_start,
            nu_fit_stop,
            fitname=None,
            fit_workplace=None,
            loglog=False,
            silent=False,
            get_conf_int=False,
            max_ev=0,
            use_facke_err=False,
            use_UL=False,
            skip_minimizer=False):

        self.silent=silent
        #print('-->nu_fit_start', nu_fit_start)

        self._prepare_fit( fit_Model, sed_data, nu_fit_start, nu_fit_stop, fitname=fitname, fit_workplace=fit_workplace,
                     loglog=loglog, silent=silent, get_conf_int=get_conf_int, use_facke_err=use_facke_err,use_UL=use_UL)


        if skip_minimizer == False:
            self.minimizer.fit(self,max_ev=max_ev,silent=silent)
        else:
            pass

        #print('-->nu_fit_start', nu_fit_start)

        return self.get_fit_results(fit_Model,nu_fit_start,nu_fit_stop,fitname,loglog=loglog,silent=silent)

    def get_fit_results(self, fit_Model, nu_fit_start, nu_fit_stop, fitname, silent=False, loglog=False):
        self.reset_to_best_fit()
        best_fit = FitResults(fitname,
                              self,
                              fit_Model.parameters,
                              self.minimizer.calls,
                              self.minimizer.mesg,
                              self.minimizer.success,
                              self.minimizer.chisq,
                              self.minimizer.dof,
                              self.minimizer.chisq_red,
                              self.minimizer.null_hyp_sig,
                              self.out_dir,
                              chisq_no_UL=self.minimizer.chisq_no_UL,
                              dof_no_UL=self.minimizer.dof_no_UL,
                              chisq_red_no_UL=self.minimizer.chisq_red_no_UL,
                              null_hyp_sig_no_UL=self.minimizer.null_hyp_sig_no_UL)

        if silent == False:
            best_fit.show_report()

        #_nu1=fit_Model.nu_min
        #_nu2=fit_Model.nu_max
        #print('-->nu_fit_start', nu_fit_start)

        fit_Model.set_nu_grid(nu_min=nu_fit_start, nu_max=nu_fit_stop)
        fit_Model.eval(fill_SED=True, loglog=loglog, phys_output=True)

        res_bestfit = self.minimizer.residuals_Fit(self.minimizer.pout,
                                                 self.fit_par_free,
                                                 self.nu_fit,
                                                 self.nuFnu_fit,
                                                 self.err_nuFnu_fit,
                                                 self.fit_Model,
                                                 self.loglog,
                                                 self.UL)

        #msk1 = self.nu_fit > nu_fit_start
        #msk2 = self.nu_fit < nu_fit_stop

        #res_bestfit=res_bestfit[msk1*msk2]

        if loglog == True:
            fit_Model.SED.fill(nu_residuals=np.power(10, self.nu_fit), residuals=res_bestfit)
        else:
            fit_Model.SED.fill(nu_residuals=self.nu_fit, residuals=res_bestfit)


        if silent == False:
            print(section_separator)

        #fit_Model.set_nu_grid(nu_min=_nu1, nu_max=_nu2)
        fit_Model.eval(fill_SED=True)
        #print('-->nu_fit_start bf A%e'%nu_fit_start)
        #print('-->nu_fit_start bf B%e'%self.fit_Model.nu_min_fit)
        return best_fit

    def reset_to_best_fit(self):
        for pi in range(len(self.fit_par_free)):
            self.fit_par_free[pi].set(val=self.minimizer.pout[pi])
            self.fit_par_free[pi].best_fit_val = self.minimizer.pout[pi]
            self.fit_par_free[pi].best_fit_err = self.minimizer.errors[pi]
            if  hasattr(self.minimizer,'asymm_errors'):
                self.fit_par_free[pi].err_p=self.minimizer.asymm_errors[pi][0]
                self.fit_par_free[pi].err_m=self.minimizer.asymm_errors[pi][1]

        self.fit_Model.eval()


class Minimizer(object):

    def __init__(self,model):
        self.model=model
        self._progress_iter = cycle(['|', '/', '-', '\\'])


    def fit(self,model,max_ev=None,use_UL=False,silent=False):
        self.use_UL = use_UL
        self.calls=0
        self.res_check=None
        self.molde=model
        self.silent=silent
        self._fit(max_ev)
        self._fit_stats()
        self._set_fit_errors()



    def _fit_stats(self):
        self.dof = len(self.model.nu_fit) - self.model.free_pars
        self.chisq = self.get_chisq()
        self.chisq_red = self.chisq / float(self.dof)
        self.null_hyp_sig = 1.0 - chi2.cdf(self.chisq, self.dof)
        if self.use_UL==True:
            self.dof_no_UL = len(self.model.nu_fit) - self.model.free_pars - self.model.UL.sum()
            self.use_UL=False
            self.chisq_no_UL = self.get_chisq()
            self.use_UL = True
            self.chisq_red_no_UL = self.chisq_no_UL / float(self.dof_no_UL)
            self.null_hyp_sig_no_UL=1.0 - chi2.cdf(self.chisq_no_UL, self.dof_no_UL)
        else:
            self.dof_no_UL = None
            self.chisq_red_no_UL = None
            self.chisq_no_UL=None
            self.null_hyp_sig_no_UL=None
        self.success = True
        self.status = 1


    def _set_fit_errors(self):
        if self.covar is None:
            print("!Warning, no covariance matrix produced")
            self.errors = np.zeros(len(self.pout))
        else:
            self.errors = [np.sqrt(np.fabs(self.covar[pi, pi]) * self.chisq_red) for pi in range(len(self.model.fit_par_free))]


    def _progess_bar(self, _res_sum, res_UL):
        if np.mod(self.calls, 10) == 0 and self.calls != 0:
            _c= ' ' * 256
            print("\r%s"%_c,end="")
            print("\r%s minim function calls=%d, chisq=%f UL part=%f" % (next(self._progress_iter),self.calls, _res_sum, -2.0*np.sum(res_UL)), end="")


    def residuals_Fit(self,p,
                      fit_par,
                      nu_data,
                      nuFnu_data,
                      err_nuFnu_data,
                      best_fit_SEDModel,
                      loglog,
                      UL,
                      chisq=False,
                      use_UL=False,
                      silent=False):

        #hide_cursor()
        _warn=False
        for pi in range(len(fit_par)):
            fit_par[pi].set(val=p[pi])
            if np.isnan(p[pi]):
                _warn=True
                print('warning nan for par',pi,' old paramter value was',self._par_check[pi])

        if _warn==True:
            best_fit_SEDModel.show_pars()
            print('res_sum',self._res_sum_chekc)
            print('res_chekc', self._res_chekc)
            print('res_UL_chekc', self._res_UL_chekc)


        model = best_fit_SEDModel.eval(nu=nu_data, fill_SED=False, get_model=True, loglog=loglog)


        _res_sum, _res, _res_UL=log_like(nuFnu_data,model,err_nuFnu_data,UL,use_UL=use_UL)

        self._res_sum_chekc=_res_sum
        self._res_chekc = _res
        self._res_UL_chekc = _res_UL
        self._par_check=p

        self.calls +=1

        if silent==False:
            self._progess_bar(_res_sum, _res_UL)
            print("\r", end="")

        if chisq==True:
            res=_res_sum
        else:
            res=_res

        #print('res', res, chisq,use_UL)
        #show_cursor()

        return res


def log_like(data,model,data_error,UL,use_UL=False):
    res_no_UL = (data[~UL] - model[~UL]) / (data_error[~UL])
    res = (data - model) / (data_error)
    res_UL = [0]
    if UL.sum() > 0 and use_UL == True:
        res_UL = UL_log_like(data[UL], model[UL], data_error[UL])

    res_sum=np.sum(res_no_UL * res_no_UL) - 2.0*np.sum(res_UL)
    return res_sum,res,res_UL

def UL_log_like(y_UL,y_model,y_err):
    y = (y_UL - y_model) / (np.sqrt(2) *y_err)
    x=0.5*(1.0+sp.special.erf(y))
    x[x==0]=1E-200
    return  np.log(x)



class LSBMinimizer(Minimizer):

    def __init__(self, model):
        super(LSBMinimizer, self).__init__(model)

    def _fit(self, max_ev,):
        bounds = [(par.fit_range_min, par.fit_range_max) for par in self.model.fit_par_free]
        max_nfev = 0 if (max_ev == 0 or max_ev == None) else max_ev
        pout, covar, info, mesg, success = leastsqbound(self.residuals_Fit,
                                                        self.model.pinit,
                                                        args=(self.model.fit_par_free,
                                                              self.model.nu_fit,
                                                              self.model.nuFnu_fit,
                                                              self.model.err_nuFnu_fit,
                                                              self.model.fit_Model,
                                                              self.model.loglog,
                                                              self.model.UL,
                                                              False,
                                                              False,
                                                              self.silent),
                                                        xtol=5.0E-8,
                                                        ftol=5.0E-8,
                                                        full_output=1,
                                                        bounds=bounds,
                                                        maxfev=max_nfev)

        self.mesg = mesg
        self.covar = covar
        self.chisq =  sum(info["fvec"] * info["fvec"])

        self.pout = pout

    def get_chisq(self):
        return self.chisq


class LSMinimizer(Minimizer):

    def __init__(self,model ):
        super(LSMinimizer, self).__init__(model)


    def _fit(self,max_ev):
        bounds = ([-np.inf if par.fit_range_min is None else par.fit_range_min for par in self.model.fit_par_free],
                  [np.inf if par.fit_range_max is None else par.fit_range_max for par in self.model.fit_par_free])

        max_nfev = None if (max_ev == 0 or max_ev is None) else max_ev
        fit = least_squares(self.residuals_Fit,
                            self.model.pinit,
                            args=(self.model.fit_par_free,
                                  self.model.nu_fit,
                                  self.model.nuFnu_fit,
                                  self.model.err_nuFnu_fit,
                                  self.model.fit_Model,
                                  self.model.loglog,
                                  self.model.UL,
                                  False,
                                  False,
                                  self.silent),

                            xtol=1.0E-8,
                            ftol=1.0E-8,
                            jac='3-point',
                            loss='linear',
                            f_scale=0.01,
                            bounds=bounds,
                            max_nfev=max_nfev,
                            verbose=0)
        #print(vars(fit))
        self.covar=np.linalg.inv(2 * np.dot(fit.jac.T, fit.jac))
        self.chisq=sum(fit.fun * fit.fun)
        self.mesg = fit.message
        self.pout = fit.x

    def get_chisq(self):
        return self.chisq




class MinutiMinimizer(Minimizer):

    def __init__(self,model):
        if minuit_installed==True:
            pass
        else:
            raise RuntimeError('iminuit non istalled')

        super(MinutiMinimizer, self).__init__(model)



    def _fit(self,max_ev=None):
        bounds = [(par.fit_range_min, par.fit_range_max) for par in self.model.fit_par_free]
        self._set_minuit_func(self.model.pinit, bounds)
        #print('=>,')
        if max_ev is None or max_ev==0:
            max_ev =10000

        self.mesg=self.minuit_fun.migrad(ncall=max_ev)
        self.pout=[self.minuit_fun.values[k] for k in self.minuit_fun.values.keys()]

    def _set_fit_errors(self):
        self.errors = [self.minuit_fun.errors[k] for k in self.minuit_fun.errors.keys()]





    def _set_minuit_func(self, p_init, bounds,p_error=None):
        #print('=>Hi')
        if p_error==None:
            p_error=[0.1]*len(p_init)

        p_names = ['par_{}'.format(_) for _ in range(len(p_init))]
        error_names = ['error_par_{}'.format(_) for _ in range(len(p_init))]
        p_bound_names = ['limit_par_{}'.format(_) for _ in range(len(p_init))]

        #This dict contains all the kw for
        #par
        kwdarg = {}
        for n, p, bn, b,en,e in zip(p_names, p_init, p_bound_names, bounds,error_names,p_error):
            kwdarg[n] = p
            kwdarg[bn] = b
            kwdarg[en] = e
            #print( kwdarg[n] ,p, kwdarg[bn] ,b, kwdarg[en] ,e)

        self.minuit_par_name_dict={}
        self.minuit_bounds_dict={}
        self.par_dict = {}
        for ID,par in enumerate(self.model.fit_par_free):
            self.minuit_par_name_dict[par.name]=p_names[ID]
            self.minuit_bounds_dict[par.name]=bounds[ID]
            self.par_dict[par.name]=par

        self.minuit_fun = iminuit.Minuit(
            fcn=self.chisq_func,
            forced_parameters=p_names,
            pedantic=False,
            errordef=1,
            **kwdarg)

    def chisq_func(self, *p):
        self.p = p
        return self.residuals_Fit(p,
                          self.model.fit_par_free,
                          self.model.nu_fit,
                          self.model.nuFnu_fit,
                          self.model.err_nuFnu_fit,
                          self.model.fit_Model,
                          self.model.loglog,
                          self.model.UL,
                          chisq=True,
                          use_UL=self.use_UL,
                          silent=self.silent)



    def get_chisq(self):
        #print ('p',self.p)
        return self.chisq_func(*self.p)

    def minos_errors(self,par=None):
        if par is not None:
            par=self.minuit_par_name_dict[par]
        self.calls = 0
        self.minuit_fun.minos(var=par)
        self.minuit_fun.print_param()
        self.asymm_errors = [(self.minuit_fun.merrors[(k,1.0)],self.minuit_fun.merrors[(k,-1.0)]) for k in self.minuit_fun.errors.keys()]
        self.model.reset_to_best_fit()


    def profile(self,par,bound=2,subtract_min=True):
        self.calls = 0
        bound=self._set_bounds(par,bound=bound)

        x, y =  self.minuit_fun.profile(self.minuit_par_name_dict[par],
                                        bound=bound,
                                        subtract_min=subtract_min)

        self.model.reset_to_best_fit()
        return x,y


    def mnprofile(self,par,bound=2,subtract_min=True):
        """

        Parameters
        ----------
        par
        bound
        subtract_min

        Returns
        -------

        """
        self.calls = 0
        bound = self._set_bounds(par, bound=bound)
        x, y,r =  self.minuit_fun.mnprofile(self.minuit_par_name_dict[par],
                                          bound=bound,
                                          subtract_min=subtract_min)
        #print('-->r',r)
        self.model.reset_to_best_fit()
        return x,y,r



    def draw_mnprofile(self,par,bound=2):
        """

        Parameters
        ----------
        par
        bound

        Returns
        -------

        """
        self.calls = 0
        bound = self._set_bounds(par, bound=bound)
        x,y,r=self.mnprofile(par,bound,subtract_min=True)
        fig,ax = self._draw_profile(x, y, par)
        self.model.reset_to_best_fit()
        return x, y, fig,ax

    def draw_profile(self,par,bound=2):
        self.calls = 0
        bound = self._set_bounds(par, bound=bound)
        x, y = self.profile(par, subtract_min=True, bound=bound)
        fig,ax=self._draw_profile(x,y,par)
        self.model.reset_to_best_fit()
        return x,y,fig,ax

    def _draw_profile(self,x,y,par):
        self.calls = 0
        x = np.array(x)
        y = np.array(y)

        fig, ax = plt.subplots()

        ax.plot(x, y)
        ax.grid(True)
        ax.set_xlabel(par)
        ax.set_ylabel('FCN')


        minpos = np.argmin(y)
        ymin = y[minpos]
        tmpy = y - ymin
        # now scan for minpos to the right until greater than one
        up = self.minuit_fun.errordef
        righty = np.power(tmpy[minpos:] - up, 2)
        right_min = np.argmin(righty)
        rightpos = right_min + minpos
        lefty = np.power((tmpy[:minpos] - up), 2)
        left_min = np.argmin(lefty)
        leftpos = left_min
        le = x[minpos] - x[leftpos]
        re = x[rightpos] - x[minpos]

        ax.axvspan(x[leftpos], x[rightpos], facecolor='g', alpha=0.5)

        plt.figtext(0.5, 0.5,
                        '%s = %7.3e ( -%7.3e , +%7.3e)' % (par, x[minpos], le, re),
                        ha='center')



        return fig,ax





    def contour(self,par_1,par_2,bound=2,bins=20,subtract_min=True):
        self.calls = 0
        if np.shape(bound)==():

            bound_1 = self._set_bounds(par_1,bound)
            bound_2 = self._set_bounds(par_2,bound)

        elif np.shape(bound)==(2,) or np.shape(bound)==(2,2):

            bound_1 = self._set_bounds(par_1, bound[0])
            bound_2 = self._set_bounds(par_2, bound[1])
        else:
            raise  RuntimeError('bound must be scalar or (2,) or (2,2)')

        bound = [bound_1, bound_2]
        x, y, z= self.minuit_fun.contour(self.minuit_par_name_dict[par_1],
                                         self.minuit_par_name_dict[par_2],
                                         subtract_min=subtract_min,
                                         bound=bound,
                                         bins=bins)
        self.model.reset_to_best_fit()
        return x,y,z



    def draw_contour(self,par_1,par_2,bound=2,levels=np.arange(5)):
        self.calls = 0
        l=self.minuit_fun.errordef*levels
        x,y,z=self.contour(par_1,par_2,bound=bound)
        fig, ax = plt.subplots()
        CS = ax.contour(x, y, z, levels,)
        ax.clabel(CS, fontsize=9, inline=1)
        ax.set_xlabel(par_1)
        ax.set_ylabel(par_2)
        self.model.reset_to_best_fit()
        return x,y,z,fig,ax

    def _set_bounds(self,par,bound):
        p = self.par_dict[par]
        if np.shape(bound)==():
            bound = [p.best_fit_val - p.best_fit_err * bound, p.best_fit_val + p.best_fit_err * bound]

        elif np.shape(bound)==(2,):
            pass
        else:
            raise  RuntimeError('bound shape',np.shape(bound),'for par',p.name,'is wrong, has to be scalare or (2,)')
        #print('bound for par', par, ' set to', bound)
        bound = self._check_bounds(par, bound)
        #print('bound for par',par,' updated to',bound)
        return bound

    def _check_bounds(self,par,bound):

        if  self.minuit_bounds_dict[par][0] is not None:
            if bound[0]<self.minuit_bounds_dict[par][0]:
                bound[0]=self.minuit_bounds_dict[par][0]

        if self.minuit_bounds_dict[par][1] is not None:
            if bound[1]>self.minuit_bounds_dict[par][1]:
                bound[1]=self.minuit_bounds_dict[par][1]

        return bound

def fit_SED(fit_Model, sed_data, nu_fit_start, nu_fit_stop, fitname=None, fit_workplace=None, loglog=False, silent=False,
            get_conf_int=False, max_ev=0, use_facke_err=False, minimizer='lsb', use_UL=False):
    mm = ModelMinimizer(minimizer)
    #print('-->nu_fit_start',nu_fit_start)
    return mm,mm.fit(fit_Model,
                  sed_data,
                  nu_fit_start,
                  nu_fit_stop,
                  fitname=fitname,
                  fit_workplace=fit_workplace,
                  loglog=loglog,
                  silent=silent,
                  get_conf_int=get_conf_int,
                  max_ev=max_ev,
                  use_facke_err=use_facke_err,
                  use_UL=use_UL)
