
__author__ = "Andrea Tramacere"

from itertools import cycle


import scipy as sp

import numpy as np

import warnings

from scipy.stats import chi2

import dill as pickle

try:
    import iminuit
    from iminuit import Minuit
    minuit_installed=True
except:
    minuit_installed=False


from tqdm.auto import tqdm


from .plot_sedfit import plt

from scipy.optimize import least_squares

from .leastsqbound.leastsqbound import leastsqbound

from .output import section_separator,WorkPlace,makedir
from .utils import JetkerneltException

from .data_loader import ObsData,Data
from .data_loader import lin_to_log





__all__ = ['FitResults','fit_SED','Minimizer','LSMinimizer','LSBMinimizer', 'MinuitMinimizer', 'ModelMinimizer']



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

        self.update_report()
        self.wd=wd
         

    def update_report(self):
        # self.model_table=self.parameters._par_table
        self.parameters._build_best_fit_par_table()
        self.parameters._build_par_table()

    def _show_report(self,):
        self.update_report()

        print('-------------------------------------------------------------------------')
        print("Fit report")
        print("")
        print("Model: %s"%self.name)

        self.parameters.show_pars()
        print("")
        print("converged=%s"%self.success)
        print("calls=%d"%self.calls)
        print("mesg=")
        if hasattr(self,'mesg'):
            try:
                from IPython.display import display
                display(self.mesg)
            except:
                print(self.mesg)

        print("dof=%d"%self.dof)
        print("chisq=%f, chisq/red=%f null hypothesis sig=%f"%(self.chisq,self.chisq_red,self.null_hyp_sig))
        if self.dof_no_UL is not None:
            print("")
            print("stats without the UL")
            print("dof  UL=%d" % self.dof_no_UL)
            print(
                "chisq=%f, chisq/red=%f null hypothesis sig=%f" % (self.chisq_no_UL, self.chisq_red_no_UL, self.null_hyp_sig_no_UL))
            print("")
        print("")
        print("best fit pars")
        self.parameters.show_best_fit_pars()

        print('-------------------------------------------------------------------------')
        print("")

    @property
    def bestfit_table(self):
        return self.parameters.best_fit_par_table

    def _update_asymm_errors(self):
        if hasattr(self.mm):
            for pi in range(len(self.mm.fit_par_free)):
                if  hasattr(self.mm.minimizer,'asymm_errors'):
                    self.mm.fit_par_free[pi].err_p=self.mm.minimizer.asymm_errors[pi][0]
                    self.mm.fit_par_free[pi].err_m=self.mm.minimizer.asymm_errors[pi][1]

    def show_report(self):
        try:
            self._show_report()
        except:
            print('problem in formatting text for report')
        #except Exception as e:
        #    raise(RuntimeWarning,'problem in formatting text for report',e)

    @classmethod
    def load_report(cls,file_name):
        c = pickle.load(open(file_name, "rb"))
        return c

    def save_report(self,name=None):
        if name is None:
            name = 'best_fit_report.pkl'
        mm =self.mm
        self.mm=None
        setattr(self, 'mm', mm)
        try:
            pickle.dump(self, open(name, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
        except:
            o = FitResults(name=self.name,
                           mm=None,
                           parameters=self.parameters,
                           calls=self.calls,
                           mesg=None,
                           success=self.success,
                           chisq=self.chisq,
                           dof=self.dof,
                           chisq_red=self.chisq_red,
                           null_hyp_sig=self.null_hyp_sig,
                           wd=self.wd,
                           chisq_no_UL=self.chisq_no_UL,
                           dof_no_UL=self.dof_no_UL,
                           chisq_red_no_UL=self.chisq_red_no_UL,
                           null_hyp_sig_no_UL=self.null_hyp_sig_no_UL)

            pickle.dump(o, open(name, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)












class ModelMinimizer(object):


    def __init__(self,minimizer_type):
        __accepted__ = ['lsb', 'minuit', 'ls','sherpa']


        if minimizer_type == 'lsb':
            self.minimizer=LSBMinimizer(self)

        elif minimizer_type=='ls':
            self.minimizer=LSMinimizer(self)

        elif minimizer_type=='minuit':
            self.minimizer=MinuitMinimizer(self)

        elif minimizer_type == 'sherpa':
            from .sherpa_plugin import  SherpaMinimizer
            self.minimizer = SherpaMinimizer(self)

        elif minimizer_type not in __accepted__:
            raise RuntimeError('minimizer ', minimizer_type, 'not accepted, please choose among', __accepted__)

        else:
            raise RuntimeError('minimizer factory failed')

        self.minimizer_type=minimizer_type
        #print('minimizer',minimizer_type)

    def save_model(self, file_name):
        _m= self.minimizer
        self.minimizer=self.minimizer_type

        pickle.dump(self, open(file_name, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
        self.minimizer =_m

    @classmethod
    def load_model(cls, file_name):

        try:
            c = pickle.load(open(file_name, "rb"))
            if isinstance(c, ModelMinimizer):
                c.__init__(c.minimizer)
                return  c
            else:
                raise RuntimeError('The model you loaded is not valid please check the file name')

        except Exception as e:
            raise RuntimeError(e)

    def _handle_sed_data(self,sed_data,loglog=False,use_fake_err=False):
        if sed_data.data['dnuFnu_data'] is None:
            sed_data.data['dnuFnu_data'] = np.ones(sed_data.data['nu_data'].size)

        if loglog is False:
            nu_fit = sed_data.data['nu_data']
            nuFnu_fit = sed_data.data['nuFnu_data']
            if use_fake_err == False:
                err_nuFnu_fit = sed_data.data['dnuFnu_data']
            else:
                err_nuFnu_fit = sed_data.data['dnuFnu_fake']
        else:
            nu_fit = sed_data.data['nu_data_log']
            nuFnu_fit = sed_data.data['nuFnu_data_log']
            err_nuFnu_fit = sed_data.data['dnuFnu_data_log']
            if use_fake_err == False:
                err_nuFnu_fit = sed_data.data['dnuFnu_data_log']
            else:
                err_nuFnu_fit = sed_data.data['dnuFnu_fake_log']

        UL = sed_data.data['UL']

        data = np.zeros(nu_fit.size,
                             dtype=[('x', nu_fit.dtype), ('y', nuFnu_fit.dtype), ('dy', err_nuFnu_fit.dtype),
                                    ('UL', UL.dtype)])
        data['x'] = nu_fit
        data['y'] = nuFnu_fit
        data['dy'] = err_nuFnu_fit
        data['UL'] = UL
        return data

    def _handle_data_xy(self, data, loglog=False, use_fake_err=False):
        #if data['dy'] is None:
        #    data['dy'] = np.ones(data['x'].size)

        if loglog is False:
            x = data.table['x']
            y = data.table['y']
            if use_fake_err == False:
                dy = data.table['dy']
            else:
                dy = data.table['dy_fake']
        else:
            x = np.log10(data.table['x'])
            y, dy = lin_to_log(val=data.table['y'], err=data.table['dy'])
            if use_fake_err == False:
                pass
            else:
                y, dy = lin_to_log(val=data.table['y'], err=data.table['dy_fake'])

        UL = data.table['UL']

        _data = np.zeros(x.size,
                        dtype=[('x', x.dtype), ('y', y.dtype), ('dy', dy.dtype),
                               ('UL', UL.dtype)])
        _data['x'] = x
        _data['y'] = y
        _data['dy'] = dy
        _data['UL'] = UL
        return _data


    def _prepare_fit(self,
                     fit_model,
                     data,
                     nu_fit_start,
                     nu_fit_stop,
                     fitname=None,
                     fit_workplace=None,
                     loglog=False,
                     silent=False,
                     get_conf_int=False,
                     use_fake_err=False,
                     use_UL=False):
        #print('--> DEBUG WorkPlace.flag',WorkPlace.flag,fit_workplace)

        if fitname is None:
            fitname = fit_model.name
            #print('--> DEBUG ')

        if fit_workplace is None:
            fit_workplace = WorkPlace()
            out_dir = fit_workplace.out_dir + '/' + fitname + '/'
        else:
            out_dir = fit_workplace.out_dir + '/' + fitname + '/'

        makedir(out_dir)

        for model in fit_model.components.components_list:
            if model.model_type == 'jet':
                model.set_path(out_dir)

        if isinstance(data,ObsData):
            self.data=self._handle_sed_data(data,use_fake_err=use_fake_err,loglog=loglog)
        elif isinstance(data,Data):
            self.data= self._handle_data_xy(data,use_fake_err=use_fake_err,loglog=loglog)
        else:
            raise  RuntimeError('data type must be Obs  Data or Data')

        #if sed_data.data['dnuFnu_data'] is None:
        #    sed_data.data['dnuFnu_data'] = np.ones(sed_data.data['nu_data'].size)

        # filter data points
        if loglog is True:
            nu_fit_start = np.log10(nu_fit_start)
            nu_fit_stop = np.log10(nu_fit_stop)

        msk1 = self.data['x'] > nu_fit_start
        msk2 = self.data['x'] < nu_fit_stop
        msk_zero_error = self.data['dy'] > 0.0
        self.use_UL=use_UL
        if use_UL == True:
            msk = msk1 * msk2 * msk_zero_error

        else:
            msk = msk1 * msk2 * msk_zero_error*~self.data['UL']

        self.data = self.data[msk]

        if silent == False:
            print("filtering data in fit range = [%e,%e]" % (nu_fit_start, nu_fit_stop))
            print("data length", self.data.size)



        for par in fit_model.parameters.par_array:
            if get_conf_int == True:
                par.set_fit_initial_value(par.best_fit_val)
            else:
                par.set_fit_initial_value(par.val)

        fit_par_free = [par for par in fit_model.parameters.par_array if par.frozen is False and (par.linked is False or par.root is True) and par._is_dependent is False]

        #remove duplicates and keeps the order
        fit_par_free = sorted(set(fit_par_free), key=fit_par_free.index)

        pinit = [par.get_fit_initial_value() for par in fit_par_free]

        # bounds
        free_pars = 0

        for pi in range(len(fit_model.parameters.par_array)):

            if fit_model.parameters.par_array[pi].frozen == False:
                free_pars += 1

        if silent == False:
            print(section_separator)
            print("*** start fit process ***")
            #print("initial pars: ")
            #fit_model.parameters.show_pars()
            print("----- ")

        self.out_dir=out_dir
        self.pinit=pinit
        self.pout=None
        self.free_pars=free_pars
        self.fit_par_free=fit_par_free
        self.fit_model=fit_model
        self.loglog=loglog

        if hasattr(self.fit_model,'nu_min_fit'):
            self.fit_model.nu_min_fit = nu_fit_start
            self.fit_model.nu_max_fit = nu_fit_stop



        if len(fit_par_free)>self.data['x'].size and isinstance(self.minimizer,LSBMinimizer):
            m='number of data points: %d is lower than number of free pars: %d'%(self.data['x'], len(fit_par_free))
            raise  JetkerneltException(message=m)

    def fit(self,
            fit_model,
            sed_data,
            nu_fit_start,
            nu_fit_stop,
            fitname=None,
            fit_workplace=None,
            loglog=False,
            silent=False,
            get_conf_int=False,
            max_ev=0,
            use_fake_err=False,
            use_UL=False,
            skip_minimizer=False,
            repeat=1):

        self.silent=silent
        #print('-->nu_fit_start', nu_fit_start)

        self._prepare_fit( fit_model, sed_data, nu_fit_start, nu_fit_stop, fitname=fitname, fit_workplace=fit_workplace,
                     loglog=loglog, silent=silent, get_conf_int=get_conf_int, use_fake_err=use_fake_err,use_UL=use_UL)

        fit_model.set_nu_grid(nu_min=nu_fit_start*0.5, nu_max=nu_fit_stop*1.5)


        for i in range(repeat):
            if skip_minimizer == False:
                if repeat>1:
                    if silent is False:
                        print('fit run:',i)
                    if i>0:
                        self.pinit = [v for v in self.minimizer.pout]
                        if silent is False:
                            print('- old chisq=%5.5e' % (self.minimizer.chisq))
                        old_chisq = self.minimizer.chisq

                self.minimizer.fit(self,max_ev=max_ev,silent=silent,use_UL=use_UL)

                self.pout = self.minimizer.pout
                self.errors = self.minimizer.errors
                if hasattr(self.minimizer, 'asymm_errors'):
                    self.asymm_errors = self.minimizer.asymm_errors


                self.reset_to_best_fit()
                self.minimizer._fit_stats()
                if silent is False:
                    print('- best chisq=%5.5e'%(self.minimizer.chisq))
                    print()

                new_chisq = self.minimizer.chisq
                if i==0:
                    old_chisq = new_chisq

                chisq_min = min(new_chisq, old_chisq)
                #if i==0 or new_chisq<=chisq_min:
                #    self.pout = self.minimizer.pout
                #    self.errors = self.minimizer.errors
                #    if  hasattr(self.minimizer,'asymm_errors'):
                #        self.asymm_errors = self.minimizer.asymm_errors
                #elif new_chisq>=chisq_min:
                #    self.pout = self.minimizer.pout
                #    self.errors = self.minimizer.errors
                #    if hasattr(self.minimizer, 'asymm_errors'):
                #        self.asymm_errors = self.minimizer.asymm_errors

                #if i<repeat-1:
                #    self.pinit = [par.val for par in self.fit_par_free]
                #    if silent is False:
                #        print()
            else:
                pass


        return self.get_fit_results(fit_model,nu_fit_start,nu_fit_stop,fitname,loglog=loglog,silent=silent)

    def get_fit_results(self, fit_model, nu_fit_start, nu_fit_stop, fitname, silent=False, loglog=False):

        self.reset_to_best_fit()
        self.minimizer._fit_stats()
        best_fit = FitResults(fitname,
                              self,
                              fit_model.parameters,
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

        fit_model.set_nu_grid(nu_min=nu_fit_start , nu_max=nu_fit_stop)
        fit_model.eval(fill_SED=True, loglog=loglog, phys_output=True)

        res_bestfit = self.minimizer.residuals_Fit(self.minimizer.pout,
                                                 self.fit_par_free,
                                                 self.data,
                                                 self.fit_model,
                                                 self.loglog,
                                                 use_UL=self.use_UL)



        if loglog == True:
            fit_model.SED.fill(nu_residuals=np.power(10, self.data['x']), residuals=res_bestfit)
        else:
            fit_model.SED.fill(nu_residuals=self.data['x'], residuals=res_bestfit)


        if silent == False:
            print(section_separator)

        fit_model.eval(fill_SED=True)
        if self.minimizer._post_fit_warnings!='':
                print('there are  fit warnings messages, use  the .show_fit_warnings() or access the member .minimizer._post_fit_warnings')

        return best_fit

    def show_fit_warnings(self):
        print(self.minimizer._post_fit_warnings)


    def reset_to_best_fit(self):

        for pi in range(len(self.fit_par_free)):
            self.fit_par_free[pi].set(val=self.pout[pi])
            self.fit_par_free[pi].best_fit_val = self.pout[pi]
            self.fit_par_free[pi].best_fit_err = self.errors[pi]
            if  hasattr(self,'asymm_errors'):
                self.fit_par_free[pi].err_p=self.asymm_errors[pi][0]
                self.fit_par_free[pi].err_m=self.asymm_errors[pi][1]

        self.fit_model.eval()


class Minimizer(object):

    def __init__(self,model):
        self.model=model
        self._progress_iter = cycle(['|', '/', '-', '\\'])
        self._post_fit_warnings=''
        self.pbar=None
            #self._progess_bar(pbar, _res_sum, _res_sum_UL)

    def fit(self,model,
            max_ev=None,
            use_UL=False,
            silent=False):
        if silent is False:
            self.pbar = tqdm(total=None)
        self.use_UL = use_UL
        self.calls=0
        self.res_check=None
        self.model=model
        self.silent=silent
        self._original_silent_state = silent
        self._fit(max_ev)
        self._fit_stats()
        self._set_fit_errors()


    def _force_silent(self):
        self.silent=True

    def _restore_silent_state(self):
        self.silent = self._original_silent_state


    def _fit_stats(self):
        self.dof = len(self.model.data['x']) - self.model.free_pars
        self.chisq = self.get_chisq()
        self.chisq_red = self.chisq / float(self.dof)
        self.null_hyp_sig = 1.0 - chi2.cdf(self.chisq, self.dof)
        if self.use_UL==True:
            self.dof_no_UL = len(self.model.data['x']) - self.model.free_pars - self.model.data['UL'].sum()
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


    def _progess_bar(self, _res_sum, res_sum_UL):
        #if self.calls==1:
        #    print("minim function calls=%d, chisq=%5.5e UL part=%f" %(self.calls, _res_sum, -2.0*np.sum(res_UL)))
        #    print()
        if (np.mod(self.calls, 10) == 0 and self.calls != 0)  :
            #_c= ' ' * 256
            #print("\r%s"%_c,end="")
            #print("\r%s minim function calls=%d, chisq=%5.5e UL part=%f" % (next(self._progress_iter),self.calls, _res_sum, res_sum_UL), end="")
            m="minim. function calls=%d, chisq=%5.5e UL part=%f" %(self.calls, _res_sum, res_sum_UL)
            self.pbar.n=self.calls
            self.pbar.set_description(m)
            self.pbar.update(1)

    def residuals_Fit(self,
                      p,
                      fit_par,
                      data,
                      best_fit_SEDModel,
                      loglog,
                      chisq=False,
                      use_UL=False,
                      silent=False):
        #print('==>p', p)
        #print('==>fit_par_fre', self.model.fit_par_free)
        _warn=False
        for pi in range(len(fit_par)):
            _old_v=fit_par[pi].val
            #print('==> set', fit_par[pi].val,fit_par[pi].val_min,p[pi])
            fit_par[pi].set(val=p[pi])
            if np.isnan(p[pi]):
                _warn=True
                s1 = str('warning nan for par %s' % fit_par[pi].name)
                s2 = ', old parameter value was %s' % str(fit_par[pi].val)
                s =s1 + s2
                warnings.warn(s)
                self._post_fit_warnings+=s+'\n'
        #if _warn==True:
        #    best_fit_SEDModel.show_pars()
        #    print('res_sum',self._res_sum_chekc)
        #    print('res_chekc', self._res_chekc)
        #    print('res_UL_chekc', self._res_UL_chekc)


        model = best_fit_SEDModel.eval(nu=data['x'], fill_SED=False, get_model=True, loglog=loglog)


        _res_sum, _res, _res_sum_UL=_eval_res(data['y'], model, data['dy'], data['UL'], use_UL=use_UL)

        self._res_sum_chekc=_res_sum
        self._res_chekc = _res
        self._res_UL_chekc = _res_sum_UL
        #self._par_check=p
        #print('--> model', model[0],_res_sum)
        self.calls +=1

        if silent==False:
            #with tqdm(total=None, file=sys.stdout) as pbar:
            if self.pbar is not None:
                self._progess_bar(_res_sum, _res_sum_UL)
            #print("\r", end="")

        if chisq==True:
            res=_res_sum
        else:
            res=_res


        return res


    def get_chisq(self):
        model = self.model.fit_model.eval(nu=self.model.data['x'],
                                fill_SED=False,
                                get_model=True,
                                loglog=self.model.loglog)

        _res_sum, _res, _res_sum_UL = _eval_res(self.model.data['y'],
                                            model,
                                            self.model.data['dy'],
                                            self.model.data['UL'],
                                            use_UL=self.use_UL)
        return _res_sum

def _eval_res(data, model, data_error, UL, use_UL=False):

    res_no_UL = (data[~UL] - model[~UL]) / (data_error[~UL])
    res = (data - model) / (data_error)
    res_UL_log = [0]
    #print('UL.sum() ',UL.sum(),use_UL)
    if UL.sum() > 0 and use_UL is True:

        res_UL_log, res_UL = _eval_res_UL(data[UL], model[UL], data_error[UL])
        #print('ress_UL', res_UL)
        res[UL]=res_UL

    res_sum_UL= -2.0*np.sum(res_UL_log)
    res_sum=np.sum(res_no_UL * res_no_UL) + res_sum_UL
    return res_sum,res,res_sum_UL

def _eval_res_UL(y_UL,y_model,y_err):
    y = (y_UL - y_model) / (np.sqrt(2) *y_err)
    res_ul_chi2 = 0.5 * (1.0 + sp.special.erf(y))
    res_ul_chi2[res_ul_chi2 < 1E-300] = 1E-300

    res_ul = np.copy(y)
    res_ul[y > 0] = 0

    return  np.log(res_ul_chi2), res_ul



class LSBMinimizer(Minimizer):

    def __init__(self, model):

        super(LSBMinimizer, self).__init__(model)
        self.xtol=5.0E-8
        self.ftol = 5.0E-8
        self.factor=100.

    def _fit(self, max_ev,):
        #if self.use_UL is True:
        #    raise  RuntimeError('lsb minimizer currently is not supporting UL')

        bounds = [(par.fit_range_min, par.fit_range_max) for par in self.model.fit_par_free]
        max_nfev = 0 if (max_ev == 0 or max_ev == None) else max_ev
        pout, covar, info, mesg, success = leastsqbound(self.residuals_Fit,
                                                        self.model.pinit,
                                                        args=(self.model.fit_par_free,
                                                              self.model.data,
                                                              self.model.fit_model,
                                                              self.model.loglog,
                                                              False,
                                                              self.use_UL,
                                                              self.silent,),
                                                        xtol=self.xtol,
                                                        ftol=self.ftol,
                                                        factor=self.factor,
                                                        full_output=1,
                                                        bounds=bounds,
                                                        maxfev=max_nfev)

        self.mesg = mesg
        self.covar = covar
        #self.chisq =  sum(info["fvec"] * info["fvec"])

        self.pout = pout

    #def get_chisq(self):
    #    return self.chisq


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
                                  self.model.data,
                                  self.model.fit_model,
                                  self.model.loglog,
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

    #def get_chisq(self):
    #    return self.chisq




class MinuitMinimizer(Minimizer):

    def __init__(self,model):
        if minuit_installed==True:
            pass
        else:
            raise RuntimeError('iminuit non installed')

        super(MinuitMinimizer, self).__init__(model)



    def _fit(self,max_ev=None):
        bounds = [(par.fit_range_min, par.fit_range_max) for par in self.model.fit_par_free]
        self._set_minuit_func(self.model.pinit, bounds)
        #print('=>,')
        if max_ev is None or max_ev==0:
            max_ev =10000

        self.mesg=self.minuit_fun.migrad(ncall=max_ev)
        #print('==> mesg',self.mesg)
        if iminuit.__version__ < "2":
            self.pout=[self.minuit_fun.values[k] for k in self.minuit_fun.values.keys()]
            self.p=[self.minuit_fun.values[k] for k in self.minuit_fun.values.keys()]
        else:
            self.pout=[None]*len(self.minuit_fun.values)
            self.p = [None] * len(self.minuit_fun.values)
            for ID,val in enumerate(self.minuit_fun.values):
                self.pout[ID] = val
                self.p[ID] = val

    def _set_fit_errors(self):
        if iminuit.__version__ < "2":
            self.errors = [self.minuit_fun.errors[k] for k in self.minuit_fun.errors]
        else:
            self.errors = [None] * len(self.minuit_fun.values)
            for ID,err in enumerate(self.minuit_fun.errors):
                self.errors[ID] = err





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
        #print('==> iminuit.__version__',iminuit.__version__)
        for ID,par in enumerate(self.model.fit_par_free):
            self.minuit_par_name_dict[par.name]=p_names[ID]
            self.minuit_bounds_dict[par.name]=bounds[ID]
            self.par_dict[par.name]=par
        if iminuit.__version__ < "2":
            self.minuit_fun = iminuit.Minuit(fcn=self.chisq_func,
                name=p_names,
                pedantic=False,
                errordef=1,
                **kwdarg)
        else:
            values=tuple(par.val for ID,par in enumerate(self.model.fit_par_free))
            #print('=> values',values)
            self.minuit_fun = iminuit.Minuit(self.chisq_func,
                                             values,
                                             name=p_names)
            self.minuit_fun.limits=[ bounds[ID] for ID,par in enumerate(self.model.fit_par_free) ]
            #print('=>  self.minuit_fun.limits',  self.minuit_fun.limits)
            self.minuit_fun.errordef = Minuit.LEAST_SQUARES
    def chisq_func(self, *p):
        if iminuit.__version__ < "2":
            self.p = p
        else:
            self.p = p[0]

        #print('==>p',self.p)
        #print('==>fit_par_fre', self.model.fit_par_free)
        return self.residuals_Fit(self.p,
                          self.model.fit_par_free,
                          self.model.data,
                          self.model.fit_model,
                          self.model.loglog,
                          chisq=True,
                          use_UL=self.use_UL,
                          silent=self.silent)



    #def get_chisq(self):
    #    #print ('p',self.p)
    #    return self.chisq_func(*self.p)

    def minos_errors(self,par=None):
        try:
            if par is not None:
                par=self.minuit_par_name_dict[par]
            self.calls = 0
            self.minuit_fun.minos(var=par)
            self.minuit_fun.print_param()
            self.asymm_errors = [(self.minuit_fun.merrors[(k,1.0)],self.minuit_fun.merrors[(k,-1.0)]) for k in self.minuit_fun.errors.keys()]
            self.model.reset_to_best_fit()
        except:
            raise RuntimeWarning('Fit quality not sufficient to run Minos')

    def profile(self,par,bound=2,subtract_min=True):
        self._force_silent()
        try:
            self.calls = 0
            bound=self._set_bounds(par,bound=bound)

            x, y =  self.minuit_fun.profile(self.minuit_par_name_dict[par],
                                            bound=bound,
                                            subtract_min=subtract_min)

            self.model.reset_to_best_fit()
            self._restore_silent_state()
            return x,y
        except:
            self._restore_silent_state()
            raise RuntimeWarning('Fit quality not sufficient to run profile')

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
        try:
            self._force_silent()
            self.calls = 0
            bound = self._set_bounds(par, bound=bound)
            x, y,r =  self.minuit_fun.mnprofile(self.minuit_par_name_dict[par],
                                              bound=bound,
                                              subtract_min=subtract_min)
            #print('-->r',r)
            self.model.reset_to_best_fit()
            self._restore_silent_state()
            return x,y,r
        except:
            self._restore_silent_state()
            raise RuntimeWarning('Fit quality not sufficient to run mnprofile')


    def draw_mnprofile(self,par,bound=2):
        """

        Parameters
        ----------
        par
        bound

        Returns
        -------

        """

        try:
            self._force_silent()
            self.calls = 0
            bound = self._set_bounds(par, bound=bound)
            x,y,r=self.mnprofile(par,bound,subtract_min=True)
            fig,ax = self._draw_profile(x, y, par)
            self.model.reset_to_best_fit()
            self._restore_silent_state()
            return x, y, fig,ax

        except:
            self._restore_silent_state()
            raise RuntimeWarning('Fit quality not sufficient to run draw_mnprofile')

    def draw_profile(self,par,bound=2):
        try:
            self._force_silent()
            self.calls = 0
            bound = self._set_bounds(par, bound=bound)
            x, y = self.profile(par, subtract_min=True, bound=bound)
            fig,ax=self._draw_profile(x,y,par)
            self.model.reset_to_best_fit()
            self._restore_silent_state()
            return x,y,fig,ax
        except:
            self._restore_silent_state()
            raise RuntimeWarning('Fit quality not sufficient to run draw_profile')

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
        try:
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
            self._force_silent()
            x, y, z= self.minuit_fun.contour(self.minuit_par_name_dict[par_1],
                                             self.minuit_par_name_dict[par_2],
                                             subtract_min=subtract_min,
                                             bound=bound,
                                             bins=bins)
            self.model.reset_to_best_fit()
            self._restore_silent_state()
            return x,y,z
        except:
            self._restore_silent_state()
            raise RuntimeWarning('Fit quality not sufficient to run contour')


    def draw_contour(self,par_1,par_2,bound=2,levels=np.arange(5)):
        try:
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
        except:
            raise RuntimeWarning('Fit quality not sufficient to run draw_contour')

    def mncontour(self,par_1,par_2,bound=2,bins=20,subtract_min=True):
        try:
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
            self._force_silent()
            x, y, z= self.minuit_fun.mncontour(self.minuit_par_name_dict[par_1],
                                             self.minuit_par_name_dict[par_2],
                                             subtract_min=subtract_min,
                                             bound=bound,
                                             bins=bins)
            self.model.reset_to_best_fit()
            self._restore_silent_state()
            return x,y,z
        except:
            self._restore_silent_state()
            raise RuntimeWarning('Fit quality not sufficient to run mcontour')


    def draw_mncontour(self,par_1,par_2,bound=2,levels=np.arange(5)):
        try:
            self.calls = 0
            l=self.minuit_fun.errordef*levels
            x,y,z=self.mncontour(par_1,par_2,bound=bound)
            fig, ax = plt.subplots()
            CS = ax.contour(x, y, z, levels,)
            ax.clabel(CS, fontsize=9, inline=1)
            ax.set_xlabel(par_1)
            ax.set_ylabel(par_2)
            self.model.reset_to_best_fit()
            return x,y,z,fig,ax
        except:
            raise RuntimeWarning('Fit quality not sufficient to run draw_mcontour')

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

def fit_SED(fit_model, sed_data, nu_fit_start, nu_fit_stop, fitname=None, fit_workplace=None, loglog=False, silent=False,
            get_conf_int=False, max_ev=0, use_fake_err=False, minimizer='lsb', use_UL=False,repeat=3):

    mm = ModelMinimizer(minimizer)
    return mm,mm.fit(fit_model,
                  sed_data,
                  nu_fit_start,
                  nu_fit_stop,
                  fitname=fitname,
                  fit_workplace=fit_workplace,
                  loglog=loglog,
                  silent=silent,
                  get_conf_int=get_conf_int,
                  max_ev=max_ev,
                  use_fake_err=use_fake_err,
                  use_UL=use_UL,
                  repeat=repeat)



def fit_XY(fit_model, data, x_fit_start,x_fit_stop,fitname=None,silent=False,get_conf_int=False,minimizer='minuit',
           use_UL=False, repeat=1):

    mm = ModelMinimizer(minimizer)
    return mm, mm.fit(fit_model,
                      data,
                      nu_fit_start=x_fit_start,
                      nu_fit_stop=x_fit_stop,
                      fitname=fitname,
                      fit_workplace=None,
                      loglog=False,
                      silent=silent,
                      get_conf_int=get_conf_int,
                      max_ev=0,
                      use_fake_err=False,
                      use_UL=use_UL,
                      repeat=repeat)