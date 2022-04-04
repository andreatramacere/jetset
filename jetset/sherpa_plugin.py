__author__ = "Andrea Tramacere"


try:
    from sherpa.models.model import ArithmeticModel, modelCacher1d, RegriddableModel1D
    from sherpa.models.parameter import Parameter
    from sherpa import data
    from sherpa.fit import Fit
    from sherpa.stats import Chi2
    from sherpa.optmethods import LevMar
    from sherpa.fit import Fit
    from sherpa.stats import Chi2
    from sherpa import data as sherpa_data

    sherpa_installed = True

except:
    raise  ImportError('to use sherpa plugin you need to install sherpa: https://sherpa.readthedocs.io/en/latest/install.html')

import  numpy as np
from .plot_sedfit import  PlotSED
from .minimizer import  Minimizer

__all__=['JetsetSherpaModel','plot_sherpa_model']


class JetsetSherpaModel(RegriddableModel1D):
    """
    authomatic sherpa model generator
    """

    def __init__(self, jetset_model,par_list=None,clone=False):
        if clone is True:
            self._jetset_model = jetset_model.clone()
        else:
            self._jetset_model=jetset_model

        self._jp_list = []
        self._jp_par_array = []
        self._jp_list_names = []
        setattr(self, '_jetset_ncalls', 0)
        keep=True
        #print('-->, ', par_list)
        for p in self._jetset_model.parameters.par_array:
            if par_list is None:
                keep = True
            else:
                keep = p in par_list
            #print('-->, ',p, p.name, keep)
            if p is not None and keep is True:

                if p.name.lower() in self._jp_list_names or p.name.upper() in self._jp_list_names:
                    name = p.name + '_sh'
                    print('jetset model name', p.name, 'renamed to ', name, 'due to sherpa internal naming convention')
                else:
                    name = p.name
                if p.fit_range_min is not None:
                    val_min=p.fit_range_min
                else:
                    val_min = p.val_min

                if p.fit_range_max is not None:
                    val_max = p.fit_range_max
                else:
                    val_max = p.val_max

                sh_p = Parameter(self._jetset_model.name, name, p.val, min=val_min, max=val_max, units=p.units)
                setattr(self, sh_p.name, sh_p)
                p._sherpa_ref = sh_p
                if np.isnan(sh_p.max):
                    sh_p.max = sh_p.hard_max
                if np.isnan(sh_p.min):
                    sh_p.min = sh_p.hard_min

                self._jp_list.append(sh_p)
                self._jp_par_array.append(p)
                self._jp_list_names.append(p.name)
        RegriddableModel1D.__init__(self, jetset_model.name,(p._sherpa_ref for p in self._jp_par_array))


    def calc(self, pars, x):
        for ID, p in enumerate(self._jp_list):
            j_p = self._jp_par_array[ID]
            j_p.val = p.val
        self._jetset_ncalls +=1
        return self._jetset_model.eval(get_model=True, nu=x)

    def plot_model(self, fit_range, model_range=[1E10, 1E30], nu_grid_size=200, plot_obj=None, sed_data=None):
        self._jetset_model.set_nu_grid(model_range[0], model_range[1], nu_grid_size)
        self._jetset_model.eval()
        plot_obj = self._jetset_model.plot_model(plot_obj=plot_obj, sed_data=sed_data)
        plot_obj.add_model_residual_plot(data=sed_data, model=self._jetset_model,
                                         fit_range=[fit_range[0], fit_range[1]])


def plot_sherpa_model(sherpa_model, fit_range=None, model_range=[1E10, 1E30], nu_grid_size=200, sed_data=None,
                      add_res=False, plot_obj=None, label=None, line_style=None):
    if fit_range is not None:
        x = np.logspace(np.log10(fit_range[0]), np.log10(fit_range[1]), nu_grid_size)
    else:
        x = np.logspace(np.log10(model_range[0]), np.log10(model_range[1]), nu_grid_size)
    y = sherpa_model(x)

    if plot_obj is None:
        plot_obj = PlotSED(frame='obs', density=False)

    if sed_data is not None:
        plot_obj.add_data_plot(sed_data=sed_data)

    plot_obj.add_xy_plot(x, y, label=label, line_style=line_style)

    if add_res is True and fit_range is not None:
        nufnu_res = sherpa_model(sed_data.data['nu_data'])
        y_res = (sed_data.data['nuFnu_data'] - nufnu_res) / sed_data.data['dnuFnu_data']
        x_res = sed_data.data['nu_data']
        plot_obj.add_xy_residual_plot(x=x_res, y=y_res, fit_range= [fit_range[0], fit_range[1]])

    return plot_obj


class SherpaMinimizer(Minimizer):

    def __init__(self, model,method=LevMar(),stat=Chi2()):
        if sherpa_installed is True:
            pass
        else:
            raise ImportError('sherpa not installed, \n to use sherpa plugin you need to install sherpa: https://sherpa.readthedocs.io/en/latest/install.html')

        super(SherpaMinimizer, self).__init__(model)
        self._method=method
        self._stat=stat
        self._sherpa_model = None
        self._sherpa_data = None
        self.pbar = None

    def _create_sherpa_model(self):
        self._sherpa_model = JetsetSherpaModel(jetset_model = self.model.fit_model, par_list=self.model.fit_par_free)

    def _create_sherpa_data(self):
        self._sherpa_data = sherpa_data.Data1D("sed", self.model.data['x'], self.model.data['y'], staterror=self.model.data['dy'])

    @property
    def sherpa_fitter(self):
        return self._sherpa_fitter

    @property
    def calls(self):
        if self._sherpa_model is not None:
            return self._sherpa_model._jetset_ncalls
        else:
            return None

    @calls.setter
    def calls(self,n):
        if self._sherpa_model is not None:
            self._sherpa_model._jetset_ncalls = n


    def _fit(self, max_ev,):
        self._create_sherpa_model()
        self._create_sherpa_data()
        self._sherpa_model._jetset_ncalls = 0


        self._sherpa_fitter=Fit(self._sherpa_data,self._sherpa_model, method=self._method,stat=self._stat)

        self.mesg  = self._sherpa_fitter.fit()
        self.covar = self.mesg.covar
        self.pout = [p for p in self.mesg.parvals]
        self.p = [p for p in self.mesg.parvals]

    def _set_fit_errors(self):
        self.errors = [np.sqrt(np.fabs(self.covar[pi, pi])) for pi in range(len(self.model.fit_par_free))]