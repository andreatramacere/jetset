import pytest
from .base_class import TestBase

class TestEmcee(TestBase):

    def integration_suite(self,fit_dict=None,sed_number=None,plot=False):
        if sed_number is not None and fit_dict is None:
            from .test_model_fit import TestModelFit
            t=TestModelFit()
            fit_dict=t.integration_suite(plot=plot,sed_number=sed_number,skip_minuit=True)
        elif fit_dict and sed_number is None:
            pass
        else:
            raise RecursionError("please provide either fit_dict or sed_number")
        self.run_emcee(fit_dict,plot=plot)

    def test(self,plot=False):
        from .test_model_fit import TestModelFit
        t=TestModelFit()
        fit_dict=t.integration_suite(plot=plot,sed_number=1,skip_minuit=True)
        self.run_emcee(fit_dict=fit_dict,plot=plot)

    def run_emcee(self,fit_dict=None,model_minimizer=None,sed_data=None,plot=False):

        if fit_dict is not None:
            model_minimizer = fit_dict['model_minimizer']
            sed_data = fit_dict['sed_data']
        elif sed_data is None or model_minimizer is None:
            raise RuntimeError("please, provide either fit_dict, or both  sed_data and model_minimizer")
        else:
            pass

        from jetset.mcmc import McmcSampler

        mcmc = McmcSampler(model_minimizer)
        labels = ['N', 'B', 'beam_obj', 's', 'gamma0_log_parab']
        model_name = 'jet_leptonic'
        use_labels_dict = {model_name: labels}
        mcmc.run_sampler(nwalkers=64, burnin=10, steps=50, bound=5.0, bound_rel=True, threads=None,
                        walker_start_bound=0.005, use_labels_dict=use_labels_dict)

        print(mcmc.acceptance_fraction)
        if plot is True:

            p = mcmc.plot_model(sed_data=sed_data, fit_range=[11., 27.4], size=50)
            p.setlim(y_min=1E-13, x_min=1E6, x_max=3E28)

        mcmc.save('mcmc_sampler.pkl')

        ms = McmcSampler.load('mcmc_sampler.pkl')

        if plot is True:
            p = ms.plot_model(sed_data=sed_data, fit_range=[11., 27.4], size=50)
            p.setlim(y_min=1E-13, x_min=1E6, x_max=3E28)