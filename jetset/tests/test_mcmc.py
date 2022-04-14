import pytest
from .base_class import TestBase

class TestEmcee(TestBase):

    def integration_suite(self,fit_dict=None,sed_number=None,plot=False):
        if sed_number is not None and fit_dict is None:
            from .test_model_fit import TestModelFit
            t=TestModelFit()
            fit_dict=t.integration_suite(plot=plot,sed_number=sed_number)
        elif fit_dict and sed_number is None:
            pass
        else:
            raise RecursionError("please provide either fit_dict or sed_number")
        self.run_emcee(fit_dict,plot=plot)


    def test_all(self,sed_number=1,plot=False):
        self.integration_suite(sed_number=sed_number,plot=plot)

    def run_emcee(self,fit_dict=None,model_minimizer=None,sed_data=None,plot=False):
        #fit_model = fit_dict['fit_model']

        if fit_dict is not None:
            model_minimizer = fit_dict['model_minimizer']
            sed_data = fit_dict['sed_data']
        elif sed_data is None or model_minimizer is None:
            raise RuntimeError("please, provide either fit_dict, or both  sed_data and model_minimizer")
        else:
            pass

        from jetset.mcmc import McmcSampler

        mcmc = McmcSampler(model_minimizer)
        mcmc.run_sampler(nwalkers=100, burnin=10, steps=20,threads=1)
        labels = ['N', 'B', 'beam_obj', 's', 'gamma0_log_parab']
        model_name = 'jet_leptonic'
        use_labels_dict = {model_name: labels}

        mcmc.run_sampler(nwalkers=128, burnin=10, steps=50, bound=5.0, bound_rel=True, threads=None,
                        walker_start_bound=0.005, use_labels_dict=use_labels_dict)

        print(mcmc.acceptance_fraction)
        if plot is True:

            p = mcmc.plot_model(sed_data=sed_data, fit_range=[11., 27.4], size=50)
            p.rescale(y_min=-13, x_min=6, x_max=28.5)

        mcmc.save('mcmc_sampler.pkl')

        ms = McmcSampler.load('mcmc_sampler.pkl')



        mcmc.run_sampler(nwalkers=50, burnin=10, steps=40, bound=5.0, bound_rel=True, threads=None,
                        walker_start_bound=0.005, use_labels_dict=use_labels_dict)

        if plot is True:
            p = ms.plot_model(sed_data=sed_data, fit_range=[11., 27.4], size=50)
            p.rescale(y_min=-13, x_min=6, x_max=28.5)