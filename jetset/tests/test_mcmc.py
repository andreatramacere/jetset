import pytest
import numpy as np
from .base_class import TestBase

class TestEmcee(TestBase):

    def integration_suite(self,fit_dict=None,sed_number=None,plot=False):
        if sed_number is not None and fit_dict is None:
            from .test_model_fit import prepare_asset
            fit_dict=prepare_asset(plot=plot,sed_number=sed_number,skip_minuit=True)
        elif fit_dict and sed_number is None:
            pass
        else:
            raise RecursionError("please provide either fit_dict or sed_number")
        self.run_emcee(fit_dict,plot=plot)

    def test(self,plot=False):
        from .test_model_fit import prepare_asset
        fit_dict=prepare_asset(plot=plot,sed_number=1,skip_minuit=True)
        self.run_emcee(fit_dict=fit_dict,plot=plot)
        self.run_ultranest(fit_dict=fit_dict,plot=plot)

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

        # Check freeze/free workflow on sampler parameters.
        free_non_linked = [
            p for p in mcmc._par_array_sampler
            if (getattr(p, '_is_dependent', False) is False and getattr(p, '_linked', False) is False)
        ]
        assert len(free_non_linked) > 0
        par_to_toggle = free_non_linked[0]
        model_name = par_to_toggle.model.name
        par_name = par_to_toggle.name
        n_free_before = len(mcmc._par_array_sampler)

        mcmc.model.parameters.freeze(model_name, par_name)
        frozen_par = mcmc.model.parameters.get_par_by_name(model_name, par_name)
        assert frozen_par.frozen is True
        assert len(mcmc._par_array_sampler) == n_free_before - 1

        mcmc.model.parameters.free(model_name, par_name)
        thawed_par = mcmc.model.parameters.get_par_by_name(model_name, par_name)
        assert thawed_par.frozen is False
        assert len(mcmc._par_array_sampler) == n_free_before

        # Redefine labels for a subset of sampled parameters.
        custom_labels = {}
        for par in mcmc._par_array_sampler[:min(5, len(mcmc._par_array_sampler))]:
            new_label = f"{par.model.name}:{par.name}"
            mcmc.set_plot_label(par_name=par.name, plot_label=new_label, comp_name=par.model.name)
            custom_labels[(par.model.name, par.name)] = new_label

        mcmc.set_bounds(bound=5.0,bound_rel=True)
        mcmc.run_sampler(nwalkers=64, burnin=50, steps=500)

        assert np.isfinite(mcmc.acceptance_fraction)
        assert 0.0 <= mcmc.acceptance_fraction <= 1.0
        assert mcmc.samples.ndim == 2
        assert mcmc.samples.shape[0] > 0
        assert mcmc.samples.shape[1] == len(mcmc._par_array_sampler)
        assert mcmc.chain.ndim == 3
        assert mcmc.chain.shape[0] == 64
        assert mcmc.chain.shape[2] == len(mcmc._par_array_sampler)
        assert mcmc.log_prob_chain.ndim == 2
        assert mcmc.log_prob_chain.shape[0] == 64
        assert np.all(np.isfinite(mcmc.samples))
        assert np.all(np.isfinite(mcmc.samples_log_prob))

        for par in mcmc._par_array_sampler:
            assert np.isfinite(par.best_fit_mcmc_val)
            assert par.q_16 is not None
            assert par.q_50 is not None
            assert par.q_84 is not None

        # Check autocorrelation handling and burn-in retuning.
        tau = np.asarray(mcmc.sampler.get_autocorr_time(tol=0), dtype=float)
        assert tau.ndim == 1
        assert tau.shape[0] == len(mcmc._par_array_sampler)
        assert np.all(np.isfinite(tau))
        assert np.all(tau > 0)

        tau_max = float(np.max(tau))
        tau_coeff = min(1.0, (mcmc.steps - 10) / max(tau_max, 1.0))
        expected_burnin = int(tau_coeff * tau_max)
        mcmc.tune_burnin()
        assert mcmc.burnin == expected_burnin
        assert mcmc.burnin < mcmc.steps
        assert mcmc.samples.shape[0] > 0
        assert mcmc.samples_log_prob.shape[0] == mcmc.samples.shape[0]
        assert np.all(np.isfinite(mcmc.samples))
        assert np.all(np.isfinite(mcmc.samples_log_prob))

        if plot is True:

            p = mcmc.plot_model(sed_data=sed_data, fit_range=[11., 27.4], size=50)
            p.setlim(y_min=1E-13, x_min=1E6, x_max=3E28)

        mcmc.save('mcmc_sampler.pkl')

        ms = McmcSampler.load('mcmc_sampler.pkl')
        assert ms.samples.shape == mcmc.samples.shape
        assert ms.chain.shape == mcmc.chain.shape
        assert ms.log_prob_chain.shape == mcmc.log_prob_chain.shape
        np.testing.assert_allclose(ms.samples, mcmc.samples)
        np.testing.assert_allclose(ms.samples_log_prob, mcmc.samples_log_prob)
        assert np.isfinite(ms.acceptance_fraction)
        assert 0.0 <= ms.acceptance_fraction <= 1.0

        for par in ms._par_array_sampler:
            assert np.isfinite(par.best_fit_mcmc_val)
            assert par.q_16 is not None
            assert par.q_50 is not None
            assert par.q_84 is not None

        for (comp_name, p_name), expected_label in custom_labels.items():
            assert ms.get_par(p_name, comp_name=comp_name).plot_label == expected_label

        if plot is True:
            p = ms.plot_model(sed_data=sed_data, fit_range=[11., 27.4], size=50)
            p.setlim(y_min=1E-13, x_min=1E6, x_max=3E28)

        
    def run_ultranest(self,fit_dict=None,model_minimizer=None,sed_data=None,plot=False):
        if fit_dict is not None:
            model_minimizer = fit_dict['model_minimizer']
            sed_data = fit_dict['sed_data']
        elif sed_data is None or model_minimizer is None:
            raise RuntimeError("please, provide either fit_dict, or both  sed_data and model_minimizer")
        else:
            pass
        from jetset.mcmc_ultranest import UltraNestSampler
        mcmc = UltraNestSampler(model_minimizer)

          # Check freeze/free workflow on sampler parameters.
        free_non_linked = [
            p for p in mcmc._par_array_sampler
            if (getattr(p, '_is_dependent', False) is False and getattr(p, '_linked', False) is False)
        ]
        assert len(free_non_linked) > 0
        par_to_toggle = free_non_linked[0]
        model_name = par_to_toggle.model.name
        par_name = par_to_toggle.name
        n_free_before = len(mcmc._par_array_sampler)

        mcmc.model.parameters.freeze(model_name, par_name)
        frozen_par = mcmc.model.parameters.get_par_by_name(model_name, par_name)
        assert frozen_par.frozen is True
        assert len(mcmc._par_array_sampler) == n_free_before - 1

        mcmc.model.parameters.free(model_name, par_name)
        thawed_par = mcmc.model.parameters.get_par_by_name(model_name, par_name)
        assert thawed_par.frozen is False
        assert len(mcmc._par_array_sampler) == n_free_before
        mcmc.model.freeze_all()
        for ID in range(3):
            mcmc.model.parameters.par_array[ID].free()

        print(mcmc.model.parameters)
        # Redefine labels for a subset of sampled parameters.
        custom_labels = {}
        for par in mcmc._par_array_sampler[:min(5, len(mcmc._par_array_sampler))]:
            new_label = f"{par.model.name}:{par.name}"
            mcmc.set_plot_label(par_name=par.name, plot_label=new_label, comp_name=par.model.name)
            custom_labels[(par.model.name, par.name)] = new_label

        mcmc.set_bounds(bound=5.0,bound_rel=True)
        mcmc.run_sampler(min_num_live_points=64,
                 dlogz=1.0,
                 dKL=np.inf,
                 frac_remain=0.5,
                 nsteps=1,
                 max_num_improvement_loops=1, 
                 min_ess=10)

        # Post-run checks.
        assert mcmc.samples.ndim == 2
        assert mcmc.samples.shape[0] > 0
        assert mcmc.samples.shape[1] == len(mcmc._par_array_sampler)
        assert mcmc.chain.ndim == 3
        assert mcmc.chain.shape[0] == 1
        assert mcmc.chain.shape[2] == len(mcmc._par_array_sampler)
        assert mcmc.log_prob_chain.ndim == 2
        assert mcmc.log_prob_chain.shape[0] == 1
        assert np.all(np.isfinite(mcmc.samples))
        assert np.all(np.isfinite(mcmc.samples_log_prob) | np.isnan(mcmc.samples_log_prob))
        assert mcmc.calls > 0
        assert mcmc.calls_OK == mcmc.calls
        assert mcmc.result is not None

        if mcmc.log_evidence is not None:
            assert np.isfinite(mcmc.log_evidence)
        if mcmc.log_evidence_err is not None:
            assert np.isfinite(mcmc.log_evidence_err)

        for par in mcmc._par_array_sampler:
            assert np.isfinite(par.best_fit_mcmc_val)
            assert par.q_16 is not None
            assert par.q_50 is not None
            assert par.q_84 is not None

        for (comp_name, p_name), expected_label in custom_labels.items():
            assert mcmc.get_par(p_name, comp_name=comp_name).plot_label == expected_label
