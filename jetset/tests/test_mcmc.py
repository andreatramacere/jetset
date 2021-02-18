import pytest
from .test_model_fit import test_model_fit

def test_emcee(plot=False):
    from jetset.mcmc import McmcSampler

    fit_model, model_minimizer,sed_data = test_model_fit(minimizer='lsb', sed_number=1)

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