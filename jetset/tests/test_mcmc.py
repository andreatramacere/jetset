import pytest

def test_emcee(jet_bf_model,model_minimizer):
    from jetset.mcmc import McmcSampler

    mcmc = McmcSampler(model_minimizer)
    mcmc.run_sampler(nwalkers=150, burnin=10, steps=50,threads=2)
