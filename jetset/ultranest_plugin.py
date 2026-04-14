"""UltraNest plugin built on top of the JetSeT MCMC sampler interface."""

__author__ = "Andrea Tramacere"

import os
import time

import numpy as np

try:
    from ultranest import ReactiveNestedSampler
    from ultranest.stepsampler import SliceSampler
    from ultranest.stepsampler import generate_mixture_random_direction
    ultranest_installed = True
except Exception:
    on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
    if on_rtd:
        ReactiveNestedSampler = object
        ultranest_installed = False
    else:
        ultranest_installed = False

from .mcmc import McmcSampler, emcee_log_like

__all__ = ['UltraNestSampler', 'UltranestSampler']


class UltraNestSampler(McmcSampler):
    """UltraNest backend that reuses most of ``McmcSampler`` post-processing.

    Notes
    -----
    The class keeps the parameter/bounds workflow from :class:`~jetset.mcmc.McmcSampler`
    and swaps the sampling engine with ``ultranest.ReactiveNestedSampler``.
    Posterior samples are stored in ``self.samples`` with the same shape
    conventions used by ``McmcSampler`` so existing corner/model plotting
    helpers continue to work.
    """

    def __init__(self, model_minimizer, build_mcmc_parameters=True):
        """Create a new `UltraNestSampler` instance.

        Parameters
        ----------
        model_minimizer : object
            Initialized model-minimizer object.
        build_mcmc_parameters : bool, optional
            If ``True``, initialize JetSeT MCMC-style parameter metadata.
        """
        if ultranest_installed is not True:
            raise ImportError(
                'to use UltraNest plugin you need to install ultranest: '
                'https://johannesbuchner.github.io/UltraNest/'
            )

        super(UltraNestSampler, self).__init__(
            model_minimizer=model_minimizer,
            build_mcmc_parameters=build_mcmc_parameters,
        )
        self.result = None
        self.log_evidence = None
        self.log_evidence_err = None
        self.posterior_weighted_points = None
        self.posterior_weights = None
        self.posterior_weighted_logl = None
        self._ultranest_param_names = None

    def _check_ultranest_bounds(self):
        missing = []
        invalid = []
        for par in self._par_array_sampler:
            low = par.mcmc_bound_min
            high = par.mcmc_bound_max
            if low is None or high is None:
                missing.append(f'{par.model.name}.{par.name}')
                continue

            try:
                low = float(low)
                high = float(high)
                if not np.isfinite(low) or not np.isfinite(high) or low >= high:
                    invalid.append(f'{par.model.name}.{par.name}: [{low}, {high}]')
            except Exception:
                invalid.append(f'{par.model.name}.{par.name}: [{low}, {high}]')

        if missing or invalid:
            err_msg = [
                'UltraNest requires finite [min, max] bounds for all sampling parameters.'
            ]
            if missing:
                err_msg.append('Missing bounds: ' + ', '.join(missing))
            if invalid:
                err_msg.append('Invalid bounds: ' + ', '.join(invalid))
            raise RuntimeError('\n'.join(err_msg))

        self._build_sampler_bounds()

    def _build_ultranest_parameter_names(self):
        names = []
        used = set()
        for idx, par in enumerate(self._par_array_sampler):
            name = par.name
            if name in used:
                name = f'{par.name}_{par.model.name}'
            if name in used:
                name = f'{name}_{idx}'
            names.append(name)
            used.add(name)
        return names

    def _prior_transform(self, cube):
        cube = np.asarray(cube, dtype=float)
        theta = np.zeros(self.ndim, dtype=float)
        for idx, bounds in enumerate(self._bounds_sampler):
            low, high = bounds
            theta[idx] = low + cube[idx] * (high - low)
        return theta

    def _extract_posterior_samples(self, posterior_samples_size=None, rnd_seed=0):
        result = self.result if isinstance(self.result, dict) else {}
        weighted = result.get('weighted_samples', {})

        points = None
        weights = None
        logl = None
        if isinstance(weighted, dict):
            points = np.asarray(weighted.get('points', []))
            weights = np.asarray(weighted.get('weights', []))
            logl = np.asarray(weighted.get('logl', []))

        if points is not None and points.ndim == 2 and points.shape[1] == self.ndim:
            self.posterior_weighted_points = points

            if weights is not None and weights.ndim == 1 and weights.shape[0] == points.shape[0]:
                weights = np.asarray(weights, dtype=float)
                if np.any(~np.isfinite(weights)) or np.sum(weights) <= 0:
                    weights = np.ones(points.shape[0], dtype=float) / float(points.shape[0])
                else:
                    weights = weights / np.sum(weights)

                self.posterior_weights = weights
                if logl is not None and logl.ndim == 1 and logl.shape[0] == points.shape[0]:
                    self.posterior_weighted_logl = logl
                else:
                    self.posterior_weighted_logl = None

                if posterior_samples_size is None:
                    posterior_samples_size = points.shape[0]
                posterior_samples_size = max(1, int(posterior_samples_size))

                rng = np.random.default_rng(rnd_seed)
                draw_idx = rng.choice(
                    np.arange(points.shape[0]),
                    size=posterior_samples_size,
                    replace=True,
                    p=weights,
                )

                samples = points[draw_idx]
                if self.posterior_weighted_logl is not None:
                    samples_logl = self.posterior_weighted_logl[draw_idx]
                else:
                    samples_logl = np.full(samples.shape[0], np.nan)

                return samples, samples_logl

            if logl is not None and logl.ndim == 1 and logl.shape[0] == points.shape[0]:
                return points, logl

            return points, np.full(points.shape[0], np.nan)

        samples = np.asarray(result.get('samples', []))
        if samples.ndim != 2 or samples.shape[1] != self.ndim:
            raise RuntimeError(
                'unable to build posterior samples from UltraNest output; '
                'please check result content'
            )
        return samples, np.full(samples.shape[0], np.nan)

    def run_sampler(
        self,
        min_num_live_points=400,
        dlogz=0.5,
        min_ess=None,
        frac_remain=0.01,
        max_ncalls=None,
        use_UL=False,
        loglog=False,
        resume='subfolder',
        log_dir=None,
        show_status=True,
        posterior_samples_size=None,
        rnd_seed=0,
        **run_kwargs,
    ):
        """Run posterior sampling with UltraNest.

        Parameters
        ----------
        min_num_live_points : int, optional
            Minimum number of live points.
        dlogz : float, optional
            Stopping criterion on remaining log-evidence.
        min_ess : int, optional
            Minimum effective sample size target.
        frac_remain : float, optional
            Fractional prior volume stopping criterion.
        max_ncalls : int, optional
            Maximum number of likelihood calls.
        use_UL : bool, optional
            If ``True``, include upper-limit terms in the likelihood.
        loglog : bool, optional
            If ``True``, evaluate model and data in log10 space.
        resume : str, optional
            UltraNest resume mode.
        log_dir : str, optional
            UltraNest output directory.
        show_status : bool, optional
            If ``True``, show UltraNest progress/status.
        posterior_samples_size : int, optional
            Number of equal-weight posterior samples to draw from weighted
            nested-sampling samples. If ``None``, use all weighted points.
        rnd_seed : int, optional
            Seed for posterior resampling.
        **run_kwargs : dict
            Additional keyword arguments passed to ``ReactiveNestedSampler.run``.
        """
        self._check_ultranest_bounds()

        self.calls = 0
        self.calls_OK = 0
        self.use_UL = use_UL
        self.ndim = len(self._par_array_sampler)

        if log_dir is None:
            log_dir = f'ultranest_{self.model.name}'

        calls_counter = {'count': 0}

        def _loglike(theta):
            calls_counter['count'] += 1
            return emcee_log_like(
                theta,
                self.model,
                self.data,
                use_UL,
                self._par_array_sampler,
                loglog,
            )

        self._ultranest_param_names = self._build_ultranest_parameter_names()
        self.sampler = ReactiveNestedSampler(
            self._ultranest_param_names,
            _loglike,
            self._prior_transform,
            log_dir=log_dir,
            resume=resume,
        )
        self.sampler.stepsampler = SliceSampler(nsteps=1 *self.ndim, generate_direction=generate_mixture_random_direction)
        run_args = {
            'min_num_live_points': min_num_live_points,
            'dlogz': dlogz,
            'frac_remain': frac_remain,
            'show_status': show_status,
        }
        if min_ess is not None:
            run_args['min_ess'] = min_ess
        if max_ncalls is not None:
            run_args['max_ncalls'] = max_ncalls
        run_args.update(run_kwargs)

        print('ultranest run starting')
        print('')
        start = time.time()
        self.result = self.sampler.run(**run_args)
        end = time.time()
        comp_time = end - start
        print('ultranest run done, took %2.2f seconds' % comp_time)

        self.samples, self.samples_log_prob = self._extract_posterior_samples(
            posterior_samples_size=posterior_samples_size,
            rnd_seed=rnd_seed,
        )
        self.chain = self.samples[np.newaxis, :, :]
        self.log_prob_chain = self.samples_log_prob[np.newaxis, :]
        self.burnin = 0

        if isinstance(self.result, dict):
            self.log_evidence = self.result.get('logz')
            self.log_evidence_err = self.result.get('logzerr')

        self.calls = calls_counter['count']
        self.calls_OK = self.calls
        self.calls_tot = self.calls if max_ncalls is None else max_ncalls
        self.acceptance_fraction = np.nan

        self.reset_to_mcmc_best_fit()

    def reset_to_mcmc_best_fit(self, verbose=True):
        """Reset model parameters to the UltraNest best-fit point."""
        theta = None
        if isinstance(self.result, dict):
            ml = self.result.get('maximum_likelihood', {})
            if isinstance(ml, dict):
                theta = ml.get('point')

        if theta is None:
            if (
                hasattr(self, 'samples_log_prob')
                and self.samples_log_prob is not None
                and np.any(np.isfinite(self.samples_log_prob))
            ):
                theta = self.samples[np.nanargmax(self.samples_log_prob)]
            else:
                theta = np.median(self.samples, axis=0)

        if theta is None:
            raise RuntimeError('unable to determine best-fit point from UltraNest output')

        theta = np.asarray(theta, dtype=float)
        if theta.shape[0] != len(self._par_array_sampler):
            raise RuntimeError(
                'best-fit point dimensionality mismatch: '
                f'{theta.shape[0]} != {len(self._par_array_sampler)}'
            )

        if verbose:
            print('----------------------------')
            print('UltraNest best fit solution')

        for idx, par in enumerate(self._par_array_sampler):
            par.val = theta[idx]
            self.model.parameters.set_par(
                model_name=par.model.name,
                par_name=par.name,
                val=par.val,
            )
            par.best_fit_mcmc_val = par.val

            quantiles = self.get_par_quantiles(
                par_name=par.name,
                comp_name=par.model.name,
                quantiles=(0.16, 0.5, 0.84),
            )
            par.q_16 = quantiles[0]
            par.q_50 = quantiles[1]
            par.q_84 = quantiles[2]

            if verbose:
                print(
                    'comp: %s par: %s ultranest best fit val: %s '
                    'quantiles(0.16,0.5,0.84): %s'
                    % (par.model.name, par.name, par.val, quantiles)
                )

        if verbose:
            if self.log_evidence is not None:
                print(
                    'logZ=%s logZerr=%s'
                    % (str(self.log_evidence), str(self.log_evidence_err))
                )
            print('----------------------------')

    @property
    def logz(self):
        """Return the estimated Bayesian log-evidence."""
        return self.log_evidence

    @property
    def logzerr(self):
        """Return the uncertainty on the Bayesian log-evidence."""
        return self.log_evidence_err


UltranestSampler = UltraNestSampler
