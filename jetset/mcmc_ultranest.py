"""UltraNest backend built on top of the JetSeT MCMC sampler interface."""

__author__ = "Andrea Tramacere"

import os
import time
import shutil
import subprocess
from pathlib import Path
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

__all__ = ['UltraNestSampler', 'run_open_mpi']


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
                'IMPORTANT: UltraNest support requires ultranest>=4.0.\n'
                'Install with pip: pip install "ultranest>=4.0"\n'
                'Install with conda: conda install -c conda-forge "ultranest>=4.0"\n'
                'More info: https://johannesbuchner.github.io/UltraNest/'
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
        ultranest_output_dir=None,
        show_status=True,
        posterior_samples_size=None,
        rnd_seed=0,
        use_stepsampler=True,
        nsteps=2,
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
        ultranest_output_dir : str, optional
            This is where UltraNest writes its run files. 
            If ``None`` ultranest_output_dir = f'ultranest_{self.model.name}'
        show_status : bool, optional
            If ``True``, show UltraNest progress/status.
        posterior_samples_size : int, optional
            Number of equal-weight posterior samples to draw from weighted
            nested-sampling samples. If ``None``, use all weighted points.
        rnd_seed : int, optional
            Seed for posterior resampling.
        use_stepsampler : bool, optional
            If ``True``, attach a ``SliceSampler`` to the UltraNest sampler.
        nsteps : int, optional
            Slice-sampler steps per parameter dimension. The effective value is
            ``nsteps * self.ndim``. Default is ``2``.
        **run_kwargs : dict
            Additional keyword arguments passed to ``ReactiveNestedSampler.run``.

        Notes
        -----
        On completion, posterior samples are stored in ``self.samples`` and
        model parameters are reset to the mcmc best-fit point.
        """
        self._check_ultranest_bounds()

        self.calls = 0
        self.calls_OK = 0
        self.use_UL = use_UL
        self.ndim = len(self._par_array_sampler)

        if ultranest_output_dir is None:
            ultranest_output_dir = f'ultranest_{self.model.name}'

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
            log_dir=ultranest_output_dir,
            resume=resume,
        )
        
        if use_stepsampler:
            self.sampler.stepsampler = SliceSampler(nsteps=nsteps *self.ndim, generate_direction=generate_mixture_random_direction)
        
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



def run_open_mpi(sampler,
                 n_proc=8,
                 ultranest_output_dir=None,
                 min_num_live_points=64,
                 nsteps=1,
                 max_num_improvement_loops=1,
                 min_ess=100,
                 num_c_threads=1,
                 extra_mpirun_args=''):

    """Run ``UltraNestSampler.run_sampler`` through ``mpirun``.

    Parameters
    ----------
    sampler : UltraNestSampler
        Configured sampler instance to serialize and execute.
    n_proc : int, optional
        Number of MPI ranks passed to ``mpirun -np``. Default is ``8``.
    ultranest_output_dir : str, optional
        Name passed to ``run_sampler``, and resolved within run_mpi/ directory
    min_num_live_points : int, optional
        Minimum live points forwarded to ``run_sampler``.
    nsteps : int, optional
        Slice-sampler step multiplier forwarded to ``run_sampler``.
    max_num_improvement_loops : int, optional
        Value forwarded to ``run_sampler``.
    min_ess : int, optional
        Minimum effective sample size target forwarded to ``run_sampler``.
    num_c_threads : int, optional
        Thread count forwarded to ``mcmc.model.set_num_c_threads`` in the
        generated MPI worker script.
    extra_mpirun_args : str, optional
        Extra command-line arguments inserted after ``mpirun`` and before
        ``-np`` (for example hostfile or binding flags).

    Returns
    -------
    UltraNestSampler
        Sampler reloaded from ``run_mpi/sampler.pkl`` after MPI completion.

    Notes
    -----
    The function deletes and recreates ``./run_mpi`` in the current working
    directory, writes a helper Python script there, and temporarily sets
    ``OMP_NUM_THREADS`` to ``n_proc`` for the parent process while the MPI
    command runs.
    """
    CURRENT_DIR = Path.cwd().resolve()
    run_dir = CURRENT_DIR / 'run_mpi'
    script_name = 'run_jetset_ultranest.py'
    OMP_NUM_THREADS_INITIAL=os.getenv("OMP_NUM_THREADS")
    os.environ["OMP_NUM_THREADS"] = "%s"%str(int(n_proc))

    command = f"mpirun {extra_mpirun_args} -np {n_proc} python {script_name}"

    if ultranest_output_dir is None:
            ultranest_output_dir = f'ultranest_{sampler.model.name}'
    try:
        if run_dir.exists():
            shutil.rmtree(run_dir)
        run_dir.mkdir(parents=True, exist_ok=True)

        sampler_path = run_dir / 'sampler.pkl'
        sampler.save(str(sampler_path))
        #(run_dir / out_dir).mkdir(parents=True, exist_ok=True)
    
        sampler_run = f"""
mcmc.model.set_num_c_threads({int(num_c_threads)})

mcmc.run_sampler(
                min_num_live_points={int(min_num_live_points)},
                nsteps={int(nsteps)},
                max_num_improvement_loops={int(max_num_improvement_loops)},
                min_ess={float(min_ess)},
                ultranest_output_dir={repr(str(ultranest_output_dir))},
            )
    """

        script = "from jetset.mcmc_ultranest import UltraNestSampler\n"
        script += "mcmc=UltraNestSampler.load('sampler.pkl')\n"
        script += f"{sampler_run}\n"
        script += "mcmc.save('sampler.pkl')\n"
        print("====== ultranest script ========")
        print(script)
        script_path = run_dir / script_name
        with open(script_path, 'w') as f:
            f.write(script)
        print("================================")

        print("executing", command)
        subprocess.run(command, shell=True, check=True, cwd=run_dir)
        mcmc= UltraNestSampler.load(str(sampler_path))
        
    except Exception as e:
        raise RuntimeError(e)

    finally:

        if OMP_NUM_THREADS_INITIAL is not None:
            os.environ["OMP_NUM_THREADS"] = "%s"%str(int(OMP_NUM_THREADS_INITIAL))
        else:
            os.unsetenv("OMP_NUM_THREADS")


    return mcmc
