import pytest
import numpy as np
from .base_class import TestBase

class TestUltranest(TestBase):

    def integration_suite(self,fit_dict=None,sed_number=None,plot=False):
        if sed_number is not None and fit_dict is None:
            from .test_model_fit import prepare_asset
            fit_dict=prepare_asset(plot=plot,sed_number=sed_number,skip_minuit=True)
        elif fit_dict and sed_number is None:
            pass
        else:
            raise RecursionError("please provide either fit_dict or sed_number")
        self.run_ultranest(fit_dict,plot=plot)

    def test(self,plot=False,run_ultranest=True):
        from .test_model_fit import prepare_asset
        fit_dict=prepare_asset(plot=plot,sed_number=1,skip_minuit=True)
        self.run_ultranest(fit_dict=fit_dict,plot=plot)

    def test_open_mpi(self, tmp_path, monkeypatch):
        """Opt-in real MPI regression: JETSET_TEST_MPI=1 pytest -k test_open_mpi."""
        import os
        import shutil
        import subprocess
        import sys
        from pathlib import Path

        if os.environ.get('JETSET_TEST_MPI') != '1':
            pytest.skip('set JETSET_TEST_MPI=1 to run the MPI integration test')
        pytest.importorskip('ultranest')
        pytest.importorskip('mpi4py.MPI')
        launcher = shutil.which('mpirun')
        if launcher is None:
            pytest.skip('mpirun is not installed')

        # Fail an explicitly enabled test if ranks cannot communicate.
        probe = subprocess.run(
            [launcher, '-np', '2', sys.executable, '-c',
             'from mpi4py import MPI; '
             'assert MPI.COMM_WORLD.Get_size() == 2'],
            capture_output=True, text=True, timeout=30,
        )
        assert probe.returncode == 0, probe.stdout + probe.stderr

        from jetset import mcmc_ultranest
        from .test_model_fit import prepare_asset

        monkeypatch.chdir(tmp_path)
        fit_dict = prepare_asset(plot=False, sed_number=1, skip_minuit=True)
        sampler = mcmc_ultranest.UltraNestSampler(fit_dict['model_minimizer'])
        sampler.model.freeze_all()
        sampler.model.parameters.par_array[0].free()
        sampler.set_bounds(bound=5.0, bound_rel=True)
        ndim = len(sampler._par_array_sampler)
        real_run = subprocess.run
        run_dirs = []

        def run_with_instrumentation(command, **kwargs):
            run_dir = Path(kwargs['cwd'])
            run_dirs.append(run_dir)
            script_path = run_dir / 'run_jetset_ultranest.py'
            # Instrument actual worker saves, without changing their rank guard.
            # Each rank has its own marker, so recording cannot hide a write race.
            instrumentation = r"""
from pathlib import Path
import numpy as np
from mpi4py import MPI
from jetset.mcmc_ultranest import UltraNestSampler
assert MPI.COMM_WORLD.Get_size() == 2
np.random.seed(42 + MPI.COMM_WORLD.Get_rank())
_original_save = UltraNestSampler.save
_original_run_sampler = UltraNestSampler.run_sampler

def _record_save(self, name):
    marker = Path('save-rank-%d.txt' % MPI.COMM_WORLD.Get_rank())
    with marker.open('a') as stream:
        stream.write(str(name) + '\n')
    return _original_save(self, name)

def _bounded_run(self, **kwargs):
    kwargs.update(max_ncalls=256, show_status=False)
    return _original_run_sampler(self, **kwargs)

UltraNestSampler.save = _record_save
UltraNestSampler.run_sampler = _bounded_run
"""
            script_path.write_text(instrumentation + script_path.read_text())
            kwargs.update(timeout=180, capture_output=True, text=True)
            try:
                return real_run(command, **kwargs)
            except subprocess.CalledProcessError as error:
                pytest.fail('MPI worker failed:\n' + error.stdout + error.stderr)

        monkeypatch.setattr(mcmc_ultranest.subprocess, 'run', run_with_instrumentation)
        environment_before = dict(os.environ)
        first_run_files = None
        for _ in range(2):
            result = mcmc_ultranest.run_open_mpi(
                sampler, n_proc=2, min_num_live_points=64, nsteps=1,
                max_num_improvement_loops=0, min_ess=10,
            )
            run_dir = run_dirs[-1]
            assert sorted(p.name for p in run_dir.glob('save-rank-*.txt')) == [
                'save-rank-0.txt'
            ]
            assert (run_dir / 'save-rank-0.txt').read_text().splitlines() == [
                'sampler.pkl'
            ]
            input_sampler = mcmc_ultranest.UltraNestSampler.load(
                str(run_dir / 'sampler_input.pkl'))
            assert input_sampler.result is None
            reloaded = mcmc_ultranest.UltraNestSampler.load(str(run_dir / 'sampler.pkl'))
            assert result.samples.ndim == 2
            assert result.samples.shape[0] > 0
            assert result.samples.shape[1] == ndim
            assert np.all(np.isfinite(result.samples))
            assert np.all(np.isfinite(result.samples_log_prob))
            assert np.isfinite(result.logz)
            assert np.isfinite(result.logzerr)
            np.testing.assert_array_equal(reloaded.samples, result.samples)
            assert reloaded.logz == result.logz
            if first_run_files is None:
                first_run_files = {
                    path: path.read_bytes() for path in (
                        run_dir / 'sampler_input.pkl', run_dir / 'sampler.pkl')
                }
            else:
                assert run_dirs[0] != run_dir
                for path, contents in first_run_files.items():
                    assert path.read_bytes() == contents
        assert dict(os.environ) == environment_before

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
