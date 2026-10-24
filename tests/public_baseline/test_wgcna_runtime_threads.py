"""Resource configuration and real R parity, without study-data computation."""
from __future__ import annotations

import os
from pathlib import Path
import sys
import threading

import numpy as np
import pandas as pd
import pytest

from katana_cosplice import r_wgcna


@pytest.fixture
def clean_budget(monkeypatch):
    for name in ("KATANA_NETWORK_THREADS", "KATANA_CPUS", "SLURM_CPUS_PER_TASK", "PBS_NP", "NSLOTS",
                 "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setattr(os, "sched_getaffinity", lambda _: set(range(256)), raising=False)


def test_local_unconfigured_does_not_claim_all_machine_cores(clean_budget):
    assert r_wgcna.allocated_wgcna_threads() == 1


def test_scheduler_budget_caps_environment_and_explicit_request(clean_budget, monkeypatch):
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    monkeypatch.setenv("OPENBLAS_NUM_THREADS", "256")
    assert r_wgcna.allocated_wgcna_threads() == 8
    assert r_wgcna.allocated_wgcna_threads(64) == 8
    assert r_wgcna.allocated_wgcna_threads(3) == 3


def test_affinity_and_user_budget_bound_threads(clean_budget, monkeypatch):
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "64")
    monkeypatch.setenv("KATANA_NETWORK_THREADS", "16")
    monkeypatch.setattr(os, "sched_getaffinity", lambda _: {2, 3, 4, 5})
    assert r_wgcna.allocated_wgcna_threads() == 4


def test_nextflow_budget_caps_user_and_inherited_blas_budget(clean_budget, monkeypatch):
    monkeypatch.setenv("KATANA_CPUS", "4")
    monkeypatch.setenv("KATANA_NETWORK_THREADS", "64")
    monkeypatch.setenv("OPENBLAS_NUM_THREADS", "256")
    assert r_wgcna.allocated_wgcna_threads() == 4


def test_local_omp_budget_and_invalid_threads(clean_budget, monkeypatch):
    monkeypatch.setenv("OMP_NUM_THREADS", "3")
    assert r_wgcna.allocated_wgcna_threads() == 3
    for value in (0, -1, 1.5, True):
        with pytest.raises(ValueError, match="positive integer"):
            r_wgcna.allocated_wgcna_threads(value)


def test_child_wgcna_and_blas_receive_same_bounded_budget():
    environment = r_wgcna._threaded_r_env(Path(sys.executable), 4)
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "ALLOW_WGCNA_THREADS"):
        assert environment[key] == "4"
    assert environment["OMP_DYNAMIC"] == "FALSE"


def test_progress_is_streamed_before_child_exits(monkeypatch):
    progress_received = threading.Event()
    original_write = sys.stderr.write

    class Log:
        def write(self, text):
            if "started-before-exit" in text:
                progress_received.set()
            return original_write(text)

        def flush(self):
            return None

    monkeypatch.setattr(r_wgcna.sys, "stderr", Log())
    result = []
    runner = threading.Thread(target=lambda: result.append(r_wgcna._run_wgcna(
        [sys.executable, "-u", "-c",
         "import time; print('started-before-exit', flush=True); time.sleep(2); print('complete', flush=True)"],
        environment=os.environ.copy(), timeout=10)))
    runner.start()
    assert progress_received.wait(1.5), "R progress must not remain buffered until process completion"
    assert runner.is_alive()
    runner.join(10)
    assert result[0][0] == 0
    assert "complete" in result[0][1]


def test_failure_tail_is_bounded(capsys):
    code, tail = r_wgcna._run_wgcna(
        [sys.executable, "-u", "-c", "import sys; [print(i) for i in range(120)]; sys.exit(7)"],
        environment=os.environ.copy(), timeout=10)
    assert code == 7
    assert len(tail.splitlines()) == 100
    assert tail.splitlines()[0] == "20"
    assert "119" in capsys.readouterr().err


def test_real_wgcna_threads_preserve_labels_and_tom(clean_budget):
    from katana_cosplice.benchmark import _rscript_with_wgcna

    executable, diagnostic = _rscript_with_wgcna()
    if executable is None:
        pytest.skip(diagnostic)
    rng = np.random.default_rng(348)
    latent = rng.normal(size=(3, 24))
    matrix = pd.DataFrame(np.repeat(latent, 20, axis=0) + rng.normal(0, .07, (60, 24)),
                          index=[f"event-{i}" for i in range(60)])
    parameters = dict(soft_power=6, min_module_size=10, merge_cut_height=.25,
                      deep_split=2, pam_respects_dendro=False, rscript=executable)
    one = r_wgcna.exact_wgcna_labels(matrix, threads=1, **parameters)
    two = r_wgcna.exact_wgcna_labels(matrix, threads=2, **parameters)
    assert np.array_equal(one, two)
    one_labels, one_tom = r_wgcna.exact_wgcna_labels(matrix, threads=1, return_tom=True, **parameters)
    two_labels, two_tom = r_wgcna.exact_wgcna_labels(matrix, threads=2, return_tom=True, **parameters)
    assert np.array_equal(one_labels, two_labels)
    np.testing.assert_allclose(one_tom, two_tom, rtol=1e-12, atol=1e-12)
