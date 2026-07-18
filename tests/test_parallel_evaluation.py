"""
Tests for parallel batch evaluation (bug report 2026-06-23).

Covers three defects that together made `processes > 1` unusable with slow
evaluators such as Fred docking:

1. `from_config()` hardcoded `processes=1`, so evaluation parallelism could
   not be enabled through the config API at all.
2. `pool.map(sampler.evaluate, ...)` pickled the sampler, and therefore the
   evaluator, and therefore any SWIG-wrapped OpenEye object it held --
   raising `TypeError: cannot pickle 'SwigPyObject' object` on the first batch.
3. Importing TACTICS loaded the OpenEye toolkits at module level, which could
   segfault the interpreter via a libexpat conflict with prompt_toolkit.

Note on start methods: multiprocessing under `spawn` (the macOS/Windows
default) re-imports the entry-point module in each worker, so these tests must
run from an importable module -- which pytest provides.
"""
import os
import pickle
import sys
import tempfile
import shutil

import numpy as np
import pytest

from TACTICS.library_enumeration import SynthesisPipeline, ReactionConfig, ReactionDef
from TACTICS.thompson_sampling import ThompsonSamplingConfig, ThompsonSampler
from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig
from TACTICS.thompson_sampling.core.parallel_evaluator import ParallelEvaluator
from TACTICS.thompson_sampling.strategies.config import GreedyConfig

REACTION_SMARTS = "[C:1](=O)[OH].[N:2]>>[C:1](=O)[N:2]"
ALL_COMBOS = [[0, 0], [0, 1], [1, 0], [1, 1]]


class _Unpicklable:
    """Stands in for a SWIG-wrapped OpenEye object (OEDock, ROCS shape)."""

    def __reduce__(self):
        raise TypeError("cannot pickle 'SwigPyObject' object")


class TestParallelEvaluationConfig:
    """Bug 1: `processes` must be reachable through the config API."""

    def setup_method(self):
        self.temp_dir = tempfile.mkdtemp()
        self.config = _make_config(self.temp_dir, processes=1)

    def teardown_method(self):
        shutil.rmtree(self.temp_dir)

    def test_config_exposes_processes_with_safe_default(self):
        assert self.config.processes == 1
        assert self.config.min_cpds_per_core == 10

    def test_processes_must_be_positive(self):
        with pytest.raises(Exception):
            _make_config(self.temp_dir, processes=0)

    def test_from_config_propagates_processes(self):
        """The regression: from_config used to hardcode processes=1."""
        config = _make_config(self.temp_dir, processes=4)
        sampler = ThompsonSampler.from_config(config)

        assert sampler.processes == 4, (
            "from_config ignored config.processes -- evaluation would run "
            "single-threaded regardless of the allocated cores"
        )
        sampler.close()

    def test_from_config_propagates_min_cpds_per_core(self):
        config = _make_config(self.temp_dir, processes=2, min_cpds_per_core=3)
        sampler = ThompsonSampler.from_config(config)

        assert sampler.min_cpds_per_core == 3
        sampler.close()

    def test_from_config_binds_evaluator_config_for_workers(self):
        """Workers need the config to rebuild their own evaluator."""
        config = _make_config(self.temp_dir, processes=2)
        sampler = ThompsonSampler.from_config(config)

        assert sampler.evaluator_config is not None
        assert sampler.parallel_evaluator._worker_context is not None
        sampler.close()


class TestSamplerPickling:
    """Bug 2, part 1: the sampler must survive being sent to a worker."""

    def setup_method(self):
        self.temp_dir = tempfile.mkdtemp()

    def teardown_method(self):
        shutil.rmtree(self.temp_dir)

    def test_sampler_pickles_without_evaluator_pool_or_logger(self):
        config = _make_config(self.temp_dir, processes=2)
        sampler = ThompsonSampler.from_config(config)

        restored = pickle.loads(pickle.dumps(sampler))

        # Dropped on the way out...
        assert restored.evaluator is None, "evaluator must not be pickled"
        assert restored.parallel_evaluator is None, "pool must not be pickled"
        # ...and the logger is rebuilt on the way in.
        assert restored.logger is not None
        # Search state that workers do need is preserved.
        assert len(restored.reagent_lists) == len(sampler.reagent_lists)
        sampler.close()

    def test_sampler_pickles_even_with_unpicklable_evaluator(self):
        """A SWIG-holding evaluator must not block sending the sampler."""
        config = _make_config(self.temp_dir, processes=2)
        sampler = ThompsonSampler.from_config(config)
        sampler.evaluator.swig_handle = _Unpicklable()

        # Would raise TypeError if __getstate__ did not drop the evaluator.
        restored = pickle.loads(pickle.dumps(sampler))

        assert restored.evaluator is None
        sampler.close()


class TestParallelEvaluationCorrectness:
    """Bug 2, part 2: parallel evaluation must run, and agree with sequential."""

    def setup_method(self):
        self.temp_dir = tempfile.mkdtemp()

    def teardown_method(self):
        shutil.rmtree(self.temp_dir)

    def _evaluate(self, processes):
        config = _make_config(self.temp_dir, processes=processes,
                              min_cpds_per_core=1)
        sampler = ThompsonSampler.from_config(config)
        try:
            return sampler.parallel_evaluator.evaluate_batch(
                sampler.evaluate, ALL_COMBOS)
        finally:
            sampler.close()

    def test_parallel_matches_sequential(self):
        sequential = self._evaluate(processes=1)
        parallel = self._evaluate(processes=2)

        assert len(parallel) == len(ALL_COMBOS)
        assert parallel == sequential, (
            "parallel evaluation returned different results than sequential"
        )

    def test_empty_batch_short_circuits(self):
        evaluator = ParallelEvaluator(processes=2)
        # Must not create a pool for an empty batch.
        assert evaluator.evaluate_batch(lambda c: c, []) == []
        assert evaluator._pool is None


class TestWorkerRebuildsOwnEvaluator:
    """Bug 2, part 3: the evaluator is constructed in the worker, not shipped.

    This is what makes OpenEye-backed evaluators usable with processes > 1:
    OEDock never crosses the pickle boundary, so its unpicklability is
    irrelevant. These tests use a plain unpicklable stand-in so they run
    without an OpenEye licence.
    """

    def setup_method(self):
        self.temp_dir = tempfile.mkdtemp()

    def teardown_method(self):
        shutil.rmtree(self.temp_dir)

    def test_unpicklable_evaluator_still_evaluates_in_parallel(self):
        config = _make_config(self.temp_dir, processes=2, min_cpds_per_core=1)
        sampler = ThompsonSampler.from_config(config)

        # Attach a SWIG stand-in. Under the old implementation pool.map would
        # pickle the bound method -> sampler -> evaluator -> this object, and
        # raise TypeError on the very first batch.
        sampler.evaluator.swig_handle = _Unpicklable()

        results = sampler.parallel_evaluator.evaluate_batch(
            sampler.evaluate, ALL_COMBOS)

        assert len(results) == len(ALL_COMBOS)
        assert all(np.isfinite(score) for _, _, score in results)
        sampler.close()

    def test_binding_a_context_resets_an_existing_pool(self):
        """Rebinding must not leave workers holding a stale evaluator."""
        evaluator = ParallelEvaluator(processes=2)
        evaluator._pool = _DummyPool()
        dummy = evaluator._pool

        evaluator.bind_worker_context(sampler=None, evaluator_config=None)

        assert dummy.closed, "existing pool should be shut down on rebind"
        assert evaluator._pool is None
        assert evaluator._worker_context is None


class TestOpenEyeImportIsDeferred:
    """Bug 3: importing TACTICS must not load the OpenEye toolkits.

    Loading OpenEye at module scope initialises libexpat; prompt_toolkit
    (via tqdm) later re-enters expat from pyexpat at module scope, and the
    conflict segfaults the interpreter during `import TACTICS`.
    """

    def test_importing_evaluators_does_not_import_openeye(self):
        # Importing the module must leave the OpenEye globals unbound. This
        # holds whether or not openeye is installed.
        import TACTICS.thompson_sampling.core.evaluators as evaluators

        source = open(evaluators.__file__).read()
        assert "from openeye import" not in source.split("def _ensure_openeye")[0], (
            "OpenEye must not be imported at module scope"
        )

    def test_ensure_openeye_is_idempotent_or_raises_cleanly(self):
        import TACTICS.thompson_sampling.core.evaluators as evaluators

        try:
            evaluators._ensure_openeye()
        except ImportError as exc:
            # No licence/install: the error must name the package and be actionable.
            assert "openeye" in str(exc).lower()
            return

        # Installed: modules are cached and a second call is a no-op.
        assert evaluators.oechem is not None
        first = evaluators.oechem
        evaluators._ensure_openeye()
        assert evaluators.oechem is first

    def test_ml_classifier_does_not_depend_on_openeye_import(self):
        """joblib used to be imported inside the OpenEye try/except."""
        import TACTICS.thompson_sampling.core.evaluators as evaluators

        source = open(evaluators.__file__).read()
        header = source.split("class Evaluator")[0]
        assert "import joblib" not in header, (
            "joblib must not be bound to the OpenEye import block"
        )


class _DummyPool:
    """Minimal stand-in recording that close/join were called."""

    def __init__(self):
        self.closed = False

    def close(self):
        self.closed = True

    def join(self):
        pass


def _make_config(temp_dir, processes=1, min_cpds_per_core=None):
    """Build a small 2x2 lookup-scored config in ``temp_dir``."""
    acids = os.path.join(temp_dir, "acids.smi")
    amines = os.path.join(temp_dir, "amines.smi")
    scores = os.path.join(temp_dir, "scores.csv")

    if not os.path.exists(acids):
        with open(acids, "w") as f:
            f.write("CC(=O)O\tacid1\nCCC(=O)O\tacid2\n")
        with open(amines, "w") as f:
            f.write("CN\tam1\nCCN\tam2\n")
        with open(scores, "w") as f:
            f.write("Product_Code,Scores\n")
            f.write("acid1_am1,0.1\nacid1_am2,0.2\n")
            f.write("acid2_am1,0.3\nacid2_am2,0.4\n")

    pipeline = SynthesisPipeline(ReactionConfig(
        reactions=[ReactionDef(reaction_smarts=REACTION_SMARTS, step_index=0)],
        reagent_file_list=[acids, amines],
    ))

    kwargs = dict(
        synthesis_pipeline=pipeline,
        num_ts_iterations=2,
        strategy_config=GreedyConfig(mode="maximize"),
        evaluator_config=LookupEvaluatorConfig(ref_filename=scores),
        processes=processes,
        batch_size=2,
        hide_progress=True,
    )
    if min_cpds_per_core is not None:
        kwargs["min_cpds_per_core"] = min_cpds_per_core

    return ThompsonSamplingConfig(**kwargs)
