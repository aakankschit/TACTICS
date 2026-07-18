"""
Parallel evaluation wrapper for Thompson Sampling.

Handles multi-process evaluation of compound batches using multiprocessing.Pool.

Why this module is not a plain ``pool.map(sampler.evaluate, ...)``
------------------------------------------------------------------
``pool.map`` serialises the task function and its arguments with ``pickle`` to
send them down a pipe to each worker. ``sampler.evaluate`` is a bound method,
so pickling it means pickling the sampler, which owns the evaluator. For
OpenEye-backed evaluators the evaluator holds a SWIG-wrapped C++ object
(``OEDock`` for Fred, the shape engine for ROCS), and SWIG objects raise
``TypeError: cannot pickle 'SwigPyObject' object`` on any pickle attempt.
The result was that ``processes > 1`` crashed on the very first batch for
exactly the slow evaluators that need parallelism most.

The fix is to never send an evaluator across the pipe. Each worker builds its
*own* evaluator instead, once, from the evaluator's picklable Pydantic config
(a receptor file path plus options). See :func:`_init_worker`.
"""
import multiprocessing
from typing import Any, Callable, List, Optional, Tuple

# Worker-process globals. Populated once per worker by _init_worker() and
# reused for every task that worker handles. They only ever exist inside a
# worker process; in the parent they stay None.
_worker_sampler = None


def _init_worker(sampler, evaluator_config):
    """Pool initializer: rebuild this worker's own evaluator.

    Runs exactly once per worker process, at pool creation. ``sampler`` arrives
    without its evaluator (see ``ThompsonSampler.__getstate__``), and this
    function constructs a fresh evaluator here, inside the worker, from the
    picklable ``evaluator_config``. The evaluator -- and any SWIG object it
    owns, such as ``OEDock`` -- is therefore created in the worker and never
    crosses the pickle boundary.

    Building it in the initializer rather than per task also means the design
    unit is read and the docking engine constructed once per worker, not once
    per molecule.

    Args:
        sampler: A ThompsonSampler pickled without its evaluator, pool, or logger.
        evaluator_config: Picklable evaluator config (a Pydantic model).
    """
    global _worker_sampler

    # Imported here rather than at module scope to avoid a circular import
    # (factories -> evaluators -> ... -> sampler -> parallel_evaluator).
    from ..factories import create_evaluator

    sampler.set_evaluator(create_evaluator(evaluator_config))
    _worker_sampler = sampler


def _worker_evaluate(choice_list):
    """Top-level task function, therefore picklable.

    Only the lightweight ``choice_list`` (a list of reagent indices) is sent
    over the pipe. The sampler and evaluator are already resident in the
    worker, installed by :func:`_init_worker`.
    """
    return _worker_sampler.evaluate(choice_list)


class ParallelEvaluator:
    """
    Wrapper for parallel evaluation of compounds with persistent pool.

    The pool is created once and reused across all evaluations to avoid
    the overhead of repeatedly spawning processes.

    Two execution paths:

    - ``processes == 1``: sequential, calling ``evaluate_fn`` directly. No
      pickling and no worker processes.
    - ``processes > 1``: a persistent pool whose workers each construct their
      own evaluator via :func:`_init_worker`. Requires a worker context to
      have been bound with :meth:`bind_worker_context`; without one, falls
      back to sending ``evaluate_fn`` itself, which works for picklable
      evaluators but not for OpenEye-backed ones.
    """

    def __init__(self, processes: int = 1):
        """
        Initialize parallel evaluator.

        Args:
            processes: Number of CPU cores to use. If 1, uses sequential evaluation.
        """
        self.processes = processes
        self._pool: Optional[Any] = None
        self._worker_context: Optional[Tuple[Any, Any]] = None

    def bind_worker_context(self, sampler, evaluator_config) -> None:
        """Record what each worker needs in order to rebuild its own evaluator.

        Called by ``ThompsonSampler.set_evaluator`` when the evaluator's config
        is known (i.e. when the sampler was built by ``from_config``, or when
        ``set_evaluator`` was given an explicit ``evaluator_config``).

        Any existing pool is shut down, because its workers were initialized
        from the previous context and would otherwise keep using a stale
        evaluator.

        Args:
            sampler: The ThompsonSampler whose ``evaluate`` the workers run.
            evaluator_config: Picklable evaluator config used to rebuild the
                evaluator inside each worker. If None, the context is cleared.
        """
        self.close()
        self._worker_context = (
            None if evaluator_config is None else (sampler, evaluator_config)
        )

    def _ensure_pool(self):
        """Create the process pool if it doesn't exist."""
        if self.processes > 1 and self._pool is None:
            if self._worker_context is not None:
                self._pool = multiprocessing.Pool(
                    self.processes,
                    initializer=_init_worker,
                    initargs=self._worker_context,
                )
            else:
                # No config to rebuild from; workers get whatever the parent
                # sends them. Fine for picklable evaluators, fatal for SWIG.
                self._pool = multiprocessing.Pool(self.processes)

    def close(self):
        """Close the process pool if it exists."""
        if self._pool is not None:
            self._pool.close()
            self._pool.join()
            self._pool = None

    def evaluate_batch(self,
                      evaluate_fn: Callable,
                      compound_list: List[Any]) -> List[Any]:
        """
        Evaluate a batch of compounds in parallel.

        Args:
            evaluate_fn: Function that takes a compound representation (e.g. list of reagent indices)
                        and returns [score, smiles, name]. Used directly when running
                        sequentially, and as the fallback task function when no worker
                        context has been bound.
            compound_list: List of compound representations to evaluate

        Returns:
            List of evaluation results [score, smiles, name] for each compound

        Raises:
            TypeError: if ``processes > 1`` with no worker context bound and the
                evaluator cannot be pickled. The message explains how to fix it.
        """
        if not compound_list:
            return []

        # Sequential evaluation for single process
        if self.processes <= 1:
            return [evaluate_fn(compound) for compound in compound_list]

        self._ensure_pool()

        # Preferred path: workers already hold their own evaluator.
        if self._worker_context is not None:
            return self._pool.map(_worker_evaluate, compound_list)

        # Fallback path: ship the bound method. Works for picklable evaluators
        # (Lookup, FP, MW, DB); fails loudly and usefully for SWIG-backed ones.
        try:
            return self._pool.map(evaluate_fn, compound_list)
        except TypeError as exc:
            if "pickle" not in str(exc).lower():
                raise
            raise TypeError(
                "Cannot run this evaluator with processes > 1: it holds an "
                "object that cannot be pickled (OpenEye OEDock and ROCS shape "
                "objects are SWIG-wrapped C++ objects and are never "
                "picklable).\n\n"
                "Fix: build the sampler with ThompsonSampler.from_config(config), "
                "or pass the evaluator's config explicitly via "
                "sampler.set_evaluator(evaluator, evaluator_config=...). Either "
                "lets each worker construct its own evaluator instead of "
                "receiving one over the pipe."
            ) from exc

    def __del__(self):
        """Cleanup: close pool when evaluator is garbage collected."""
        self.close()
