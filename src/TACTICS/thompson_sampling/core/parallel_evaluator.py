"""
Parallel evaluation wrapper for Thompson Sampling.

Handles multi-process evaluation of compound batches using multiprocessing.Pool.
"""
from multiprocessing import Pool
from typing import List, Callable, Any, Optional


class ParallelEvaluator:
    """
    Wrapper for parallel evaluation of compounds with persistent pool.

    The pool is created once and reused across all evaluations to avoid
    the overhead of repeatedly spawning processes.
    """

    def __init__(self, processes: int = 1):
        """
        Initialize parallel evaluator.

        Args:
            processes: Number of CPU cores to use. If 1, uses sequential evaluation.
        """
        self.processes = processes
        self._pool: Optional[Pool] = None

    def _ensure_pool(self):
        """Create the process pool if it doesn't exist."""
        if self.processes > 1 and self._pool is None:
            self._pool = Pool(self.processes)

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
            evaluate_fn: Function that takes a compound representation (e.g., list of reagent indices)
                        and returns [score, smiles, name]
            compound_list: List of compound representations to evaluate

        Returns:
            List of evaluation results [score, smiles, name] for each compound
        """
        if not compound_list:
            return []

        # Use parallel evaluation if processes > 1
        if self.processes > 1:
            self._ensure_pool()
            results = self._pool.map(evaluate_fn, compound_list)
        else:
            # Sequential evaluation for single process
            results = [evaluate_fn(compound) for compound in compound_list]

        return results

    def __del__(self):
        """Cleanup: close pool when evaluator is garbage collected."""
        self.close()
