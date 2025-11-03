import random
from typing import List, Optional, Tuple, TYPE_CHECKING
import math
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm.auto import tqdm

from ..strategies.base_strategy import SelectionStrategy
from ..legacy.disallow_tracker import DisallowTracker
from .reagent import Reagent
from ..utils.ts_logger import get_logger
from ..utils.ts_utils import read_reagents
from .evaluators import DBEvaluator, LookupEvaluator
from .parallel_evaluator import ParallelEvaluator
from ..warmup import WarmupStrategy, StandardWarmup

if TYPE_CHECKING:
    from ..config import ThompsonSamplingConfig


class ThompsonSampler:
    """
    Unified Thompson Sampler that accepts any selection strategy.

    Parameters:
    -----------
    selection_strategy : SelectionStrategy
        The selection strategy to use (GreedySelection, RouletteWheelSelection, etc.)

    batch_size : int, default=1
        Number of compounds to SAMPLE per cycle from the strategy.
        - batch_size=1: Sample one compound per cycle (standard Thompson Sampling)
        - batch_size>1: Sample multiple compounds per cycle (batch Thompson Sampling)
        Note: This is independent of parallel evaluation settings.

    processes : int, default=1
        Number of CPU cores to use for parallel evaluation.
        - processes=1: Sequential evaluation (no multiprocessing overhead)
        - processes>1: Parallel evaluation using multiprocessing.Pool
        Recommendation: Use processes=1 for fast evaluators (LookupEvaluator, DBEvaluator)
        and processes>1 for slow evaluators (ROCSEvaluator, FredEvaluator, ML models).

    min_cpds_per_core : int, default=10
        Minimum compounds to accumulate per CPU core before triggering parallel evaluation.
        Evaluation threshold = processes * min_cpds_per_core.
        - Higher values: Less frequent evaluation, lower overhead, but more memory
        - Lower values: More frequent evaluation, higher overhead, but less memory
        Example: processes=4, min_cpds_per_core=10 → evaluate every 40 compounds
    """

    def __init__(self,
                 selection_strategy: SelectionStrategy,
                 warmup_strategy: WarmupStrategy = None,
                 log_filename: str = None,
                 batch_size: int = 1,
                 max_resamples: int = None,
                 processes: int = 1,
                 min_cpds_per_core: int = 10):
        self.selection_strategy = selection_strategy
        self.warmup_strategy = warmup_strategy or StandardWarmup()
        self.reagent_lists = []
        self.reaction = None
        self.evaluator = None
        self.logger = get_logger(__name__, filename=log_filename)
        self._disallow_tracker = None
        self.batch_size = batch_size
        self.max_resamples = max_resamples
        self.hide_progress = False
        self.num_prods = 0
        self.processes = processes
        self.min_cpds_per_core = min_cpds_per_core
        self.parallel_evaluator = ParallelEvaluator(processes=processes)

        # Log multiprocessing configuration
        if self.processes > 1:
            self.logger.info(f"Multiprocessing enabled: {self.processes} processes, "
                           f"min_cpds_per_core={self.min_cpds_per_core}, "
                           f"batch_threshold={self.processes * self.min_cpds_per_core}")

    @classmethod
    def from_config(cls, config: 'ThompsonSamplingConfig') -> 'ThompsonSampler':
        """
        Create a ThompsonSampler from a Pydantic configuration.

        This factory method supports both legacy and modern config formats.

        Args:
            config: ThompsonSamplingConfig with either:
                - Modern: strategy_config, warmup_config, evaluator_config
                - Legacy: selection_strategy, evaluator_class_name, etc.

        Returns:
            ThompsonSampler: Configured sampler instance

        Example (Modern):
            >>> from TACTICS.thompson_sampling import ThompsonSamplingConfig
            >>> from TACTICS.thompson_sampling.strategies.config import RouletteWheelConfig
            >>> from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig
            >>>
            >>> config = ThompsonSamplingConfig(
            ...     reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
            ...     reagent_file_list=["acids.smi", "amines.smi"],
            ...     num_ts_iterations=1000,
            ...     strategy_config=RouletteWheelConfig(mode="maximize", alpha=0.1),
            ...     evaluator_config=LookupEvaluatorConfig(ref_filename="scores.csv")
            ... )
            >>> sampler = ThompsonSampler.from_config(config)
        """
        from ..factories import create_strategy, create_warmup, create_evaluator
        from ..config import ThompsonSamplingConfig

        # Determine if using modern or legacy config
        use_modern = config.strategy_config is not None

        if use_modern:
            # Modern approach: Use factories
            strategy = create_strategy(config.strategy_config)
            warmup = create_warmup(config.warmup_config) if config.warmup_config else StandardWarmup()
            evaluator = create_evaluator(config.evaluator_config)
        else:
            # Legacy approach: Import and instantiate from strings
            import importlib
            import json
            from ..strategies import GreedySelection, RouletteWheelSelection, UCBSelection, EpsilonGreedySelection

            # Create strategy from string
            if config.selection_strategy == "greedy":
                strategy = GreedySelection(mode=config.mode)
            elif config.selection_strategy == "roulette_wheel":
                params = config.strategy_params or {}
                strategy = RouletteWheelSelection(mode=config.mode, **params)
            elif config.selection_strategy == "ucb":
                params = config.strategy_params or {}
                strategy = UCBSelection(mode=config.mode, **params)
            elif config.selection_strategy == "epsilon_greedy":
                params = config.strategy_params or {}
                strategy = EpsilonGreedySelection(mode=config.mode, **params)
            else:
                raise ValueError(f"Unknown selection_strategy: {config.selection_strategy}")

            # Create evaluator from class name
            module = importlib.import_module("TACTICS.thompson_sampling.core.evaluators")
            evaluator_class = getattr(module, config.evaluator_class_name)
            evaluator_arg = config.evaluator_arg
            if isinstance(evaluator_arg, dict):
                evaluator_arg = json.dumps(evaluator_arg)
            evaluator = evaluator_class(evaluator_arg)

            # Default warmup
            warmup = StandardWarmup()

        # Create sampler instance
        sampler = cls(
            selection_strategy=strategy,
            warmup_strategy=warmup,
            log_filename=config.log_filename,
            batch_size=config.batch_size,
            max_resamples=config.max_resamples,
            processes=config.processes,
            min_cpds_per_core=config.min_cpds_per_core
        )

        # Set up sampler
        sampler.set_evaluator(evaluator)
        sampler.read_reagents(config.reagent_file_list)
        sampler.set_hide_progress(config.hide_progress)

        return sampler

    def set_hide_progress(self, hide_progress: bool) -> None:
        """Hide the progress bars"""
        self.hide_progress = hide_progress

    def close(self) -> None:
        """
        Close the parallel evaluator and clean up resources.

        Call this when done with the sampler to properly shut down
        the multiprocessing pool.
        """
        if self.parallel_evaluator:
            self.parallel_evaluator.close()

    def __del__(self):
        """Cleanup: close parallel evaluator when sampler is garbage collected."""
        self.close()

    def read_reagents(self, reagent_file_list, num_to_select: Optional[int] = None):
        """Read reagents from file list"""
        self.reagent_lists = read_reagents(reagent_file_list, num_to_select=num_to_select)
        self.num_prods = math.prod([len(x) for x in self.reagent_lists])
        self.logger.info(f"{self.num_prods:.2e} possible products")
        self._disallow_tracker = DisallowTracker([len(x) for x in self.reagent_lists])

    def get_num_prods(self) -> int:
        """Get the total number of possible products"""
        return self.num_prods

    def set_evaluator(self, evaluator):
        """
        Define the evaluator.

        Automatically disables multiprocessing for fast evaluators (LookupEvaluator, DBEvaluator)
        where pickle overhead exceeds evaluation time.
        """
        self.evaluator = evaluator

        # Auto-detect fast evaluators and warn about multiprocessing inefficiency
        if self.processes > 1:
            fast_evaluators = (LookupEvaluator, DBEvaluator)
            if isinstance(evaluator, fast_evaluators):
                evaluator_name = type(evaluator).__name__
                self.logger.warning(
                    f"⚠️  Multiprocessing with {evaluator_name} may be slower than sequential! "
                    f"These evaluators perform fast lookups where pickle overhead >> evaluation time. "
                    f"Consider setting processes=1 for better performance."
                )

    def set_reaction(self, rxn_smarts):
        """Define the reaction"""
        self.reaction = AllChem.ReactionFromSmarts(rxn_smarts)

    def evaluate(self, choice_list: List[int]) -> Tuple[str, str, float]:
        """
        Evaluate a single set of reagents.

        NOTE: This method does NOT update reagent scores. Score updates must be done
        by the caller after evaluation to ensure compatibility with multiprocessing.

        Args:
            choice_list: List of reagent indices for each component

        Returns:
            Tuple of (product_smiles, product_name, score)
        """
        selected_reagents = []
        for idx, choice in enumerate(choice_list):
            component_reagent_list = self.reagent_lists[idx]
            selected_reagents.append(component_reagent_list[choice])
        prod = self.reaction.RunReactants([reagent.mol for reagent in selected_reagents])
        product_name = "_".join([reagent.reagent_name for reagent in selected_reagents])
        res = np.nan
        product_smiles = "FAIL"
        if prod:
            prod_mol = prod[0][0]
            Chem.SanitizeMol(prod_mol)
            product_smiles = Chem.MolToSmiles(prod_mol, isomericSmiles=True)
            if isinstance(self.evaluator, DBEvaluator):
                res = self.evaluator.evaluate(product_name)
                res = float(res)
            elif isinstance(self.evaluator, LookupEvaluator):
                res = self.evaluator.evaluate(product_name)
            else:
                res = self.evaluator.evaluate(prod_mol)
        return product_smiles, product_name, res

    def evaluate_batch(self, choice_lists: List[List[int]]) -> List[Tuple[str, str, float]]:
        """
        Evaluate a batch of reagent combinations in parallel.

        Args:
            choice_lists: List of choice_lists, where each choice_list is reagent indices

        Returns:
            List of tuples (product_smiles, product_name, score) for each combination
        """
        return self.parallel_evaluator.evaluate_batch(self.evaluate, choice_lists)

    def warm_up(self, num_warmup_trials=3):
        """
        Warm-up phase using configured warmup strategy.

        The warmup strategy determines how reagent combinations are generated
        to initialize reagent posteriors before the main search begins.

        Args:
            num_warmup_trials: Number of trials per reagent

        Returns:
            pl.DataFrame: Warmup results with columns ["score", "SMILES", "Name"]
        """
        warmup_results = []

        # Log warmup strategy information
        strategy_name = self.warmup_strategy.get_name()
        expected_evals = self.warmup_strategy.get_expected_evaluations(
            self.reagent_lists,
            num_warmup_trials
        )
        self.logger.info(
            f"Warmup strategy: {strategy_name}, "
            f"num_trials={num_warmup_trials}, "
            f"expected_evaluations={expected_evals}"
        )

        # Generate warmup combinations using strategy
        warmup_combinations = self.warmup_strategy.generate_warmup_combinations(
            self.reagent_lists,
            num_warmup_trials,
            self._disallow_tracker
        )

        self.logger.info(f"Generated {len(warmup_combinations)} warmup combinations")

        # Evaluate all warmup combinations (in parallel if processes > 1)
        if warmup_combinations:
            if self.processes > 1:
                self.logger.info(
                    f"Evaluating warmup combinations using {self.processes} processes..."
                )

            results = self.evaluate_batch(warmup_combinations)

            # Update reagent scores in main process after parallel evaluation
            for combination, (product_smiles, product_name, score) in zip(warmup_combinations, results):
                if np.isfinite(score):
                    warmup_results.append([score, product_smiles, product_name])
                    # Add scores to reagents used in this combination
                    for component_idx, reagent_idx in enumerate(combination):
                        self.reagent_lists[component_idx][reagent_idx].add_score(score)

        # Calculate warmup statistics
        warmup_scores = [ws[0] for ws in warmup_results]

        if not warmup_scores:
            raise RuntimeError("No valid warmup evaluations! Cannot initialize priors.")

        self.logger.info(
            f"Warmup score stats: "
            f"cnt={len(warmup_scores)}, "
            f"mean={np.mean(warmup_scores):0.4f}, "
            f"std={np.std(warmup_scores):0.4f}, "
            f"min={np.min(warmup_scores):0.4f}, "
            f"max={np.max(warmup_scores):0.4f}"
        )

        # Initialize priors for all reagents
        prior_mean = np.mean(warmup_scores)
        prior_std = np.std(warmup_scores)

        for i in range(0, len(self.reagent_lists)):
            for j in range(0, len(self.reagent_lists[i])):
                reagent = self.reagent_lists[i][j]
                try:
                    reagent.init_prior(prior_mean=prior_mean, prior_std=prior_std)
                except ValueError:
                    self.logger.info(
                        f"Skipping reagent {reagent.reagent_name} - "
                        f"no successful evaluations during warmup"
                    )
                    self._disallow_tracker.retire_one_synthon(i, j)

        # Report best score based on mode
        if self.selection_strategy.mode in ["maximize", "maximize_boltzmann"]:
            best_warmup_score = max(warmup_scores)
        else:
            best_warmup_score = min(warmup_scores)

        self.logger.info(f"Top score found during warmup: {best_warmup_score:.3f}")

        # Convert to polars DataFrame
        warmup_df = pl.DataFrame(warmup_results, schema=["score", "SMILES", "Name"], orient="row")
        return warmup_df

    def search(self, num_cycles=100):
        """
        Generic search loop that works with any strategy.

        Supports both single-compound (batch_size=1) and batch mode (batch_size>1).
        """
        if self.batch_size > 1:
            return self._search_batch(num_cycles)
        else:
            return self._search_single(num_cycles)

    def _search_single(self, num_cycles):
        """
        Single compound selection per cycle.

        Now supports thermal cycling component rotation, adaptive temperature control,
        and uniqueness tracking for RouletteWheelSelection strategy.

        Returns:
            pl.DataFrame: Search results with columns ["score", "SMILES", "Name"]
        """
        out_list = []
        rng = np.random.default_rng()
        n_components = len(self.reagent_lists)
        unique_compounds = set()
        n_resamples = 0

        for i in tqdm(range(num_cycles), desc="Search", disable=self.hide_progress):
            selected_reagents = [DisallowTracker.Empty] * len(self.reagent_lists)

            # Select reagents for each component
            for component_idx in random.sample(range(len(self.reagent_lists)),
                                              len(self.reagent_lists)):
                reagent_list = self.reagent_lists[component_idx]
                selected_reagents[component_idx] = DisallowTracker.To_Fill

                disallow_mask = self._disallow_tracker.get_disallowed_selection_mask(selected_reagents)

                # Use strategy to select reagent
                selected_idx = self.selection_strategy.select_reagent(
                    reagent_list=reagent_list,
                    disallow_mask=disallow_mask,
                    rng=rng,
                    component_idx=component_idx,
                    iteration=i
                )
                selected_reagents[component_idx] = selected_idx

            self._disallow_tracker.update(selected_reagents)

            # Check for uniqueness (for RouletteWheel strategy)
            comb_key = '_'.join(str(idx) for idx in selected_reagents)
            is_unique = comb_key not in unique_compounds

            if is_unique:
                unique_compounds.add(comb_key)
                n_resamples = 0
            else:
                n_resamples += 1

            # Update adaptive temperature for RouletteWheel strategy
            # Use batch_size=1 for single mode
            if hasattr(self.selection_strategy, 'update_temperature'):
                self.selection_strategy.update_temperature(
                    n_unique=1 if is_unique else 0,
                    batch_size=1
                )

            # Rotate thermal cycling component for RouletteWheel strategy
            if hasattr(self.selection_strategy, 'rotate_component'):
                self.selection_strategy.rotate_component(n_components)

            # Check stopping criteria (if max_resamples is set)
            if self.max_resamples and n_resamples >= self.max_resamples:
                self.logger.info(f"Stopping: {n_resamples} consecutive resamples")
                break

            # Evaluate and update (only if unique or no uniqueness tracking)
            smiles, name, score = self.evaluate(selected_reagents)
            if np.isfinite(score):
                out_list.append([score, smiles, name])
                # Update reagent posteriors
                for comp_idx, reagent_idx in enumerate(selected_reagents):
                    self.reagent_lists[comp_idx][reagent_idx].add_score(score)

            if i % 100 == 0 and out_list:
                scores = [x[0] for x in out_list]
                if self.selection_strategy.mode in ["maximize", "maximize_boltzmann"]:
                    best_score = max(scores)
                else:
                    best_score = min(scores)
                self.logger.info(f"Iteration {i}: Best score = {best_score:.3f}")

        # Convert to polars DataFrame
        search_df = pl.DataFrame(out_list, schema=["score", "SMILES", "Name"], orient="row")
        return search_df

    def _search_batch(self, num_cycles):
        """
        Batch selection mode: sample multiple compounds per cycle.

        Implements uniqueness tracking, adaptive temperature control,
        and parallel evaluation for efficiency.

        Returns:
            pl.DataFrame: Search results with columns ["score", "SMILES", "Name"]
        """
        out_list = []
        rng = np.random.default_rng()
        unique_compounds = set()
        n_resamples = 0
        n_components = len(self.reagent_lists)

        # Accumulator for compounds to evaluate in parallel
        compounds_to_evaluate = []
        min_cpds_per_batch = self.processes * self.min_cpds_per_core

        pbar = tqdm(total=num_cycles, desc="Search", disable=self.hide_progress)

        cycle = 0
        while cycle < num_cycles:
            # Build matrix of selected reagents (batch_size x n_components)
            matrix = []

            for component_idx in range(n_components):
                reagent_list = self.reagent_lists[component_idx]

                # Use batch selection if available
                if hasattr(self.selection_strategy, 'select_batch'):
                    selected_indices = self.selection_strategy.select_batch(
                        reagent_list=reagent_list,
                        batch_size=self.batch_size,
                        disallow_mask=None,  # Handled after batch generation
                        rng=rng,
                        component_idx=component_idx,
                        iteration=cycle
                    )
                else:
                    # Fallback to multiple single selections
                    selected_indices = np.array([
                        self.selection_strategy.select_reagent(
                            reagent_list=reagent_list,
                            disallow_mask=None,
                            rng=rng,
                            component_idx=component_idx,
                            iteration=cycle
                        )
                        for _ in range(self.batch_size)
                    ])

                matrix.append(selected_indices)

            # Transpose to get list of compound combinations
            combinations = np.array(matrix).T

            # Filter for unique combinations
            n_unique = 0

            for comb in combinations:
                comb_key = '_'.join(str(idx) for idx in comb)

                if comb_key not in unique_compounds:
                    compounds_to_evaluate.append(comb)
                    unique_compounds.add(comb_key)
                    n_resamples = 0
                    n_unique += 1
                else:
                    n_resamples += 1

            # Update adaptive temperature for RouletteWheel strategy
            if hasattr(self.selection_strategy, 'update_temperature'):
                self.selection_strategy.update_temperature(n_unique, self.batch_size)

            # Rotate thermal cycling component
            if hasattr(self.selection_strategy, 'rotate_component'):
                self.selection_strategy.rotate_component(n_components)

            # Check stopping criteria
            if self.max_resamples and n_resamples >= self.max_resamples:
                self.logger.info(f"Stopping: {n_resamples} consecutive resamples")
                break

            # Trigger evaluation when we have enough compounds OR at end of cycles
            should_evaluate = (
                len(compounds_to_evaluate) >= min_cpds_per_batch or
                cycle == num_cycles - 1
            )

            if should_evaluate and compounds_to_evaluate:
                # Convert to list of lists for evaluation
                choice_lists = [comb.tolist() for comb in compounds_to_evaluate]

                # Parallel evaluation
                if self.processes > 1 and cycle % 100 == 0:
                    self.logger.info(f"Evaluating batch of {len(choice_lists)} compounds "
                                   f"using {self.processes} processes")
                results = self.evaluate_batch(choice_lists)

                # Process results
                for comb, (smiles, name, score) in zip(compounds_to_evaluate, results):
                    if np.isfinite(score):
                        out_list.append([score, smiles, name])

                        # Update reagent posteriors
                        for comp_idx, reagent_idx in enumerate(comb):
                            self.reagent_lists[comp_idx][reagent_idx].add_score(score)

                # Clear accumulator
                compounds_to_evaluate = []

            # Logging
            if cycle % 100 == 0 and out_list:
                best_score = max([x[0] for x in out_list])
                self.logger.info(
                    f"Cycle {cycle}: Best={best_score:.3f}, "
                    f"Unique this batch={n_unique}/{self.batch_size}, "
                    f"Processes={self.processes}"
                )

            pbar.update(1)
            cycle += 1

        pbar.close()
        # Convert to polars DataFrame
        search_df = pl.DataFrame(out_list, schema=["score", "SMILES", "Name"], orient="row")
        return search_df
