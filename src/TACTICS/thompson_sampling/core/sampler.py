from typing import Any, Dict, List, Optional, Tuple, TYPE_CHECKING
import math
import numpy as np
import polars as pl
from rdkit import Chem
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
    from ...library_enumeration import SynthesisPipeline


class ThompsonSampler:
    """
    Unified Thompson Sampler that accepts any selection strategy.

    Parameters:
    -----------
    synthesis_pipeline : SynthesisPipeline
        The synthesis pipeline containing reaction configuration and reagent files.
        This is the single source of truth for compound generation.

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

    def __init__(
        self,
        synthesis_pipeline: "SynthesisPipeline",
        selection_strategy: SelectionStrategy,
        warmup_strategy: WarmupStrategy = None,
        log_filename: str = None,
        batch_size: int = 1,
        max_resamples: int = None,
        processes: int = 1,
        min_cpds_per_core: int = 10,
        product_library_file: Optional[str] = None,
        cats_manager=None,
        use_boltzmann_weighting: bool = False,
        seed: Optional[int] = None,
        track_diagnostics: bool = False,
    ):
        self.synthesis_pipeline = synthesis_pipeline
        self.selection_strategy = selection_strategy
        self.warmup_strategy = warmup_strategy or StandardWarmup()
        self.reagent_lists = []
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
        self.product_smiles_dict = None
        self.cats_manager = cats_manager  # Optional CATS integration
        self.use_boltzmann_weighting = use_boltzmann_weighting

        # Master RNG for reproducibility
        self._seed = seed
        self._rng = np.random.default_rng(seed)

        # Diagnostics
        self.track_diagnostics = track_diagnostics
        self._diagnostics_records: list = []

        # Load product library if provided
        if product_library_file:
            self.load_product_library(product_library_file)

        # Log multiprocessing configuration
        if self.processes > 1:
            self.logger.info(
                f"Multiprocessing enabled: {self.processes} processes, "
                f"min_cpds_per_core={self.min_cpds_per_core}, "
                f"batch_threshold={self.processes * self.min_cpds_per_core}"
            )

        # Log Boltzmann weighting status
        if self.use_boltzmann_weighting:
            self.logger.info(
                "Using Boltzmann-weighted Bayesian updates (legacy RWS algorithm)"
            )

        # Log pipeline info
        self.logger.info(
            f"Synthesis pipeline: {synthesis_pipeline.num_steps} step(s), "
            f"{synthesis_pipeline.num_components} component(s)"
        )

    @classmethod
    def from_config(cls, config: "ThompsonSamplingConfig") -> "ThompsonSampler":
        """
        Create a ThompsonSampler from a Pydantic configuration.

        Args:
            config: ThompsonSamplingConfig with synthesis_pipeline, strategy_config,
                    warmup_config, and evaluator_config

        Returns:
            ThompsonSampler: Configured sampler instance

        Example:
            >>> from TACTICS.library_enumeration import SynthesisPipeline
            >>> from TACTICS.library_enumeration.smarts_toolkit import ReactionDef, ReactionConfig
            >>> from TACTICS.thompson_sampling import ThompsonSamplingConfig
            >>> from TACTICS.thompson_sampling.strategies.config import GreedyConfig
            >>> from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig
            >>>
            >>> # Create pipeline
            >>> rxn_config = ReactionConfig(
            ...     reactions=[ReactionDef(reaction_smarts="...", step_index=0)],
            ...     reagent_file_list=["acids.smi", "amines.smi"]
            ... )
            >>> pipeline = SynthesisPipeline(rxn_config)
            >>>
            >>> # Create Thompson Sampling config
            >>> ts_config = ThompsonSamplingConfig(
            ...     synthesis_pipeline=pipeline,
            ...     num_ts_iterations=1000,
            ...     strategy_config=GreedyConfig(mode="maximize"),
            ...     evaluator_config=LookupEvaluatorConfig(ref_filename="scores.csv")
            ... )
            >>> sampler = ThompsonSampler.from_config(ts_config)
        """
        from ..factories import create_strategy, create_warmup, create_evaluator

        # Create components from config
        strategy = create_strategy(config.strategy_config)
        warmup = (
            create_warmup(config.warmup_config)
            if config.warmup_config
            else StandardWarmup()
        )
        evaluator = create_evaluator(config.evaluator_config)

        # Get pipeline from config (single source of truth)
        pipeline = config.synthesis_pipeline

        # Wire seed to warmup strategy if not already set
        seed = getattr(config, "seed", None)
        if seed is not None and hasattr(warmup, "seed") and warmup.seed is None:
            warmup.seed = seed

        # Create sampler instance
        sampler = cls(
            synthesis_pipeline=pipeline,
            selection_strategy=strategy,
            warmup_strategy=warmup,
            log_filename=config.log_filename,
            batch_size=config.batch_size,
            max_resamples=config.max_resamples,
            processes=1,  # Default to 1, compound generation is fast
            min_cpds_per_core=10,
            product_library_file=config.product_library_file,
            use_boltzmann_weighting=config.use_boltzmann_weighting,
            seed=seed,
            track_diagnostics=config.track_diagnostics,
        )

        # Set up sampler
        sampler.set_evaluator(evaluator)
        sampler.read_reagents(pipeline.reagent_file_list)

        # Auto-detect SMARTS compatibility if enabled
        if config.auto_detect_smarts_compatibility:
            pipeline.auto_detect_compatibility(
                reagent_lists=sampler.reagent_lists,
                deprotect=config.deprotect_for_compatibility,
                desalt=config.desalt_for_compatibility,
            )
            sampler.logger.info(
                "Auto-detected SMARTS compatibility for synthesis pipeline"
            )

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

    def load_product_library(self, library_file: str) -> None:
        """
        Load pre-enumerated product library for testing mode.

        When a product library is provided, the sampler will skip reaction synthesis
        and directly lookup product SMILES from the library using product codes.
        This is useful for testing on pre-enumerated libraries where synthesis is redundant.

        Args:
            library_file: Path to CSV file with 'Product_Code' and 'SMILES' columns

        Raises:
            FileNotFoundError: If library file doesn't exist
            ValueError: If required columns are missing
        """
        import os

        if not os.path.exists(library_file):
            raise FileNotFoundError(f"Product library file not found: {library_file}")

        df = pl.read_csv(library_file)

        # Check for required columns
        if "Product_Code" not in df.columns:
            raise ValueError("Product library must have 'Product_Code' column")
        if "SMILES" not in df.columns:
            raise ValueError("Product library must have 'SMILES' column")

        # Create lookup dictionary
        self.product_smiles_dict = dict(
            zip(df["Product_Code"].to_list(), df["SMILES"].to_list())
        )
        self.logger.info(
            f"Loaded pre-enumerated product library with {len(self.product_smiles_dict):,} products"
        )
        self.logger.info(
            "Product synthesis will be skipped; using pre-enumerated SMILES"
        )

    def read_reagents(self, reagent_file_list, num_to_select: Optional[int] = None):
        """Read reagents from file list with optional Boltzmann weighting and mode"""
        self.reagent_lists = read_reagents(
            reagent_file_list,
            num_to_select=num_to_select,
            use_boltzmann_weighting=self.use_boltzmann_weighting,
            mode=self.selection_strategy.mode,
        )
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

        product_name = "_".join([reagent.reagent_name for reagent in selected_reagents])
        res = np.nan
        product_smiles = "FAIL"
        prod_mol = None

        # For LookupEvaluator and DBEvaluator, skip molecule generation entirely
        # They only need product_name for lookup
        if isinstance(self.evaluator, (LookupEvaluator, DBEvaluator)):
            # Evaluate directly by product code
            if isinstance(self.evaluator, DBEvaluator):
                res = self.evaluator.evaluate(product_name)
                res = float(res)
            else:  # LookupEvaluator
                res = self.evaluator.evaluate(product_name)
            return product_smiles, product_name, res

        # For other evaluators, we need to generate the molecule
        # Try to get product molecule: first from library, then from synthesis
        use_synthesis = False

        # Check if using pre-enumerated product library
        if self.product_smiles_dict is not None:
            # Look up pre-computed SMILES from library
            product_smiles = self.product_smiles_dict.get(product_name)

            if product_smiles:
                # Convert SMILES to mol object
                prod_mol = Chem.MolFromSmiles(product_smiles)
                if prod_mol:
                    try:
                        Chem.SanitizeMol(prod_mol)
                    except Exception as e:
                        self.logger.warning(
                            f"Failed to sanitize product {product_name}: {e}"
                        )
                        prod_mol = None
                else:
                    self.logger.warning(
                        f"Failed to parse SMILES for product {product_name}"
                    )
            else:
                # Product not in library - fall back to synthesis
                self.logger.debug(
                    f"Product {product_name} not in library, falling back to synthesis"
                )
                use_synthesis = True
        else:
            # No product library - must use synthesis
            use_synthesis = True

        # Synthesize using the pipeline
        if use_synthesis:
            reagent_mols = [r.mol for r in selected_reagents]
            reagent_keys = [r.reagent_name for r in selected_reagents]

            result = self.synthesis_pipeline.enumerate_single(
                reagent_mols, reagent_keys
            )

            if result.success:
                prod_mol = result.product
                product_smiles = result.product_smiles
            else:
                # Log enumeration error for debugging
                self.logger.debug(
                    f"Enumeration failed for {product_name}: {result.error}"
                )

        # Evaluate if we have a valid molecule
        if prod_mol:
            res = self.evaluator.evaluate(prod_mol)

        return product_smiles, product_name, res

    def evaluate_batch(
        self, choice_lists: List[List[int]]
    ) -> List[Tuple[str, str, float]]:
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
            self.reagent_lists, num_warmup_trials
        )
        self.logger.info(
            f"Warmup strategy: {strategy_name}, "
            f"num_trials={num_warmup_trials}, "
            f"expected_evaluations={expected_evals}"
        )

        # Generate warmup combinations using strategy
        warmup_combinations = self.warmup_strategy.generate_warmup_combinations(
            self.reagent_lists, num_warmup_trials, self._disallow_tracker
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
            # NOTE: Scores are NOT scaled here - strategies handle mode logic themselves
            for combination, (product_smiles, product_name, score) in zip(
                warmup_combinations, results
            ):
                if np.isfinite(score):
                    warmup_results.append([score, product_smiles, product_name])
                    # Add scores to reagents WITHOUT scaling
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

        # Initialize priors for all reagents WITHOUT SCALING
        # Strategies handle mode logic themselves
        prior_mean = np.mean(warmup_scores)  # No scaling
        prior_std = np.std(warmup_scores)  # No scaling to std

        # Check if using per-reagent variance (from BalancedWarmup)
        use_per_reagent_variance = getattr(
            self.warmup_strategy, "use_per_reagent_variance", False
        )
        shrinkage_strength = getattr(self.warmup_strategy, "shrinkage_strength", 3.0)

        if use_per_reagent_variance:
            self.logger.info(
                f"Using per-reagent variance with James-Stein shrinkage "
                f"(shrinkage_strength={shrinkage_strength})"
            )

        for i in range(0, len(self.reagent_lists)):
            for j in range(0, len(self.reagent_lists[i])):
                reagent = self.reagent_lists[i][j]
                try:
                    if use_per_reagent_variance:
                        reagent.init_prior_per_reagent(
                            global_mean=prior_mean,
                            global_std=prior_std,
                            shrinkage_strength=shrinkage_strength,
                        )
                    else:
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
        warmup_df = pl.DataFrame(
            warmup_results, schema=["score", "SMILES", "Name"], orient="row"
        )
        return warmup_df

    def search(self, num_cycles=100, max_evaluations=None):
        """
        Unified search loop that works with any batch_size.

        Supports batch_size=1 (single compound per cycle) or batch_size>1 (multiple compounds per cycle).

        Args:
            num_cycles: Maximum number of sampling cycles to run
            max_evaluations: Maximum number of unique compounds to evaluate (optional)
                            If specified, search stops after evaluating this many unique compounds

        Returns:
            pl.DataFrame: Search results with columns ["score", "SMILES", "Name"]
        """
        # Validation: warn if max_evaluations doesn't align with batch_size
        if max_evaluations is not None and max_evaluations % self.batch_size != 0:
            import warnings

            warnings.warn(
                f"max_evaluations ({max_evaluations}) is not evenly divisible by batch_size ({self.batch_size}). "
                f"This may cause slightly more evaluations than expected. "
                f"Recommend using max_evaluations as a multiple of batch_size.",
                UserWarning,
                stacklevel=2,
            )

        # Calculate total cycles for CATS progressive weighting
        if max_evaluations is not None:
            total_cycles = max_evaluations // self.batch_size
        else:
            total_cycles = num_cycles

        out_list = []
        rng = self._rng
        n_resamples = 0
        n_components = len(self.reagent_lists)

        # Accumulator for compounds to evaluate in parallel
        compounds_to_evaluate = []
        min_cpds_per_batch = self.processes * self.min_cpds_per_core

        # Use max_evaluations as progress bar total if specified
        pbar_total = max_evaluations if max_evaluations else num_cycles
        pbar = tqdm(total=pbar_total, desc="Search", disable=self.hide_progress)

        cycle = 0
        while cycle < num_cycles:
            # Check if we've reached the evaluation limit
            if max_evaluations is not None and len(out_list) >= max_evaluations:
                self.logger.info(
                    f"Reached max_evaluations limit: {max_evaluations} (evaluated {len(out_list)} compounds)"
                )
                break

            # Generate batch_size unique combinations using DisallowTracker
            combinations = []
            n_unique = 0

            for _ in range(self.batch_size):
                # Build one combination iteratively, respecting DisallowTracker
                selected_reagents = [DisallowTracker.Empty] * n_components

                # Randomize component selection order to avoid bias
                selection_order = list(range(n_components))
                rng.shuffle(selection_order)

                for component_idx in selection_order:
                    reagent_list = self.reagent_lists[component_idx]
                    selected_reagents[component_idx] = DisallowTracker.To_Fill

                    # Get disallow mask from tracker
                    disallow_mask = (
                        self._disallow_tracker.get_disallowed_selection_mask(
                            selected_reagents
                        )
                    )

                    # Select reagent with disallow constraint
                    # Pass CATS context for progressive weighting
                    selected_idx = self.selection_strategy.select_reagent(
                        reagent_list=reagent_list,
                        disallow_mask=disallow_mask,
                        rng=rng,
                        component_idx=component_idx,
                        iteration=cycle,
                        current_cycle=cycle,
                        total_cycles=total_cycles,
                    )
                    selected_reagents[component_idx] = selected_idx

                # Update DisallowTracker with this combination
                self._disallow_tracker.update(selected_reagents)
                combinations.append(selected_reagents)
                compounds_to_evaluate.append(selected_reagents)
                n_unique += 1
                n_resamples = 0

            # Rotate thermal cycling component (weighted by criticality if available)
            if hasattr(self.selection_strategy, "rotate_component_weighted"):
                criticalities = [
                    self.selection_strategy.get_component_criticality(rl) or 0.5
                    for rl in self.reagent_lists
                ]
                self.selection_strategy.rotate_component_weighted(
                    n_components, criticalities, rng=rng
                )
            elif hasattr(self.selection_strategy, "rotate_component"):
                self.selection_strategy.rotate_component(n_components)

            # Adaptive temperature: signal efficiency to strategy
            if hasattr(self.selection_strategy, "adapt_temperatures"):
                # With DisallowTracker all compounds are unique, so use
                # library coverage as a proxy for sampling efficiency.
                # As coverage grows, the strategy should increase temperatures.
                _total_evaluated = len(out_list) + len(compounds_to_evaluate)
                _coverage = _total_evaluated / max(self.num_prods, 1)
                # Map coverage to efficiency: high coverage → low efficiency
                _proxy_attempted = self.batch_size
                _proxy_unique = max(1, int(self.batch_size * (1.0 - _coverage)))
                self.selection_strategy.adapt_temperatures(_proxy_unique, _proxy_attempted)

            # Check stopping criteria
            if self.max_resamples and n_resamples >= self.max_resamples:
                self.logger.info(f"Stopping: {n_resamples} consecutive resamples")
                break

            # Trigger evaluation when we have enough compounds OR at end of cycles
            should_evaluate = (
                len(compounds_to_evaluate) >= min_cpds_per_batch
                or cycle == num_cycles - 1
            )

            if should_evaluate and compounds_to_evaluate:
                # compounds_to_evaluate already contains lists
                choice_lists = compounds_to_evaluate

                # Parallel evaluation
                if self.processes > 1 and cycle % 100 == 0:
                    self.logger.info(
                        f"Evaluating batch of {len(choice_lists)} compounds "
                        f"using {self.processes} processes"
                    )
                results = self.evaluate_batch(choice_lists)

                # Process results WITHOUT scaling
                # Strategies handle mode logic themselves
                for comb, (smiles, name, score) in zip(compounds_to_evaluate, results):
                    if np.isfinite(score):
                        out_list.append([score, smiles, name])

                        # Update reagent posteriors WITHOUT scaling
                        for comp_idx, reagent_idx in enumerate(comb):
                            self.reagent_lists[comp_idx][reagent_idx].add_score(score)

                # Collect diagnostics after posterior updates
                if self.track_diagnostics:
                    for comp_idx, reagent_list in enumerate(self.reagent_lists):
                        state = self.selection_strategy.get_component_state(
                            reagent_list, comp_idx, cycle, total_cycles
                        )
                        if state is not None:
                            self._diagnostics_records.append(state)
                        else:
                            # Fallback to legacy 3-column schema
                            crit = self.selection_strategy.get_component_criticality(
                                reagent_list
                            )
                            if crit is not None:
                                self._diagnostics_records.append({
                                    "cycle": cycle,
                                    "component_idx": comp_idx,
                                    "criticality": crit,
                                })

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

            # Update progress bar based on mode
            if max_evaluations is not None:
                # Track evaluations progress
                pbar.n = len(out_list)
                pbar.refresh()
            else:
                # Track cycles progress
                pbar.update(1)

            cycle += 1

        pbar.close()
        # Convert to polars DataFrame
        search_df = pl.DataFrame(
            out_list, schema=["score", "SMILES", "Name"], orient="row"
        )
        return search_df

    # Schema for the enhanced 18-column diagnostics DataFrame
    _ENHANCED_DIAGNOSTICS_SCHEMA = {
        "component_idx": pl.Int64,
        "current_cycle": pl.Int64,
        "total_cycles": pl.Int64,
        "criticality": pl.Float64,
        "snr": pl.Float64,
        "imbalance_strength": pl.Float64,
        "normalized_entropy": pl.Float64,
        "n_active_reagents": pl.Int64,
        "participation_ratio": pl.Float64,
        "effective_n": pl.Float64,
        "sharpening_factor": pl.Float64,
        "base_temp": pl.Float64,
        "is_heated": pl.Boolean,
        "criticality_weight": pl.Float64,
        "decay": pl.Float64,
        "cats_multiplier": pl.Float64,
        "effective_multiplier": pl.Float64,
        "final_temperature": pl.Float64,
    }

    _LEGACY_DIAGNOSTICS_SCHEMA = {
        "cycle": pl.Int64,
        "component_idx": pl.Int64,
        "criticality": pl.Float64,
    }

    def get_diagnostics(self) -> pl.DataFrame:
        """Return per-cycle diagnostics trajectory as a DataFrame.

        Requires ``track_diagnostics=True`` in config.

        For CATS-aware strategies (RouletteWheelSelection, BayesUCBSelection),
        returns the enhanced 15-column schema with full intermediate values.
        The ``current_cycle`` column serves the same role as ``cycle`` in the
        legacy schema.

        For non-CATS strategies that only support ``get_component_criticality()``,
        returns the legacy 3-column schema (cycle, component_idx, criticality).

        Returns an empty DataFrame (with the correct schema) if diagnostics
        were not tracked or the strategy doesn't support criticality.
        """
        if not self._diagnostics_records:
            return pl.DataFrame(schema=self._LEGACY_DIAGNOSTICS_SCHEMA)

        # Detect schema from first record
        first = self._diagnostics_records[0]
        if "current_cycle" in first:
            # Enhanced schema
            return pl.DataFrame(
                self._diagnostics_records,
                schema=self._ENHANCED_DIAGNOSTICS_SCHEMA,
            )
        else:
            # Legacy 3-column schema
            return pl.DataFrame(
                self._diagnostics_records,
                schema=self._LEGACY_DIAGNOSTICS_SCHEMA,
            )

    def get_posterior_landscape(self) -> pl.DataFrame:
        """Return per-reagent posterior state for all components.

        Each row contains:

        - **component_idx** – component index (0-based)
        - **reagent_name** – reagent identifier
        - **mean** – posterior mean
        - **std** – posterior standard deviation
        - **n_samples** – number of observations

        Works regardless of ``track_diagnostics`` setting.
        """
        records = []
        for comp_idx, reagent_list in enumerate(self.reagent_lists):
            for reagent in reagent_list:
                records.append({
                    "component_idx": comp_idx,
                    "reagent_name": reagent.reagent_name,
                    "mean": reagent.mean,
                    "std": reagent.std,
                    "n_samples": reagent.n_samples,
                })
        return pl.DataFrame(records)

    def get_sar_summary(self) -> Dict[str, Any]:
        """Strategy-agnostic post-hoc SAR assessment from final posteriors.

        Works for ANY strategy — computes entropy-based concentration from
        the current posterior state.  Uses the same z-score softmax approach
        as CATS ``_calculate_criticality()`` so that the post-hoc concentration
        and CATS criticality are directly comparable.

        When ``track_diagnostics=True`` and the strategy supports
        ``get_component_state()`` (CATS strategies), convergence dynamics
        are also included.

        Returns:
            Dict with keys:

            - ``components`` – list of per-component dicts
            - ``landscape_type`` – "structured_SAR" / "diffuse_SAR" /
              "mixed" / "insufficient_data"
            - ``summary_text`` – human-readable 2-3 sentence assessment
            - ``convergence_dynamics`` – (only with CATS + track_diagnostics)
              per-component convergence info
        """
        mode = self.selection_strategy.mode
        components = []

        for comp_idx, reagent_list in enumerate(self.reagent_lists):
            comp_info = self._compute_component_sar(reagent_list, comp_idx, mode)
            components.append(comp_info)

        # Classify landscape
        assessments = [c["convergence_assessment"] for c in components]
        if all(a == "insufficient_data" for a in assessments):
            landscape_type = "insufficient_data"
        elif all(a == "structured" for a in assessments):
            landscape_type = "structured_SAR"
        elif all(a == "diffuse" for a in assessments):
            landscape_type = "diffuse_SAR"
        else:
            landscape_type = "mixed"

        # Build summary text
        n_structured = sum(1 for a in assessments if a == "structured")
        n_diffuse = sum(1 for a in assessments if a == "diffuse")
        n_insufficient = sum(1 for a in assessments if a == "insufficient_data")
        n_total = len(assessments)

        if landscape_type == "structured_SAR":
            summary_text = (
                f"All {n_total} components show structured SAR — the search has "
                f"identified dominant reagents with concentrated posterior mass."
            )
        elif landscape_type == "diffuse_SAR":
            summary_text = (
                f"All {n_total} components show diffuse SAR — posterior mass is "
                f"spread across many reagents with no clear winners."
            )
        elif landscape_type == "insufficient_data":
            summary_text = (
                f"Insufficient data across all {n_total} components to assess "
                f"SAR structure. More observations are needed."
            )
        else:
            summary_text = (
                f"Mixed SAR landscape: {n_structured}/{n_total} structured, "
                f"{n_diffuse}/{n_total} diffuse"
                + (f", {n_insufficient}/{n_total} insufficient data" if n_insufficient else "")
                + ". Components differ in how concentrated the activity is."
            )

        result: Dict[str, Any] = {
            "components": components,
            "landscape_type": landscape_type,
            "summary_text": summary_text,
        }

        # Add convergence dynamics if CATS diagnostics available
        if self.track_diagnostics and self._diagnostics_records:
            first = self._diagnostics_records[0]
            if "current_cycle" in first:
                diag_df = self.get_diagnostics()
                dynamics = self._compute_convergence_dynamics(diag_df)
                result["convergence_dynamics"] = dynamics

        return result

    def _compute_component_sar(
        self, reagent_list: List, comp_idx: int, mode: str
    ) -> Dict[str, Any]:
        """Compute SAR metrics for a single component from posteriors."""
        active_reagents = [r for r in reagent_list if r.n_samples > 0]
        n_active = len(active_reagents)

        if n_active < 2:
            return {
                "component_idx": comp_idx,
                "posterior_entropy": float("nan"),
                "normalized_entropy": float("nan"),
                "concentration": float("nan"),
                "gini_coefficient": float("nan"),
                "top1_share": float("nan"),
                "top5_share": float("nan"),
                "snr": float("nan"),
                "effective_n": float("nan"),
                "participation_ratio": float("nan"),
                "dominant_reagent": None,
                "convergence_assessment": "insufficient_data",
            }

        means = np.array([r.mean for r in active_reagents])
        names = [r.reagent_name for r in active_reagents]

        if np.std(means) < 1e-10:
            return {
                "component_idx": comp_idx,
                "posterior_entropy": np.log(n_active),
                "normalized_entropy": 1.0,
                "concentration": 0.0,
                "gini_coefficient": 0.0,
                "top1_share": 1.0 / n_active,
                "top5_share": min(5, n_active) / n_active,
                "snr": 0.0,
                "effective_n": float(n_active),
                "participation_ratio": 1.0 / n_active,
                "dominant_reagent": names[0],
                "convergence_assessment": "diffuse",
            }

        if mode == "minimize":
            means = -means

        # Z-score softmax (same as CATS _calculate_criticality)
        mean_std = np.std(means)

        # SNR
        se_squared = np.array([
            r.std ** 2 / max(r.n_samples, 1) for r in active_reagents
        ])
        noise_std = np.sqrt(np.mean(se_squared))
        snr = float(mean_std / max(noise_std, 1e-10))

        imbalance_strength = 1.0 - np.exp(-max(snr - 1.0, 0.0))

        z_scores = (means - np.mean(means)) / mean_std
        z_scores *= imbalance_strength

        # Read criticality metric from strategy if available
        criticality_metric = getattr(self.selection_strategy, "criticality_metric", "ipr")
        n_adaptive_sharpening = getattr(self.selection_strategy, "n_adaptive_sharpening", True)

        if criticality_metric == "ipr" and n_adaptive_sharpening and n_active > 2:
            sharpening = max(1.0, np.sqrt(np.log(n_active)))
            z_scores *= sharpening

        exp_means = np.exp(z_scores - z_scores.max())  # numerical stability
        probabilities = exp_means / exp_means.sum()

        # Shannon entropy
        entropy = -np.sum(probabilities * np.log(probabilities + 1e-10))
        max_entropy = np.log(n_active)
        normalized_entropy = entropy / max_entropy if max_entropy > 1e-10 else 1.0

        # IPR-based concentration
        ipr = float(np.sum(probabilities ** 2))
        effective_n = 1.0 / ipr

        if criticality_metric == "ipr":
            concentration = 1.0 - (effective_n / n_active)
        else:
            concentration = 1.0 - normalized_entropy

        # Gini coefficient
        sorted_probs = np.sort(probabilities)
        n = len(sorted_probs)
        index = np.arange(1, n + 1)
        gini = float((2 * np.sum(index * sorted_probs) - (n + 1) * np.sum(sorted_probs)) / (n * np.sum(sorted_probs)))

        # Top-k shares
        sorted_desc = np.sort(probabilities)[::-1]
        top1_share = float(sorted_desc[0])
        top5_share = float(np.sum(sorted_desc[:min(5, n)]))

        # Dominant reagent
        best_idx = int(np.argmax(probabilities))
        dominant_reagent = names[best_idx]

        # Assessment
        if n_active < 5:
            assessment = "insufficient_data"
        elif concentration > 0.3:
            assessment = "structured"
        else:
            assessment = "diffuse"

        return {
            "component_idx": comp_idx,
            "posterior_entropy": float(entropy),
            "normalized_entropy": float(normalized_entropy),
            "concentration": float(concentration),
            "gini_coefficient": gini,
            "top1_share": top1_share,
            "top5_share": top5_share,
            "snr": snr,
            "effective_n": float(effective_n),
            "participation_ratio": ipr,
            "dominant_reagent": dominant_reagent,
            "convergence_assessment": assessment,
        }

    @staticmethod
    def _compute_convergence_dynamics(diag_df: pl.DataFrame) -> List[Dict[str, Any]]:
        """Extract convergence dynamics from enhanced diagnostics DataFrame."""
        dynamics = []
        for comp_idx in diag_df["component_idx"].unique().sort().to_list():
            comp_df = diag_df.filter(pl.col("component_idx") == comp_idx).sort("current_cycle")

            crits = comp_df["criticality"].to_numpy()
            cycles = comp_df["current_cycle"].to_numpy()

            # First cycle where criticality > 0.3
            structured_mask = crits > 0.3
            cycle_first = int(cycles[structured_mask][0]) if structured_mask.any() else None

            # Stable: first cycle after which criticality stays > 0.3
            cycle_stable = None
            if structured_mask.any():
                for i in range(len(crits)):
                    if np.all(crits[i:] > 0.3):
                        cycle_stable = int(cycles[i])
                        break

            # SNR trajectory slope (linear regression)
            snr_values = comp_df["snr"].to_numpy()
            valid_snr = np.isfinite(snr_values)
            if valid_snr.sum() >= 2:
                valid_cycles = cycles[valid_snr].astype(float)
                valid_snr_vals = snr_values[valid_snr]
                slope = float(np.polyfit(valid_cycles, valid_snr_vals, 1)[0])
            else:
                slope = float("nan")

            dynamics.append({
                "component_idx": comp_idx,
                "cycle_first_structured": cycle_first,
                "cycle_stable": cycle_stable,
                "snr_trajectory_slope": slope,
            })

        return dynamics
