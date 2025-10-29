from pydantic import BaseModel, Field
from typing import Literal, Optional

class ThompsonSamplingConfig(BaseModel):
    """Unified configuration for all TS variants"""

    # Core parameters
    evaluator_class_name: str
    evaluator_arg: str | dict
    reaction_smarts: str
    reagent_file_list: list[str]
    num_ts_iterations: int
    num_warmup_trials: int = 3

    # Strategy selection
    selection_strategy: Literal[
        "greedy",
        "roulette_wheel",
        "ucb",
        "epsilon_greedy",
        "boltzmann"
    ] = "greedy"

    mode: Literal["maximize", "minimize"] = "maximize"

    # Strategy-specific parameters
    strategy_params: Optional[dict] = Field(default_factory=dict)
    # For roulette_wheel: {
    #     "alpha": 0.1,
    #     "beta": 0.1,
    #     "scaling": 1.0,
    #     "alpha_increment": 0.01,
    #     "beta_increment": 0.001,
    #     "efficiency_threshold": 0.1
    # }
    # For ucb: {"c": 2.0}
    # For epsilon_greedy: {"epsilon": 0.1, "decay": 0.995}

    # Batch sampling parameters
    batch_size: int = 1  # Number of compounds to sample per cycle (1 = single mode, >1 = batch mode)
    max_resamples: Optional[int] = None  # Stop after N consecutive duplicate resamples

    # Output
    results_filename: Optional[str] = "results.csv"
    log_filename: Optional[str] = None

    # Performance / Multiprocessing
    # IMPORTANT: Set processes=1 for fast evaluators (LookupEvaluator, DBEvaluator)
    #            Multiprocessing overhead >> evaluation time for dictionary lookups!
    #            Use processes>1 ONLY for slow evaluators (ROCS, Fred, ML models)
    processes: int = 1  # Number of CPU cores for parallel evaluation
    min_cpds_per_core: int = 10  # Compounds per core before evaluation (threshold = processes * min_cpds_per_core)
    hide_progress: bool = False


class RandomBaselineConfig(BaseModel):
    """Configuration for random baseline sampling"""

    evaluator_class_name: str
    evaluator_arg: str | dict
    reaction_smarts: str
    reagent_file_list: list[str]
    num_trials: int
    num_to_save: int
    ascending_output: bool = False
    outfile_name: Optional[str] = None
    log_filename: Optional[str] = None