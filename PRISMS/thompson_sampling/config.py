from pydantic import BaseModel, Field
from typing import List, Optional, Literal, Union, Dict, Any

class StandardSamplerConfig(BaseModel):
    sampler_type: Literal["standard"]
    ts_mode: Literal["maximize", "minimize", "maximize_boltzmann", "minimize_boltzmann"]
    evaluator_class_name: str
    evaluator_arg: Union[str, Dict[str, Any]]  # Can be string or dict
    reaction_smarts: str
    num_ts_iterations: int
    reagent_file_list: List[str]
    num_warmup_trials: int
    results_filename: Optional[str] = None
    log_filename: Optional[str] = None

class EnhancedSamplerConfig(BaseModel):
    sampler_type: Literal["enhanced"]
    processes: int = Field(..., gt=0)
    scaling: float
    percent_of_library: float = Field(..., gt=0, le=1)
    minimum_no_of_compounds_per_core: int = Field(..., gt=0)
    stopping_criteria: int = Field(..., gt=0)
    evaluator_class_name: str
    evaluator_arg: Union[str, Dict[str, Any]]  # Can be string or dict
    reaction_smarts: str
    num_ts_iterations: int
    reagent_file_list: List[str]
    num_warmup_trials: int
    results_filename: Optional[str] = None
    log_filename: Optional[str] = None

class RandomBaselineConfig(BaseModel):
    """Configuration for Random Baseline sampling."""
    evaluator_class_name: str
    evaluator_arg: Union[str, Dict[str, Any]]  # Can be string or dict
    reaction_smarts: str
    reagent_file_list: List[str]
    num_trials: int = Field(..., gt=0)
    num_to_save: int = Field(..., gt=0)
    ascending_output: bool = False
    outfile_name: Optional[str] = "random_scores.csv"
    log_filename: Optional[str] = None
    results_filename: Optional[str] = None

TSConfig = Union[StandardSamplerConfig, EnhancedSamplerConfig] 