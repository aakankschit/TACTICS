"""
Direct tests for Pydantic config model validation.

Tests the new synthesis_pipeline-based API.
"""

import pytest
import tempfile
import os
from pydantic import ValidationError

from TACTICS.thompson_sampling import ThompsonSamplingConfig, RandomBaselineConfig
from TACTICS.thompson_sampling.strategies.config import (
    GreedyConfig,
    RouletteWheelConfig,
    UCBConfig,
    EpsilonGreedyConfig,
    BoltzmannConfig,
    BayesUCBConfig,
)
from TACTICS.thompson_sampling.core.evaluator_config import (
    LookupEvaluatorConfig,
    DBEvaluatorConfig,
)
from TACTICS.library_enumeration import SynthesisPipeline
from TACTICS.library_enumeration.smarts_toolkit import ReactionConfig, ReactionDef


@pytest.fixture
def temp_reagent_files():
    """Create temporary reagent files for testing."""
    temp_dir = tempfile.mkdtemp()
    reagent_file1 = os.path.join(temp_dir, "reagents1.smi")
    reagent_file2 = os.path.join(temp_dir, "reagents2.smi")

    with open(reagent_file1, "w") as f:
        f.write("CCO\tethanol\nCCCO\tpropanol\n")
    with open(reagent_file2, "w") as f:
        f.write("CC(=O)O\tacetic_acid\nCC(C)(C)O\ttert_butanol\n")

    yield reagent_file1, reagent_file2

    # Cleanup
    import shutil

    shutil.rmtree(temp_dir)


@pytest.fixture
def pipeline(temp_reagent_files):
    """Create a SynthesisPipeline for testing."""
    reagent_file1, reagent_file2 = temp_reagent_files
    config = ReactionConfig(
        reactions=[
            ReactionDef(
                reaction_smarts="[C:1](=O)[OH].[O:2]>>[C:1](=O)[O:2]",
                step_index=0,
            )
        ],
        reagent_file_list=[reagent_file1, reagent_file2],
    )
    return SynthesisPipeline(config)


def test_thompson_sampling_config_valid(pipeline):
    """
    Test that a valid ThompsonSamplingConfig can be created.
    """
    config = ThompsonSamplingConfig(
        synthesis_pipeline=pipeline,
        num_ts_iterations=10,
        num_warmup_trials=2,
        strategy_config=GreedyConfig(mode="maximize"),
        evaluator_config=DBEvaluatorConfig(db_filename="test.db"),
    )
    assert config.num_ts_iterations == 10
    assert isinstance(config.strategy_config, GreedyConfig)
    assert config.strategy_config.mode == "maximize"


def test_thompson_sampling_config_roulette_wheel(pipeline):
    """
    Test config with roulette_wheel strategy.
    """
    config = ThompsonSamplingConfig(
        synthesis_pipeline=pipeline,
        num_ts_iterations=10,
        num_warmup_trials=2,
        strategy_config=RouletteWheelConfig(
            mode="maximize",
            alpha=0.1,
            beta=0.1,
        ),
        evaluator_config=DBEvaluatorConfig(db_filename="test.db"),
    )
    assert isinstance(config.strategy_config, RouletteWheelConfig)
    assert config.strategy_config.alpha == 0.1


def test_thompson_sampling_config_missing_field():
    """
    Test that missing required field raises ValidationError.
    """
    with pytest.raises(ValidationError):
        ThompsonSamplingConfig(
            # synthesis_pipeline missing
            num_ts_iterations=10,
            num_warmup_trials=2,
            strategy_config=GreedyConfig(mode="maximize"),
            evaluator_config=DBEvaluatorConfig(db_filename="test.db"),
        )


def test_thompson_sampling_config_invalid_strategy():
    """
    Test that invalid mode raises ValidationError.
    """
    with pytest.raises(ValidationError):
        GreedyConfig(mode="invalid_mode")


def test_config_field_types(pipeline):
    """
    Test that config models enforce correct field types.
    """
    with pytest.raises(ValidationError):
        ThompsonSamplingConfig(
            synthesis_pipeline=pipeline,
            num_ts_iterations="not_an_integer",  # invalid type
            num_warmup_trials=2,
            strategy_config=GreedyConfig(mode="maximize"),
            evaluator_config=DBEvaluatorConfig(db_filename="test.db"),
        )


def test_config_optional_fields(pipeline):
    """
    Test that optional fields work correctly with defaults.
    """
    config = ThompsonSamplingConfig(
        synthesis_pipeline=pipeline,
        num_ts_iterations=10,
        num_warmup_trials=2,
        strategy_config=GreedyConfig(mode="maximize"),
        evaluator_config=DBEvaluatorConfig(db_filename="test.db"),
    )
    assert config.results_filename == "results.csv"  # default value
    assert config.log_filename is None


def test_random_baseline_config_valid(pipeline):
    """
    Test that RandomBaselineConfig can be created with valid data.
    """
    config = RandomBaselineConfig(
        synthesis_pipeline=pipeline,
        evaluator_config=LookupEvaluatorConfig(ref_filename="scores.csv"),
        num_trials=100,
        num_to_save=10,
    )
    assert config.num_trials == 100
    assert config.num_to_save == 10


def test_all_selection_strategies(pipeline):
    """
    Test that all valid selection strategies can be configured.
    """
    strategy_configs = [
        GreedyConfig(mode="maximize"),
        RouletteWheelConfig(mode="maximize"),
        UCBConfig(mode="maximize"),
        EpsilonGreedyConfig(mode="maximize"),
        BoltzmannConfig(mode="maximize_boltzmann"),
        BayesUCBConfig(mode="maximize"),
    ]

    for strategy_config in strategy_configs:
        config = ThompsonSamplingConfig(
            synthesis_pipeline=pipeline,
            num_ts_iterations=10,
            num_warmup_trials=3,
            strategy_config=strategy_config,
            evaluator_config=DBEvaluatorConfig(db_filename="test.db"),
        )
        assert config.strategy_config is strategy_config


def test_dict_evaluator_arg(pipeline):
    """
    Test that evaluator config with parameters works.
    """
    config = ThompsonSamplingConfig(
        synthesis_pipeline=pipeline,
        num_ts_iterations=10,
        num_warmup_trials=3,
        strategy_config=GreedyConfig(mode="maximize"),
        evaluator_config=LookupEvaluatorConfig(
            ref_filename="scores.csv",
            score_col="Score",
        ),
    )
    assert isinstance(config.evaluator_config, LookupEvaluatorConfig)
    assert config.evaluator_config.ref_filename == "scores.csv"
