"""Tests for Pydantic configuration models."""

import pytest
import tempfile
import os
from pydantic import ValidationError

from TACTICS.thompson_sampling.strategies.config import (
    GreedyConfig,
    RouletteWheelConfig,
    UCBConfig,
    EpsilonGreedyConfig,
    BoltzmannConfig,
)
from TACTICS.thompson_sampling.warmup.config import (
    StandardWarmupConfig,
    EnhancedWarmupConfig,
    BalancedWarmupConfig,
)
from TACTICS.thompson_sampling.core.evaluator_config import (
    LookupEvaluatorConfig,
    DBEvaluatorConfig,
    FPEvaluatorConfig,
    MWEvaluatorConfig,
    ROCSEvaluatorConfig,
    FredEvaluatorConfig,
    MLClassifierEvaluatorConfig,
)
from TACTICS.thompson_sampling.factories import (
    create_strategy,
    create_warmup,
    create_evaluator,
)
from TACTICS.thompson_sampling.config import ThompsonSamplingConfig
from TACTICS.thompson_sampling.presets import ConfigPresets, get_preset
from TACTICS.library_enumeration import SynthesisPipeline, ReactionConfig, ReactionDef


class TestStrategyConfigs:
    """Tests for selection strategy configuration models."""

    def test_greedy_config_creation(self):
        """Test GreedyConfig can be created with valid data."""
        config = GreedyConfig(mode="maximize")
        assert config.strategy_type == "greedy"
        assert config.mode == "maximize"

    def test_roulette_wheel_config_creation(self):
        """Test RouletteWheelConfig with all parameters."""
        config = RouletteWheelConfig(
            mode="minimize",
            alpha=0.2,
            beta=0.15,
        )
        assert config.strategy_type == "roulette_wheel"
        assert config.alpha == 0.2
        assert config.beta == 0.15

    def test_roulette_wheel_validation(self):
        """Test RouletteWheelConfig validates positive values."""
        with pytest.raises(ValidationError):
            RouletteWheelConfig(mode="maximize", alpha=-0.1)  # Negative alpha

    def test_ucb_config_creation(self):
        """Test UCBConfig creation."""
        config = UCBConfig(mode="maximize", c=3.0)
        assert config.strategy_type == "ucb"
        assert config.c == 3.0

    def test_epsilon_greedy_config(self):
        """Test EpsilonGreedyConfig creation and validation."""
        config = EpsilonGreedyConfig(epsilon=0.3, decay=0.99)
        assert config.epsilon == 0.3
        assert config.decay == 0.99

        with pytest.raises(ValidationError):
            EpsilonGreedyConfig(epsilon=1.5)  # > 1


class TestWarmupConfigs:
    """Tests for warmup strategy configuration models."""

    def test_standard_warmup_config(self):
        """Test StandardWarmupConfig creation."""
        config = StandardWarmupConfig()
        assert config.warmup_type == "standard"

    def test_balanced_warmup_config(self):
        """Test BalancedWarmupConfig creation."""
        config = BalancedWarmupConfig()
        assert config.warmup_type == "balanced"
        assert config.observations_per_reagent == 5  # default
        assert config.use_per_reagent_variance is True  # default

    def test_balanced_warmup_config_custom(self):
        """Test BalancedWarmupConfig with custom parameters."""
        config = BalancedWarmupConfig(
            observations_per_reagent=10,
            seed=42,
            use_per_reagent_variance=False,
            shrinkage_strength=5.0,
        )
        assert config.observations_per_reagent == 10
        assert config.seed == 42
        assert config.use_per_reagent_variance is False
        assert config.shrinkage_strength == 5.0


class TestEvaluatorConfigs:
    """Tests for evaluator configuration models."""

    def test_lookup_evaluator_config(self):
        """Test LookupEvaluatorConfig creation."""
        config = LookupEvaluatorConfig(ref_filename="scores.csv", score_col="Score")
        assert config.evaluator_type == "lookup"
        assert config.ref_filename == "scores.csv"
        assert config.score_col == "Score"

    def test_db_evaluator_config(self):
        """Test DBEvaluatorConfig creation."""
        config = DBEvaluatorConfig(db_filename="scores.db", db_prefix="test_")
        assert config.evaluator_type == "db"
        assert config.db_filename == "scores.db"

    def test_rocs_evaluator_config(self):
        """Test ROCSEvaluatorConfig creation and validation."""
        config = ROCSEvaluatorConfig(query_molfile="ref.sdf", max_confs=100)
        assert config.evaluator_type == "rocs"
        assert config.max_confs == 100

        with pytest.raises(ValidationError):
            ROCSEvaluatorConfig(query_molfile="ref.sdf", max_confs=-10)


class TestFactories:
    """Tests for factory functions."""

    def test_create_strategy_greedy(self):
        """Test creating GreedySelection from config."""
        config = GreedyConfig(mode="maximize")
        strategy = create_strategy(config)
        from TACTICS.thompson_sampling.strategies import GreedySelection

        assert isinstance(strategy, GreedySelection)
        assert strategy.mode == "maximize"

    def test_create_strategy_roulette_wheel(self):
        """Test creating RouletteWheelSelection from config."""
        config = RouletteWheelConfig(mode="minimize", alpha=0.15, beta=0.1)
        strategy = create_strategy(config)
        from TACTICS.thompson_sampling.strategies import RouletteWheelSelection

        assert isinstance(strategy, RouletteWheelSelection)
        assert strategy.alpha == 0.15

    def test_create_warmup_balanced(self):
        """Test creating BalancedWarmup from config."""
        config = BalancedWarmupConfig(observations_per_reagent=5)
        warmup = create_warmup(config)
        from TACTICS.thompson_sampling.warmup import BalancedWarmup

        assert isinstance(warmup, BalancedWarmup)
        assert warmup.observations_per_reagent == 5

    def test_create_evaluator_mw(self):
        """Test creating MWEvaluator from config."""
        config = MWEvaluatorConfig()
        evaluator = create_evaluator(config)
        from TACTICS.thompson_sampling.core.evaluators import MWEvaluator

        assert isinstance(evaluator, MWEvaluator)


class TestThompsonSamplingConfig:
    """Tests for unified ThompsonSamplingConfig."""

    def setup_method(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.reagent_file1 = os.path.join(self.temp_dir, "reagents1.smi")
        self.reagent_file2 = os.path.join(self.temp_dir, "reagents2.smi")

        with open(self.reagent_file1, "w") as f:
            f.write("CCO\tethanol\n")
        with open(self.reagent_file2, "w") as f:
            f.write("CC(=O)O\tacetic_acid\n")

        # Create a pipeline for testing
        self.reaction_smarts = "[C:1](=O)[OH].[O:2]>>[C:1](=O)[O:2]"
        reaction_config = ReactionConfig(
            reactions=[ReactionDef(reaction_smarts=self.reaction_smarts, step_index=0)],
            reagent_file_list=[self.reagent_file1, self.reagent_file2],
        )
        self.pipeline = SynthesisPipeline(reaction_config)

    def teardown_method(self):
        """Clean up test fixtures."""
        import shutil

        shutil.rmtree(self.temp_dir)

    def test_modern_config_creation(self):
        """Test creating config with modern nested configs."""
        config = ThompsonSamplingConfig(
            synthesis_pipeline=self.pipeline,
            num_ts_iterations=100,
            strategy_config=GreedyConfig(mode="maximize"),
            evaluator_config=MWEvaluatorConfig(),
        )
        assert config.strategy_config is not None
        assert config.num_ts_iterations == 100
        assert config.synthesis_pipeline is self.pipeline

    def test_legacy_config_creation(self):
        """Test that old-style config without synthesis_pipeline fails."""
        # The new API requires synthesis_pipeline
        with pytest.raises(ValidationError):
            ThompsonSamplingConfig(
                # Missing synthesis_pipeline
                num_ts_iterations=100,
                strategy_config=GreedyConfig(mode="maximize"),
                evaluator_config=MWEvaluatorConfig(),
            )

    def test_mixed_config_error(self):
        """Test that config validates required fields."""
        # With synthesis_pipeline, config should work
        config = ThompsonSamplingConfig(
            synthesis_pipeline=self.pipeline,
            num_ts_iterations=100,
            strategy_config=GreedyConfig(mode="maximize"),
            evaluator_config=MWEvaluatorConfig(),
        )
        assert config is not None

    def test_batch_size_validation(self):
        """Test batch_size validation."""
        with pytest.raises(ValidationError):
            ThompsonSamplingConfig(
                synthesis_pipeline=self.pipeline,
                num_ts_iterations=100,
                batch_size=0,  # Invalid
                strategy_config=GreedyConfig(),
                evaluator_config=MWEvaluatorConfig(),
            )

    def test_max_resamples_validation(self):
        """Test max_resamples can be set."""
        config = ThompsonSamplingConfig(
            synthesis_pipeline=self.pipeline,
            num_ts_iterations=100,
            max_resamples=500,
            strategy_config=GreedyConfig(),
            evaluator_config=MWEvaluatorConfig(),
        )
        assert config.max_resamples == 500


class TestPresets:
    """Tests for configuration presets."""

    def setup_method(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.reagent_file = os.path.join(self.temp_dir, "reagents.smi")
        with open(self.reagent_file, "w") as f:
            f.write("CCO\tethanol\n")
            f.write("CC(=O)O\tacetic_acid\n")

        # Create a pipeline for testing
        self.reaction_smarts = "[C:1](=O)[OH].[O:2]>>[C:1](=O)[O:2]"
        reaction_config = ReactionConfig(
            reactions=[ReactionDef(reaction_smarts=self.reaction_smarts, step_index=0)],
            reagent_file_list=[self.reagent_file, self.reagent_file],
        )
        self.pipeline = SynthesisPipeline(reaction_config)

    def teardown_method(self):
        """Clean up."""
        import shutil

        shutil.rmtree(self.temp_dir)

    def test_fast_exploration_preset(self):
        """Test fast_exploration preset."""
        config = ConfigPresets.fast_exploration(
            synthesis_pipeline=self.pipeline,
            evaluator_config=MWEvaluatorConfig(),
            num_iterations=500,
        )
        assert config.num_ts_iterations == 500
        assert isinstance(config.strategy_config, EpsilonGreedyConfig)
        assert config.strategy_config.mode == "maximize"
        assert config.batch_size == 1

    def test_fast_exploration_minimize_mode(self):
        """Test fast_exploration preset with minimize mode."""
        config = ConfigPresets.fast_exploration(
            synthesis_pipeline=self.pipeline,
            evaluator_config=MWEvaluatorConfig(),
            mode="minimize",
        )
        assert config.strategy_config.mode == "minimize"
        assert isinstance(config.strategy_config, EpsilonGreedyConfig)

    def test_parallel_batch_preset(self):
        """Test parallel_batch preset."""
        config = ConfigPresets.parallel_batch(
            synthesis_pipeline=self.pipeline,
            evaluator_config=MWEvaluatorConfig(),
            batch_size=50,
        )
        assert config.batch_size == 50
        assert isinstance(config.strategy_config, RouletteWheelConfig)
        assert config.strategy_config.mode == "maximize"

    def test_parallel_batch_minimize_mode(self):
        """Test parallel_batch preset with minimize mode for docking."""
        config = ConfigPresets.parallel_batch(
            synthesis_pipeline=self.pipeline,
            evaluator_config=MWEvaluatorConfig(),
            mode="minimize",
            batch_size=100,
        )
        assert config.strategy_config.mode == "minimize"
        assert config.batch_size == 100

    def test_all_presets_support_mode(self):
        """Test that all presets support mode parameter."""
        presets = [
            "fast_exploration",
            "parallel_batch",
            "conservative_exploit",
            "balanced_sampling",
            "diverse_coverage",
        ]

        for preset_name in presets:
            # Test maximize
            config = get_preset(
                preset_name,
                synthesis_pipeline=self.pipeline,
                evaluator_config=MWEvaluatorConfig(),
                mode="maximize",
            )
            assert config.strategy_config.mode == "maximize"

            # Test minimize
            config = get_preset(
                preset_name,
                synthesis_pipeline=self.pipeline,
                evaluator_config=MWEvaluatorConfig(),
                mode="minimize",
            )
            assert config.strategy_config.mode == "minimize"

    def test_get_preset_function(self):
        """Test get_preset convenience function."""
        config = get_preset(
            "balanced_sampling",
            synthesis_pipeline=self.pipeline,
            evaluator_config=MWEvaluatorConfig(),
        )
        assert isinstance(config.strategy_config, UCBConfig)

    def test_get_preset_invalid_name(self):
        """Test get_preset with invalid name."""
        with pytest.raises(ValueError, match="Unknown preset"):
            get_preset(
                "nonexistent_preset",
                synthesis_pipeline=self.pipeline,
                evaluator_config=MWEvaluatorConfig(),
            )

    def test_minimize_mode_removed(self):
        """Test that minimize_mode preset is removed (no longer needed)."""
        with pytest.raises(ValueError, match="Unknown preset"):
            get_preset(
                "minimize_mode",
                synthesis_pipeline=self.pipeline,
                evaluator_config=MWEvaluatorConfig(),
            )
