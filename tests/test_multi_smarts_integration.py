"""
Integration tests for multi-SMARTS pattern support in Thompson Sampling.

These tests demonstrate the full workflow from configuration to execution.
"""

import pytest
import tempfile
import os
from pathlib import Path

from TACTICS.thompson_sampling import ThompsonSamplingConfig
from TACTICS.thompson_sampling.core.sampler import ThompsonSampler
from TACTICS.thompson_sampling.strategies.config import GreedyConfig
from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig
from TACTICS.library_enumeration.smarts_toolkit import (
    MultiStepSynthesisConfig,
    AlternativeSMARTSConfig,
    SMARTSPatternConfig
)


class TestMultiSMARTSIntegration:
    """Integration tests for multi-SMARTS Thompson Sampling."""

    @pytest.fixture
    def temp_dir(self):
        """Create temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def reagent_files(self, temp_dir):
        """Create test reagent files."""
        # Acids
        acids_file = temp_dir / "acids.smi"
        with open(acids_file, 'w') as f:
            f.write("CC(=O)O acetic_acid\n")
            f.write("CCC(=O)O propionic_acid\n")
            f.write("c1ccccc1C(=O)O benzoic_acid\n")

        # Amines (mix of primary and secondary)
        amines_file = temp_dir / "amines.smi"
        with open(amines_file, 'w') as f:
            f.write("CCN ethylamine\n")  # Primary
            f.write("CCCN propylamine\n")  # Primary
            f.write("CCNCC diethylamine\n")  # Secondary

        return [str(acids_file), str(amines_file)]

    @pytest.fixture
    def scores_file(self, temp_dir):
        """Create lookup scores file."""
        scores_file = temp_dir / "scores.csv"
        with open(scores_file, 'w') as f:
            f.write("Product_Code,Scores\n")  # Use correct column names
            # Primary amine products (higher scores)
            f.write("acetic_acid_ethylamine,10.0\n")
            f.write("acetic_acid_propylamine,9.5\n")
            f.write("propionic_acid_ethylamine,9.0\n")
            f.write("propionic_acid_propylamine,8.5\n")
            f.write("benzoic_acid_ethylamine,8.0\n")
            f.write("benzoic_acid_propylamine,7.5\n")
            # Secondary amine products (lower scores)
            f.write("acetic_acid_diethylamine,6.0\n")
            f.write("propionic_acid_diethylamine,5.5\n")
            f.write("benzoic_acid_diethylamine,5.0\n")

        return str(scores_file)

    def test_backward_compatibility_simple_smarts(self, reagent_files, scores_file):
        """Test that simple reaction_smarts still works (backward compatible)."""
        config = ThompsonSamplingConfig(
            reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
            reagent_file_list=reagent_files,
            num_ts_iterations=5,
            num_warmup_trials=2,
            strategy_config=GreedyConfig(mode="maximize"),
            evaluator_config=LookupEvaluatorConfig(ref_filename=scores_file)
        )

        sampler = ThompsonSampler.from_config(config)
        results = sampler.search(num_cycles=5)

        assert len(results) > 0
        assert "score" in results.columns
        # Should only get primary amine products (secondary won't react with NH2 pattern)
        assert all("ethylamine" in name or "propylamine" in name
                  for name in results["Name"] if name != "FAIL")

    def test_single_step_with_alternative_smarts(self, reagent_files, scores_file):
        """Test single-step with alternative SMARTS patterns."""
        # Create config with alternatives for secondary amines
        reaction_config = MultiStepSynthesisConfig(
            alternative_smarts=AlternativeSMARTSConfig(
                primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
                alternatives=[
                    SMARTSPatternConfig(
                        pattern_id="secondary_amine",
                        reaction_smarts="[C:1](=O)[OH].[NH:2]>>[C:1](=O)[N:2]",
                        description="For secondary amines"
                    )
                ]
            ),
            reagent_file_list=reagent_files
        )

        config = ThompsonSamplingConfig(
            reaction_sequence_config=reaction_config,
            reagent_file_list=reagent_files,
            num_ts_iterations=5,
            num_warmup_trials=3,
            strategy_config=GreedyConfig(mode="maximize"),
            evaluator_config=LookupEvaluatorConfig(ref_filename=scores_file)
        )

        sampler = ThompsonSampler.from_config(config)

        # Don't manually set compatibilities - let it use primary pattern for all
        # This is a simpler test - just verify the router infrastructure works

        results = sampler.search(num_cycles=5)

        assert len(results) > 0
        assert "score" in results.columns
        # Should get products successfully
        product_names = set(results["Name"])
        product_names.discard("FAIL")

        # Check we got products
        assert len(product_names) > 0

    def test_reaction_sequence_with_thompson_sampler(self, reagent_files, scores_file):
        """Test that ReactionSequence integrates correctly with ThompsonSampler."""
        reaction_config = MultiStepSynthesisConfig(
            reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
            reagent_file_list=reagent_files
        )

        config = ThompsonSamplingConfig(
            reaction_sequence_config=reaction_config,
            reagent_file_list=reagent_files,
            num_ts_iterations=8,
            num_warmup_trials=3,
            strategy_config=GreedyConfig(mode="maximize"),
            evaluator_config=LookupEvaluatorConfig(ref_filename=scores_file)
        )

        sampler = ThompsonSampler.from_config(config)
        assert sampler.reaction_sequence is not None
        assert sampler.reaction_sequence.num_steps == 1

        results = sampler.search(num_cycles=8)

        assert len(results) > 0
        assert "score" in results.columns
        # Should find high-scoring products
        max_score = results["score"].max()
        assert max_score >= 8.0  # Should find some good products

    def test_warmup_with_reaction_sequence(self, reagent_files, scores_file):
        """Test warmup works correctly with ReactionSequence."""
        reaction_config = MultiStepSynthesisConfig(
            reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
            reagent_file_list=reagent_files
        )

        config = ThompsonSamplingConfig(
            reaction_sequence_config=reaction_config,
            reagent_file_list=reagent_files,
            num_ts_iterations=10,
            num_warmup_trials=5,
            strategy_config=GreedyConfig(mode="maximize"),
            evaluator_config=LookupEvaluatorConfig(ref_filename=scores_file)
        )

        sampler = ThompsonSampler.from_config(config)

        # Run warmup
        warmup_df = sampler.warm_up(num_warmup_trials=5)

        # Verify warmup ran successfully
        assert len(warmup_df) > 0
        assert "score" in warmup_df.columns

        # Verify we found some good molecules during warmup
        max_score = warmup_df["score"].max()
        assert max_score >= 8.0  # Should find at least one good product

    def test_config_validation_requires_reaction(self, reagent_files, scores_file):
        """Test that config validation requires either reaction_smarts or reaction_sequence_config."""
        from pydantic_core import ValidationError

        with pytest.raises(ValueError, match="Must provide either reaction_smarts"):
            # Missing both reaction_smarts and reaction_sequence_config
            config = ThompsonSamplingConfig(
                # reaction_smarts missing
                # reaction_sequence_config missing
                reagent_file_list=reagent_files,
                num_ts_iterations=10,
                strategy_config=GreedyConfig(mode="maximize"),
                evaluator_config=LookupEvaluatorConfig(ref_filename=scores_file)
            )

    def test_multiple_patterns_with_routing(self, temp_dir, reagent_files):
        """Test SMARTS routing with multiple patterns."""
        # Create simple scores
        scores_file = temp_dir / "scores.csv"
        with open(scores_file, 'w') as f:
            f.write("Product_Code,Scores\n")  # Use correct column names
            f.write("acetic_acid_ethylamine,10.0\n")
            f.write("acetic_acid_diethylamine,7.0\n")

        reaction_config = MultiStepSynthesisConfig(
            alternative_smarts=AlternativeSMARTSConfig(
                primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
                alternatives=[
                    SMARTSPatternConfig(
                        pattern_id="alt_pattern",
                        reaction_smarts="[C:1](=O)[OH].[NH:2]>>[C:1](=O)[N:2]"
                    )
                ]
            ),
            reagent_file_list=reagent_files
        )

        config = ThompsonSamplingConfig(
            reaction_sequence_config=reaction_config,
            reagent_file_list=reagent_files,
            num_ts_iterations=5,
            strategy_config=GreedyConfig(mode="maximize"),
            evaluator_config=LookupEvaluatorConfig(ref_filename=str(scores_file))
        )

        sampler = ThompsonSampler.from_config(config)

        # Verify router was created
        assert sampler.reaction_sequence is not None
        assert sampler.reaction_sequence.num_steps >= 1

        # Simple test that it runs
        results = sampler.search(num_cycles=5)
        assert len(results) > 0


class TestBackwardCompatibility:
    """Ensure new features don't break existing functionality."""

    @pytest.fixture
    def simple_config(self, tmp_path):
        """Create simple config files."""
        acids = tmp_path / "acids.smi"
        acids.write_text("CC(=O)O acetic\nCCC(=O)O prop\n")

        amines = tmp_path / "amines.smi"
        amines.write_text("CCN eth\nCCCN prop\n")

        scores = tmp_path / "scores.csv"
        scores.write_text("Product_Code,Scores\nacetic_eth,10\nacetic_prop,9\nprop_eth,8\nprop_prop,7\n")

        return {
            'reagents': [str(acids), str(amines)],
            'scores': str(scores)
        }

    def test_old_style_config_still_works(self, simple_config):
        """Verify old-style configs without reaction_sequence_config work."""
        config = ThompsonSamplingConfig(
            reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
            reagent_file_list=simple_config['reagents'],
            num_ts_iterations=3,
            strategy_config=GreedyConfig(mode="maximize"),
            evaluator_config=LookupEvaluatorConfig(
                ref_filename=simple_config['scores']
            )
        )

        sampler = ThompsonSampler.from_config(config)

        # Old behavior: reaction set, no reaction_sequence
        assert sampler.reaction is not None
        assert sampler.reaction_sequence is None

        results = sampler.search(num_cycles=3)
        assert len(results) == 3

    def test_new_style_config_works(self, simple_config):
        """Verify new-style configs with reaction_sequence_config work."""
        reaction_config = MultiStepSynthesisConfig(
            reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
            reagent_file_list=simple_config['reagents']
        )

        config = ThompsonSamplingConfig(
            reaction_sequence_config=reaction_config,
            reagent_file_list=simple_config['reagents'],
            num_ts_iterations=3,
            strategy_config=GreedyConfig(mode="maximize"),
            evaluator_config=LookupEvaluatorConfig(
                ref_filename=simple_config['scores']
            )
        )

        sampler = ThompsonSampler.from_config(config)

        # New behavior: reaction_sequence set
        assert sampler.reaction_sequence is not None

        results = sampler.search(num_cycles=3)
        assert len(results) == 3
