"""
Simple Multi-SMARTS Pattern Support Demo
=========================================

This script demonstrates the multi-SMARTS pattern functionality using
configurations that closely match the integration tests.

Author: TACTICS Development Team
Date: November 2024
"""

import tempfile
from pathlib import Path
from TACTICS.thompson_sampling import ThompsonSamplingConfig
from TACTICS.thompson_sampling.core.sampler import ThompsonSampler
from TACTICS.thompson_sampling.strategies.config import GreedyConfig
from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig
from TACTICS.library_enumeration.smarts_toolkit import (
    MultiStepSynthesisConfig,
    AlternativeSMARTSConfig,
    SMARTSPatternConfig,
)


def main():
    """Demonstrate multi-SMARTS functionality."""
    print("\n" + "=" * 80)
    print("TACTICS Multi-SMARTS Pattern Support Demo")
    print("=" * 80)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Create test reagent files
        acids = tmpdir / "acids.smi"
        acids.write_text(
            "CC(=O)O acetic_acid\n"
            "CCC(=O)O propionic_acid\n"
            "c1ccccc1C(=O)O benzoic_acid\n"
        )

        amines = tmpdir / "amines.smi"
        amines.write_text(
            "CCN ethylamine\n"
            "CCCN propylamine\n"
            "CCNCC diethylamine\n"
        )

        # Create scores file
        scores = tmpdir / "scores.csv"
        scores.write_text(
            "Product_Code,Scores\n"
            "acetic_acid_ethylamine,10.0\n"
            "acetic_acid_propylamine,9.5\n"
            "propionic_acid_ethylamine,9.0\n"
            "propionic_acid_propylamine,8.5\n"
            "benzoic_acid_ethylamine,8.0\n"
            "benzoic_acid_propylamine,7.5\n"
            "acetic_acid_diethylamine,6.0\n"
            "propionic_acid_diethylamine,5.5\n"
            "benzoic_acid_diethylamine,5.0\n"
        )

        print("\nExample 1: Backward Compatible - Simple SMARTS")
        print("-" * 80)

        # Old-style configuration
        config = ThompsonSamplingConfig(
            reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
            reagent_file_list=[str(acids), str(amines)],
            num_ts_iterations=5,
            num_warmup_trials=2,
            strategy_config=GreedyConfig(mode="maximize"),
            evaluator_config=LookupEvaluatorConfig(ref_filename=str(scores))
        )

        sampler = ThompsonSampler.from_config(config)
        results = sampler.search(num_cycles=5)

        print(f"Found {len(results)} products using simple SMARTS")
        print(f"Best score: {results['score'].max():.2f}")
        print(f"\nTop 3 products:")
        top_3 = results.sort("score", descending=True).head(3)
        print(top_3.select(["Name", "score"]))

        print("\n\nExample 2: Multi-SMARTS with Alternative Patterns")
        print("-" * 80)

        # New-style configuration with alternative SMARTS
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
            reagent_file_list=[str(acids), str(amines)]
        )

        config2 = ThompsonSamplingConfig(
            reaction_sequence_config=reaction_config,
            reagent_file_list=[str(acids), str(amines)],
            num_ts_iterations=5,
            num_warmup_trials=3,
            strategy_config=GreedyConfig(mode="maximize"),
            evaluator_config=LookupEvaluatorConfig(ref_filename=str(scores))
        )

        sampler2 = ThompsonSampler.from_config(config2)
        results2 = sampler2.search(num_cycles=5)

        print(f"Found {len(results2)} products using alternative SMARTS")
        print(f"Best score: {results2['score'].max():.2f}")
        print(f"\nAll products:")
        print(results2.select(["Name", "score"]).sort("score", descending=True))

        print("\n" + "=" * 80)
        print("Key Features Demonstrated:")
        print("=" * 80)
        print("✓ Backward compatibility - old reaction_smarts parameter still works")
        print("✓ New reaction_sequence_config parameter for advanced features")
        print("✓ Alternative SMARTS patterns for different reagent types")
        print("✓ SMARTS router automatically selects appropriate pattern")
        print("\n" + "=" * 80)
        print("Demo completed successfully!")
        print("=" * 80)


if __name__ == "__main__":
    main()
