"""
Integration tests for selection strategies with ThompsonSampler.

These tests ensure that any selection strategy (existing or new) can be properly
integrated with the ThompsonSampler class.
"""

import pytest
import os
import tempfile
import shutil
import numpy as np

from TACTICS.thompson_sampling.core.sampler import ThompsonSampler
from TACTICS.thompson_sampling.core.evaluators import LookupEvaluator, MWEvaluator
from TACTICS.thompson_sampling.strategies.base_strategy import SelectionStrategy
from TACTICS.thompson_sampling.strategies.greedy_selection import GreedySelection
from TACTICS.thompson_sampling.strategies.roulette_wheel import RouletteWheelSelection
from TACTICS.thompson_sampling.strategies.ucb_selection import UCBSelection
from TACTICS.thompson_sampling.strategies.epsilon_greedy import EpsilonGreedySelection
from TACTICS.thompson_sampling.strategies.bayes_ucb_selection import BayesUCBSelection
from TACTICS.library_enumeration import SynthesisPipeline, ReactionConfig, ReactionDef


# Amide coupling reaction
REACTION_SMARTS = "[C:1](=O)[OH].[N:2]>>[C:1](=O)[N:2]"


class TestStrategyIntegration:
    """Test that selection strategies integrate properly with ThompsonSampler"""

    def setup_method(self):
        """Set up test fixtures with temporary files."""
        self.temp_dir = tempfile.mkdtemp()

        # Create reagent files
        self.reagent_file1 = os.path.join(self.temp_dir, "acids.smi")
        self.reagent_file2 = os.path.join(self.temp_dir, "amines.smi")

        with open(self.reagent_file1, "w") as f:
            f.write("CC(=O)O\tacetic_acid\n")
            f.write("CCC(=O)O\tpropionic_acid\n")
            f.write("CCCC(=O)O\tbutyric_acid\n")
            f.write("CCCCC(=O)O\tvaleric_acid\n")

        with open(self.reagent_file2, "w") as f:
            f.write("CN\tmethylamine\n")
            f.write("CCN\tethylamine\n")
            f.write("CCCN\tpropylamine\n")
            f.write("CCCCN\tbutylamine\n")

        # Create lookup file with scores
        self.lookup_file = os.path.join(self.temp_dir, "scores.csv")
        with open(self.lookup_file, "w") as f:
            f.write("Product_Code,Scores\n")
            for acid in [
                "acetic_acid",
                "propionic_acid",
                "butyric_acid",
                "valeric_acid",
            ]:
                for amine in ["methylamine", "ethylamine", "propylamine", "butylamine"]:
                    score = np.random.uniform(0.5, 1.0)
                    f.write(f"{acid}_{amine},{score:.2f}\n")

        # Create pipeline
        self.reaction_config = ReactionConfig(
            reactions=[ReactionDef(reaction_smarts=REACTION_SMARTS, step_index=0)],
            reagent_file_list=[self.reagent_file1, self.reagent_file2],
        )
        self.pipeline = SynthesisPipeline(self.reaction_config)

    def teardown_method(self):
        """Clean up temporary files."""
        shutil.rmtree(self.temp_dir)

    def test_greedy_selection_integration(self):
        """Test that GreedySelection integrates with ThompsonSampler"""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(self.pipeline, selection_strategy=strategy)
        sampler.read_reagents(self.pipeline.reagent_file_list)
        sampler.set_evaluator(MWEvaluator())

        # Directly run search without warmup
        search_results = sampler.search(num_cycles=5)

        assert len(search_results) > 0
        sampler.close()

    def test_roulette_wheel_integration(self):
        """Test that RouletteWheelSelection integrates with ThompsonSampler"""
        strategy = RouletteWheelSelection(mode="maximize", alpha=0.1, beta=0.1)
        sampler = ThompsonSampler(self.pipeline, selection_strategy=strategy)
        sampler.read_reagents(self.pipeline.reagent_file_list)
        sampler.set_evaluator(MWEvaluator())

        # Run search
        search_results = sampler.search(num_cycles=5)

        assert len(search_results) > 0
        sampler.close()

    def test_ucb_selection_integration(self):
        """Test that UCBSelection integrates with ThompsonSampler"""
        strategy = UCBSelection(mode="maximize", c=1.0)
        sampler = ThompsonSampler(self.pipeline, selection_strategy=strategy)
        sampler.read_reagents(self.pipeline.reagent_file_list)
        sampler.set_evaluator(MWEvaluator())

        # Run search
        search_results = sampler.search(num_cycles=5)

        assert len(search_results) > 0
        sampler.close()

    def test_epsilon_greedy_integration(self):
        """Test that EpsilonGreedySelection integrates with ThompsonSampler"""
        strategy = EpsilonGreedySelection(mode="maximize", epsilon=0.1)
        sampler = ThompsonSampler(self.pipeline, selection_strategy=strategy)
        sampler.read_reagents(self.pipeline.reagent_file_list)
        sampler.set_evaluator(MWEvaluator())

        # Run search
        search_results = sampler.search(num_cycles=5)

        assert len(search_results) > 0
        sampler.close()

    def test_strategy_mode_maximize(self):
        """Test that strategies work correctly in maximize mode"""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(self.pipeline, selection_strategy=strategy)
        sampler.read_reagents(self.pipeline.reagent_file_list)
        sampler.set_evaluator(MWEvaluator())

        # Run search
        results_df = sampler.search(num_cycles=10)

        # Verify results contain scores
        assert len(results_df) > 0
        scores = results_df["score"].to_list()
        assert all(isinstance(s, (int, float)) for s in scores)
        assert all(np.isfinite(s) for s in scores)

        sampler.close()

    def test_strategy_mode_minimize(self):
        """Test that strategies work correctly in minimize mode"""
        strategy = GreedySelection(mode="minimize")
        sampler = ThompsonSampler(self.pipeline, selection_strategy=strategy)
        sampler.read_reagents(self.pipeline.reagent_file_list)
        sampler.set_evaluator(MWEvaluator())

        # Run search
        results_df = sampler.search(num_cycles=10)

        # Verify results contain scores
        assert len(results_df) > 0
        scores = results_df["score"].to_list()
        assert all(isinstance(s, (int, float)) for s in scores)
        assert all(np.isfinite(s) for s in scores)

        sampler.close()

    def test_batch_selection_single_mode(self):
        """Test batch_size=1 (single selection per cycle)"""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(
            self.pipeline, selection_strategy=strategy, batch_size=1
        )
        sampler.read_reagents(self.pipeline.reagent_file_list)
        sampler.set_evaluator(MWEvaluator())

        results = sampler.search(num_cycles=5)

        assert len(results) > 0
        sampler.close()

    def test_batch_selection_batch_mode(self):
        """Test batch_size>1 (multiple selections per cycle)"""
        strategy = RouletteWheelSelection(mode="maximize")
        sampler = ThompsonSampler(
            self.pipeline, selection_strategy=strategy, batch_size=3
        )
        sampler.read_reagents(self.pipeline.reagent_file_list)
        sampler.set_evaluator(MWEvaluator())

        results = sampler.search(num_cycles=5)

        assert len(results) > 0
        sampler.close()

    def test_strategy_with_disallow_mask(self):
        """Test that strategies respect disallow mask"""
        strategy = GreedySelection(mode="maximize")
        sampler = ThompsonSampler(self.pipeline, selection_strategy=strategy)
        sampler.read_reagents(self.pipeline.reagent_file_list)
        sampler.set_evaluator(MWEvaluator())

        # Initialize reagents by running at least one evaluation
        sampler.evaluate([0, 0])

        # Test select_reagent with disallow mask
        rng = np.random.default_rng(seed=42)
        disallow_mask = {0}  # Disallow first reagent
        selected = strategy.select_reagent(
            reagent_list=sampler.reagent_lists[0], disallow_mask=disallow_mask, rng=rng
        )

        # Should not select disallowed index
        assert selected not in disallow_mask

        sampler.close()

    def test_roulette_wheel_temperature_adaptation(self):
        """Test RouletteWheelSelection adaptive temperature control"""
        strategy = RouletteWheelSelection(mode="maximize", alpha=0.1, beta=0.1)

        # Test temperature reset
        strategy.alpha = 0.5
        strategy.reset_temperature()
        assert strategy.alpha == strategy.initial_alpha
        assert strategy.beta == strategy.initial_beta

    def test_roulette_wheel_component_rotation(self):
        """Test RouletteWheelSelection component rotation for thermal cycling"""
        strategy = RouletteWheelSelection(mode="maximize")

        # Test component rotation
        assert strategy.current_component_idx == 0
        strategy.rotate_component(n_components=3)
        assert strategy.current_component_idx == 1
        strategy.rotate_component(n_components=3)
        assert strategy.current_component_idx == 2
        strategy.rotate_component(n_components=3)
        assert strategy.current_component_idx == 0  # Should wrap around

    def test_custom_strategy_integration(self):
        """Test that a custom strategy can be integrated"""

        class CustomStrategy(SelectionStrategy):
            """Simple custom strategy that selects randomly"""

            def select_reagent(self, reagent_list, disallow_mask=None, **kwargs):
                rng = kwargs.get("rng", np.random.default_rng())
                available = list(range(len(reagent_list)))

                if disallow_mask:
                    available = [i for i in available if i not in disallow_mask]

                return rng.choice(available) if available else 0

        # Test custom strategy
        strategy = CustomStrategy(mode="maximize")
        sampler = ThompsonSampler(self.pipeline, selection_strategy=strategy)
        sampler.read_reagents(self.pipeline.reagent_file_list)
        sampler.set_evaluator(MWEvaluator())

        search_results = sampler.search(num_cycles=5)

        assert len(search_results) > 0
        sampler.close()

    def test_strategy_select_batch_method(self):
        """Test that select_batch method works for strategies that implement it"""
        strategy = RouletteWheelSelection(mode="maximize")
        sampler = ThompsonSampler(
            self.pipeline, selection_strategy=strategy, batch_size=5
        )
        sampler.read_reagents(self.pipeline.reagent_file_list)
        sampler.set_evaluator(MWEvaluator())

        # Initialize reagents
        sampler.evaluate([0, 0])

        # Test select_batch
        rng = np.random.default_rng(seed=42)
        batch = strategy.select_batch(
            reagent_list=sampler.reagent_lists[0], batch_size=5, rng=rng
        )

        assert len(batch) == 5
        assert all(0 <= idx < len(sampler.reagent_lists[0]) for idx in batch)

        sampler.close()

    def test_bayes_ucb_selection_integration(self):
        """Test that BayesUCBSelection integrates with ThompsonSampler"""
        strategy = BayesUCBSelection(
            mode="maximize", initial_p_high=0.90, initial_p_low=0.60
        )
        sampler = ThompsonSampler(self.pipeline, selection_strategy=strategy)
        sampler.read_reagents(self.pipeline.reagent_file_list)
        sampler.set_evaluator(MWEvaluator())

        # Run search
        search_results = sampler.search(num_cycles=5)

        assert len(search_results) > 0
        sampler.close()

    def test_bayes_ucb_percentile_adaptation(self):
        """Test BayesUCBSelection percentile reset"""
        strategy = BayesUCBSelection(
            mode="maximize",
            initial_p_high=0.90,
            initial_p_low=0.60,
        )

        # Test percentile modification and reset
        strategy.p_high = 0.95
        strategy.p_low = 0.50

        # Test percentile reset
        strategy.reset_percentiles()
        assert strategy.p_high == 0.90
        assert strategy.p_low == 0.60

    def test_bayes_ucb_component_rotation(self):
        """Test BayesUCBSelection component rotation for thermal cycling"""
        strategy = BayesUCBSelection(mode="maximize")

        # Test component rotation
        assert strategy.current_component_idx == 0
        strategy.rotate_component(n_components=3)
        assert strategy.current_component_idx == 1
        strategy.rotate_component(n_components=3)
        assert strategy.current_component_idx == 2
        strategy.rotate_component(n_components=3)
        assert strategy.current_component_idx == 0  # Should wrap around

    def test_bayes_ucb_with_thermal_cycling(self):
        """Test BayesUCB with thermal cycling in real sampling scenario"""
        strategy = BayesUCBSelection(
            mode="maximize",
            initial_p_high=0.90,
            initial_p_low=0.60,
        )
        sampler = ThompsonSampler(
            self.pipeline, selection_strategy=strategy, batch_size=3
        )
        sampler.read_reagents(self.pipeline.reagent_file_list)
        sampler.set_evaluator(MWEvaluator())

        # Run search
        results = sampler.search(num_cycles=10)

        assert len(results) > 0
        sampler.close()

    def test_roulette_wheel_weighted_rotation(self):
        """Test that weighted rotation heats flexible components more often."""
        from TACTICS.thompson_sampling.core.reagent import Reagent

        strategy = RouletteWheelSelection(mode="maximize")
        rng = np.random.default_rng(42)
        counts = np.zeros(3)

        # Component 0: critical (one dominant reagent)
        critical = []
        for i, m in enumerate([1.0, 1.0, 10.0]):
            r = Reagent(f"r0_{i}", "C"); r.mean = m; r.std = 0.5; r.n_samples = 10
            critical.append(r)
        # Component 1: flexible (all similar)
        flexible = []
        for i, m in enumerate([5.0, 5.0, 5.0]):
            r = Reagent(f"r1_{i}", "C"); r.mean = m; r.std = 1.0; r.n_samples = 10
            flexible.append(r)
        # Component 2: moderate
        moderate = []
        for i, m in enumerate([3.0, 5.0, 4.0]):
            r = Reagent(f"r2_{i}", "C"); r.mean = m; r.std = 0.5; r.n_samples = 10
            moderate.append(r)

        reagent_lists = [critical, flexible, moderate]

        for _ in range(1000):
            strategy.rotate_component_weighted(3, reagent_lists, rng=rng)
            counts[strategy.current_component_idx] += 1

        # Flexible component (idx=1) should be heated most often
        assert counts[1] > counts[0], (
            f"Flexible component (idx=1) should be heated more than "
            f"critical component (idx=0): {counts[1]} vs {counts[0]}"
        )


class TestRouletteWheelSoftmaxNumerics:
    """Regression tests for exp() overflow in the CATS Boltzmann softmax.

    The Boltzmann weight is exp(z) where z = (score - mean) / std / temperature.
    Because z is standardized, it is O(1) for a well-behaved score distribution.
    But on a heavy-tailed, zero-inflated landscape (e.g. DEL read counts, where
    most reagents score 0 and a few score very high) a single extreme outlier
    can push max(z) to ~sqrt(n), and a small CATS temperature divides that
    further. The naive exp(z) then overflows to inf, and inf/sum(inf) = nan,
    which makes rng.choice raise "probabilities contain NaN".

    Subtracting max(z) before exponentiating is shift-invariant -- the
    normalized probabilities are mathematically identical -- but keeps every
    exponent <= 0, so the largest weight is exp(0) = 1 and overflow is
    impossible.
    """

    def _heavy_tailed_reagents(self, n=500):
        """Zero-inflated reagents: most near zero, one extreme outlier."""
        from TACTICS.thompson_sampling.core.reagent import Reagent

        reagents = []
        for i in range(n):
            r = Reagent(f"r_{i}", "C")
            r.mean = 0.0
            r.std = 1e-6  # near-deterministic, so sampled scores ~= means
            r.n_samples = 10
            reagents.append(r)
        # One dominant reagent far out in the tail.
        reagents[0].mean = 1e6
        return reagents

    def test_select_reagent_survives_extreme_z_scores(self, monkeypatch):
        """A tiny temperature on a heavy tail must not produce NaN probabilities."""
        strategy = RouletteWheelSelection(mode="maximize")
        reagents = self._heavy_tailed_reagents()

        # Force a very small temperature; the naive exp() would overflow here.
        monkeypatch.setattr(
            strategy, "_get_component_temperature", lambda *a, **k: 1e-3
        )

        rng = np.random.default_rng(42)
        idx = strategy.select_reagent(reagents, rng=rng, component_idx=0)

        assert 0 <= idx < len(reagents)

    def test_select_batch_survives_extreme_z_scores(self, monkeypatch):
        """Same guarantee for the batch selection path."""
        strategy = RouletteWheelSelection(mode="maximize")
        reagents = self._heavy_tailed_reagents()

        monkeypatch.setattr(
            strategy, "_get_component_temperature", lambda *a, **k: 1e-3
        )

        rng = np.random.default_rng(42)
        indices = strategy.select_batch(reagents, batch_size=10, rng=rng,
                                        component_idx=0)

        assert len(indices) == 10
        assert all(0 <= i < len(reagents) for i in indices)

    def test_max_subtraction_is_shift_invariant(self):
        """The stable form gives identical probabilities where naive doesn't overflow."""
        rng = np.random.default_rng(0)
        scores = rng.normal(size=200) * 3.0 + 10.0
        z = (scores - np.mean(scores)) / np.std(scores) / 0.5

        naive = np.exp(z)
        naive_probs = naive / np.sum(naive)

        stable = np.exp(z - np.max(z))
        stable_probs = stable / np.sum(stable)

        assert np.isfinite(naive_probs).all()  # guard: naive is valid here
        np.testing.assert_allclose(stable_probs, naive_probs, rtol=1e-12)

    def test_naive_softmax_would_overflow_on_this_input(self):
        """Documents the failure mode the max-subtraction prevents."""
        z = np.zeros(500)
        z[0] = 1e4  # the standardized outlier / small-temperature case

        with np.errstate(over="ignore"):
            naive = np.exp(z)
            naive_probs = naive / np.sum(naive)

        # Naive overflows to inf and yields NaN after normalization.
        assert np.isinf(naive).any()
        assert np.isnan(naive_probs).any()

        # Stable form stays finite and normalized.
        stable = np.exp(z - np.max(z))
        stable_probs = stable / np.sum(stable)
        assert np.isfinite(stable_probs).all()
        np.testing.assert_allclose(np.sum(stable_probs), 1.0)
