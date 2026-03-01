"""
Tests for CATS diagnostics: criticality trajectory, posterior landscape,
enhanced diagnostics, SAR summary, and analysis functions.

Verifies:
1. get_posterior_landscape() returns correct shape and columns after a short run
2. get_diagnostics() returns empty DataFrame when track_diagnostics=False
3. get_diagnostics() returns populated DataFrame with correct columns when
   track_diagnostics=True using RWS strategy
4. get_diagnostics() returns empty DataFrame (no criticality) with GreedySelection
5. Criticality values are in [0, 1] range
6. Enhanced 15-column diagnostics schema for CATS strategies
7. get_component_state() interface and implementations
8. get_sar_summary() works for CATS and non-CATS strategies
9. Pure analysis functions in diagnostics.py
"""

import pytest
import math
import numpy as np
import importlib.resources
import polars as pl

from TACTICS.thompson_sampling.core.sampler import ThompsonSampler
from TACTICS.thompson_sampling.core.evaluators import LookupEvaluator
from TACTICS.thompson_sampling.strategies.greedy_selection import GreedySelection
from TACTICS.thompson_sampling.strategies.roulette_wheel import RouletteWheelSelection
from TACTICS.thompson_sampling.strategies.bayes_ucb_selection import BayesUCBSelection
from TACTICS.thompson_sampling.diagnostics import (
    compute_posterior_entropy,
    compute_convergence_point,
    compare_trajectory_vs_snapshot,
    format_sar_report,
)
from TACTICS.library_enumeration import SynthesisPipeline
from TACTICS.library_enumeration.smarts_toolkit import ReactionConfig, ReactionDef

AMIDE_COUPLING_SMARTS = (
    "[#6:1](=[O:2])[OH]."
    "[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]);!$(N[#6]~[!#6;!#16]):3]"
    ">>[#6:1](=[O:2])[#7:3]"
)


@pytest.fixture
def thrombin_paths():
    """Return paths for the bundled thrombin dataset."""
    data = importlib.resources.files("TACTICS.data.thrombin")
    return {
        "acids": str(data / "acids.smi"),
        "amines": str(data / "coupled_aa_sub.smi"),
        "scores": str(data / "product_scores.csv"),
    }


@pytest.fixture
def pipeline(thrombin_paths):
    """Create SynthesisPipeline with Thrombin data."""
    config = ReactionConfig(
        reactions=[
            ReactionDef(
                reaction_smarts=AMIDE_COUPLING_SMARTS,
                step_index=0,
            )
        ],
        reagent_file_list=[thrombin_paths["acids"], thrombin_paths["amines"]],
    )
    return SynthesisPipeline(config)


def _make_sampler(pipeline, thrombin_paths, strategy, track_diagnostics, seed=42):
    """Helper to create a configured sampler."""
    sampler = ThompsonSampler(
        synthesis_pipeline=pipeline,
        selection_strategy=strategy,
        seed=seed,
        track_diagnostics=track_diagnostics,
    )
    sampler.read_reagents(
        [thrombin_paths["acids"], thrombin_paths["amines"]]
    )
    evaluator = LookupEvaluator({"ref_filename": thrombin_paths["scores"]})
    sampler.set_evaluator(evaluator)
    sampler.set_hide_progress(True)
    return sampler


class TestPosteriorLandscape:
    """Tests for get_posterior_landscape()."""

    def test_columns_and_shape(self, pipeline, thrombin_paths):
        """Posterior landscape returns expected columns after warmup + search."""
        strategy = GreedySelection(mode="minimize")
        sampler = _make_sampler(pipeline, thrombin_paths, strategy, track_diagnostics=False)

        sampler.warm_up(num_warmup_trials=3)
        sampler.search(num_cycles=10)

        landscape = sampler.get_posterior_landscape()

        assert isinstance(landscape, pl.DataFrame)
        assert set(landscape.columns) == {
            "component_idx", "reagent_name", "mean", "std", "n_samples",
        }

        # Should have one row per reagent across all components
        total_reagents = sum(len(rl) for rl in sampler.reagent_lists)
        assert len(landscape) == total_reagents

        # component_idx should only contain 0 and 1 (two-component reaction)
        assert set(landscape["component_idx"].to_list()) == {0, 1}

        sampler.close()

    def test_works_without_diagnostics_flag(self, pipeline, thrombin_paths):
        """Posterior landscape works regardless of track_diagnostics setting."""
        strategy = RouletteWheelSelection(mode="minimize")
        sampler = _make_sampler(pipeline, thrombin_paths, strategy, track_diagnostics=False)

        sampler.warm_up(num_warmup_trials=3)
        sampler.search(num_cycles=10)

        landscape = sampler.get_posterior_landscape()
        assert len(landscape) > 0

        sampler.close()


class TestDiagnosticsOff:
    """Tests when track_diagnostics=False (default)."""

    def test_returns_empty_dataframe(self, pipeline, thrombin_paths):
        """get_diagnostics() returns empty DataFrame with correct schema."""
        strategy = RouletteWheelSelection(mode="minimize")
        sampler = _make_sampler(pipeline, thrombin_paths, strategy, track_diagnostics=False)

        sampler.warm_up(num_warmup_trials=3)
        sampler.search(num_cycles=10)

        diag = sampler.get_diagnostics()

        assert isinstance(diag, pl.DataFrame)
        assert len(diag) == 0
        assert set(diag.columns) == {"cycle", "component_idx", "criticality"}

        sampler.close()


class TestDiagnosticsRWS:
    """Tests with track_diagnostics=True using RouletteWheelSelection."""

    def test_populated_dataframe(self, pipeline, thrombin_paths):
        """get_diagnostics() returns populated DataFrame with RWS strategy."""
        strategy = RouletteWheelSelection(mode="minimize")
        sampler = _make_sampler(pipeline, thrombin_paths, strategy, track_diagnostics=True)

        sampler.warm_up(num_warmup_trials=3)
        sampler.search(num_cycles=50)

        diag = sampler.get_diagnostics()

        assert isinstance(diag, pl.DataFrame)
        assert len(diag) > 0
        # Should have enhanced schema now
        assert "current_cycle" in diag.columns
        assert "criticality" in diag.columns
        assert "component_idx" in diag.columns

        sampler.close()

    def test_criticality_in_valid_range(self, pipeline, thrombin_paths):
        """All criticality values should be in [0, 1]."""
        strategy = RouletteWheelSelection(mode="minimize")
        sampler = _make_sampler(pipeline, thrombin_paths, strategy, track_diagnostics=True)

        sampler.warm_up(num_warmup_trials=3)
        sampler.search(num_cycles=50)

        diag = sampler.get_diagnostics()
        assert diag["criticality"].min() >= 0.0
        assert diag["criticality"].max() <= 1.0

        sampler.close()

    def test_both_components_tracked(self, pipeline, thrombin_paths):
        """Both component indices should appear in diagnostics."""
        strategy = RouletteWheelSelection(mode="minimize")
        sampler = _make_sampler(pipeline, thrombin_paths, strategy, track_diagnostics=True)

        sampler.warm_up(num_warmup_trials=3)
        sampler.search(num_cycles=50)

        diag = sampler.get_diagnostics()
        assert set(diag["component_idx"].unique().to_list()) == {0, 1}

        sampler.close()


class TestDiagnosticsBayesUCB:
    """Tests with track_diagnostics=True using BayesUCBSelection."""

    def test_populated_dataframe(self, pipeline, thrombin_paths):
        """get_diagnostics() returns populated DataFrame with BayesUCB strategy."""
        strategy = BayesUCBSelection(mode="minimize")
        sampler = _make_sampler(pipeline, thrombin_paths, strategy, track_diagnostics=True)

        sampler.warm_up(num_warmup_trials=3)
        sampler.search(num_cycles=50)

        diag = sampler.get_diagnostics()

        assert isinstance(diag, pl.DataFrame)
        assert len(diag) > 0
        assert "current_cycle" in diag.columns
        assert "criticality" in diag.columns
        assert diag["criticality"].min() >= 0.0
        assert diag["criticality"].max() <= 1.0

        sampler.close()


class TestDiagnosticsGreedy:
    """Tests with GreedySelection (no criticality support)."""

    def test_returns_empty_dataframe(self, pipeline, thrombin_paths):
        """get_diagnostics() returns empty DataFrame with Greedy (no criticality)."""
        strategy = GreedySelection(mode="minimize")
        sampler = _make_sampler(pipeline, thrombin_paths, strategy, track_diagnostics=True)

        sampler.warm_up(num_warmup_trials=3)
        sampler.search(num_cycles=50)

        diag = sampler.get_diagnostics()

        assert isinstance(diag, pl.DataFrame)
        assert len(diag) == 0
        assert set(diag.columns) == {"cycle", "component_idx", "criticality"}

        sampler.close()


class TestGetComponentCriticality:
    """Unit tests for the get_component_criticality interface."""

    def test_base_returns_none(self):
        """Base SelectionStrategy.get_component_criticality returns None."""
        strategy = GreedySelection(mode="maximize")
        assert strategy.get_component_criticality([]) is None

    def test_rws_returns_float(self, pipeline, thrombin_paths):
        """RWS get_component_criticality returns a float after warmup."""
        strategy = RouletteWheelSelection(mode="minimize")
        sampler = _make_sampler(pipeline, thrombin_paths, strategy, track_diagnostics=False)
        sampler.warm_up(num_warmup_trials=5)

        crit = strategy.get_component_criticality(sampler.reagent_lists[0])
        assert isinstance(crit, float)
        assert 0.0 <= crit <= 1.0

        sampler.close()

    def test_bayes_ucb_returns_float(self, pipeline, thrombin_paths):
        """BayesUCB get_component_criticality returns a float after warmup."""
        strategy = BayesUCBSelection(mode="minimize")
        sampler = _make_sampler(pipeline, thrombin_paths, strategy, track_diagnostics=False)
        sampler.warm_up(num_warmup_trials=5)

        crit = strategy.get_component_criticality(sampler.reagent_lists[0])
        assert isinstance(crit, float)
        assert 0.0 <= crit <= 1.0

        sampler.close()


# =====================================================================
# New test classes for enhanced diagnostics
# =====================================================================


class TestEnhancedDiagnosticsRWS:
    """Verify 18-column enhanced diagnostics schema for RouletteWheelSelection."""

    EXPECTED_COLUMNS = {
        "component_idx", "current_cycle", "total_cycles",
        "criticality", "snr", "imbalance_strength",
        "normalized_entropy", "n_active_reagents",
        "participation_ratio", "effective_n", "sharpening_factor",
        "base_temp", "is_heated", "criticality_weight", "decay",
        "cats_multiplier", "effective_multiplier", "final_temperature",
    }

    def test_schema_columns(self, pipeline, thrombin_paths):
        """Enhanced diagnostics has all 18 columns."""
        strategy = RouletteWheelSelection(mode="minimize")
        sampler = _make_sampler(pipeline, thrombin_paths, strategy, track_diagnostics=True)
        sampler.warm_up(num_warmup_trials=3)
        sampler.search(num_cycles=50)

        diag = sampler.get_diagnostics()
        assert set(diag.columns) == self.EXPECTED_COLUMNS
        sampler.close()

    def test_snr_non_negative(self, pipeline, thrombin_paths):
        """SNR values should be >= 0 (or NaN for insufficient data)."""
        strategy = RouletteWheelSelection(mode="minimize")
        sampler = _make_sampler(pipeline, thrombin_paths, strategy, track_diagnostics=True)
        sampler.warm_up(num_warmup_trials=3)
        sampler.search(num_cycles=50)

        diag = sampler.get_diagnostics()
        snr_vals = diag["snr"].to_numpy()
        finite = np.isfinite(snr_vals)
        assert np.all(snr_vals[finite] >= 0.0)
        sampler.close()

    def test_final_temperature_positive(self, pipeline, thrombin_paths):
        """Final temperature should always be positive."""
        strategy = RouletteWheelSelection(mode="minimize")
        sampler = _make_sampler(pipeline, thrombin_paths, strategy, track_diagnostics=True)
        sampler.warm_up(num_warmup_trials=3)
        sampler.search(num_cycles=50)

        diag = sampler.get_diagnostics()
        assert diag["final_temperature"].min() > 0.0
        sampler.close()

    def test_criticality_weight_bounded(self, pipeline, thrombin_paths):
        """Criticality weight should be in [0, 1]."""
        strategy = RouletteWheelSelection(mode="minimize")
        sampler = _make_sampler(pipeline, thrombin_paths, strategy, track_diagnostics=True)
        sampler.warm_up(num_warmup_trials=3)
        sampler.search(num_cycles=50)

        diag = sampler.get_diagnostics()
        assert diag["criticality_weight"].min() >= 0.0
        assert diag["criticality_weight"].max() <= 1.0
        sampler.close()


class TestEnhancedDiagnosticsBayesUCB:
    """Verify enhanced diagnostics for BayesUCBSelection."""

    def test_snr_is_finite(self, pipeline, thrombin_paths):
        """BayesUCB now uses SNR dampening (aligned with RWS) — SNR should be finite."""
        strategy = BayesUCBSelection(mode="minimize")
        sampler = _make_sampler(pipeline, thrombin_paths, strategy, track_diagnostics=True)
        sampler.warm_up(num_warmup_trials=3)
        sampler.search(num_cycles=50)

        diag = sampler.get_diagnostics()
        snr_vals = diag["snr"].to_numpy()
        # Some early records may still be NaN (insufficient data), but not all
        finite = np.isfinite(snr_vals)
        assert finite.any(), "At least some SNR values should be finite"
        assert np.all(snr_vals[finite] >= 0.0)
        sampler.close()

    def test_same_schema_as_rws(self, pipeline, thrombin_paths):
        """BayesUCB uses the same 18-column schema as RWS."""
        strategy = BayesUCBSelection(mode="minimize")
        sampler = _make_sampler(pipeline, thrombin_paths, strategy, track_diagnostics=True)
        sampler.warm_up(num_warmup_trials=3)
        sampler.search(num_cycles=50)

        diag = sampler.get_diagnostics()
        assert set(diag.columns) == TestEnhancedDiagnosticsRWS.EXPECTED_COLUMNS
        sampler.close()


class TestGetComponentState:
    """Tests for the get_component_state() interface."""

    def test_base_returns_none(self):
        """Base SelectionStrategy.get_component_state returns None."""
        strategy = GreedySelection(mode="maximize")
        assert strategy.get_component_state([], 0, 0, 100) is None

    def test_rws_returns_dict_with_all_keys(self, pipeline, thrombin_paths):
        """RWS get_component_state returns dict with all expected keys."""
        strategy = RouletteWheelSelection(mode="minimize")
        sampler = _make_sampler(pipeline, thrombin_paths, strategy, track_diagnostics=False)
        sampler.warm_up(num_warmup_trials=5)

        state = strategy.get_component_state(
            sampler.reagent_lists[0], 0, 10, 100
        )
        assert isinstance(state, dict)
        expected_keys = {
            "component_idx", "current_cycle", "total_cycles",
            "criticality", "snr", "imbalance_strength",
            "normalized_entropy", "n_active_reagents",
            "participation_ratio", "effective_n", "sharpening_factor",
            "base_temp", "is_heated", "criticality_weight", "decay",
            "cats_multiplier", "effective_multiplier", "final_temperature",
        }
        assert set(state.keys()) == expected_keys
        assert state["component_idx"] == 0
        assert state["current_cycle"] == 10
        assert state["total_cycles"] == 100
        assert 0.0 <= state["criticality"] <= 1.0
        assert state["final_temperature"] > 0
        sampler.close()

    def test_bayes_ucb_returns_dict(self, pipeline, thrombin_paths):
        """BayesUCB get_component_state returns dict with same keys as RWS."""
        strategy = BayesUCBSelection(mode="minimize")
        sampler = _make_sampler(pipeline, thrombin_paths, strategy, track_diagnostics=False)
        sampler.warm_up(num_warmup_trials=5)

        state = strategy.get_component_state(
            sampler.reagent_lists[0], 0, 10, 100
        )
        assert isinstance(state, dict)
        # BayesUCB now has SNR dampening (aligned with RWS)
        assert isinstance(state["snr"], float)
        assert isinstance(state["imbalance_strength"], float)
        assert 0.0 <= state["criticality"] <= 1.0
        # IPR fields should be present
        assert "participation_ratio" in state
        assert "effective_n" in state
        assert "sharpening_factor" in state
        sampler.close()

    def test_rws_caches_state(self, pipeline, thrombin_paths):
        """RWS caches state in _last_component_states."""
        strategy = RouletteWheelSelection(mode="minimize")
        sampler = _make_sampler(pipeline, thrombin_paths, strategy, track_diagnostics=False)
        sampler.warm_up(num_warmup_trials=5)

        state = strategy.get_component_state(
            sampler.reagent_lists[0], 0, 10, 100
        )
        assert 0 in strategy._last_component_states
        assert strategy._last_component_states[0] is state
        sampler.close()


class TestSARSummary:
    """Tests for get_sar_summary()."""

    def test_works_for_cats_strategy(self, pipeline, thrombin_paths):
        """get_sar_summary() works with CATS (RWS) strategy."""
        strategy = RouletteWheelSelection(mode="minimize")
        sampler = _make_sampler(pipeline, thrombin_paths, strategy, track_diagnostics=True)
        sampler.warm_up(num_warmup_trials=3)
        sampler.search(num_cycles=50)

        summary = sampler.get_sar_summary()

        assert "components" in summary
        assert "landscape_type" in summary
        assert "summary_text" in summary
        assert len(summary["components"]) == 2  # Two-component reaction

        # Each component should have expected keys
        for comp in summary["components"]:
            assert "component_idx" in comp
            assert "concentration" in comp
            assert "convergence_assessment" in comp
            assert comp["convergence_assessment"] in {"structured", "diffuse", "insufficient_data"}

        # Landscape type should be valid
        assert summary["landscape_type"] in {
            "structured_SAR", "diffuse_SAR", "mixed", "insufficient_data"
        }

        # With track_diagnostics=True and CATS, convergence_dynamics should exist
        assert "convergence_dynamics" in summary
        assert len(summary["convergence_dynamics"]) == 2

        sampler.close()

    def test_works_for_non_cats_strategy(self, pipeline, thrombin_paths):
        """get_sar_summary() works with non-CATS (Greedy) strategy."""
        strategy = GreedySelection(mode="minimize")
        sampler = _make_sampler(pipeline, thrombin_paths, strategy, track_diagnostics=False)
        sampler.warm_up(num_warmup_trials=3)
        sampler.search(num_cycles=50)

        summary = sampler.get_sar_summary()

        assert "components" in summary
        assert "landscape_type" in summary
        assert len(summary["components"]) == 2

        # No convergence_dynamics without CATS diagnostics
        assert "convergence_dynamics" not in summary

        sampler.close()

    def test_correct_landscape_type(self, pipeline, thrombin_paths):
        """Landscape type should be consistent with component assessments."""
        strategy = RouletteWheelSelection(mode="minimize")
        sampler = _make_sampler(pipeline, thrombin_paths, strategy, track_diagnostics=True)
        sampler.warm_up(num_warmup_trials=3)
        sampler.search(num_cycles=100)

        summary = sampler.get_sar_summary()
        assessments = [c["convergence_assessment"] for c in summary["components"]]

        if all(a == "structured" for a in assessments):
            assert summary["landscape_type"] == "structured_SAR"
        elif all(a == "diffuse" for a in assessments):
            assert summary["landscape_type"] == "diffuse_SAR"
        elif all(a == "insufficient_data" for a in assessments):
            assert summary["landscape_type"] == "insufficient_data"
        else:
            assert summary["landscape_type"] == "mixed"

        sampler.close()

    def test_summary_text_is_nonempty(self, pipeline, thrombin_paths):
        """Summary text should be a non-empty string."""
        strategy = RouletteWheelSelection(mode="minimize")
        sampler = _make_sampler(pipeline, thrombin_paths, strategy, track_diagnostics=False)
        sampler.warm_up(num_warmup_trials=3)
        sampler.search(num_cycles=50)

        summary = sampler.get_sar_summary()
        assert isinstance(summary["summary_text"], str)
        assert len(summary["summary_text"]) > 10

        sampler.close()


class TestDiagnosticAnalysisFunctions:
    """Tests for pure analysis functions in diagnostics.py."""

    def test_compute_posterior_entropy_schema(self, pipeline, thrombin_paths):
        """compute_posterior_entropy returns correct schema."""
        strategy = RouletteWheelSelection(mode="minimize")
        sampler = _make_sampler(pipeline, thrombin_paths, strategy, track_diagnostics=False)
        sampler.warm_up(num_warmup_trials=3)
        sampler.search(num_cycles=50)

        landscape = sampler.get_posterior_landscape()
        entropy_df = compute_posterior_entropy(landscape, mode="minimize")

        assert isinstance(entropy_df, pl.DataFrame)
        assert set(entropy_df.columns) == {
            "component_idx", "n_active", "entropy",
            "normalized_entropy", "concentration", "snr",
        }
        assert len(entropy_df) == 2  # Two components

        # Concentration should be in [0, 1]
        concs = entropy_df["concentration"].to_numpy()
        finite = np.isfinite(concs)
        assert np.all(concs[finite] >= 0.0)
        assert np.all(concs[finite] <= 1.0)

        sampler.close()

    def test_compute_convergence_point(self, pipeline, thrombin_paths):
        """compute_convergence_point returns correct schema."""
        strategy = RouletteWheelSelection(mode="minimize")
        sampler = _make_sampler(pipeline, thrombin_paths, strategy, track_diagnostics=True)
        sampler.warm_up(num_warmup_trials=3)
        sampler.search(num_cycles=50)

        diag = sampler.get_diagnostics()
        conv = compute_convergence_point(diag)

        assert isinstance(conv, pl.DataFrame)
        assert set(conv.columns) == {
            "component_idx", "cycle_first_above", "cycle_stable",
        }
        assert len(conv) == 2

        sampler.close()

    def test_compare_trajectory_vs_snapshot(self, pipeline, thrombin_paths):
        """compare_trajectory_vs_snapshot returns correct schema."""
        strategy = RouletteWheelSelection(mode="minimize")
        sampler = _make_sampler(pipeline, thrombin_paths, strategy, track_diagnostics=True)
        sampler.warm_up(num_warmup_trials=3)
        sampler.search(num_cycles=50)

        diag = sampler.get_diagnostics()
        landscape = sampler.get_posterior_landscape()
        comp = compare_trajectory_vs_snapshot(diag, landscape, mode="minimize")

        assert isinstance(comp, pl.DataFrame)
        assert set(comp.columns) == {
            "component_idx", "snapshot_concentration",
            "trajectory_final_criticality",
            "cycle_first_above", "cycle_stable",
            "trajectory_adds_info",
        }
        assert len(comp) == 2

        sampler.close()

    def test_format_sar_report(self, pipeline, thrombin_paths):
        """format_sar_report produces a non-empty string."""
        strategy = RouletteWheelSelection(mode="minimize")
        sampler = _make_sampler(pipeline, thrombin_paths, strategy, track_diagnostics=True)
        sampler.warm_up(num_warmup_trials=3)
        sampler.search(num_cycles=50)

        summary = sampler.get_sar_summary()
        report = format_sar_report(summary)

        assert isinstance(report, str)
        assert "SAR Landscape Assessment" in report
        assert "Component 0" in report
        assert "Component 1" in report

        sampler.close()


class TestIPRCriticality:
    """Tests for IPR-based criticality metric."""

    def test_ipr_higher_than_shannon_for_peaked(self):
        """IPR should produce higher criticality than Shannon for peaked distributions."""
        from TACTICS.thompson_sampling.core.reagent import Reagent

        # Create peaked distribution: one dominant reagent
        reagents = []
        for i in range(20):
            r = Reagent(f"reagent_{i}", "CC")
            r.mean = 0.5 if i < 19 else 2.0  # One much better
            r.std = 0.3
            r.n_samples = 10
            reagents.append(r)

        strategy_ipr = RouletteWheelSelection(
            mode="maximize", criticality_metric="ipr", n_adaptive_sharpening=True
        )
        strategy_shannon = RouletteWheelSelection(
            mode="maximize", criticality_metric="shannon", n_adaptive_sharpening=False
        )

        crit_ipr = strategy_ipr._calculate_criticality(reagents)
        crit_shannon = strategy_shannon._calculate_criticality(reagents)

        assert crit_ipr > crit_shannon, (
            f"IPR criticality ({crit_ipr:.4f}) should exceed Shannon ({crit_shannon:.4f}) "
            f"for peaked distributions"
        )

    def test_nonzero_criticality_at_large_n(self):
        """IPR should produce nonzero criticality even at N=100."""
        from TACTICS.thompson_sampling.core.reagent import Reagent

        rng = np.random.default_rng(42)
        reagents = []
        for i in range(100):
            r = Reagent(f"reagent_{i}", "CC")
            # Create some structure: means drawn from bimodal distribution
            r.mean = rng.choice([0.3, 0.7]) + rng.normal(0, 0.05)
            r.std = 0.2
            r.n_samples = 10
            reagents.append(r)

        strategy = RouletteWheelSelection(
            mode="maximize", criticality_metric="ipr", n_adaptive_sharpening=True
        )

        crit = strategy._calculate_criticality(reagents)
        assert crit > 0.01, (
            f"IPR criticality at N=100 should be nonzero, got {crit:.6f}"
        )

    def test_sharpening_increases_criticality(self):
        """N-adaptive sharpening should increase criticality for large N."""
        from TACTICS.thompson_sampling.core.reagent import Reagent

        rng = np.random.default_rng(42)
        reagents = []
        for i in range(50):
            r = Reagent(f"reagent_{i}", "CC")
            r.mean = rng.normal(0.5, 0.1)
            r.std = 0.2
            r.n_samples = 10
            reagents.append(r)

        strategy_sharp = RouletteWheelSelection(
            mode="maximize", criticality_metric="ipr", n_adaptive_sharpening=True
        )
        strategy_no_sharp = RouletteWheelSelection(
            mode="maximize", criticality_metric="ipr", n_adaptive_sharpening=False
        )

        crit_sharp = strategy_sharp._calculate_criticality(reagents)
        crit_no_sharp = strategy_no_sharp._calculate_criticality(reagents)

        assert crit_sharp >= crit_no_sharp, (
            f"Sharpening ({crit_sharp:.4f}) should >= no sharpening ({crit_no_sharp:.4f})"
        )

    def test_shannon_backward_compat(self):
        """Shannon metric should produce same results as legacy behavior."""
        from TACTICS.thompson_sampling.core.reagent import Reagent

        reagents = []
        for i in range(5):
            r = Reagent(f"reagent_{i}", "CC")
            r.mean = float(i) * 0.5
            r.std = 0.5
            r.n_samples = 10
            reagents.append(r)

        strategy = RouletteWheelSelection(
            mode="maximize", criticality_metric="shannon", n_adaptive_sharpening=False
        )

        crit = strategy._calculate_criticality(reagents)
        assert 0.0 <= crit <= 1.0
        assert isinstance(crit, float)
