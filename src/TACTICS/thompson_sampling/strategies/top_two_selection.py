"""Top-Two Thompson Sampling (TT-TS) selection strategy.

Implements the Top-Two TS algorithm (Russo 2020) for best-arm identification.
Standard TS optimizes cumulative regret, but recovery is a top-K identification
problem. TT-TS redirects ~β fraction of the sampling budget toward reagents
where there is genuine uncertainty about top-K membership.

Supports asymmetric thermal cycling: heated components get inflated posterior
std (more disagreement → more challenger exploration), cooled components get
deflated std (more agreement → more exploitation). Component criticality via
GMIC is used strictly for weighted rotation — it does NOT modulate temperature.

Reference:
    Russo, D. (2020). Simple Bayesian Algorithms for Best-Arm Identification.
    Operations Research, 68(6), 1625-1647.
"""

import numpy as np
from typing import Any, Dict, List, Optional

from .base_strategy import SelectionStrategy


class TopTwoSelection(SelectionStrategy):
    """Top-Two Thompson Sampling with asymmetric thermal cycling.

    Two orthogonal mechanisms:
      1. **TT-TS selection** — draws two posterior samples; when they disagree
         on the best reagent, explores the challenger with probability β.
      2. **Thermal cycling** — scales posterior std asymmetrically across
         components. The heated component gets inflated uncertainty (more
         disagreement → more exploration), cooled components get deflated
         uncertainty (more agreement → more exploitation).

    Component criticality is computed via GMIC (Gaussian Mutual Information
    Criticality) and used ONLY for weighted rotation — deciding which component
    gets heated. It does not produce a temperature multiplier.

    Args:
        mode: "maximize" or "minimize" optimization direction.
        beta: Probability of selecting the challenger when samples disagree.
            β=0.5 gives equal exploration/exploitation balance.
        heated_scale: Multiplier on posterior std for the heated component.
            >1 inflates uncertainty → more TT-TS disagreement → exploration.
        cooled_scale: Multiplier on posterior std for cooled components.
            <1 deflates uncertainty → more TT-TS agreement → exploitation.
        min_observations: Minimum observations before trusting GMIC.
    """

    def __init__(
        self,
        mode: str = "maximize",
        beta: float = 0.5,
        heated_scale: float = 1.5,
        cooled_scale: float = 0.75,
        min_observations: int = 5,
        adaptive_temperature: bool = False,
        scale_increment: float = 0.01,
        cooled_scale_increment: float = 0.001,
        efficiency_threshold: float = 0.10,
        heated_scale_max: float = 5.0,
        adaptive_disagreement: bool = True,
        disagreement_window: int = 200,
        disagreement_high_threshold: float = 0.8,
        disagreement_low_threshold: float = 0.3,
        disagreement_decay_rate: float = 0.95,
        ema_alpha: float = 0.02,
        heated_scale_min: float = 1.0,
        gmic_convergence_gate: Optional[float] = None,
        max_growth_per_step: Optional[float] = None,
    ):
        super().__init__(mode)
        self.beta = beta
        self.initial_heated_scale = heated_scale
        self.initial_cooled_scale = cooled_scale
        self.heated_scale = heated_scale
        self.cooled_scale = cooled_scale
        self.min_observations = min_observations

        # Adaptive thermal cycling parameters (legacy efficiency-based)
        self.adaptive_temperature = adaptive_temperature
        self.scale_increment = scale_increment
        self.cooled_scale_increment = cooled_scale_increment
        self.efficiency_threshold = efficiency_threshold
        self.heated_scale_max = heated_scale_max

        # Bidirectional disagreement-rate adaptation
        self.adaptive_disagreement = adaptive_disagreement
        self.disagreement_high_threshold = disagreement_high_threshold
        self.disagreement_low_threshold = disagreement_low_threshold
        self.disagreement_decay_rate = disagreement_decay_rate
        self.ema_alpha = ema_alpha
        self.heated_scale_min = heated_scale_min

        # Growth limiters for adaptive disagreement
        self.gmic_convergence_gate = gmic_convergence_gate
        self.max_growth_per_step = max_growth_per_step

        # Per-component state
        self._heated_scale_per_component: Dict[int, float] = {}
        self._component_disagreement_ema: Dict[int, float] = {}
        self._component_disagreement_counts: Dict[int, int] = {}
        self._component_gmic: Dict[int, float] = {}

        # Global disagreement tracking (for diagnostics / backward compat)
        self.disagreement_window = disagreement_window
        self._disagreement_buffer: list = []
        self._disagreement_rate: float = 1.0
        self._total_selections: int = 0

        # Thermal cycling state
        self.current_component_idx = 0

    def select_reagent(self, reagent_list: List, disallow_mask=None, **kwargs) -> int:
        """Select a reagent using Top-Two Thompson Sampling.

        1. Scale posterior std by thermal cycling (heated/cooled).
        2. Draw two independent posterior samples with scaled std.
        3. Find the best reagent under each sample.
        4. If they disagree, select the challenger with probability β.
        5. If they agree, exploit the believed-best.
        """
        rng = kwargs.get("rng", np.random.default_rng())
        component_idx = kwargs.get("component_idx", 0)

        stds = np.array([r.std for r in reagent_list])
        mu = np.array([r.mean for r in reagent_list])
        n = len(reagent_list)

        # Thermal cycling: scale posterior std (per-component heated_scale)
        is_heated = component_idx == self.current_component_idx
        if is_heated:
            scale = self._heated_scale_per_component.get(
                component_idx, self.heated_scale
            )
        else:
            scale = self.cooled_scale
        effective_stds = stds * scale

        # Draw two independent posterior samples
        sample_1 = rng.normal(size=n) * effective_stds + mu
        sample_2 = rng.normal(size=n) * effective_stds + mu

        # Apply disallow mask
        if disallow_mask:
            mask_arr = np.array(list(disallow_mask))
            sample_1[mask_arr] = np.nan
            sample_2[mask_arr] = np.nan

        # Find best under each sample
        if self.mode in ("maximize", "maximize_boltzmann"):
            best_1 = int(np.nanargmax(sample_1))
            best_2 = int(np.nanargmax(sample_2))
        else:
            best_1 = int(np.nanargmin(sample_1))
            best_2 = int(np.nanargmin(sample_2))

        # Track disagreement for adaptive cycling (per-component)
        disagreed = best_1 != best_2
        self._record_disagreement(disagreed, component_idx)

        # Top-Two decision
        if disagreed:
            return best_2 if rng.random() < self.beta else best_1
        else:
            return best_1

    # ===== Disagreement-rate adaptive thermal cycling =====

    def _record_disagreement(self, disagreed: bool, component_idx: int = 0) -> None:
        """Record a disagreement event and update per-component EMA.

        Args:
            disagreed: Whether the two posterior samples disagreed.
            component_idx: Which component this selection was for.
        """
        self._total_selections += 1

        # Per-component EMA update
        if component_idx not in self._component_disagreement_ema:
            self._component_disagreement_ema[component_idx] = 1.0
            self._component_disagreement_counts[component_idx] = 0

        self._component_disagreement_counts[component_idx] += 1
        old_ema = self._component_disagreement_ema[component_idx]
        self._component_disagreement_ema[component_idx] = (
            self.ema_alpha * float(disagreed) + (1.0 - self.ema_alpha) * old_ema
        )

        # Global rolling window (for diagnostics / backward compat)
        self._disagreement_buffer.append(disagreed)
        if len(self._disagreement_buffer) > self.disagreement_window:
            self._disagreement_buffer.pop(0)
        if len(self._disagreement_buffer) >= self.disagreement_window:
            self._disagreement_rate = (
                sum(self._disagreement_buffer) / len(self._disagreement_buffer)
            )

    def adapt_heated_scale(self) -> bool:
        """Adapt per-component heated_scale based on disagreement rate.

        Called once per iteration (after rotation). Checks the EMA
        disagreement rate for the currently-heated component and adjusts
        its per-component heated_scale:

        - If disagreement > high_threshold: decay toward 1.0
        - If disagreement < low_threshold: grow away from 1.0
        - Otherwise: stable, no change

        Returns:
            True if heated_scale was adjusted, False otherwise.
        """
        if not self.adaptive_disagreement:
            return False

        comp_idx = self.current_component_idx
        if comp_idx not in self._component_disagreement_ema:
            return False

        # Require minimum observations before adapting
        min_obs = max(20, int(1.0 / self.ema_alpha))
        if self._component_disagreement_counts.get(comp_idx, 0) < min_obs:
            return False

        rate = self._component_disagreement_ema[comp_idx]
        current_scale = self._heated_scale_per_component.get(
            comp_idx, self.heated_scale
        )

        if rate > self.disagreement_high_threshold:
            # Too much disagreement — reduce toward 1.0
            new_scale = 1.0 + (current_scale - 1.0) * self.disagreement_decay_rate
            self._heated_scale_per_component[comp_idx] = max(
                self.heated_scale_min, new_scale
            )
            return True

        elif rate < self.disagreement_low_threshold:
            # GMIC convergence gate: skip inflation if component is already solved
            if self.gmic_convergence_gate is not None:
                gmic = self._component_gmic.get(comp_idx, 0.0)
                if gmic > self.gmic_convergence_gate:
                    return False

            # Too little disagreement — increase away from 1.0
            new_scale = 1.0 + (current_scale - 1.0) / self.disagreement_decay_rate

            # Optional per-step growth cap
            if self.max_growth_per_step is not None:
                new_scale = min(new_scale, current_scale + self.max_growth_per_step)

            self._heated_scale_per_component[comp_idx] = min(
                self.heated_scale_max, new_scale
            )
            return True

        return False

    @property
    def disagreement_rate(self) -> float:
        """Current global rolling disagreement rate."""
        return self._disagreement_rate

    @property
    def component_disagreement_rates(self) -> Dict[int, float]:
        """Per-component EMA disagreement rates."""
        return dict(self._component_disagreement_ema)

    @property
    def effective_heated_scale(self) -> float:
        """Current heated_scale for the heated component."""
        return self._heated_scale_per_component.get(
            self.current_component_idx, self.heated_scale
        )

    # ===== GMIC criticality (for rotation weighting only) =====

    def _calculate_gmic(self, reagent_list: List) -> float:
        """Calculate Gaussian Mutual Information Criticality for a component.

        GMIC = 0.5 * log(1 + var(means) / mean(noise_vars))

        High GMIC = critical component (clear winners among reagents).
        Low GMIC = flexible component (all reagents similar).

        Used strictly for weighted rotation — does NOT affect temperature
        or selection.

        Returns:
            GMIC value >= 0. Returns 0.0 if insufficient data.
        """
        active = [r for r in reagent_list if r.n_samples > 0]
        if len(active) < 2:
            return 0.0

        if min(r.n_samples for r in active) < self.min_observations:
            return 0.0

        means = np.array([r.mean for r in active])
        signal = np.var(means)
        noise = np.mean([r.std ** 2 for r in active])
        return float(0.5 * np.log1p(signal / max(noise, 1e-10)))

    def get_component_criticality(self, reagent_list: List) -> float:
        """Return GMIC criticality for a component (used by sampler for rotation)."""
        return self._calculate_gmic(reagent_list)

    # ===== Thermal cycling rotation =====

    def rotate_component(self, n_components: int):
        """Round-robin component rotation (fallback)."""
        self.current_component_idx = (self.current_component_idx + 1) % n_components

    def adapt_temperatures(self, n_unique, n_attempted):
        """Adapt thermal cycling scales based on sampling efficiency.

        Mirrors RouletteWheelSelection.adapt_temperatures but operates on
        posterior std scales instead of Boltzmann temperatures. When posteriors
        tighten, TT-TS samples agree more often (less challenger exploration),
        causing recovery to stall. Increasing heated_scale counteracts this by
        inflating uncertainty on the heated component, generating more
        disagreement and thus more exploration.

        Args:
            n_unique: Number of unique compounds generated in this batch.
            n_attempted: Number of compounds attempted in this batch.

        Returns:
            True if scales were adjusted, False otherwise.
        """
        if not self.adaptive_temperature:
            return False

        adjusted = False
        efficiency = n_unique / max(n_attempted, 1)

        if efficiency < self.efficiency_threshold and self.heated_scale < self.heated_scale_max:
            self.heated_scale += self.scale_increment
            adjusted = True

        if n_unique == 0:
            self.cooled_scale = max(0.01, self.cooled_scale - self.cooled_scale_increment)
            adjusted = True

        return adjusted

    def reset_temperature(self):
        """Reset thermal cycling scales and disagreement state to initial values."""
        self.heated_scale = self.initial_heated_scale
        self.cooled_scale = self.initial_cooled_scale
        self._heated_scale_per_component = {}
        self._component_disagreement_ema = {}
        self._component_disagreement_counts = {}
        self._component_gmic = {}
        self._disagreement_buffer = []
        self._disagreement_rate = 1.0
        self._total_selections = 0

    def get_component_state(
        self,
        reagent_list: List,
        component_idx: int,
        current_cycle: int,
        total_cycles: int,
    ) -> Optional[Dict[str, Any]]:
        """Return TT-TS diagnostic state for a component."""
        gmic = self._calculate_gmic(reagent_list)
        is_heated = component_idx == self.current_component_idx
        comp_scale = self._heated_scale_per_component.get(
            component_idx, self.heated_scale
        )
        ema_rate = self._component_disagreement_ema.get(
            component_idx, float("nan")
        )

        return {
            "component_idx": component_idx,
            "current_cycle": current_cycle,
            "total_cycles": total_cycles,
            "gmic": gmic,
            "criticality": gmic,
            "is_heated": is_heated,
            "heated_scale": comp_scale,
            "cooled_scale": self.cooled_scale,
            "effective_scale": comp_scale if is_heated else self.cooled_scale,
            "disagreement_ema": ema_rate,
            "disagreement_global": self._disagreement_rate,
            "n_active_reagents": sum(1 for r in reagent_list if r.n_samples > 0),
        }

    def rotate_component_weighted(self, n_components: int, reagent_lists, rng=None):
        """Rotate to next heated component using GMIC-weighted probabilities.

        Flexible components (low GMIC) get heated more often because they
        benefit most from the inflated uncertainty that TT-TS uses to
        generate challenger candidates.

        Weight = 1 / (1 + gmic): high GMIC → low weight → less heating.
        """
        if rng is None:
            rng = np.random.default_rng()
        gmics = np.array([self._calculate_gmic(rl) for rl in reagent_lists], dtype=float)

        # Cache GMIC values for adaptive disagreement gate
        for i, g in enumerate(gmics):
            self._component_gmic[i] = float(g)

        flexibility = 1.0 / (1.0 + gmics)
        heat_probs = flexibility / flexibility.sum()

        self.current_component_idx = int(rng.choice(n_components, p=heat_probs))
