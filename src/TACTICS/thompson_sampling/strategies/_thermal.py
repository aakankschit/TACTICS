"""Shared thermal-cycling and GMIC machinery for component-aware strategies.

Private module. ``RouletteWheelSelection``, ``TopTwoSelection`` and
``BayesUCBSelection`` each heat one reaction component at a time and pick the
next component to heat with a flexibility-weighted draw. The two mixins here
hold that common code so it exists once:

* :class:`ThermalCyclingMixin` — the heated-component index and the two
  rotation methods the sampler dispatches on.
* :class:`GMICCriticalityMixin` — Gaussian Mutual Information Criticality and
  the GMIC-derived rotation weights used by RWS and TT-TS.

Strategies that do not cycle components (Greedy, UCB, EpsilonGreedy) do not
use these; ``SelectionStrategy`` stays free of thermal concepts so the
sampler's ``hasattr(strategy, "rotate_component_weighted")`` check keeps its
meaning.
"""

from typing import Dict, List

import numpy as np


class ThermalCyclingMixin:
    """Heated-component bookkeeping and rotation.

    Subclasses call :meth:`_init_thermal_cycling` from ``__init__`` and
    implement :meth:`_rotation_flexibility`.
    """

    current_component_idx: int

    def _init_thermal_cycling(self) -> None:
        self.current_component_idx = 0

    def rotate_component(self, n_components: int) -> None:
        """Round-robin rotation to the next component."""
        self.current_component_idx = (self.current_component_idx + 1) % n_components

    def rotate_component_weighted(self, n_components: int, reagent_lists, rng=None) -> None:
        """Rotate to the next heated component with flexibility-weighted probabilities.

        Draws exactly one sample from ``rng`` after deterministic weight
        computation, so seeded runs are reproducible across strategies.

        Parameters:
            n_components: Number of reaction components.
            reagent_lists: One list of Reagent objects per component.
            rng: ``numpy.random.Generator``; a fresh default generator if None.
        """
        if rng is None:
            rng = np.random.default_rng()
        flexibility = self._rotation_flexibility(reagent_lists)
        heat_probs = flexibility / flexibility.sum()
        self.current_component_idx = int(rng.choice(n_components, p=heat_probs))

    def _rotation_flexibility(self, reagent_lists) -> np.ndarray:
        """Per-component heating weights (higher = heated more often)."""
        raise NotImplementedError


class GMICCriticalityMixin(ThermalCyclingMixin):
    """GMIC criticality plus GMIC-weighted rotation.

    GMIC = 0.5 * log(1 + var(posterior means) / mean(posterior variances)).
    High GMIC = critical component (clear winners among reagents); low GMIC =
    flexible component (reagents look alike). Flexible components are heated
    more often: weight = 1 / (1 + GMIC).

    Subclasses call :meth:`_init_gmic_state` from ``__init__``.
    """

    _cached_gmics: Dict[int, float]

    def _init_gmic_state(self) -> None:
        self._init_thermal_cycling()
        # Per-component GMIC, refreshed on every weighted rotation.
        self._cached_gmics = {}

    def _calculate_gmic_details(self, reagent_list: List) -> tuple:
        """Return ``(gmic, details)`` for one component.

        ``details`` holds ``signal_var``, ``mean_noise_var`` and
        ``n_active_reagents``. GMIC is 0.0 when fewer than two reagents have
        been observed.

        There is deliberately no minimum-observation gate here. One used to
        exist in TopTwoSelection (2026-06): it zeroed a component's GMIC
        whenever its least-observed active reagent had fewer than N samples.
        On large components (adenine's 688 isocyanides) a single straggler
        pinned GMIC to zero every cycle — on 25/28 benchmark libraries — which
        over-weighted that component in rotation. Removing the gate lifted
        adenine TT-TS top-100 recovery 87.1 -> 93.2 and halved its replicate
        variance (sd 18.4 -> 10.6), with no change on the libraries where the
        gate never fired.
        """
        active = [r for r in reagent_list if r.n_samples > 0]
        details = {
            "signal_var": float("nan"),
            "mean_noise_var": float("nan"),
            "n_active_reagents": len(active),
        }
        if len(active) < 2:
            return 0.0, details
        means = np.array([r.mean for r in active])
        signal = float(np.var(means))
        noise = float(np.mean([r.std ** 2 for r in active]))
        details["signal_var"] = signal
        details["mean_noise_var"] = noise
        gmic = 0.5 * np.log1p(signal / max(noise, 1e-10))
        return float(gmic), details

    def _calculate_gmic(self, reagent_list: List) -> float:
        """Gaussian Mutual Information Criticality for one component."""
        return self._calculate_gmic_details(reagent_list)[0]

    def get_component_criticality(self, reagent_list: List) -> float:
        """GMIC criticality (used by the sampler for rotation)."""
        return self._calculate_gmic(reagent_list)

    def _rotation_flexibility(self, reagent_lists) -> np.ndarray:
        gmics = np.array([self._calculate_gmic(rl) for rl in reagent_lists], dtype=float)
        for i, g in enumerate(gmics):
            self._cached_gmics[i] = float(g)
        return 1.0 / (1.0 + gmics)
