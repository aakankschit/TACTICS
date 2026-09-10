"""Extending — a minimal component-aware strategy built on the GMIC mixin."""
# requires: none
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _thrombin import thrombin_pipeline, thrombin_scores

# [start:strategy]
import numpy as np
from TACTICS.thompson_sampling.strategies.base_strategy import SelectionStrategy
from TACTICS.thompson_sampling.strategies._thermal import GMICCriticalityMixin


class HeatedArgmax(GMICCriticalityMixin, SelectionStrategy):
    """Argmax of posterior draws, with the heated component's std inflated."""

    def __init__(self, mode="maximize", heat=2.0):
        super().__init__(mode)
        self._init_gmic_state()      # current_component_idx, _cached_gmics
        self.heat = heat

    def select_reagent(self, reagent_list, disallow_mask=None, **kwargs):
        rng = kwargs["rng"]                                   # never the global RNG
        heated = kwargs.get("component_idx", 0) == self.current_component_idx
        scale = self.heat if heated else 1.0
        draws = np.array([r.mean + scale * r.std * rng.standard_normal() for r in reagent_list])
        if disallow_mask:
            draws[list(disallow_mask)] = -np.inf if self.mode == "maximize" else np.inf
        return int(np.argmax(draws) if self.mode == "maximize" else np.argmin(draws))
    # rotate_component_weighted (GMIC-weighted) and get_component_criticality
    # are inherited from the mixin; the sampler finds them with hasattr.
# [end:strategy]


def main():
    from TACTICS import ThompsonSampler, LookupEvaluator
    pipeline = thrombin_pipeline()
    sampler = ThompsonSampler(pipeline, selection_strategy=HeatedArgmax(mode="minimize"),
                              batch_size=50, seed=3, track_diagnostics=True)
    sampler.read_reagents(pipeline.reagent_file_list)
    sampler.set_evaluator(LookupEvaluator({"ref_filename": thrombin_scores()}))
    sampler.warm_up(num_warmup_trials=2)
    results = sampler.search(num_cycles=8)
    sampler.close()
    assert len(results) > 0
    # the mixin's rotation ran and cached a GMIC per component
    assert set(sampler.selection_strategy._cached_gmics) == {0, 1}


if __name__ == "__main__":
    main()
