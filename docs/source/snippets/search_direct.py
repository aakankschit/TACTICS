"""Block 3 — drive the sampler directly, without a config object."""
# requires: none
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _thrombin import thrombin_pipeline, thrombin_scores


def main():
    pipeline = thrombin_pipeline()
    # [start:direct]
    from TACTICS import ThompsonSampler, TopTwoSelection, LookupEvaluator

    sampler = ThompsonSampler(
        pipeline,
        selection_strategy=TopTwoSelection(mode="minimize"),
        batch_size=50,
        seed=7,
    )                                                   # warmup defaults to EnhancedWarmup()
    sampler.read_reagents(pipeline.reagent_file_list)
    sampler.set_evaluator(LookupEvaluator({"ref_filename": thrombin_scores()}))

    sampler.warm_up(num_warmup_trials=3)
    results = sampler.search(num_cycles=10, max_evaluations=400)   # stop after 400 scored
    sampler.close()
    # [end:direct]

    assert 0 < len(results) <= 400 + 50 * 3   # warmup rows are not counted by max_evaluations


if __name__ == "__main__":
    main()
