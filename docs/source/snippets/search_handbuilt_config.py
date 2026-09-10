"""Block 3 — build the config yourself instead of taking a preset."""
# requires: none
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _thrombin import thrombin_pipeline, thrombin_scores


def main():
    pipeline = thrombin_pipeline()
    # [start:config]
    from TACTICS import ThompsonSampler
    from TACTICS.thompson_sampling import (
        ThompsonSamplingConfig, TopTwoConfig, EnhancedWarmupConfig, LookupEvaluatorConfig,
    )

    config = ThompsonSamplingConfig(
        synthesis_pipeline=pipeline,
        evaluator_config=LookupEvaluatorConfig(ref_filename=thrombin_scores()),
        strategy_config=TopTwoConfig(mode="minimize", beta=0.5, heated_scale=2.0),
        warmup_config=EnhancedWarmupConfig(),
        num_warmup_trials=3,
        num_ts_iterations=15,
        batch_size=50,
        use_boltzmann_weighting=True,   # the posterior update the presets use
        seed=42,                        # reproducible selection and rotation
    )
    sampler = ThompsonSampler.from_config(config)
    sampler.warm_up(num_warmup_trials=config.num_warmup_trials)
    results = sampler.search(num_cycles=config.num_ts_iterations)
    sampler.close()

    results.write_parquet("run_seed42.parquet")   # keep what you found
    # [end:config]

    assert len(results) > 0 and os.path.exists("run_seed42.parquet")


if __name__ == "__main__":
    main()
