"""Block 5 — what the search learned about each component."""
# requires: none
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _thrombin import thrombin_pipeline, thrombin_scores


def main():
    pipeline = thrombin_pipeline()
    from TACTICS import ThompsonSampler, get_preset
    from TACTICS.thompson_sampling import LookupEvaluatorConfig
    config = get_preset(synthesis_pipeline=pipeline,
                        evaluator_config=LookupEvaluatorConfig(ref_filename=thrombin_scores()),
                        mode="minimize", num_iterations=30, batch_size=50)
    # [start:track]
    config.track_diagnostics = True          # record per-cycle component state
    config.seed = 1

    sampler = ThompsonSampler.from_config(config)
    sampler.warm_up(num_warmup_trials=config.num_warmup_trials)
    results = sampler.search(num_cycles=config.num_ts_iterations)

    diagnostics = sampler.get_diagnostics()          # one row per (cycle, component)
    landscape = sampler.get_posterior_landscape()    # one row per reagent: mean, std, n_samples
    summary = sampler.get_sar_summary()              # dict: per-component structure + prose
    sampler.close()

    print(diagnostics.columns)
    print(summary["landscape_type"])
    # [end:track]

    # [start:analyse]
    from TACTICS.thompson_sampling.diagnostics import (
        compute_posterior_entropy, compute_convergence_point, format_sar_report,
    )

    print(format_sar_report(summary))
    entropy = compute_posterior_entropy(landscape, mode="minimize")
    converged = compute_convergence_point(diagnostics, threshold=0.3)
    print(entropy.join(converged, on="component_idx"))

    diagnostics.write_parquet("diagnostics.parquet")   # for post-hoc analysis and plots
    landscape.write_parquet("landscape.parquet")
    # [end:analyse]

    assert {"current_cycle", "component_idx", "criticality", "gmic"} <= set(diagnostics.columns)
    assert landscape.height == sum(len(rl) for rl in sampler.reagent_lists)
    assert summary["landscape_type"] in {"structured_SAR", "diffuse_SAR", "mixed", "insufficient_data"}


if __name__ == "__main__":
    main()
