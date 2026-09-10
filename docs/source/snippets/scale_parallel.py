"""Block 4 — evaluate in parallel with a structure-based scoring function."""
# requires: none
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _thrombin import thrombin_pipeline


def main():
    pipeline = thrombin_pipeline(small=True)
    # [start:parallel]
    from TACTICS import ThompsonSampler, get_preset
    from TACTICS.thompson_sampling import FPEvaluatorConfig

    config = get_preset(
        synthesis_pipeline=pipeline,
        # similarity to a query needs the product structure, so every score
        # costs a synthesis + fingerprint: worth parallelising
        evaluator_config=FPEvaluatorConfig(query_smiles="CC(=O)Nc1ccc(O)cc1"),
        mode="maximize",
        num_iterations=6,
        batch_size=40,
    )
    config.processes = 2            # worker processes
    config.min_cpds_per_core = 10   # evaluate once 2 x 10 compounds have accumulated

    sampler = ThompsonSampler.from_config(config)   # workers rebuild the evaluator from its config
    sampler.warm_up(num_warmup_trials=2)
    results = sampler.search(num_cycles=config.num_ts_iterations)
    sampler.close()                 # shuts the worker pool down
    # [end:parallel]

    assert len(results) > 0


if __name__ == "__main__":   # required: spawn-based multiprocessing re-imports this file
    main()
