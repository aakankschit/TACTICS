"""Block 6 — plot a TT-TS diagnostic trajectory."""
# requires: viz
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _thrombin import thrombin_pipeline, thrombin_scores


def main():
    from TACTICS import ThompsonSampler, get_preset
    from TACTICS.thompson_sampling import LookupEvaluatorConfig
    config = get_preset(synthesis_pipeline=thrombin_pipeline(),
                        evaluator_config=LookupEvaluatorConfig(ref_filename=thrombin_scores()),
                        mode="minimize", num_iterations=30, batch_size=50)
    config.track_diagnostics = True
    s = ThompsonSampler.from_config(config); s.set_hide_progress(True)
    s.warm_up(num_warmup_trials=config.num_warmup_trials); s.search(num_cycles=30)
    diagnostics = s.get_diagnostics(); s.close()

    # [start:plot]
    import polars as pl
    from TACTICS.library_analysis.diagnostic_plots import plot_ttts_diagnostic

    # the combined plots average over replicates, so they expect a `replicate`
    # column; a single run is replicate 0
    diag = diagnostics.with_columns(pl.lit(0).alias("replicate"))

    fig = plot_ttts_diagnostic(diag, title="thrombin — TT-TS")
    fig.savefig("ttts_diagnostic.png", dpi=150, bbox_inches="tight")
    # [end:plot]

    assert os.path.exists("ttts_diagnostic.png")


if __name__ == "__main__":
    main()
