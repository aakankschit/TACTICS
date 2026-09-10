"""Block 6 — compare two methods' recovery of the known top compounds."""
# requires: viz
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _thrombin import thrombin_pipeline, thrombin_scores


def run(preset: str, seed: int):
    from TACTICS import ThompsonSampler, get_preset
    from TACTICS.thompson_sampling import LookupEvaluatorConfig
    config = get_preset(preset, synthesis_pipeline=thrombin_pipeline(),
                        evaluator_config=LookupEvaluatorConfig(ref_filename=thrombin_scores()),
                        mode="minimize", num_iterations=20, batch_size=50)
    config.seed = seed
    s = ThompsonSampler.from_config(config); s.set_hide_progress(True)
    s.warm_up(num_warmup_trials=config.num_warmup_trials)
    df = s.search(num_cycles=config.num_ts_iterations); s.close()
    return df


def main():
    # [start:benchmarks]
    import polars as pl
    from TACTICS.library_analysis import TS_Benchmarks

    # one search() DataFrame per replicate, per method
    runs = {
        "TT-TS (recommended)": [run("recommended", seed) for seed in (1, 2)],
        "RWS / CATS":          [run("recommended_rws", seed) for seed in (1, 2)],
    }
    reference = (
        pl.read_parquet(thrombin_scores())
        .rename({"Product_Code": "Name", "Scores": "score"})
        .select(["score", "Name"])              # same columns, same order, as search() output
    )

    bench = TS_Benchmarks(
        no_of_cycles=2,                    # number of replicates per method (not search cycles)
        methods_list=list(runs),
        TS_runs_data=runs,
        reference_data=reference,
        top_n=100,                         # "hit" = one of the 100 best reference compounds
        sort_type="minimize",
    )
    bench.plot_barplot_TS_results(save_path="recovery_bar.html", show_plot=False)
    bench.plot_line_performance_with_error_bars(save_path="recovery_curve.html", show_plot=False)
    # [end:benchmarks]

    assert os.path.exists("recovery_bar.html") and os.path.exists("recovery_curve.html")


if __name__ == "__main__":
    main()
