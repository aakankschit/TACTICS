"""Block 3 — the smallest complete search, on the bundled thrombin library."""
# requires: none
from importlib.resources import files


def main():
    # [start:search]
    from TACTICS import ThompsonSampler, get_preset
    from TACTICS.library_enumeration import SynthesisPipeline, ReactionConfig, ReactionDef
    from TACTICS.thompson_sampling import LookupEvaluatorConfig

    data = files("TACTICS.data.thrombin")  # bundled example: 130 acids x 3844 amines

    # 1. Describe the library: one reaction, one reagent file per component
    pipeline = SynthesisPipeline(ReactionConfig(
        reactions=[ReactionDef(
            reaction_smarts="[#6:1](=[O:2])[OH].[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]):3]"
                            ">>[#6:1](=[O:2])[#7:3]",
            step_index=0,
        )],
        reagent_file_list=[str(data / "acids.smi"), str(data / "coupled_aa_sub.smi")],
    ))

    # 2. Describe how a product is scored (here: a precomputed docking table)
    evaluator = LookupEvaluatorConfig(ref_filename=str(data / "product_scores.parquet"))

    # 3. Take the tuned preset, run, and read the results
    config = get_preset(
        synthesis_pipeline=pipeline,
        evaluator_config=evaluator,
        mode="minimize",        # docking scores: lower is better
        num_iterations=20,      # cycles; 1000+ for a real screen
        batch_size=50,          # compounds per cycle
    )
    sampler = ThompsonSampler.from_config(config)
    sampler.warm_up(num_warmup_trials=config.num_warmup_trials)
    results = sampler.search(num_cycles=config.num_ts_iterations)
    sampler.close()

    print(results.sort("score").head(5))
    # [end:search]

    assert results.columns == ["score", "SMILES", "Name"]
    assert len(results) > 0
    return results


if __name__ == "__main__":
    main()
