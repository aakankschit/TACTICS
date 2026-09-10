"""Block 2 — score products with your own Python function."""
# requires: none
from importlib.resources import files


# [start:function]
from rdkit.Chem import Descriptors, QED


def drug_likeness(mol) -> float:
    """Any callable taking an RDKit Mol and returning a float will do."""
    return QED.qed(mol) - 0.01 * max(0.0, Descriptors.MolWt(mol) - 500)
# [end:function]


def main():
    from TACTICS import ThompsonSampler, get_preset
    from TACTICS.library_enumeration import SynthesisPipeline, ReactionConfig, ReactionDef
    data = files("TACTICS.data.thrombin")
    pipeline = SynthesisPipeline(ReactionConfig(
        reactions=[ReactionDef(
            reaction_smarts="[#6:1](=[O:2])[OH].[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]):3]"
                            ">>[#6:1](=[O:2])[#7:3]",
            step_index=0,
        )],
        reagent_file_list=[str(data / "acids.smi"), str(data / "amino_acids_no_fmoc.smi")],
    ))

    # [start:run]
    from TACTICS.thompson_sampling import CustomEvaluatorConfig

    config = get_preset(
        synthesis_pipeline=pipeline,
        evaluator_config=CustomEvaluatorConfig(scoring_function=drug_likeness),
        mode="maximize",
        num_iterations=10,
        batch_size=20,
    )
    sampler = ThompsonSampler.from_config(config)
    sampler.warm_up(num_warmup_trials=config.num_warmup_trials)
    results = sampler.search(num_cycles=config.num_ts_iterations)
    sampler.close()

    print(results.sort("score", descending=True).head(3))
    # [end:run]

    assert len(results) > 0 and results["SMILES"][0] != "FAIL"


if __name__ == "__main__":
    main()
