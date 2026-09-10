"""Block 1 — define a two-component amide library and make one product."""
# requires: none
from importlib.resources import files


def main():
    # [start:define]
    from TACTICS.library_enumeration import SynthesisPipeline, ReactionConfig, ReactionDef

    data = files("TACTICS.data.thrombin")

    config = ReactionConfig(
        reactions=[ReactionDef(
            # reactant 1: carboxylic acid, reactant 2: primary/secondary amine
            reaction_smarts="[#6:1](=[O:2])[OH].[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]):3]"
                            ">>[#6:1](=[O:2])[#7:3]",
            step_index=0,
            description="Amide coupling",
        )],
        reagent_file_list=[str(data / "acids.smi"), str(data / "coupled_aa_sub.smi")],
    )
    pipeline = SynthesisPipeline(config)

    print(pipeline.num_components, "components,", pipeline.num_steps, "step")
    # [end:define]

    # [start:single]
    from rdkit import Chem
    from TACTICS.library_enumeration import read_reagent_file

    acids = read_reagent_file(str(data / "acids.smi"))      # [(smiles, name), ...]
    amines = read_reagent_file(str(data / "coupled_aa_sub.smi"))

    result = pipeline.enumerate_single(
        [Chem.MolFromSmiles(acids[0][0]), Chem.MolFromSmiles(amines[0][0])],
        reagent_keys=[acids[0][1], amines[0][1]],
    )
    print(result.product_name, "->", result.product_smiles)
    # [end:single]

    # [start:validate]
    check = config.reactions[0].validate_reaction(
        reagent_files=config.reagent_file_list,
    )
    for position, pct in check.coverage_stats.items():
        print(f"reactant {position}: {pct:.0f}% of reagents match the template")
    # [end:validate]

    assert pipeline.num_components == 2 and pipeline.num_steps == 1
    assert result.success and result.product_smiles
    assert all(v > 90 for v in check.coverage_stats.values())


if __name__ == "__main__":
    main()
