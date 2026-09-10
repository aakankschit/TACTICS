"""Block 1 — a two-step synthesis: step 1's product feeds step 2."""
# requires: none


def main():
    # [start:multistep]
    from TACTICS.library_enumeration import (
        ReactionConfig, ReactionDef, StepInput, InputSource, DeprotectionSpec,
    )

    config = ReactionConfig(
        reactions=[
            # step 0: Boc-protected amine + acid -> amide (the Boc survives)
            ReactionDef(
                reaction_smarts="[#6:1](=[O:2])[OH].[#7X3;H1,H2:3]>>[#6:1](=[O:2])[#7:3]",
                step_index=0,
                # remove Boc from the *product* so step 1 sees a free amine
                deprotections=[DeprotectionSpec(target="product", group="Boc")],
            ),
            # step 1: that amine + a second acid
            ReactionDef(
                reaction_smarts="[#6:1](=[O:2])[OH].[#7X3;H1,H2:3]>>[#6:1](=[O:2])[#7:3]",
                step_index=1,
            ),
        ],
        reagent_file_list=["acids_A.smi", "boc_amines.smi", "acids_B.smi"],
        step_inputs={
            0: [StepInput(source=InputSource.REAGENT_FILE, file_index=0),
                StepInput(source=InputSource.REAGENT_FILE, file_index=1)],
            1: [StepInput(source=InputSource.REAGENT_FILE, file_index=2),
                StepInput(source=InputSource.PREVIOUS_STEP, step_index=0)],
        },
    )
    print(config.num_steps, "steps; multi-step:", config.is_multi_step)
    # [end:multistep]

    assert config.num_steps == 2 and config.is_multi_step
    assert config.get_inputs_for_step(1)[1].source == InputSource.PREVIOUS_STEP


if __name__ == "__main__":
    main()
