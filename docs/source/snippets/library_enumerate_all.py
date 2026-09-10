"""Block 1 — enumerate every product without Thompson Sampling."""
# requires: none
from importlib.resources import files


def main():
    from TACTICS.library_enumeration import SynthesisPipeline, ReactionConfig, ReactionDef
    data = files("TACTICS.data.thrombin")
    # a small library so this runs in a second: the acids x the 62 amino acids
    pipeline = SynthesisPipeline(ReactionConfig(
        reactions=[ReactionDef(
            reaction_smarts="[#6:1](=[O:2])[OH].[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]):3]"
                            ">>[#6:1](=[O:2])[#7:3]",
            step_index=0,
        )],
        reagent_file_list=[str(data / "acids.smi"), str(data / "amino_acids_no_fmoc.smi")],
    ))

    # [start:enumerate]
    from TACTICS.library_enumeration import (
        results_to_dataframe, summarize_failures, write_enumerated_library,
    )

    results = pipeline.enumerate_library(show_progress=False)   # list of EnumerationResult
    products = results_to_dataframe(results)                    # Polars: product_name, product_smiles, ...
    print(len(products), "products;", summarize_failures(results)["failures"], "failures")

    write_enumerated_library(results, "library.smi", format="smi")
    # [end:enumerate]

    assert len(products) > 1000 and "product_smiles" in products.columns


if __name__ == "__main__":
    main()
