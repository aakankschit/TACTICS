"""Shared setup for the snippets: the bundled thrombin amide library."""
from importlib.resources import files

AMIDE = ("[#6:1](=[O:2])[OH].[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]):3]"
         ">>[#6:1](=[O:2])[#7:3]")


def thrombin_pipeline(small: bool = False):
    from TACTICS.library_enumeration import SynthesisPipeline, ReactionConfig, ReactionDef
    data = files("TACTICS.data.thrombin")
    amines = "amino_acids_no_fmoc.smi" if small else "coupled_aa_sub.smi"
    return SynthesisPipeline(ReactionConfig(
        reactions=[ReactionDef(reaction_smarts=AMIDE, step_index=0)],
        reagent_file_list=[str(data / "acids.smi"), str(data / amines)],
    ))


def thrombin_scores() -> str:
    return str(files("TACTICS.data.thrombin") / "product_scores.parquet")
