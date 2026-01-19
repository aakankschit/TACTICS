"""Thrombin inhibitor tutorial dataset.

This dataset contains:
- acids.smi: 130 carboxylic acids
- amino_acids_no_fmoc.smi: 62 amino acids (Fmoc deprotected)
- coupled_aa_sub.smi: 3,844 coupled amino acid derivatives
- product_scores.csv: ~500K products with docking scores

The data can be accessed using importlib.resources:

    import importlib.resources

    with importlib.resources.files("TACTICS.data.thrombin") as data_dir:
        acids_file = data_dir / "acids.smi"
        scores_file = data_dir / "product_scores.csv"
"""
