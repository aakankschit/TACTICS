# TACTICS Data Directory

Input data for TACTICS experiments (NOT distributed with package, local only).

## Structure

- **`reagents/`**: Reagent SMILES files organized by dataset
  - `thrombin/` - Thrombin inhibitor building blocks
  - `seh/` - Soluble epoxide hydrolase (DEL) building blocks
  - `adenine/` - Adenine receptor dataset
  - `quinazoline/` - Quinazoline dataset

- **`scores/`**: Pre-computed score lookup tables by dataset
  - `thrombin/` - Docking scores for thrombin library

- **`libraries/`**: Enumerated full libraries (gitignored, regenerable)

## Usage

All data files are local only and not tracked in git. They are used for testing and experiments.

```python
# Example: Load reagent files
reagent_files = [
    "data/reagents/thrombin/acids.smi",
    "data/reagents/thrombin/coupled_aa_sub.smi"
]

# Example: Load scores
score_file = "data/scores/thrombin/product_scores.csv"
```

## Notes

- Datasets are **not** included in the pip package
- All files in `data/` are **gitignored** (local testing only)
- Each dataset subdirectory contains a README with dataset details
