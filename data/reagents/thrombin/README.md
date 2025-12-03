# Thrombin Dataset

## Source

- **Origin**: Internal screening library
- **Purpose**: Thrombin inhibitor optimization via amide coupling

## Reaction

- **SMARTS**: `[#6:1](=[O:2])[OH].[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]);!$(N[#6]~[!#6;!#16]):3]>>[#6:1](=[O:2])[#7:3]`
- **Type**: Amide coupling (carboxylic acid + amine)
- **Components**: 2 (acids + amines)

## Components

| File | Description | Count |
|------|-------------|-------|
| `acids.smi` | Carboxylic acids | ~150 |
| `amino_acids.smi` | Amino acids (with Fmoc protection) | ~20 |
| `amino_acids_deprotected.smi` | Amino acids (deprotected) | ~20 |
| `amino_acids_no_fmoc.smi` | Amino acids (no Fmoc) | ~20 |
| `coupled_aa_sub.smi` | Coupled amino acid derivatives | ~200 |

## Library Size

- **Total combinations**: ~30,000 (150 acids × 200 amines)
- **Enumerated**: Not pre-enumerated (generated on-the-fly)

## Scores

- **Available**: Yes
- **Type**: FRED docking scores (minimization)
- **Location**: `data/scores/thrombin/`
  - `product_scores.csv` - Consolidated scores
  - `smiles_scores.csv` - SMILES-indexed scores

## Usage

```python
from TACTICS.thompson_sampling import run_ts

reagent_files = [
    "data/reagents/thrombin/acids.smi",
    "data/reagents/thrombin/coupled_aa_sub.smi"
]

score_file = "data/scores/thrombin/product_scores.csv"
```
