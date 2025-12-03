# sEH (Soluble Epoxide Hydrolase) Dataset

## Source

- **Origin**: DNA-Encoded Library (DEL) synthetic data
- **Reference**: DEL-based combinatorial library screening
- **Purpose**: sEH inhibitor discovery

## Reaction

- **Type**: Multi-step DEL synthesis
- **Components**: 3 building block sets (BB1, BB2, BB3)

## Components

| File | Description | Count |
|------|-------------|-------|
| `DEL_seH/bb1_list.csv` | Building block set 1 | Variable |
| `DEL_seH/bb2_list.csv` | Building block set 2 | Variable |
| `DEL_seH/bb3_list.csv` | Building block set 3 | Variable |
| `DEL_seH/unique_binder_bbs/` | Unique binder building blocks | Subset |
| `DEL_seH/del_nonbinders.csv` | Non-binding compounds (1.5GB, local only) | Large |
| `DEL_seH/total_compounds.csv` | All compounds (917MB, local only) | Large |

## Library Size

- **Total combinations**: Large DEL library (millions+)
- **Enumerated**: Partially (large files local only)

## Scores

- **Available**: Embedded in DEL data
- **Type**: DEL binding affinity data
- **Location**: Within dataset files

## Notes

- Large files (`del_nonbinders.csv`, `total_compounds.csv`) are **local only**
- Not tracked in git due to size (> 2GB total)
- Building block lists are the essential files for library generation

## Usage

```python
# Use building block lists for library generation
bb_files = [
    "data/reagents/seh/DEL_seH/bb1_list.csv",
    "data/reagents/seh/DEL_seH/bb2_list.csv",
    "data/reagents/seh/DEL_seH/bb3_list.csv"
]
```
