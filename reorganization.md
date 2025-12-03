# TACTICS Repository Reorganization Plan

## Proposed Directory Structure

```
TACTICS/
│
│ ══════════════════════════════════════════════════════════════════
│ PUBLISHABLE PACKAGE (included in pip install)
│ ══════════════════════════════════════════════════════════════════
├── src/TACTICS/                        # Core package
│   ├── __init__.py
│   ├── thompson_sampling/
│   ├── library_enumeration/
│   └── library_analysis/
│
│ ══════════════════════════════════════════════════════════════════
│ RESEARCH INFRASTRUCTURE (NOT distributed - personal use only)
│ ══════════════════════════════════════════════════════════════════
├── experiment_utils.py                 # Experiment management (at root, NOT in package)
│
├── data/                               # Input data (git tracked)
│   ├── reagents/                       # Reagent SMILES files by dataset
│   │   ├── thrombin/
│   │   │   ├── acids.smi
│   │   │   ├── amines.smi
│   │   │   └── README.md               # Dataset description, source, SMARTS
│   │   ├── seh/
│   │   ├── adenine/
│   │   └── quinazoline/
│   │
│   ├── scores/                         # Pre-computed score lookup tables
│   │   ├── thrombin/
│   │   │   ├── docking_scores.csv      # Consolidated from products_*_score.txt
│   │   │   └── README.md               # How scores were generated
│   │   └── ...
│   │
│   └── libraries/                      # Enumerated full libraries (gitignored - large)
│       ├── .gitkeep
│       └── README.md                   # Instructions for regenerating
│
├── experiments/                        # Research experiment scripts (git tracked)
│   ├── README.md                       # How to run experiments, conventions
│   │
│   ├── strategy_benchmarks/            # Comparing selection strategies
│   │   ├── rws_vs_bayesucb.py
│   │   └── README.md
│   │
│   ├── cats_investigation/             # CATS-specific research
│   │   ├── warmup_cycle_sweep.py
│   │   └── README.md
│   │
│   ├── warmup_analysis/                # Warmup strategy research
│   │   └── README.md
│   │
│   └── _templates/                     # Template scripts for new experiments
│       └── experiment_template.py
│
├── outputs/                            # Generated outputs (GITIGNORED)
│   ├── .gitignore                      # Contains: *\n!.gitignore\n!README.md
│   ├── README.md                       # Explains output structure
│   └── runs/
│       └── {YYYY-MM-DD}_{experiment_name}/
│           ├── config.yaml
│           ├── logs/
│           └── figures/
│
│ ══════════════════════════════════════════════════════════════════
│ USER-FACING (distributed with package)
│ ══════════════════════════════════════════════════════════════════
├── examples/                           # Clean, documented examples for USERS
│   ├── README.md
│   ├── quickstart_thrombin.py
│   └── notebooks/
│
├── tutorials/                          # In-depth tutorials
│
├── tests/                              # Unit/integration tests
│
├── docs/                               # Documentation
│
├── CLAUDE.md                           # AI assistant guide (NOT distributed)
├── README.md
├── pyproject.toml
└── LICENSE
```

---

## File Migration Plan

### Phase 1: Create New Directories

```bash
mkdir -p data/{reagents,scores,libraries}/{thrombin,seh,adenine,quinazoline}
mkdir -p experiments/{strategy_benchmarks,cats_investigation,warmup_analysis,_templates}
mkdir -p outputs/runs
mkdir -p examples/notebooks
mkdir -p scripts
```

### Phase 2: Move Input Data

| Current Location | New Location |
|-----------------|--------------|
| `examples/input_files/acids.smi` | `data/reagents/thrombin/acids.smi` |
| `examples/input_files/amino_acids*.smi` | `data/reagents/thrombin/` |
| `examples/input_files/adenine/` | `data/reagents/adenine/` |
| `examples/input_files/quinazoline/` | `data/reagents/quinazoline/` |
| `examples/input_files/DEL_seH/` | `data/reagents/seh/` |
| `examples/docking_scores/` | `data/scores/thrombin/` (consolidate first) |

### Phase 3: Move Experiment Scripts

| Current Location | New Location |
|-----------------|--------------|
| `examples/rws_bayesucb_cats_rocs.py` | `experiments/strategy_benchmarks/` |
| `examples/strategy_recovery_comparison.py` | `experiments/strategy_benchmarks/` |
| `examples/thrombin_strategy_comparison.py` | `experiments/strategy_benchmarks/` |
| `examples/balanced_warmup_variance_*.py` | `experiments/warmup_analysis/` |
| `examples/legacy_vs_current_*.py` | `experiments/validation/` |
| Root `debug_*.py` files | `experiments/_debug/` or DELETE |
| Root `test_fixed_*.py`, `test_new_*.py` | `experiments/_debug/` or DELETE |

### Phase 4: Move/Delete Output Files

| Current Location | Action |
|-----------------|--------|
| `*.png`, `*.pdf`, `*.csv` at root | DELETE (regenerate with new structure) |
| `examples/*.png`, `*.pdf`, `*.csv` | DELETE (regenerate) |
| `build/`, `dist/` | Keep (add to .gitignore if not already) |

### Phase 5: Clean Examples Folder

Keep only user-facing examples:
- `examples/thrombin_quickstart.py` → `examples/quickstart_thrombin.py`
- `examples/adenine_dataset.ipynb` → `examples/notebooks/`
- `examples/quinazoline_dataset.ipynb` → `examples/notebooks/`
- `examples/sEH_analysis.ipynb` → `examples/notebooks/`

DELETE or move to experiments:
- All comparison scripts
- All output files
- `test.py`, `test_marimo.py`

---

## Dataset README Template

Each dataset in `data/reagents/{dataset}/` should have a README.md:

```markdown
# {Dataset Name} Dataset

## Source
- Origin: [Paper/Database/Internal]
- Reference: [Citation if applicable]

## Reaction
- SMARTS: `[reaction SMARTS pattern]`
- Type: [Amide coupling / Suzuki / etc.]

## Components
| File | Description | Count |
|------|-------------|-------|
| acids.smi | Carboxylic acids | 150 |
| amines.smi | Primary amines | 200 |

## Library Size
- Total combinations: 30,000
- Enumerated: Yes/No (see data/libraries/{dataset}/)

## Scores
- Available: Yes/No
- Type: [Docking / ROCS / FEP / ChemProp]
- Location: `data/scores/{dataset}/`
```

---

## .gitignore Additions

```gitignore
# Outputs (regenerable)
outputs/runs/

# Large enumerated libraries (regenerable)
data/libraries/**/*.csv
data/libraries/**/*.smi
!data/libraries/.gitkeep
!data/libraries/README.md

# Build artifacts
build/
dist/
*.egg-info/
```

---

## Migration Commands

```bash
# Run from TACTICS root

# 1. Create structure
mkdir -p data/{reagents,scores,libraries}/{thrombin,seh,adenine,quinazoline}
mkdir -p experiments/{strategy_benchmarks,cats_investigation,warmup_analysis,_templates,_debug}
mkdir -p outputs/runs
mkdir -p examples/notebooks
mkdir -p scripts

# 2. Move reagent files
mv examples/input_files/acids.smi data/reagents/thrombin/
mv examples/input_files/amino_acids*.smi data/reagents/thrombin/
mv examples/input_files/coupled_aa_sub.smi data/reagents/thrombin/
cp -r examples/input_files/adenine/* data/reagents/adenine/
cp -r examples/input_files/quinazoline/* data/reagents/quinazoline/
cp -r examples/input_files/DEL_seH/* data/reagents/seh/

# 3. Move score files (after consolidation)
mv examples/docking_scores/ data/scores/thrombin/

# 4. Move experiment scripts
mv examples/rws_bayesucb_*.py experiments/strategy_benchmarks/
mv examples/strategy_*.py experiments/strategy_benchmarks/
mv examples/thrombin_strategy_comparison.py experiments/strategy_benchmarks/
mv examples/balanced_warmup_*.py experiments/warmup_analysis/
mv examples/legacy_vs_current*.py experiments/validation/

# 5. Move debug scripts
mv debug_*.py experiments/_debug/
mv test_fixed_*.py experiments/_debug/
mv test_new_*.py experiments/_debug/

# 6. Clean root outputs
rm -f *.png *.pdf *.csv  # Only after confirming these are regenerable

# 7. Move notebooks
mv examples/*.ipynb examples/notebooks/

# 8. Create gitignore for outputs
echo -e "*\n!.gitignore\n!README.md" > outputs/.gitignore
```