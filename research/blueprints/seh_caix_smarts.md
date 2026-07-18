# SMARTS discovery for sEH and CAIX DEL benchmark libraries

Empirically derived reaction SMARTS for two DEL combinatorial libraries that have reagents
but no SMARTS entry in `data/reactions.csv`. Validated against ground-truth products where
available (sEH). All scratch artifacts in `outputs/debug/seh_smarts/` (gitignored).

Method: read `smarts_toolkit/{config,_validator}.py` + `synthesis_pipeline.py` for the real
API, censused functional groups in every BB set, deduced bond-forming chemistry by diffing
reagents against product substructure, then validated by reproducing known canonical
products (stereo-stripped on BOTH sides) and by per-position coverage.

Sampling caveat (important): `total_compounds.csv` is sorted, so any head slice contains a
tiny number of distinct bb3. Always sample with `.gather_every(N)` over a large `n_rows`,
never a head slice.

---

## sEH — VERDICT: enumerable in TACTICS = YES

**97.4% of ground-truth products reproduced exactly** with an alternative-SMARTS amide +
N-methyl-amide + methylamine-cap scheme (3898/4000, broad sample, stereo-insensitive;
`outputs/debug/seh_smarts/unified.json`). The 2.6% residual is almost entirely
amino-acid-bb1 cyclization/lactam edge cases (`(bb1=both, bb2=acid, bb3=amine)` = 91 of the
102 misses); the 270 non-amine/non-acid bb3 (sulfonyl chlorides, isocyanates, aryl/heteroaryl
halides) are a separate, real chemistry gap that an amide-only library simply excludes.

### Library chemistry (deduced + ground-truth-confirmed)

sEH is a **3-cycle amide-coupling DEL**: a central di-acid scaffold (bb2) bearing two
building blocks (bb1, bb3) as amides, with leftover acids capped by methylamine.

Full distinct-BB functional-group census, streamed over all 4,433,801 rows
(`outputs/debug/seh_smarts/real_bb_census.txt`):

- **bb2 (338) = central scaffold, ALL carboxylic acid**: 124 acid-only, 214 amino-acid, **110 di-carboxylic acids**. The di-acids form the bis-amide core.
- **bb1 (617) = ALL have a primary/secondary amine**: 291 are amino acids (amine + acid), 326 amine-only, 0 acid-only. bb1 always contributes its amine.
- **bb3 (4529)**: 2363 amine, 1897 acid, 1 both, **270 "neither"**. The 270 "neither" are **sulfonyl chlorides, isocyanates, and aryl/heteroaryl halides** (163 carry a C–halide; the rest are mostly ArSO2Cl / R–N=C=O), confirmed in `outputs/debug/seh_smarts/bb3_neither.json` — NOT Boc-protected amines. These react by sulfonamide formation / urea / SNAr, i.e. chemistry outside the amide scheme, and are the main coverage gap.

Bonds and the N-methyl: bb1's amine and bb3's amine each form an amide to a bb2 acid (when
bb3 is itself an acid, its acid couples to a free amine instead). The N-methyl seen in ~73%
of products has **two distinct origins**, both required to reach 95.8%:
1. **N-methylation of a coupling amide nitrogen** → tertiary `CN(R)C(=O)` (the common case when bb1/bb3 are simple amines).
2. **Methylamine capping of a leftover free carboxylic acid** → `C(=O)NHCH3` (e.g. an amino-acid bb1's own acid, or an unconsumed bb2 acid).

A unified assembler with amide + N-methyl-amide + methylamine-cap (both coupling polarities,
0–2 caps) reproduces **97.4% exactly** (3898/4000). Failure taxonomy of the 2.6% residual
(`outputs/debug/seh_smarts/unified.json`): 91 `(bb1=both, bb2=acid, bb3=amine)` — amino-acid
bb1 cyclization/lactam edge cases — and 11 `(bb1=amine, bb2=acid, bb3=amine)` di-acid stereo
edge cases. The 270 non-amine/non-acid bb3 (ArSO2Cl / R-N=C=O / Ar-halide) are excluded by an
amide-only scheme and are NOT counted here (they appear as `bb3=neither`, a separate gap).

### Best reaction SMARTS

Amine filter used throughout (excludes amides, hydrazines, N–heteroatom amines):
`#7X3;H2,H1;!$(N=*);!$(N-[!#6]);!$(NC=O)`

| role | SMARTS |
|------|--------|
| amide coupling | `[CX3:1](=[OX1])[OX2H1].[#7X3;H2,H1;!$(N=*);!$(N-[!#6]);!$(NC=O):2]>>[CX3:1](=O)[#7:2]` |
| N-methyl amide coupling | `[CX3:1](=[OX1])[OX2H1].[#7X3;H2,H1;!$(N=*);!$(N-[!#6]);!$(NC=O):2]>>[CX3:1](=O)[#7:2][CH3]` |
| methylamine cap of free acid | `[CX3:1](=[OX1])[OX2H1]>>[CX3:1](=O)NC` |

(All four SMARTS in this report — these three plus the CAIX reductive-amination and Fmoc —
parse cleanly under `AllChem.ReactionFromSmarts`; verified.)

### Per-position coverage (ValidationResult.coverage_stats), single amide config

- bb2 acid template: **100.0%** (all 338 scaffolds have a couplable acid)
- bb1 amine template: **100.0%** (all 617 bb1 have an amine)
- bb3 amine template: **55.0%** — EXPECTED, not a defect: only 2363/4529 bb3 are amines; 1897 are acids (couple via the acid template) and 270 are sulfonyl-Cl/halide (need different chemistry).

Full bb3 coverage needs **alternative SMARTS** per step (amine-coupling OR acid-coupling),
via `step_modes={step: "alternative"}`. The 270 non-amine/non-acid bb3 are out of scope for
an amide-only library (drop them, or add sulfonamide/urea/SNAr patterns to include them).

### Valid vs dropped reagents/products

- Valid reagents: bb2 100% (338/338, all acid); bb1 100% (617/617, all amine); bb3 ~94% directly couplable (2363 amine + 1897 acid); 270 (~6%) bb3 are non-amide chemistry (sulfonyl-Cl / halide). 0 invalid SMILES in any list.
- Products reproduced exactly: 97.4% with the unified amide/N-methyl-amide/cap scheme; the 2.6% miss is amino-acid-bb1 lactam/cyclization + di-acid stereo edge cases. (The 270 non-amide bb3 are a separate, deliberately-excluded class, not part of this 2.6%.)
- A naive SINGLE fixed-polarity 3-step config (`run_single_reaction` emits only `products[0][0]`) reproduces ~62% exact while still producing a *valid* member for ~95% of combinations (`outputs/debug/seh_smarts/{modelB_clean,definitive}.json`); the gap is di-acid regiochemistry that the alternative-SMARTS config closes.

### Copy-pasteable ReactionConfig (practical 3-step; ~95% valid products)

```python
from TACTICS.library_enumeration import SynthesisPipeline
from TACTICS.library_enumeration.smarts_toolkit import (
    ReactionConfig, ReactionDef, StepInput, InputSource,
)

AMINE = "#7X3;H2,H1;!$(N=*);!$(N-[!#6]);!$(NC=O)"
AMIDE = f"[CX3:1](=[OX1])[OX2H1].[{AMINE}:2]>>[CX3:1](=O)[#7:2]"
CAP   = "[CX3:1](=[OX1])[OX2H1]>>[CX3:1](=O)NC"   # methylamine cap

# Reagent order: [bb2 scaffold, bb1, bb3]
config = ReactionConfig(
    reactions=[
        ReactionDef(reaction_smarts=AMIDE, step_index=0),   # bb2 acid + bb1 amine (acid-template-first; inputs are [bb2 acid, bb1 amine] — OK)
        # step 1 inputs are [PREVIOUS_STEP product (carries the free acid), bb3 amine] -> acid-first matches:
        ReactionDef(reaction_smarts=AMIDE, step_index=1),   # remaining acid + bb3 amine
        ReactionDef(reaction_smarts=CAP,   step_index=2),   # cap any leftover free acid
    ],
    reagent_file_list=[
        "data/reagents/seh/DEL_seH/bb2_list.csv",
        "data/reagents/seh/DEL_seH/bb1_list.csv",
        "data/reagents/seh/DEL_seH/bb3_list.csv",
    ],
    step_inputs={
        0: [StepInput(source=InputSource.REAGENT_FILE, file_index=0),
            StepInput(source=InputSource.REAGENT_FILE, file_index=1)],
        1: [StepInput(source=InputSource.PREVIOUS_STEP, step_index=0),
            StepInput(source=InputSource.REAGENT_FILE, file_index=2)],
        2: [StepInput(source=InputSource.PREVIOUS_STEP, step_index=1)],
    },
)
pipeline = SynthesisPipeline(config)
```

To approach the 97.4% exact ceiling, add alternative SMARTS at steps 0/1 so each step tries
amide AND N-methyl-amide in both polarities (`step_modes={0:"alternative", 1:"alternative"}`).
The N-methyl-amide pattern is the same AMIDE with `>>[CX3:1](=O)[#7:2][CH3]`.

### `reactions.csv` entry

The library is multi-step + cap, so it does not collapse to one row. If a single
representative SMARTS is required for the lookup table, use the amide:
```
seh,"[CX3:1](=[OX1])[OX2H1].[#7X3;H2,H1;!$(N=*);!$(N-[!#6]);!$(NC=O):2]>>[CX3:1](=O)[#7:2]"
```
Otherwise register the multi-step config above directly in the benchmark runner.

### Reagent-file note

Use **`bb1_list.csv` / `bb2_list.csv` / `bb3_list.csv`** (column `SMILES`, auto-detected by
the validator). Their counts (617 / 338 / 4529) and full functional-group census match the
distinct `bb1/bb2/bb3` columns of `total_compounds.csv` exactly — they are the true reagent
sets.

### Expected enumerated product count

Full combinatorial size = 617 × 338 × 4529 ≈ **9.44e8** (~944M). Ground-truth
`total_compounds.csv` has 4,433,801 measured members. Enumerate on demand — Thompson
Sampling never materializes the full space.

---

## CAIX — VERDICT: enumerable in TACTICS = YES (2-step enum 100% = 400/400 verified; Fmoc removal post-enumeration)

No ground-truth product file, so confidence rests on coverage + chemical reasoning, but the
BB roles are unambiguous from the functional-group census.

### Library chemistry (deduced)

3-component reductive amination + amide, with Fmoc removal:

- **AL (67) = aromatic aldehydes** (100% have `[CX3H1]=O`; no acids/amines)
- **AM (52) = primary amines** (100% primary NH2)
- **AA (62) = Fmoc-protected amino acids** (100% carboxylic acid, 100% Fmoc carbamate)

Sequence: (1) reductive amination AL + AM → secondary amine; (2) amide coupling of AA's acid
onto that secondary amine; (3) Fmoc removal to expose the warhead amine (the CAIX
zinc-binding handle).

### Best SMARTS

| role | SMARTS |
|------|--------|
| reductive amination (aldehyde + amine → CH2–N) | `[#6:4][CX3H1:1]=[OX1].[#7X3;H2,H1;!$(N=*);!$(N-[!#6]);!$(NC=O):3]>>[#6:4][CH2:1][#7:3]` |
| amide coupling (acid + amine → amide) | `[CX3:1](=[OX1])[OX2H1].[#7X3;H2,H1;!$(N=*);!$(N-[!#6]);!$(NC=O):2]>>[CX3:1](=O)[#7:2]` |
| Fmoc removal (carbamate → free amine) | `[N:1]C(=O)OCC1c2ccccc2-c2ccccc21>>[N:1]` |

### Coverage + enumeration

- reductive amination: aldehyde template 100.0% (AL), amine template 100.0% (AM)
- amide: AA acid template 100.0%
- **2-step pipeline enum success (reductive amination + amide, no Fmoc step): 100% = 400/400** (`outputs/debug/seh_smarts/caix_fixed.json`)
- Fmoc removal SMARTS verified standalone on a real product: cleanly converts the carbamate → free amine (`outputs/debug/seh_smarts/{caix_defmoc_test,caix_fixed}.json`).

### CRITICAL config gotcha — reactant-template order MUST match step_inputs order

`enumerate_single` feeds the inputs of a step positionally to the reaction's reactant
templates (input[0] → template[0], input[1] → template[1]). It does NOT reorder by
functional group. So if a step's `step_inputs` are `[PREVIOUS_STEP (an amine), REAGENT_FILE
(an acid)]`, the step-1 SMARTS must be written **amine-template-first**:
`[{AMINE}:2].[CX3:1](=[OX1])[OX2H1]>>[CX3:1](=O)[#7:2]`. Writing it acid-first silently
yields zero products (verified: acid-first = 0/300, amine-first = 400/400;
`outputs/debug/seh_smarts/{caix_debug2,caix_fixed}.json`). The same rule applies to the sEH
config above — keep each step's template order aligned with its `step_inputs`.

### Fmoc removal — do it as a POST-enumeration step

In-pipeline Fmoc removal via `DeprotectionSpec(group="Fmoc", target="product")` with a custom
`ProtectingGroupInfo` did not fire in testing. **Working path:** run the 2-step config
(reductive amination + amide → 100% products), then apply the Fmoc-removal reaction SMARTS
`[N:1]C(=O)OCC1c2ccccc2-c2ccccc21>>[N:1]` to each product as a one-line post-processing pass
(verified: deprotects a real product to the free amine, `caix_fixed.json`). If leaving Fmoc
on is acceptable for benchmark scoring, the post-step is optional. Either way CAIX enumerates
fully.

### Reagent-file note

CAIX `.smi` files have a 3-column `BB,cas,SMILES` header. The validator's CSV loader
auto-detects the `SMILES` column, but the **`.smi` loader splits on whitespace and would
take `BB` as the SMILES**. Either pass them as `.csv` (rename/symlink), or pre-write plain
1-column `.smi` (as done: `outputs/debug/seh_smarts/caix_{AL,AM,AA}.smi`).

### Expected enumerated product count

67 × 52 × 62 = **216,008** products — small enough to enumerate exhaustively if desired.

### Copy-pasteable ReactionConfig (CAIX, 2-step; apply Fmoc removal post-hoc)

```python
AMINE = "#7X3;H2,H1;!$(N=*);!$(N-[!#6]);!$(NC=O)"
config = ReactionConfig(
    reactions=[
        ReactionDef(reaction_smarts=f"[#6:4][CX3H1:1]=[OX1].[{AMINE}:3]>>[#6:4][CH2:1][#7:3]", step_index=0),
        # step-1 amide written AMINE-FIRST to match step_inputs order [PREVIOUS_STEP amine, AA acid]:
        ReactionDef(reaction_smarts=f"[{AMINE}:2].[CX3:1](=[OX1])[OX2H1]>>[CX3:1](=O)[#7:2]", step_index=1),
    ],
    reagent_file_list=[
        "data/reagents/CAIX/AL_SEL2.smi",   # aldehyde  (pass as .csv OR clean .smi)
        "data/reagents/CAIX/AM_SEL2.smi",   # amine
        "data/reagents/CAIX/AA_SEL2.smi",   # Fmoc amino acid
    ],
    step_inputs={
        0: [StepInput(source=InputSource.REAGENT_FILE, file_index=0),
            StepInput(source=InputSource.REAGENT_FILE, file_index=1)],
        1: [StepInput(source=InputSource.PREVIOUS_STEP, step_index=0),
            StepInput(source=InputSource.REAGENT_FILE, file_index=2)],
    },
)
# Post-enumeration Fmoc removal (per product mol):
#   from rdkit.Chem import AllChem
#   defmoc = AllChem.ReactionFromSmarts("[N:1]C(=O)OCC1c2ccccc2-c2ccccc21>>[N:1]")
#   ps = defmoc.RunReactants((product_mol,)); free = ps[0][0] if ps else product_mol
```

---

## DEL objective note (read_count vs docking)

Both libraries carry an experimental **DEL read_count** enrichment signal (sEH:
`total_compounds.csv` `read_count`, plus `del_binders.txt` / `del_nonbinders.csv`; 4.43M
measured members). For a TACTICS benchmark the natural objective is the **DEL read_count
itself** — a real, noisy, experimentally-measured target that exercises Thompson Sampling on
authentic assay signal and is the cheapest LookupEvaluator path. Docking the enumerated
products is a viable alternative (and lets sEH/CAIX join the existing docking-benchmark set),
but read_count is the more faithful DEL task and needs no new computation. Recommendation:
read_count as primary objective; docking as an optional second objective for cross-signal
studies.
