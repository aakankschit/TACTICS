# Handoff — Docking-Benchmark Expansion (sEH + KinDEL + JNK3-ROCS)

**Date:** 2026-05-30
**Purpose:** Self-contained execution plan so this work can resume after a context clear.
**Scope decision (user, 2026-05-30):** Expand the Paper-1 docking benchmark beyond the
current 4 base chemistries (adenine, amide, quinazoline, thrombin) by adding new
combinatorial docking datasets. Headline new case study = **sEH (soluble epoxide
hydrolase)**, our own lab's DEL data, with an experimental-binder overlay. Move fast.

---

## 1. The three new dataset streams

| stream | chemistry / comp | size target | target/score | status |
|---|---|---|---|---|
| **sEH** (centerpiece) | 3-comp amide-coupling DEL | **~60M** | dock vs sEH receptor (PDB **4OD0**, TBD-confirm) + DEL read_count overlay | SMARTS recovered; needs reagent reconciliation + enumeration + docking |
| **KinDEL MAPK14 + DDR1** | 3-comp trisynthon | subset to 1–10M (81M full) | insitro docked poses / Kd + (optional) our docking | usable under license; FETCH not vendor |
| **ROCS-libs-on-JNK3** | reuse enumerable ROCS libs (1×2-comp + 1×3-comp) | native | dock vs JNK3 2ZDT (already prepped) | lowest priority; SMARTS already in reactions.csv |

User framing for JNK3 reuse (important for manuscript): the paper's claim is **"find top
hits without screening the whole library,"** treating reference scores as ground truth *by
definition of the benchmark* — NOT claiming docking is physical truth (docking ranks; FEP
would be truth). So generating docking data on JNK3 for ROCS libraries to get varied
reagent-scoring landscapes is legitimate; we are not asserting JNK3 biological relevance.

---

## 2. STANDARDIZED FORMAT (all new datasets MUST match this exactly)

Confirmed by inspecting `data/reagents/adenine/` + `data/scores/adenine.parquet` on 2026-05-30.

### 2a. Reagent files
- Path: `data/reagents/<lib_id>/<name_key>_reagent_<i>.smi` for component i = 0,1,2,…
- Format: **whitespace-delimited, NO header**, one reagent per line: `<SMILES> <reagent_name>`
  - e.g. `C1=CC(=C(N=C1)N)CO amidine_001`
- `<reagent_name>` must be a **stable canonical ID** (the product code is built from these).

### 2b. Score / lookup parquet
- Path: `data/scores/<lib_id>.parquet`
- Columns: **`Product_Code` (str)`, `Scores` (f64)`** — exactly these two names.
- `Product_Code` = `"_".join(reagent_name for each component in order)`
  - e.g. `amidine_003_isocyanide_db_472_aldehyde_335`
- This is the LookupEvaluator key. **CRITICAL:** the sampler keys on reagent *names*
  joined by `_`; it never re-derives structure during search (verified in
  `core/sampler.py:322` + `core/evaluators.py:159`). So the parquet `Product_Code` MUST
  equal `_`.join of the reagent names in the `.smi` files, same component order. Any
  mismatch silently yields all-NaN lookups.
- Optimization mode: docking = **minimize** (FRED/HYBRID Chemgauss4, lower = better);
  DEL read_count = **maximize**.

### 2c. reactions.csv
- Path: `data/reactions.csv`, header `rxn_id,SMARTS`.
- Add one row keyed by the SMARTS lib key (see Local Norm #11 below). The SMARTS must
  **parse under `AllChem.ReactionFromSmarts`** but — for a pure lookup benchmark — is never
  applied during search (it only needs to be syntactically valid). Use a real SMARTS if the
  library is also enumerated in-pipeline.

### 2d. Triple-config registration (Local Norm #12 — edit ALL THREE or recovery silently breaks)
For every new **docking (minimize)** library, add the lib_id to:
1. `examples/_benchmark_config.py` → `MINIMIZE_LIBRARIES` (set, line ~25)
2. `examples/run_posthoc_analysis.py` → `MINIMIZE_LIBRARIES` (set, line ~52)
3. `examples/_benchmark_config.py` → `SCORE_TO_REAGENT_FOLDER` (line ~31) **if** folder name ≠ lib_id
4. `scripts/stage_benchmark_data.sh` → `REAGENT_FOLDER` map (line ~42) **if** folder name ≠ lib_id

Symptom of a miss: posthoc recovery ≈ 0 (lib treated as MAXIMIZE → recovery vs worst
compounds), or a Slurm stage failure for missing reagent files.

### 2e. Name-key resolution (Local Norm #11)
`run_benchmark.py`: `folder_name = SCORE_TO_REAGENT_FOLDER.get(lib_id, lib_id)`;
`name_key = folder_name if lib_id.endswith("_hybrid") else lib_id`. Reagent files at
`data/reagents/{folder_name}/{name_key}_reagent_*.smi`; SMARTS key in reactions.csv =
`name_key`. Do NOT add ad-hoc fallbacks; preserve the `endswith("_hybrid")` discriminator.

---

## 3. sEH — the centerpiece (full recipe)

### 3a. Source & provenance (no license friction — our own lab's data)
- Paper: Zhang et al., *J. Chem. Inf. Model.* 2023, 63(16):5120, DOI 10.1021/acs.jcim.3c00588
  ("Building Block-Based Binding Predictions for DNA-Encoded Libraries"). Mobley co-author.
  Open access (PMC10466377). Code/data: github.com/MobleyLab/DEL_analysis.
- Local files (already on disk): `data/reagents/seh/DEL_seH/`
  - `bb1_list.csv`, `bb2_list.csv`, `bb3_list.csv` (columns `SMILES,iso_SMILES`; 617/338/4529)
  - `del_binders.txt` (columns `bb1,bb2,bb3,structure,read_count`; **116,666 binders**, read_count 81–4712)
  - `del_nonbinders.csv` (1.6 GB, local only), `total_compounds.csv` (962 MB, 4.43M rows, local only)

### 3b. Binder definition (the experimental ground truth)
Binary binder vs non-binder: NGS read counts statistically ≠ 0 at 95% (Poisson). Min
read_count for a binder = 81. We will overlay these as "experimental hits" to show TACTICS
recovers them. **Binder labels live in `del_binders.txt`.**

### 3c. Chemistry / SMARTS (RECOVERED — see `research/blueprints/seh_caix_smarts.md`)
3-cycle amide-coupling DEL: central di-acid scaffold (bb2) + amine bb1 + amine/acid bb3,
leftover acids methylamine-capped. The recovered multi-step `ReactionConfig` reproduces
**97.4% of ground-truth products exactly**. The ~270 non-amide bb3 (sulfonyl-Cl /
isocyanate / aryl-halide) are out of scope for the amide scheme (a deliberate, accepted
drop). Blueprint has copy-pasteable `ReactionConfig` + per-position coverage + the
reactant-template-order gotcha. **NOTE:** for a pure docking-lookup benchmark we may not
even need to enumerate via TACTICS — we can enumerate products directly (RDKit) for docking,
then score by Product_Code. SMARTS only needed if enumerating in-pipeline.

### 3d. Binder-retention vs library-size (computed 2026-05-30 from del_binders.txt)
Amide-compatible binders = **109,663 (94.0%** of 116,666; the 5.9% lost all have a
non-amide bb3 — uniform across potency, only 51 of top-1000 lost).
Reagent union of amide-compatible binders (canonical): **bb1=606, bb2=388, bb3=1111** →
full Cartesian to contain ALL = ~261M (too big).

**DECISION: 60M.** Recipe = keep all binder bb1 × bb2, keep **top-K bb3 by #binders**:

| lib size | K bb3 | binders kept | % amide-compat | % all 116,666 |
|---|---|---|---|---|
| 40M | 170 | 99,753 | 91% | 85.5% |
| **60M** | **~255** | **104,207** | **95%** | **89%** |
| 80M | 340 | 106,160 | 97% | 91% |
| 261M | 1111 | 109,663 | 100% | 94% |

606×388×255 = **59.96M ≈ 60M, retains 95% of amide-compatible binders.** Can co-optimize
(trim bb1/bb2 too) to squeeze a few % more at fixed size if desired.

### 3e. ⚠️ OPEN RISK — name reconciliation (do this FIRST, blocks everything)
The binder file's BB columns do NOT match the master reagent lists 1:1:
- binders use **665 unique bb1** strings but `bb1_list.csv` has 617; **bb2: 666 vs 338**;
  canonical: binder bb2 union = 388 but master = 338.
This is canonicalization/deprotection-form drift (master lists are post-deprotection; binder
columns may be a different form). **Implication:** some binder reagents are NOT in the master
lists, so "include the binders" requires building the reagent set from the binder requirement
+ fill, NOT just subsetting the master lists. Step 0 = reconcile to one canonical ID scheme,
then define library reagent sets = (amide-compatible binder-used reagents) ∪ (top master
reagents to reach 60M). This guarantees binder inclusion and a consistent Product_Code space.

### 3f. Receptor (CONFIRM before docking)
No single dominant sEH docking structure. Recommended for OpenEye FRED/HYBRID (needs a holo
structure with a drug-like inhibitor to define the design unit — apo 1S8O won't do):
- **4OD0** (lead pick) — human sEH + urea inhibitor, canonical sEH pharmacophore (catalytic
  Tyr/Tyr/Asp triad), clean single bound ligand. RCSB 4OD0.
- alts: **3PDC** (used for SBVS of combinatorial libs vs sEH), **3I28** (lead-opt costructure),
  **3WKE** (t-AUCB tool inhibitor).
User has no in-house sEH receptor → use a PDB. **Pending: confirm 4OD0 preps cleanly** (pull
ligand/pocket, build .oedu design unit). Mobley lab may have a preferred structure — ask.

### 3g. Objective choice (design flag, doesn't change construction)
- Docking score objective → TACTICS recovers top *docking* hits; binder overlay then measures
  docking's enrichment, not TACTICS skill.
- **DEL read_count objective → the sharpest "TACTICS finds experimental binders in few
  iterations" claim, no docking needed.**
- Recommend doing BOTH: dock for the docking-recovery benchmark; read_count variant for the
  binder claim. (`maximize` mode for read_count.)

---

## 4. KinDEL (MAPK14 + DDR1) — license VERIFIED, FETCH not vendor

- Verified by reading `raw.githubusercontent.com/insitro/kindel/main/LICENSE.md` (2026-05-30):
  custom Insitro license **incorporating CC BY-NC-SA 4.0**. **"benchmarking" and "performance
  evaluation" are explicitly permitted, for-profit or not.** Prohibits commercial use.
  Grant is **non-transferable + ShareAlike + no patent license**.
- **Plan (user-approved):** academic paper, open-source code, data **NOT in package
  distribution**. Datasets go in the separate manuscript-reproduction repo, fetched via a
  download script pointing at the public S3 `s3://kin-del-2024/data` + citation (same pattern
  as `download_scores.py`). **Do NOT re-host their data** (non-transferable/SA) — fetch from
  origin. No commercial use.
- Data ships per-cycle reacted synthons (378/1128/191 → 81M) + prepared receptor PDBs
  (MAPK14: 3KQ7/3S3I/5WJJ/5XYY/6SFI; DDR1: 6FEX) + docked poses + off-DNA Kd held-out sets.
  Enumerable via fragment-stitch on `[Fr]` markers or sample the provided enumerated parquet;
  subset 81M → 1–10M.
- Manuscript angle (user): "deliver benchmarking docking datasets for developing active-learning
  methods."

---

## 5. Recovered sub-agent artifacts (don't redo)
- `research/blueprints/seh_caix_smarts.md` — full sEH (97.4%) + CAIX SMARTS recipes, coverage,
  ReactionConfig snippets, gotchas. (The SMARTS sub-agent hit the rate limit and returned an
  empty final message, but its report was written to disk — recovered.)
- `outputs/debug/seh_smarts/` — scratch JSONs/SMI (gitignored). `FINAL_SUMMARY.txt` has the
  headline numbers.
- **CAIX is too small** (216K products, below the 1M floor) — not a benchmark candidate; skip.
- Sulfonamide / carbonic-anhydrase candidate **DROPPED** (catalytic-Zn docking too noisy).

---

## 6. ORDERED NEXT STEPS (resume here)

0. **(sEH) Reconcile names** — canonicalize bb1/bb2/bb3 across `del_binders.txt` ↔ master
   `bb*_list.csv`; pick a stable reagent-ID scheme; resolve the 665-vs-617 / 388-vs-338 drift.
   Output: a mapping {canonical reagent → stable ID} per component.
1. **(sEH) Select the 60M reagent set** — all binder bb1×bb2 (amide-compatible) + top-~255 bb3
   by binder usage, ∪ any binder-required reagents missing from master. Verify final
   |bb1|×|bb2|×|bb3| ≈ 60M and that 95% of amide-compatible binders are inside.
2. **(sEH) Write standardized reagent files** `data/reagents/seh/seh_reagent_{0,1,2}.smi`
   (SMILES + stable ID, no header). Emit a `binder_map` (Product_Code → read_count) for the
   recovery overlay.
3. **(sEH) Confirm receptor 4OD0** (or lab structure) → build `.oedu` design unit.
4. **(sEH) Enumerate + dock on HPC3** (omega + FRED/HYBRID) → `data/scores/seh.parquet`
   (`Product_Code`, `Scores`), Product_Code matching the reagent IDs exactly.
5. **(sEH) Register** in reactions.csv + triple-config (§2c–2d). Add to benchmark runner.
6. **(sEH) Run benchmark + overlay** — report docking recovery AND experimental-binder hit
   rate (and/or the read_count-objective variant).
7. **(KinDEL)** stand up download script + enumerable subset (1–10M) in the separate repo.
8. **(JNK3-ROCS, lowest priority)** dock 1×2-comp + 1×3-comp ROCS libs vs JNK3 2ZDT.

**Do not** kick off enumeration/docking/HPC jobs without explicit user go-ahead.

---

## 7. This-session figure work (already done, for completeness)
`tutorials/manuscript_plots_docking.py` §6 budget-sensitivity panel: bottom row converted to a
per-library **stacked decomposition** (per-chemistry segments sum to the bin mean), bar colors
moved to the **unused Tableau-10 members** (adenine green #59A14F, quinazoline pink #FF9DA7,
amide brown #9C755F, thrombin gold #EDC949) so they don't collide with the method-line colors;
legend title dropped (was the white-gap cause). Standalone per-library-exemplars figure removed
per user. Compiles; pixel-verified.
