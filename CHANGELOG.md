# Changelog

All notable changes to this project are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2026-09-08

A streamlining release. The package is the same algorithm with the dead
weight removed: `import TACTICS` drops from 7.1 s to 0.04 s, six runtime
dependencies are gone (about 200 MB off a fresh install), the wheel shrinks
from 13.5 MB to under 5 MB, and roughly 3,700 lines of unreachable or
duplicated code leave `src/`. Search behaviour of every preset is unchanged;
the one default that moves is documented under Changed.

### Removed

- **`thompson_sampling.legacy`** — the pre-1.0 Thompson Sampling and RWS
  implementations (17 modules). Nothing in the package or tests imported them.
  To reproduce Zhao et al. (2025) with that code, install `chem-tactics==1.2.0`.
- **`legacy_rws` preset.** Same reason. The Boltzmann-weighted posterior update
  it used is unchanged and remains the update rule of `recommended` and
  `recommended_rws`.
- **`StandardWarmup` / `StandardWarmupConfig`** — random-partner warmup, the
  weakest of the three and kept only as a comparison arm.
- **`baseline.py`** (`run_random_baseline`, `run_exhaustive_baseline`,
  `RandomBaselineConfig`) — unreachable: it read config fields that did not
  exist and passed a keyword `create_reagents` never accepted.
- **Inert config fields** that were stored but never read:
  `RouletteWheelConfig.{exploration_phase_end, transition_phase_end,
  min_observations, cats_exploration_fraction}`,
  `BayesUCBConfig.{exploration_phase_end, transition_phase_end}`,
  `TopTwoConfig.min_observations`, and
  `ThompsonSamplingConfig.max_resamples` (its early-stop branch compared
  against a counter that was never incremented).
- **`library_analysis.LibraryAnalysis`, `LibraryVisualization`,
  `compile_product_scores`, `compile_product_smiles`** and eight
  `diagnostic_plots` functions with no callers. `TS_Benchmarks` and the eight
  plot functions used by the tutorials are unchanged.
- **`library_enumeration.conformer_gen`** (never imported; unguarded OpenEye
  import), **`LibraryEnumerator`**, **`initializer`**, and a handful of
  unreferenced helpers.
- **Dependencies:** `pandas`, `dill`, `useful_rdkit_utils`, `seaborn`.
  `matplotlib` and `altair` move to the new optional `[viz]` extra.

### Changed

- **Direct `ThompsonSampler(...)` construction with no `warmup_strategy` now
  defaults to `EnhancedWarmup()`** instead of `StandardWarmup()`. This matches
  what `from_config()` and every preset already did, so preset users see no
  change.
- **Polars only.** `LookupEvaluator` reads its table with Polars; the SMARTS
  validator reads CSV reagent files with Polars. No public API accepted or
  returned pandas objects, so signatures are unchanged.
- **Bundled thrombin scores ship as Parquet** (`product_scores.parquet`,
  4.7 MB) instead of a 12 MB CSV.
- **Strategy and warmup config models reject unknown fields** (`extra="forbid"`),
  so scripts still passing a removed knob fail with a `ValidationError`
  instead of silently ignoring it.
- **`get_diagnostics()` on a strategy that records no component state returns
  an empty frame with columns `current_cycle`, `component_idx`,
  `criticality`** — the three columns every strategy-specific schema shares.
  The old 3-column `cycle`-named fallback schema was reachable by no strategy.
- `RouletteWheelSelection`, `TopTwoSelection` and `BayesUCBSelection` share
  their thermal-cycling and GMIC code through mixins in
  `strategies/_thermal.py`. Verified RNG-identical on seeded runs; the only
  visible change is that `TopTwoSelection._component_gmic` is now
  `_cached_gmics`, the name RWS already used.
- `RouletteWheelSelection.select_batch` and `BayesUCBSelection.select_batch`
  are removed; both fall through to the `SelectionStrategy` default, which the
  sampler never called anyway.

### Performance

- **Lazy package re-exports.** `TACTICS`, `TACTICS.thompson_sampling`,
  `TACTICS.thompson_sampling.core` and `TACTICS.library_analysis` resolve
  their names on first access (PEP 562). Every existing import path keeps
  working; a config-only import no longer loads RDKit, scipy or sqlitedict.
- `scipy.stats` (BayesUCB) and `matplotlib.pyplot` (diagnostic plots) are
  imported inside the function that needs them.
- `useful_rdkit_utils` — 3.85 s of the old import and ~149 MB via
  umap/pynndescent/numba/llvmlite — is replaced by direct RDKit calls that
  are bit-identical.
- `tests/test_import_time.py` pins these guarantees in fresh subprocesses.

### Fixed

- `diagnostic_plots` used `plt.cm.get_cmap`, removed in matplotlib 3.9.
- `tutorials/thompson_sampling_tutorial.py` imported a `BoltzmannConfig` that
  never existed and could not be opened.
- Three tests skipped as "warmup edge case with small test data" pass with the
  Enhanced default and are un-skipped.

## [1.2.0] - 2026-07-18

Parallel evaluation with slow evaluators (Fred docking, ROCS, ML models) was
non-functional in 1.1.0. Three defects, reported by Donald van Pinxteren on
2026-06-23, combined to make `processes > 1` unusable; all three are fixed.

### Fixed

- **Evaluation parallelism could not be enabled at all.**
  `ThompsonSampler.from_config()` hardcoded `processes=1`, so setting
  `processes` anywhere in a user script had no effect and there was no
  API-level way to turn on parallel evaluation. A 1000-iteration Fred docking
  screen ran single-threaded on a 128-core allocation with no error or
  warning. `ThompsonSamplingConfig` now exposes `processes` (and
  `min_cpds_per_core`), and `from_config()` passes them through.

- **`processes > 1` always crashed with OpenEye evaluators.**
  `ParallelEvaluator` called `pool.map(sampler.evaluate, ...)`, which pickles
  the bound method and therefore the sampler, the evaluator, and any
  SWIG-wrapped C++ object it holds. `OEDock` (Fred) and the ROCS shape engine
  raise `TypeError: cannot pickle 'SwigPyObject' object`, so the run died on
  the first batch before any docking happened.

  Workers now build their own evaluator. The pool is created with an
  `initializer` that constructs the evaluator inside each worker from its
  picklable Pydantic config, once per worker rather than once per molecule.
  The OpenEye object is never pickled and never crosses the pipe. This works
  under both `fork` and `spawn`, so it is correct on Linux, macOS and Windows
  rather than only where `fork` is the default.

  Verified end to end with a real `OEDock` under `spawn`: the evaluator is
  confirmed unpicklable, and all products still dock correctly across workers.

- **Segfault on `import TACTICS` when OpenEye and prompt_toolkit are both installed.**
  `evaluators.py` imported the OpenEye toolkits at module level, which
  initialises the global C-level `libexpat` parser. `tqdm` pulls in
  `prompt_toolkit`, whose progress-bar formatter calls
  `xml.dom.minidom.parseString()` at module level and re-enters `libexpat`
  through `pyexpat`, conflicting with OpenEye's initialisation and killing the
  interpreter (exit 139) before user code ran. OpenEye is now imported lazily,
  on first construction of an OpenEye-backed evaluator.

- **`MLClassifierEvaluator` failed when OpenEye was absent.** `joblib` was
  imported inside the OpenEye `try/except`, leaving it undefined in
  environments without OpenEye despite being unrelated to it.

- **`exp()` overflow in the CATS Boltzmann softmax.** On heavy-tailed,
  zero-inflated score landscapes a single outlier combined with a small CATS
  temperature overflowed `exp()` to `inf`, producing NaN probabilities and
  aborting the search with "probabilities contain NaN". The softmax now
  subtracts `max(z)` before exponentiating, which is shift-invariant and
  leaves the resulting probabilities unchanged.

### Added

- `ThompsonSamplingConfig.processes` and `.min_cpds_per_core` — evaluation
  parallelism, defaulting to `1` (unchanged behaviour).
- `ThompsonSampler.set_evaluator(evaluator, evaluator_config=...)` — supplies
  the picklable recipe workers use to rebuild the evaluator. `from_config()`
  provides it automatically; a warning is emitted if `processes > 1` without it.
- `LookupEvaluatorConfig.default_score` — score for product codes absent from
  the lookup table. Defaults to `None` (existing NaN behaviour); set to `0.0`
  for sparse libraries such as DEL read counts, where an unmeasured
  combination is a true non-binder rather than missing data.
- Layer-1 and layer-2 search diagnostics in `library_analysis.diagnostic_plots`
  (`plot_gmic_directed_exploration`, `plot_adaptive_intensity`,
  `plot_reagent_usage_action_panel`, `plot_gmic_vs_oracle`).

### Changed

- **TT-TS GMIC min-observation gate removed.** `TopTwoSelection._calculate_gmic`
  returned `0.0` whenever the least-observed active reagent fell below
  `min_observations`; `RouletteWheelSelection` never had this gate. On large
  components a single under-observed reagent pinned the whole component's GMIC
  to zero, distorting the rotation. Removing it raises adenine TT-TS top-100
  recovery from 87.1 to 93.2 and roughly halves replicate variance
  (sd 18.4 → 10.6). `min_observations` is still accepted for backward
  compatibility but is now inert in `TopTwoSelection`.

### Notes for users

Set `processes` to the number of cores you have allocated when using a slow
evaluator; leave it at `1` for `LookupEvaluator`/`DBEvaluator`, where process
overhead exceeds lookup cost. On macOS and Windows the default start method is
`spawn`, so entry-point scripts must be guarded with
`if __name__ == "__main__":`. Under `spawn`, and on clusters generally, the
design unit (`.oedu`) must be readable from every node — automatic on shared
filesystems (NFS, GPFS).

## [1.1.0] - earlier

See git history for releases prior to this changelog.
