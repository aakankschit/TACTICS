# Changelog

All notable changes to this project are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
