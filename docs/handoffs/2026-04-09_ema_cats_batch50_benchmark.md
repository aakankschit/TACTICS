# EMA CATS + Batch Size Benchmark Results (2026-04-07 to 2026-04-09)

## Summary

Three full benchmarks were run to evaluate the rotation-evaluation ordering fix, EMA-smoothed CATS multiplier, and batch_size=50 hypothesis. The EMA fix improves TT-TS reliability on hard docking libraries (variance reduction) without affecting aggregate recovery. Batch_size=50 is a negative result — GMIC is batch-size invariant.

## Benchmarks Run

| Benchmark | Config | Output |
|-----------|--------|--------|
| `full_benchmark` (baseline) | batch=100, no EMA | `outputs/full_benchmark/` |
| `full_benchmark_ema_cats` | batch=100, EMA decay=0.1, rotation fix | `outputs/full_benchmark_ema_cats/` |
| `full_benchmark_batch50` | batch=50, EMA decay=0.1, rotation fix | `outputs/full_benchmark_batch50/` |

All benchmarks: 218,400 trials, 24 libraries, 5 methods, 20 replicates, 10 top-N thresholds.

## Aggregate Recovery (Top-100)

| Method | batch=100 (original) | batch=100 + EMA | batch=50 + EMA |
|--------|---------------------|-----------------|----------------|
| TACTICS Enhanced-TT-TS (GMIC) | 86.5 +/- 15.6 | 86.4 +/- 15.7 | 86.3 +/- 15.8 |
| TACTICS Enhanced-RWS (GMIC) | 85.8 +/- 15.7 | 85.9 +/- 15.5 | 85.7 +/- 15.9 |
| Legacy Enhanced-RWS | 82.9 +/- 15.6 | 82.9 +/- 15.6 | 82.8 +/- 15.7 |
| TACTICS Balanced-Greedy | 82.0 +/- 19.0 | 82.0 +/- 19.0 | 82.1 +/- 18.9 |
| Legacy Standard-Greedy | 81.2 +/- 19.6 | 81.2 +/- 19.7 | 81.1 +/- 19.7 |

**Conclusion**: No meaningful change in aggregate recovery across any configuration.

## Docking Library Results (where EMA fix had visible effect)

### Adenine (hardest SAR, 3-comp, 17% docking failures)

| Method | Original | + EMA | + batch=50 |
|--------|----------|-------|------------|
| TT-TS | 93.9% +/- 13.9% | **96.0% +/- 4.0%** | **96.8% +/- 1.1%** |
| RWS | 71.8% +/- 27.6% | 78.5% +/- 24.0% | 72.8% +/- 27.0% |
| Legacy RWS | 40.1% +/- 23.8% | 45.1% +/- 27.1% | 45.9% +/- 27.0% |
| Balanced-Greedy | 31.6% +/- 0.8% | 31.5% +/- 0.9% | 31.6% +/- 1.1% |
| Legacy Greedy | 31.8% +/- 1.1% | 31.4% +/- 1.5% | 31.0% +/- 2.4% |

TT-TS std progression: **13.9% -> 4.0% -> 1.1%** across the three benchmarks.

### Thrombin (imbalanced 130x3844, 3-comp)

| Method | Original | + EMA | + batch=50 |
|--------|----------|-------|------------|
| TT-TS | 81.0% +/- 2.2% | **81.8% +/- 1.4%** | 81.8% +/- 1.8% |
| RWS | 80.5% +/- 2.4% | 80.2% +/- 2.9% | **81.2% +/- 2.2%** |

### Amide and Quinazoline

Flat across all configurations — near ceiling (amide ~97%) or greedy-reachable (quinazoline ~85%).

## ROCS Libraries

No change. All 20 ROCS libraries show deltas within +/-0.6% recovery and +/-0.5% std across all three benchmarks. The EMA fix and batch size change are neutral on ROCS.

## GMIC Stability Analysis

Diagnostic data confirms GMIC is stable from cycle 0:

| Library | Stabilization cycle | Mean cycle-to-cycle delta |
|---------|-------------------|--------------------------|
| adenine | 0-1 | 0.001-0.003 |
| mannich | 0-1 | 0.002-0.021 |
| orru | 0-13 | 0.003-0.029 |
| rxn206 | 0-6 | 0.001-0.018 |
| quinazoline | 0-10 | 0.0002-0.002 |

GMIC is computed from cumulative posteriors, not the current batch. Batch size does not affect signal quality.

## TT-TS Disagreement Analysis

Most components spend majority of time in saturated region (EMA > 0.8):

| Library | Component | % Saturated (>0.8) | % Neutral | % Resolved (<0.3) |
|---------|-----------|--------------------|-----------|--------------------|
| adenine | comp 0 | 22% | 56% | 22% |
| adenine | comp 1 | 95% | 5% | 0% |
| quinazoline | comp 0 | 94% | 6% | 1% |
| orru | comp 2 | 0-6% | 40-66% | 34-56% |

Only orru comp 2 consistently enters the resolved region, triggering heated_scale inflation (3.5-4.0x). The adaptive mechanism is primarily unidirectional (decaying heated_scale from saturation).

## EMA Relative GMIC (RWS Diagnostic)

With `cats_ema_decay=0.1`, the EMA relative GMIC provides stable directional signal for CATS multiplier:

| Library | Component | Raw GMIC | EMA Relative GMIC | CATS Multiplier |
|---------|-----------|----------|-------------------|-----------------|
| adenine | comp 0 (critical) | 2.00 | 1.376 | 0.860 (cooled) |
| adenine | comp 1 (flexible) | 0.90 | 0.618 | 1.164 (heated) |
| adenine | comp 2 (neutral) | 1.46 | 1.005 | 0.997 (neutral) |
| quinazoline | all | 2.0-2.4 | 0.91-1.06 | ~1.0 (flat) |

The EMA smooths out per-cycle noise, giving the CATS multiplier a consistent directional signal. Previously (without EMA), the multiplier oscillated randomly.

## Statistical Significance

TACTICS methods significantly outperform Legacy on 16/24 libraries (all pairwise comparisons p < 0.05). On the remaining 8, differences are negligible effect size — these are libraries where greedy methods already perform well (easy SAR). Full Tukey HSD results in `outputs/full_benchmark_batch50/analysis/tukey_hsd_results.parquet`.

## Dead End: batch_size=50

**Hypothesis**: Smaller batches -> more frequent GMIC rotations -> better adaptation -> improved recovery.

**Result**: No improvement. GMIC is batch-size invariant because it's computed from cumulative posteriors. Smaller batches double the number of rotation decisions but provide the same information per decision.

**What this rules out**: Any batch_size reduction as a tuning lever for GMIC-based methods. The cumulative posterior computation makes GMIC insensitive to batch granularity. Halving, quartering, or any other reduction will not change the rotation signal quality. The only effect is increased per-cycle overhead.

## Infrastructure Changes

- Moved 21 `.sub` files from repo root to `jobs/`, created `jobs/MANIFEST.md`
- Fixed merge job time limit: 1h -> 8h (was timing out on large libraries)
- Fixed `pl.concat()` schema mismatch in merge: added `how="diagonal"` for diagnostics tier (TT-TS has 16 cols, CATS has 27)
- Fixed posthoc job: GPU partition -> standard (CPU-only work)
- Fixed diagnostic benchmark: added missing `cats_ema_decay=0.1` in `RouletteWheelConfig`
- Updated diagnostic plots notebook to save PDFs to `plots/` directory

## Recommendation

- **Keep batch_size=100** as default for TACTICS methods
- **Keep `cats_ema_decay=0.1`** as default for RWS — provides smoother CATS multiplier
- **Keep rotation-evaluation ordering fix** — improves reliability on hard libraries
- **Revert `_benchmark_config.py`** back to batch_size=100 before next benchmark

## Heptabase Cards Needed

- [DEAD END] batch_size=50 does not improve recovery (marked PENDING in STATUS.md)
- [RESULT] EMA CATS + rotation fix reliability improvement (optional — finding is in STATUS.md)
