# TACTICS

<p align="center">
  <img src="https://raw.githubusercontent.com/aakankschit/TACTICS/main/docs/source/_static/images/TACTICS_logo.png" alt="TACTICS Logo" width="600">
</p>

**Find the best compounds in a combinatorial library without enumerating all of them.**

TACTICS (Thompson Sampling-Assisted Chemical Targeting and Iterative Compound
Selection) uses Thompson Sampling to search ultra-large combinatorial libraries,
scoring only a small fraction of the products while still recovering most of the
top hits. A typical run evaluates **1–2% of the library** and recovers **~85–95%
of the true top-100**.

📖 [Documentation](https://aakankschit.github.io/TACTICS/) · 🧪 [Tutorials](#learn-more) · 📋 [Changelog](CHANGELOG.md)

---

## Is this for you?

TACTICS fits if you have:

- A **combinatorial library** defined by a reaction plus lists of reagents
  (2 or 3 components), and
- A **scoring function** that is too slow to run on every product — docking,
  shape similarity, or an ML model — or a big table of precomputed scores.

If your library is small enough to score exhaustively, you don't need TACTICS.

## Install

```bash
pip install chem-tactics
```

Requires Python 3.11+.

<details>
<summary>Optional extras and development install</summary>

```bash
pip install chem-tactics[tutorials]   # interactive marimo notebooks
pip install chem-tactics[test]        # test dependencies

# development install
git clone https://github.com/aakankschit/TACTICS.git
cd TACTICS
pip install -e ".[test]"
```
</details>

## Quickstart

Screen a two-component amide library against a table of precomputed scores:

```python
from TACTICS.library_enumeration import SynthesisPipeline, ReactionConfig, ReactionDef
from TACTICS.thompson_sampling import ThompsonSampler, get_preset
from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

# 1. Describe the library: one reaction + one reagent file per component
pipeline = SynthesisPipeline(ReactionConfig(
    reactions=[ReactionDef(
        reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
        step_index=0,
    )],
    reagent_file_list=["acids.smi", "amines.smi"],
))

# 2. Describe how to score a product
evaluator = LookupEvaluatorConfig(ref_filename="scores.csv")

# 3. Take a tuned preset and run
config = get_preset(synthesis_pipeline=pipeline, evaluator_config=evaluator)
sampler = ThompsonSampler.from_config(config)

sampler.warm_up(num_warmup_trials=config.num_warmup_trials)
results = sampler.search(num_cycles=config.num_ts_iterations)
sampler.close()

print(results.sort("score", descending=True).head(10))
```

`results` is a Polars DataFrame of every product that was evaluated, with
columns `score`, `SMILES` and `Name`.

> With `LookupEvaluator` and `DBEvaluator`, TACTICS skips molecule generation
> entirely and scores by product code, so `SMILES` reads `FAIL` and products
> are identified by `Name` (the underscore-joined reagent names). Evaluators
> that need a molecule — docking, ROCS, fingerprints — populate `SMILES`.

## Slow scoring functions: use more cores

Docking and ROCS take seconds per molecule, so evaluation dominates runtime.
Set `processes` to the number of cores you have — each worker builds its own
docking session, so OpenEye evaluators parallelize correctly:

```python
from TACTICS.thompson_sampling.core.evaluator_config import FredEvaluatorConfig

config = get_preset(
    synthesis_pipeline=pipeline,
    evaluator_config=FredEvaluatorConfig(design_unit_file="receptor.oedu"),
    mode="minimize",   # docking scores: lower is better
    batch_size=100,
)
config.processes = 32   # or os.cpu_count()

sampler = ThompsonSampler.from_config(config)
```

> **Leave `processes=1` for `LookupEvaluator` and `DBEvaluator`.** Their lookups
> are faster than the cost of shipping work to a worker, so parallelism makes
> those runs *slower*.

> **On macOS/Windows**, guard your entry point with
> `if __name__ == "__main__":` — the default `spawn` start method re-imports
> your script in every worker. On a cluster, make sure the `.oedu` receptor is
> readable from every node (automatic on NFS/GPFS).

## What can I plug in?

**Scoring functions** (`evaluator_config`):

| Evaluator | Use it for | Speed |
|---|---|---|
| `LookupEvaluatorConfig` | A CSV/Parquet table of precomputed scores | instant |
| `DBEvaluatorConfig` | Scores in a SQLite database | instant |
| `FPEvaluatorConfig` | Fingerprint similarity to a reference ligand | fast |
| `ROCSEvaluatorConfig` | 3D shape/colour overlay (OpenEye) | slow |
| `FredEvaluatorConfig` | Docking into a receptor (OpenEye) | slow |
| `MLClassifierEvaluatorConfig` | A trained scikit-learn model | varies |

**Presets** (`get_preset(name, ...)`):

| Preset | What it does |
|---|---|
| `recommended` *(default)* | Top-Two Thompson Sampling — best overall |
| `recommended_rws` | Roulette wheel with criticality-aware thermal cycling |
| `baseline` | Balanced warmup + greedy, for measuring what the search adds |
| `legacy_rws` | Reproduces Zhao et al. 2025 |

Start with `recommended`. If you are benchmarking, run `recommended` and
`recommended_rws` and take the better result — they favour different library
structures.

<details>
<summary>Building a config by hand instead of using a preset</summary>

Presets are just `ThompsonSamplingConfig` objects, so you can construct one
directly to control the strategy, warmup and search budget:

```python
from TACTICS.thompson_sampling import ThompsonSamplingConfig, ThompsonSampler
from TACTICS.thompson_sampling.strategies.config import TopTwoConfig
from TACTICS.thompson_sampling.warmup.config import EnhancedWarmupConfig

config = ThompsonSamplingConfig(
    synthesis_pipeline=pipeline,
    evaluator_config=evaluator,
    strategy_config=TopTwoConfig(mode="maximize", beta=0.5),
    warmup_config=EnhancedWarmupConfig(),
    num_warmup_trials=3,
    num_ts_iterations=1000,
    batch_size=100,
    processes=1,
    seed=42,
)

sampler = ThompsonSampler.from_config(config)
```

Configs are Pydantic models, so invalid values raise `ValidationError` at
construction rather than failing mid-run. Other selection strategies —
`GreedyConfig`, `RouletteWheelConfig`, `UCBConfig`, `EpsilonGreedyConfig`,
`BayesUCBConfig` — are available for baselines and comparisons.

See the [configuration docs](https://aakankschit.github.io/TACTICS/) for the
full reference.
</details>

## Learn more

Interactive [marimo](https://marimo.io) notebooks in `tutorials/` (they default
to a bundled thrombin dataset, so they run with no setup):

```bash
marimo run tutorials/thompson_sampling_tutorial.py    # compare strategies
marimo edit tutorials/library_enumeration_tutorial.py # build a library
```

| Tutorial | What it covers |
|---|---|
| `thompson_sampling_tutorial.py` | Comparing selection strategies |
| `library_enumeration_tutorial.py` | `SynthesisPipeline` and enumeration |
| `reaction_config_builder.py` | Building and validating reaction SMARTS |
| `custom_evaluator_tester.py` | Writing your own evaluator |

- **API reference and guides**: [TACTICS Documentation](https://aakankschit.github.io/TACTICS/)
- **Runnable scripts**: see `examples/`
- **Build docs locally**: `cd docs && make html`

## Contributing

Fork, branch, add tests for new behaviour, make sure `pytest tests/` passes,
then open a pull request.

## Citation

```bibtex
@software{tactics,
    title={TACTICS: Thompson Sampling-Assisted Chemical Targeting and
           Iterative Compound Selection for Drug Discovery},
    author={Aakankschit Nandkeolyar},
    year={2025},
    url={https://github.com/aakankschit/TACTICS}
}
```

## Support

Open an issue on GitHub, or contact anandkeo@uci.edu.

Licensed under the MIT License — see [LICENSE](LICENSE).

---

**This work is based on [previous work](https://pubs.acs.org/doi/10.1021/acs.jcim.3c01790) by Patrick Walters.**
**This project is a collaboration between the University of California Irvine, Leiden University and Groningen University.**
