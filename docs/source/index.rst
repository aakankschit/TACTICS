.. TACTICS documentation master file, created by
   sphinx-quickstart on Thu Apr 24 15:14:03 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to TACTICS documentation!
=================================

TACTICS (Thompson Sampling-Assisted Chemical Targeting and Iterative Compound Selection for Drug Discovery) is a comprehensive Python package for Thompson Sampling-based optimization of chemical combinatorial libraries with a unified, flexible architecture.

Key Features
------------

* **Unified Thompson Sampler**: Single sampler class that accepts different selection strategies
* **Multiple Selection Strategies**: Greedy, Roulette Wheel (thermal cycling), UCB, Epsilon-Greedy
* **Flexible Warmup Strategies**: Standard, Stratified, Enhanced, Latin Hypercube sampling
* **Multiple Evaluators**: Lookup, Database, Fingerprint, ML models, ROCS, FRED docking
* **Parallel Evaluation**: Built-in multiprocessing support for expensive evaluators
* **Extensible Design**: Easy integration of custom strategies, warmup methods, and evaluators
* **Library Enumeration**: Efficient generation of combinatorial reaction products
* **Library Analysis**: Comprehensive analysis and visualization tools

Quick Start
-----------

.. code-block:: python

    from TACTICS.thompson_sampling.core.sampler import ThompsonSampler
    from TACTICS.thompson_sampling.strategies import GreedySelection
    from TACTICS.thompson_sampling.core.evaluators import LookupEvaluator

    # Create selection strategy
    strategy = GreedySelection(mode="maximize")

    # Create sampler
    sampler = ThompsonSampler(selection_strategy=strategy)

    # Setup
    sampler.read_reagents(["reagents1.smi", "reagents2.smi"])
    sampler.set_reaction("[C:1]=[O:2]>>[C:1][O:2]")
    sampler.set_evaluator(LookupEvaluator({"ref_filename": "scores.csv"}))

    # Run optimization
    warmup_results = sampler.warm_up(num_warmup_trials=3)
    search_results = sampler.search(num_cycles=100)

    # Cleanup
    sampler.close()

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   thompson_sampling
   api_reference
   configuration
   library_enumeration
   library_analysis

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

