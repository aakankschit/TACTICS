Search
======

Block 3 — configure and run. See the :doc:`guide </guides/03_search>`.

Configuration
-------------

.. autoclass:: TACTICS.thompson_sampling.config.ThompsonSamplingConfig
   :members: reagent_file_list, num_components, num_steps

Presets
-------

.. autofunction:: TACTICS.thompson_sampling.presets.get_preset

.. autoclass:: TACTICS.thompson_sampling.presets.ConfigPresets
   :members:

Sampler
-------

.. autoclass:: TACTICS.thompson_sampling.core.sampler.ThompsonSampler
   :members: from_config, warm_up, search, close, read_reagents, set_evaluator, set_hide_progress, evaluate, evaluate_batch, load_product_library, get_num_prods

Selection strategies
--------------------

Configs
~~~~~~~

.. automodule:: TACTICS.thompson_sampling.strategies.config
   :members: TopTwoConfig, RouletteWheelConfig, GreedyConfig, UCBConfig, EpsilonGreedyConfig, BayesUCBConfig

Classes
~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.strategies.top_two_selection.TopTwoSelection
   :members: select_reagent, adapt_heated_scale, adapt_temperatures, reset_temperature, get_component_state

.. autoclass:: TACTICS.thompson_sampling.strategies.roulette_wheel.RouletteWheelSelection
   :members: select_reagent, adapt_temperatures, reset_temperature, get_component_state

.. autoclass:: TACTICS.thompson_sampling.strategies.greedy_selection.GreedySelection
   :members: select_reagent

.. autoclass:: TACTICS.thompson_sampling.strategies.ucb_selection.UCBSelection
   :members: select_reagent

.. autoclass:: TACTICS.thompson_sampling.strategies.epsilon_greedy.EpsilonGreedySelection
   :members: select_reagent

.. autoclass:: TACTICS.thompson_sampling.strategies.bayes_ucb_selection.BayesUCBSelection
   :members: select_reagent, get_component_criticality, get_component_state, reset_percentiles

Warmup strategies
-----------------

.. automodule:: TACTICS.thompson_sampling.warmup.config
   :members: EnhancedWarmupConfig, BalancedWarmupConfig

.. autoclass:: TACTICS.thompson_sampling.warmup.enhanced.EnhancedWarmup
   :members: generate_warmup_combinations, get_expected_evaluations, get_name

.. autoclass:: TACTICS.thompson_sampling.warmup.balanced.BalancedWarmup
   :members: generate_warmup_combinations, get_expected_evaluations, get_name

Reagent
-------

.. autoclass:: TACTICS.thompson_sampling.core.reagent.Reagent
   :members: add_score, sample, init_prior, mean, std, n_samples

Factories
---------

.. automodule:: TACTICS.thompson_sampling.factories
   :members: create_strategy, create_warmup, create_evaluator
