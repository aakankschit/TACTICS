Scoring
=======

Block 2 — turn a product into a number. See the :doc:`guide </guides/02_scoring>`.

Each evaluator has a paired Pydantic config; pass the config to
:class:`~TACTICS.thompson_sampling.config.ThompsonSamplingConfig` and the
sampler builds the evaluator (and rebuilds it inside each worker when
``processes > 1``).

Evaluator configs
-----------------

.. automodule:: TACTICS.thompson_sampling.core.evaluator_config
   :members:

Evaluator classes
-----------------

.. automodule:: TACTICS.thompson_sampling.core.evaluators
   :members: Evaluator, LookupEvaluator, CustomEvaluator, DBEvaluator, FPEvaluator, MWEvaluator, ROCSEvaluator, FredEvaluator, MLClassifierEvaluator
