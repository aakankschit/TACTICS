Scale
=====

Block 4 — parallel evaluation and pre-enumerated libraries. See the :doc:`guide </guides/04_scale>`.

The user-facing knobs are fields on
:class:`~TACTICS.thompson_sampling.config.ThompsonSamplingConfig`
(``processes``, ``min_cpds_per_core``, ``product_library_file``) and the
sampler methods :meth:`~TACTICS.thompson_sampling.core.sampler.ThompsonSampler.set_evaluator`
and :meth:`~TACTICS.thompson_sampling.core.sampler.ThompsonSampler.load_product_library`.
The machinery behind ``processes > 1``:

.. autoclass:: TACTICS.thompson_sampling.core.parallel_evaluator.ParallelEvaluator
   :members: bind_worker_context, evaluate_batch, close
