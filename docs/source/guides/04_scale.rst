.. _guide-scale:

4. Scale
========

**What this block is.** The search loop spends almost all its time in the
scorer. When the scorer is slow — docking, shape overlay, a model that
takes a second per molecule — this block spreads that cost over worker
processes, or skips the synthesis step by handing the sampler
pre-enumerated products. Nothing else about the run changes.

Core
----

Two fields on the config control parallel scoring:

.. literalinclude:: /snippets/scale_parallel.py
   :language: python
   :start-after: # [start:parallel]
   :end-before: # [end:parallel]
   :dedent:

- ``processes`` — worker processes for evaluation. Selection, posterior
  updates and synthesis bookkeeping stay in the parent.
- ``min_cpds_per_core`` — the sampler accumulates ``processes ×
  min_cpds_per_core`` products before dispatching a batch to the pool (and
  flushes whatever is left on the last cycle). Bigger batches, fewer
  round-trips; the default 10 is fine.

**How workers get their scorer.** Evaluators such as FRED and ROCS hold
C++ objects that cannot be pickled, so the evaluator itself never crosses
the process boundary. Each worker is started once with the *config* and
builds its own evaluator from it. That is why
``ThompsonSampler.from_config`` is the parallel-safe path: it hands the
config through. It is also why a ``CustomEvaluatorConfig`` function must be
importable by name (module level).

**When it pays.** Only when a single score costs more than the ~ms it takes
to pickle a reagent tuple and collect the result. For ``LookupEvaluator``
and ``DBEvaluator`` the sampler skips synthesis entirely and the lookup is
microseconds; ``processes > 1`` there is slower, and the sampler logs a
warning saying so.

**What it produces:** the same results DataFrame as Block 3. ``close()``
now matters — it shuts the pool down.

Build on it
-----------

Pre-enumerated products
~~~~~~~~~~~~~~~~~~~~~~~

If the library has already been enumerated (Block 1's
``enumerate_library`` → ``write_enumerated_library(..., format="csv")``,
or another tool), give the sampler the file and it will look products up
instead of running the reaction:

.. code-block:: python

   config.product_library_file = "library.csv"   # columns: Product_Code, SMILES

The lookup is by product name; a miss falls back to synthesis, so a partial
file is fine. This removes the RDKit reaction cost from every evaluation —
worth it for structure-based scorers on large libraries, irrelevant for
lookup scorers (which never synthesise anyway).

Budgeting a run
~~~~~~~~~~~~~~~

- ``sampler.search(num_cycles, max_evaluations=N)`` stops after ``N`` scored
  products even if cycles remain — the natural knob when the budget is a
  number of docking runs.
- Warmup is separate and not capped: Enhanced warmup costs
  ``max(component sizes) × num_warmup_trials`` evaluations before the
  search starts. On 3,844 amines × 5 trials that is ~19 k docking runs;
  drop ``num_warmup_trials`` to 2–3 for expensive scorers, or use
  Balanced warmup (``sum(sizes) × K``) when one component is very large.
- ``batch_size`` × ``num_ts_iterations`` is the search budget. With
  ``processes = 32`` and ``min_cpds_per_core = 10``, a ``batch_size`` of
  320 or more keeps every worker busy each cycle.

On a cluster
~~~~~~~~~~~~

- Set ``config.processes`` to the cores you were allocated (e.g. from
  ``SLURM_CPUS_PER_TASK``); the pool is a plain ``multiprocessing`` pool.
- ``config.hide_progress = True`` keeps tqdm out of the job log;
  ``get_preset(..., output_dir=...)`` or ``config.log_filename`` writes the
  run log to a file.
- Write results yourself at the end (``results.write_parquet(...)``); the
  sampler holds them in memory only.
- OpenEye licences are checked per process — make sure workers can see
  ``OE_LICENSE``.
- Put the run in ``if __name__ == "__main__":``. With the ``spawn`` start
  method (macOS, Windows) each worker re-imports the script; an unguarded
  script starts workers recursively.

.. admonition:: Gotchas
   :class: warning

   - **Direct construction and** ``processes > 1``: call
     ``sampler.set_evaluator(evaluator, evaluator_config=cfg)`` — without the
     config, workers have nothing to rebuild from, and an unpicklable
     evaluator raises a ``TypeError`` that says so.
   - **A lambda scoring function fails with** ``processes > 1`` even though
     it works with ``processes = 1``.
   - **Parallelism does not change the search.** Batches are dispatched in
     evaluation order and posteriors update once per batch, so ``seed``
     still reproduces the run — but posteriors update *less often* than
     with ``batch_size = 1``, which is the trade you make for throughput.
   - ``processes`` counts *evaluation* workers; RDKit reaction enumeration
     for ``enumerate_library`` has its own ``n_jobs``.

Reference
---------

:doc:`/reference/scale` — :class:`~TACTICS.thompson_sampling.core.parallel_evaluator.ParallelEvaluator`;
the ``processes``, ``min_cpds_per_core`` and ``product_library_file`` fields
of :class:`~TACTICS.thompson_sampling.config.ThompsonSamplingConfig`;
:meth:`~TACTICS.thompson_sampling.core.sampler.ThompsonSampler.set_evaluator`.

**Next:** :doc:`05_inspect` — what the search learned.
