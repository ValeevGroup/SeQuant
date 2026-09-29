Batched Evaluation Architecture
===================================

This page documents the runtime machinery behind :doc:`the user-facing batching guide </user/guide/batching>`: how a
:class:`sequant::BatchPolicy` decision is actually carried out when numerically evaluating a forest of equations. It assumes the
vocabulary established there (mode, batch, slice, contracted/external role, persistent/volatile, peak memory) and the cost-model
architecture in :doc:`cost_model`; neither is redefined here.

Two execution strategies
---------------------------

Numeric evaluation is reached through :func:`sequant::evaluate`, dispatched by ``BatchPolicy::scheduler``
(:class:`sequant::BatchScheduler`) onto one of two strategies:

- **``forest_descent``** (the default): each equation in a forest is evaluated independently, one tree at a time, by the original
  tree-walking evaluator. Batching is realized by a *custom evaluator* (:func:`sequant::make_batched_custom_evaluator`, or the
  ``BatchPolicy``-driven adapter :func:`sequant::make_evaluator`) installed on the node cache: when the walk reaches a node whose
  contracted or external mode is batchable, the custom evaluator re-enters the subtree once per slice and accumulates/scatters the
  partial results, instead of letting a single one-shot contraction run.
- **``ordered``**: the entire forest is fused into one `DAG <https://en.wikipedia.org/wiki/Directed_acyclic_graph>`_ and driven by a single, explicitly-scheduled, table-driven walk
  (``SeQuant/core/eval/ordered_executor.hpp``). Its advantage over per-tree ``forest_descent`` is cross-equation sharing: a value that
  is structurally identical across several equations in the forest — a common intermediate contraction, or a batch loop over the same
  index — is identified once and built once, rather than rebuilt independently inside every tree that needs it.

Both strategies produce the same numeric result for the same ``BatchPolicy``; ``ordered`` exists to recover the sharing that
per-tree evaluation cannot see, at the cost of a more involved scheduling pipeline, described below.

The ``ordered`` pipeline
----------------------------

Evaluating a forest under ``BatchScheduler::ordered`` proceeds through four stages, each with its own SeQuant header:

1. **Loop identity** (``SeQuant/core/eval/dag_scope.hpp``, ``peak_profile.hpp``, ``lifetime_mask.hpp``). Before anything is scheduled,
   every physical batch loop implied by the cost model's decision needs a stable identity distinct from where it happens to sit in any
   one tree — the same auxiliary-index loop may appear inside several different equations. ``LoopKey`` provides that identity;
   ``compute_dag_boulevard`` walks the whole forest once to assign it, and from it derives *value identity*: two occurrences (one
   use-site of a node, in one equation) are recognized as the same shareable value exactly when they agree on which loop instances
   they are sliced by.
2. **Schedule construction** (``ordered_schedule.hpp``). ``build_ordered_schedule`` lowers the loop-identity facts into an
   ``OrderedSchedule``: a tree of nested ``ScopeBlock``s (one per physical batch loop), each holding the build/assemble steps that
   belong at that nesting level. Two decisions happen here that a hand-written batched loop would otherwise have to get right by hand:
   *escape placement* — deciding, for each value, whether it is purely local to one loop iteration (built fresh every time, needing no
   special handling) or must be accumulated/scattered across iterations into something that outlives the loop — and *pass splitting* —
   when a producer/consumer ordering conflict would otherwise make a value depend on a not-yet-finished batch of itself, the affected
   values are separated into ordered sub-passes within the same loop nest rather than left to race.
3. **The cell table** (``cell_table.hpp``, ``cell_table_builder.hpp``). Rather than let the executor *infer* what to read from where
   while it runs, the schedule is first lowered into an explicit table: one entry per distinct value-form actually resident at some
   scope, and one entry per read of it, naming exactly which scope produces it, how long it lives, and who consumes it. A validator
   checks this table for internal consistency before any tensor operation runs, so a scheduling bug is caught as a fast, deterministic
   failure rather than a wrong number arrived at after an expensive run.
4. **The executor** (``ordered_executor.hpp``). Walks the cell table exactly as written: every operand of every operation is resolved
   through a table lookup, not by re-deriving it from the tree. This is a deliberate design choice — the executor performs no
   scheduling decisions of its own, only realizes ones the earlier stages already made, which keeps its own logic small relative to the
   scheduling machinery it depends on.

A correctness caveat for sparse backends
---------------------------------------------

Some numeric backends (TiledArray in particular) screen individual result tiles against a norm bound (typically a `Cauchy–Schwarz
<https://en.wikipedia.org/wiki/Cauchy%E2%80%93Schwarz_inequality>`_ bound) computed for the *full* contraction. Evaluating only a fraction of a
batched reduction and comparing each partial result against that same full-contraction bound can wrongly discard tiles that would have
survived the complete sum — silently dropping a real contribution. The executor therefore accepts a *scope guard* factory (instantiated once
per batch loop, released via `RAII <https://en.wikipedia.org/wiki/Resource_acquisition_is_initialization>`_ as nested loops close) whose job is to
scale such a screening threshold down for the duration of each batched block, so a partial sum is screened against a bound appropriate
to its own, smaller size rather than the full one. Any backend/executor integration that batches over a screened or otherwise
sparsity-aware backend must supply a correct scope-guard factory; omitting it does not fail loudly — it silently returns a
slightly-wrong numeric result.

Current scope
----------------

As of this writing, the ``ordered`` pipeline does not: choose the fusion of loop instances optimally where producer/consumer
connectivity allows more than one valid choice (a valid, but not necessarily optimal, choice is recorded); reorder loop nesting across
index spaces for cost (the realized nesting order follows how the loops were discovered, not a separate cost search); make a
cost-aware choice between materializing and recomputing a value across a forced pass split (a pass split always materializes); or
share a transient value across multiple reads (an unrecorded intermediate is recomputed inside every production tree that needs it).
It also supports scattering into at most one position of one value per loop instance; a schedule that would require more throws rather
than silently doing the wrong thing. None of this affects correctness — every one of these is a possible future efficiency
improvement, not a currently-missing correctness guarantee — but they bound how much benefit ``ordered`` scheduling can currently
realize relative to a hypothetical, fully cost-optimal fusion.
