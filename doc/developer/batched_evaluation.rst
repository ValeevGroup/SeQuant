Batched Evaluation Architecture
===================================

This page documents the runtime machinery behind :doc:`the user-facing batching guide </user/guide/batching>`: how a
:class:`sequant::BatchPolicy` decision is actually carried out when numerically evaluating a *forest*: the evaluation
trees of several equations (e.g. all residual equations of a coupled-cluster method) passed to :func:`sequant::evaluate`
together. It assumes the vocabulary established there (mode, batch, slice, contracted/external role,
persistent/volatile, peak memory) and the cost-model architecture in :doc:`cost_model`; neither is redefined here.

Two execution strategies
---------------------------

Numeric evaluation is reached through :func:`sequant::evaluate`, dispatched by ``BatchPolicy::scheduler``
(:enum:`sequant::BatchScheduler`) onto one of two strategies:

- **``forest_descent``** (the default): each equation in a forest is evaluated independently, one tree at a time, by the original
  tree-walking evaluator. Batching is realized by a *custom evaluator* (:func:`sequant::make_batched_custom_evaluator`, or the
  ``BatchPolicy``-driven adapter :func:`sequant::make_evaluator`) installed on the node cache: when the walk reaches a node whose
  contracted or external mode is batchable, the custom evaluator re-enters the subtree once per slice and accumulates/scatters the
  partial results, instead of letting a single one-shot contraction run.
- **``ordered``**: the entire forest is fused into one `DAG <https://en.wikipedia.org/wiki/Directed_acyclic_graph>`_ and driven by a single, explicitly-scheduled, table-driven walk
  (``SeQuant/core/eval/ordered_executor.hpp``). Loop identity and the cell table make compatible value forms and batch loops shareable
  across equations, with one production per shared form at its scheduled scope.

Forest descent can already reuse ordinary cached intermediates across equations. The ordered scheduler additionally fuses compatible
batch loops across the forest; the per-tree custom evaluator does not perform that fusion. Both strategies produce the same
numeric result for the same ``BatchPolicy``; the ordered scheduling pipeline is described below.

``forest_descent`` builds no schedule or cell table. Its custom evaluator hoists values to the enclosing scope where
they are reused and slices each fetched value to the current batch.

Array layout and mode frames
--------------------------------

The runtime's array-mode list is ``EvalExpr::canon_indices()``. ``indices_annot()`` presents proto-free outer modes
first, then proto-carrying inner composites, separated by a semicolon. A leaf evaluator must return an array in the
requesting node's canonical layout.

Labels pair operand modes within one contraction. Decisions across occurrences use positions in the canonical layout
instead, because two occurrences can have different labels or permutations. A composite index such as ``a<i,j>``
contributes the array mode ``a`` over a proto-dependent range; slicing a plain ``i`` mode does not also slice that
composite mode.

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
   ``OrderedSchedule``: a tree of nested ``ScopeBlock`` objects (one per physical batch loop), each holding the build/assemble steps that
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

.. _batched-evaluation-loop-identity:

Loop identity and value identity
------------------------------------

A *value* is an identity in the fused DAG, an *occurrence* is a use of it by a consumer at a particular operand
position, and a *cell* is one form of a value resident at one scope. Structural node identity alone cannot identify a
value under batching: identical expressions sliced by different physical loops must remain distinct.

``LoopKey`` identifies a physical loop. An index space only constrains which loops can match: two loops over the same
space can be distinct and have different combine roles.

``compute_dag_boulevard`` identifies loops before values. It merges occurrence positions that must share a loop, both
along parent/child edges and across trees, and refuses a merge that would join two positions of one occurrence, mix
contracted and external loops, put a reader of a complete value inside the reduction it awaits, or violate nesting
constraints. A value's key then combines its structural identity with the loops slicing or reducing it and its
operands' keys; without slicing, it is just the structural identity.

The optimizer records on each node the modes it slices and the loops it opens; these are shared by every occurrence of
the node. Two occurrences can still be sliced differently in different equations, so the ordered scheduler also
records each occurrence's *home*: the enclosing loops that slice its own result modes (``occurrence_home()``). Value
identity tells apart occurrences with different homes. Forest descent instead uses the intersection over
occurrences (``sliced_modes()``).

Placement, escapes, and pass splits
----------------------------------------

The builder realizes a chain per loop instance, honoring the recorded outer/inner constraints. Loops that never
co-occur on a home-sliced value form independent nests. Dependency edges from occurrences topologically order builds
and nested blocks; external requirements of an inner block propagate outward to the level that can satisfy them.

The ordered driver's optional ``mode_order`` ranks index spaces from outermost to innermost. It only breaks ties left
by the recorded nesting constraints; there is no cost search over nest orders.

Legality classifies a value relative to each loop as ``LoopLocal``, ``Reduction``, or ``LoopCarried``. A wholly local
value is a ``BuildStep`` in its home block. A reduction or loop-carried value instead escapes through the block's
outputs: ``AccumulateSum`` for reduction, ``AccumulateScatter`` for carried result modes. Multiple escapes form a
chain, with raw production at the deepest site and forwarding at shallower sites. An escape can skip an interposed
loop on which the value is invariant.

A consumer requiring a completed loop-carried value, or a completed reduction while it is inside that reduction's
loop, moves to a later pass. A nest containing multiple passes becomes ordered sibling blocks, run in pass order; a
later pass reads the earlier pass's assembled value. Values needed across a split are materialized, not recomputed.
Unsupported cases raise ``sequant::Exception``.

.. _batched-evaluation-cell-table:

Cell forms, reads, and lifetimes
------------------------------------

A ``TableCell`` is one form of a value: the scope it lives at, the positions it is sliced on, the loops over which it
is still a partial sum, and how it is produced (leaf, build, or sum/scatter assembly). A consumer reads the deepest form
whose residency encloses it.

Each operand occurrence is a separate ``Read``, so a consumer using one value twice reads it twice. A read names its
source cell and how to slice it. A cell's lifetime is the number of reads it will serve, each weighted by how many times
the loops enclosing the consumer but not the source replay it; the executor frees the cell when the count reaches zero.
A cell invariant on an enclosing loop is built on the first batch and reused until a loop it depends on advances. Only
unbound, non-volatile values at the persistent-to-volatile frontier persist across evaluations.

``validate_cell_table`` checks the table for consistency before any tensor operation runs: every read has a visible
source of the right form, assemblies close the right partial or scattered form, lifetimes match the weighted reads,
and no value is built twice at one scope. It does not check array ranges; the dry-run backend does.

Executing the table
-----------------------

The executor resolves each declared read against the cell registry, spends one unit of the source's lifetime, and
applies the read's slices; the last non-persistent read moves the value out, which makes in-place accumulation safe.
Each batch clears the cells bound to that loop, executes the block's steps, and folds outputs into their assembly
cells (summing partials, or scattering blocks into the destination). When the loop closes, the assembled cell becomes
visible to the next link of its escape chain. Persistent cells are published to the cache's ``PersistentValueStore``
and restored from it on later evaluations, in which case their productions are skipped.

A correctness caveat for sparse backends
---------------------------------------------

Some numeric backends (TiledArray in particular) screen individual result tiles against a norm bound (typically a `Cauchy-Schwarz
<https://en.wikipedia.org/wiki/Cauchy%E2%80%93Schwarz_inequality>`_ bound) computed for the *full* contraction. Evaluating only a fraction of a
batched reduction and comparing each partial result against that same full-contraction bound can wrongly discard tiles that would have
survived the complete sum — silently dropping a real contribution. The executor therefore accepts a *scope guard* factory (instantiated once
per batch loop, released via `RAII <https://en.wikipedia.org/wiki/Resource_acquisition_is_initialization>`_ as nested loops close) whose job is to
scale such a screening threshold down for the duration of each batched block, so a partial sum is screened against a bound appropriate
to its own, smaller size rather than the full one. Any backend/executor integration that batches over a screened or otherwise
sparsity-aware backend must supply a correct scope-guard factory; omitting it does not fail loudly — it silently returns a
slightly-wrong numeric result.

.. _batched-evaluation-dry-runs:

Dry runs and diagnostics
-----------------------------

The sizing backend in ``SeQuant/core/eval/backends/dryrun/`` stores descriptors instead of tensors and checks
shared-label ranges on every operation. :func:`sequant::eval::dryrun::meter` drives the same ``evaluate`` entry point
as a real run and reports the realized peak and work, checking the runtime's scheduling and accounting rather than
the optimizer's per-product estimate. It assumes real double-precision data.

``SeQuant/core/eval/ordered_dump.hpp`` provides environment-gated dumps of each pipeline stage;
:doc:`/user/guide/evaluation` explains how to read evaluation traces.

Current scope
----------------

The ``ordered`` pipeline does not: choose the fusion of loop instances optimally where producer/consumer
connectivity allows more than one valid choice (a valid, but not necessarily optimal, choice is recorded); reorder loop nesting across
index spaces for cost (the realized nesting honors recorded constraints and the caller's ``mode_order`` ranking); make a
cost-aware choice between materializing and recomputing a value across a forced pass split (a pass split always materializes); or
share a transient value across multiple reads (an unrecorded intermediate is recomputed inside every production tree that needs it).
A scatter assembly writes along one position of its output; one needing several positions throws.
None of this affects correctness — every one of these is a possible future efficiency
improvement, not a currently-missing correctness guarantee — but they bound how much benefit ``ordered`` scheduling can currently
realize relative to a hypothetical, fully cost-optimal fusion.
