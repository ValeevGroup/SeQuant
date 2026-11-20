Cost Model and Single-Term Optimization
==========================================

:doc:`sequant::optimize() </user/guide/optimize>` picks a contraction order for every product it sees (a *term*, in the
code's naming, hence ``run_single_term_opt``) by solving, for each term independently, a subset `dynamic program
<https://en.wikipedia.org/wiki/Dynamic_programming>`_ (DP) over its tensors: which pairwise contraction to form first,
second, and so on, so as to
minimize some notion of cost. This page documents the architecture behind that DP — the extension point for anyone implementing a new
cost objective — and the roofline and connectivity refinements shipped on top of it. It complements the :doc:`user-facing optimize() guide
</user/guide/optimize>` and the :doc:`batched-evaluation architecture page <batched_evaluation>`, neither of which go into this level of
detail.

The ``CostModel`` concept and driver
----------------------------------------

The DP itself — the subset lattice, the bipartition enumeration, and the bookkeeping to turn a set of per-subset decisions back into a
contraction sequence — is implemented once, in :func:`sequant::opt::detail::solve_single_term` and
:func:`sequant::opt::detail::run_single_term_opt`. What varies between cost objectives is only the recurrence: how expensive a given
contraction is, and how to combine per-subset results. That variation point is factored out as a C++20 concept,
:concept:`sequant::opt::detail::CostModel`, requiring a type to provide:

- two associated types, ``State`` (a DP cell) and ``Context`` (precomputed tables and mutable scratch), and
- six methods: ``build_context``, ``leaf``, ``init``, ``relax`` (the recurrence step, called once per candidate bipartition of a
  subset), ``finalize`` (a post-subset hook), and ``reconstruct`` (turns the solved table into an ``EvalSequence``).

The driver visits subsets in an order that guarantees both children of a bipartition are solved before ``relax`` runs.

The driver is templated on the concept rather than dispatching through a virtual base class because ``State`` differs
per model: a cost record for ``AdditiveModel``, a Pareto frontier for ``PeakModel``, and one frontier per batch context
for ``PeakBatchedModel``. A common virtual interface would have to type-erase it.

A type satisfying this concept can be passed directly to ``run_single_term_opt`` to obtain a contraction order under an arbitrary,
user-defined cost function — this is the intended extension point for a custom objective, rather than modifying the driver itself.

Three model families satisfy the concept, implementing :enum:`sequant::ObjectiveFunction`'s objectives. Only
``DenseTimeSpace`` and ``DenseTimeSpaceBatched`` are intended for production; ``DenseFLOPs`` remains the default and the
baseline that matches the literature's operation-count metric. ``DenseSize``, ``DenseSpaceTime`` and
``DenseSpaceTimeBatched`` are deprecated: do not use them, and do not build on them.

- ``AdditiveModel`` (parameterized by a cost functor) implements ``ObjectiveFunction::DenseFLOPs``: a plain additive
  DP over flop count, with no notion of memory residency. With subnet CSE enabled, each canonically distinct
  subnetwork's cost is counted once, however many times it appears in the term.
- ``PeakModel`` implements ``DenseTimeSpace``: an "all-co-resident" pebble-game DP that tracks, at each step, the total size of every
  tensor simultaneously resident (the currently-forming result plus whatever its sibling subtree's inputs still occupy), maintaining a
  `Pareto frontier <https://en.wikipedia.org/wiki/Pareto_front>`_ of ``(peak, flops)`` points per subset so the lexicographic optimum can be
  read off at the end.
- ``PeakBatchedModel`` implements ``DenseTimeSpaceBatched``: the same pebble-game idea, extended with a second
  dimension — for each subset, a DP cell per enclosing *ordered batch nest* over the batchable indices (see :doc:`the user-facing batching guide
  </user/guide/batching>` for what "batchable," "contracted," and "peak" mean) — so slicing a mode is just one more move the DP can make
  to trade flops for a lower peak, on the same Pareto-frontier footing as a plain reordering.

.. note::
   For a product of exactly one or two tensors there is nothing to optimize (only one possible contraction sequence exists), so
   ``run_single_term_opt``/``run_single_term_opt_axes`` short-circuit before ever building a ``Context`` or invoking the DP. In
   particular, a bare two-tensor product never receives a batching annotation, even under a batched objective and a policy that would
   otherwise batch it — batching only becomes a real DP decision once there are at least three factors to choose an order over.

The peak recurrence and frontier
------------------------------------

For a subset :math:`n` of tensors, let :math:`S_n` be the footprint of its result, :math:`L_n` the sum of its
input-leaf footprints, and :math:`P_n` a candidate schedule's peak, all in elements. For a bipartition :math:`n=l\cup
r`, evaluating the left child first costs

.. math::
   P_n^{l\,\mathrm{first}} = \max\!\left(L_r + P_l,\ S_l + P_r,\ S_l + S_r + S_n\right).

The right-first alternative exchanges :math:`l` and :math:`r`. The driver considers both orders for every
child-frontier pair; a leaf starts with :math:`P_n=S_n` and zero operation cost. The :math:`L_r` term counts the
untouched sibling's resident inputs while the left subtree is computed. Dropping it would describe a different memory
model.

Peak and performance cost are combined only at the root, so a single best point per subset is insufficient: a child's
slightly larger peak can be hidden below its parent's unavoidable peak while its lower operation cost still improves
the whole schedule. ``PeakModel`` therefore retains all non-dominated ``(peak, performance cost)`` points.
``DenseTimeSpace`` selects by performance cost, then peak.

.. _cost-model-batched-search:

Ordered batch contexts and costing
--------------------------------------

A batched DP cell is keyed by a subset and a *context*: the ordered nest of sliced modes enclosing it, outermost first.
A sliced mode's footprint is its block size (``batch_target_size``, capped by the full extent). Nest order matters for
recomputation, because an intermediate can be hoisted above loops over modes it does not carry. Nest depth is bounded;
when enumerating ordered nests would be too expensive, the model falls back to unordered sets of sliced modes and loses
that distinction.

At each contraction the DP may open loops over two kinds of mode: a contracted mode summed at that node, or an external
mode of the final result carried by the node. Contracted opens honor ``persistent_only``; external opens require
``batch_spectator_indices``. Where to open an external loop is part of the same search, not a separate pass after
optimization.

Contracted and external modes bound memory at different scopes. A contracted mode is summed at one node, so slicing it
there bounds that node's intermediate while the partials accumulate. An external mode is never summed: it is carried
onto the result of the node and of every ancestor up to the root. Slicing it at a node bounds only the subtree below
that node; the node's assembled result and every ancestor are still materialized whole. An external loop is therefore
useful only when it encloses the subtree producing the large values, which is why the search chooses the opening node
instead of opening external modes where they first appear.

When a loop opens at a node, the peak recurrence above uses the children's sliced footprints and also counts the
node's accumulator (or scatter destination), which stays resident while the children are evaluated batch by batch.

Two charges prevent slicing from appearing free in the performance objective: a node that cannot be hoisted out of
an enclosing loop is recomputed once per batch of that loop; and a persistent value sliced by a loop that a volatile
node closes must be rebuilt every iteration, so it is weighted by ``volatile_weight`` like a volatile one.

``DenseTimeSpaceBatched`` also prefers fewer sliced modes among otherwise equal schedules; root selection under the
budget is described in :doc:`/user/guide/batching`. The chosen nests are recorded on the expression for the runtime
(see :ref:`batched-evaluation-loop-identity`).

.. _cost-model-roofline:

Roofline performance cost
--------------------------

:class:`sequant::RooflineParams` (the ``OptimizeOptions::roofline`` field) switches the peak objectives' performance cost
from flop count to a wall-time `roofline <https://en.wikipedia.org/wiki/Roofline_model>`_ proxy that distinguishes
compute-bound from bandwidth-bound contractions, which raw flop count cannot do on its own:

.. math::
   \mathrm{cost} = \max(\mathrm{flops},\ \beta \cdot Q), \qquad
   Q = \max\!\left(\mathrm{traffic},\ \kappa \cdot \frac{\mathrm{flops}}{\sqrt{M / c_0}}\right)

where :math:`\beta` (``machine_balance``) is the machine's FLOPs-per-element-of-traffic ratio, :math:`\mathrm{traffic}` is the contraction's
operand-plus-result footprint (elements), :math:`M` (``fast_mem_elems``) is the capacity of the relevant fast-memory level, :math:`c_0`
(``block_tiles``) is the number of resident tiles a blocked implementation needs (roughly 3, for the two operands and the result of a
blocked GEMM), and :math:`\kappa` (``block_prefactor``) folds in backend-specific constants (FMA width, packing overhead, ...). With
:math:`\beta \le 0` (the default) this reduces exactly to :math:`\mathrm{cost} = \mathrm{flops}`, the plain tie-break.

For double-precision data, convert a compute rate :math:`F` (operations/second) and bandwidth :math:`B` (bytes/second)
to :math:`\beta=8F/B`; convert the fast-memory capacity from bytes to elements too. Use the same operation-count
convention as the cost counter. For a distributed calculation, the same formula can be instantiated one memory level
up: node memory becomes the fast level (``fast_mem_elems`` is the per-node resident element count) and the
interconnect the slow level, so :math:`\beta=8F_\mathrm{agg}/B_\mathrm{net}` from the aggregate compute rate and network
bandwidth. Traffic is computed from unsliced operands, so the proxy does not credit slicing with a smaller working
set.

This formula naturally covers three regimes: a contraction with little data movement relative to its flop count is compute-bound and
the cost stays close to :math:`\mathrm{flops}`; one with a great deal of traffic per flop (e.g. contracting a single shared index over
otherwise-large operands) is bandwidth-bound and the :math:`\beta Q` term dominates; and one whose working set does not fit in the fast
memory level pays the additional re-read penalty captured by the second term inside :math:`Q`. To calibrate :math:`\beta`/:math:`\kappa`
for a given machine, time a couple of representative contractions spanning both regimes and fit the two constants to match observed
wall time.

Outer-product DP pruning
----------------------------

The subset DP as described above considers every bipartition of every subset of tensors, including ones that would form a disconnected
subnetwork — an *outer product* of two pieces that share no summed index. For expressions that never contain a genuine outer product as
a top-level term (true of every coupled-cluster residual summand, by the `linked-cluster theorem
<https://en.wikipedia.org/wiki/Linked-cluster_theorem>`_), this wastes a large fraction of the
search: :func:`sequant::opt::detail::outer_product_connectivity` restricts the DP's subset lattice to only the subsets whose induced
subgraph — under "tensors A and B are adjacent iff they share a contracted (non-target) index" — is connected, collapsing the search
from exponential-in-tensor-count down to roughly the number of connected sub-networks.

Two subtleties keep this restriction correct rather than merely fast:

- **Hyperedges are kept, not pruned.** An index shared by three or more tensors (or one that is only contracted at a later step) still
  creates an edge between every pair of its carriers, so a legitimate Hadamard-type intermediate is never wrongly excluded by this
  adjacency test — it only rules out subsets that share *no* contracted index at all.
- **A genuine product term falls back to the unpruned search automatically.** If the connectivity check finds the *entire* network
  disconnected (i.e. the term itself is an honest outer product, not a residual summand), pruning is disabled for that term rather than
  incorrectly excluding the one bipartition that would actually need to be considered.

Pruning is controlled by ``OptimizeOptions::prune_outer_products`` (default ``true``) and can be force-disabled at
runtime via the ``SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING`` environment variable (any non-empty value, including ``0``,
disables it), which is useful when validating a suspected pruning-related regression without recompiling. Pruning
assumes an optimal contraction order can be formed entirely from connected subnetworks. This is checked empirically
(parity tests comparing pruned and unpruned search results), not proven in general. For inputs where that assumption may
not hold, disable pruning to retain the exhaustive search.

Protoindices do not create edges: two composites sharing only their protoindices are not adjacent.

The connectivity mask also restricts the per-subset tables built before the search, not just the search itself; in the
batched model those tables are replicated per sliced context and would otherwise dominate the run time.
