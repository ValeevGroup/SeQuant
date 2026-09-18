Cost Model and Single-Term Optimization
==========================================

:doc:`sequant::optimize() </user/guide/optimize>` picks a contraction order for every product it sees by solving, for each term
independently, a subset `dynamic program <https://en.wikipedia.org/wiki/Dynamic_programming>`_ (DP) over its tensors: which pairwise
contraction to form first, second, and so on, so as to
minimize some notion of cost. This page documents the architecture behind that DP — the extension point for anyone implementing a new
cost objective — and two self-contained refinements shipped on top of it. It complements the :doc:`user-facing optimize() guide
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

A type satisfying this concept can be passed directly to ``run_single_term_opt`` to obtain a contraction order under an arbitrary,
user-defined cost function — this is the intended extension point for a custom objective, rather than modifying the driver itself.

Four built-in models satisfy the concept, corresponding to :class:`sequant::ObjectiveFunction`'s non-batched/batched pairs:

- ``AdditiveModel`` (parameterized by a cost functor) implements ``ObjectiveFunction::DenseFLOPs`` and ``DenseSize``: a plain additive
  DP over flop count or intermediate storage, with no notion of memory residency.
- ``PeakModel`` implements ``DenseSpaceTime``: an "all-co-resident" pebble-game DP that tracks, at each step, the total size of every
  tensor simultaneously resident (the currently-forming result plus whatever its sibling subtree's inputs still occupy), maintaining a
  `Pareto frontier <https://en.wikipedia.org/wiki/Pareto_front>`_ of ``(peak, flops)`` points per subset so the lexicographic optimum can be
  read off at the end.
- ``PeakBatchedModel`` implements ``DenseSpaceTimeBatched``/``DenseTimeSpaceBatched``: the same pebble-game idea, extended with a second
  dimension — for each subset, a DP cell per *sliced-set context* over the batchable indices (see :doc:`the user-facing batching guide
  </user/guide/batching>` for what "batchable," "contracted," and "peak" mean) — so slicing a mode is just one more move the DP can make
  to trade flops for a lower peak, on the same Pareto-frontier footing as a plain reordering.

.. note::
   For a product of exactly one or two tensors there is nothing to optimize (only one possible contraction sequence exists), so
   ``run_single_term_opt``/``run_single_term_opt_axes`` short-circuit before ever building a ``Context`` or invoking the DP. In
   particular, a bare two-tensor product never receives a batching annotation, even under a batched objective and a policy that would
   otherwise batch it — batching only becomes a real DP decision once there are at least three factors to choose an order over.

The roofline tie-break
--------------------------

The peak-first objectives (``DenseSpaceTime``, ``DenseSpaceTimeBatched``) rank schedules primarily by peak memory; among schedules tied
on peak, a secondary cost breaks the tie. By default that secondary cost is just flop count, but :class:`sequant::RooflineParams` (an
``OptimizeOptions::roofline`` field) switches it to a wall-time `roofline <https://en.wikipedia.org/wiki/Roofline_model>`_ proxy that correctly distinguishes compute-bound from
bandwidth-bound contractions, which raw flop count cannot do on its own:

.. math::
   \mathrm{cost} = \max(\mathrm{flops},\ \beta \cdot Q), \qquad
   Q = \max\!\left(\mathrm{traffic},\ \kappa \cdot \frac{\mathrm{flops}}{\sqrt{M / c_0}}\right)

where :math:`\beta` (``machine_balance``) is the machine's FLOPs-per-byte-of-traffic ratio, :math:`\mathrm{traffic}` is the contraction's
operand-plus-result footprint (elements), :math:`M` (``fast_mem_elems``) is the capacity of the relevant fast-memory level, :math:`c_0`
(``block_tiles``) is the number of resident tiles a blocked implementation needs (roughly 3, for the two operands and the result of a
blocked GEMM), and :math:`\kappa` (``block_prefactor``) folds in backend-specific constants (FMA width, packing overhead, ...). With
:math:`\beta \le 0` (the default) this reduces exactly to :math:`\mathrm{cost} = \mathrm{flops}`, the plain tie-break.

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

Pruning is controlled by ``OptimizeOptions::prune_outer_products`` (default ``true``) and can be force-disabled at runtime via the
``SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING`` environment variable, which is useful when validating a suspected pruning-related regression
without recompiling. The restriction's soundness is checked empirically (a parity test comparing pruned and unpruned search results),
not proven in general — treat it, as SeQuant's own tests do, as an assumption believed to hold for the class of expressions SeQuant
generates rather than as a theorem covering every conceivable input product.
