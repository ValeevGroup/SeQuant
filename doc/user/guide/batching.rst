Batched (Memory-Bounded) Evaluation
======================================

A tensor contraction such as a `density-fitted <https://en.wikipedia.org/wiki/Resolution_of_the_identity>`_ 2-electron integral,
:math:`g_{ij}^{ab} \approx \sum_K B^a_{iK} B^b_{jK}`,
can involve one mode (here :math:`K`, an auxiliary fitting-basis index — SeQuant's convention labels it :code:`Κ`, see
:func:`sequant::mbpt::add_df_spaces`) that is far larger than the others. Materializing the full,
unsliced intermediate for such a contraction can dominate peak memory even though the final result is comparatively small. **Batching**
addresses this: SeQuant can slice a large mode into blocks and accumulate (or scatter) the contraction block by block, trading a
controlled amount of extra bookkeeping/recompute for a lower peak-memory footprint. Which modes get batched, and how
aggressively, is decided automatically by a cost-model search — the same :doc:`optimize() <optimize>` machinery that chooses contraction
order — rather than hand-coded.

Vocabulary
------------

A few terms are used precisely throughout this page and the API:

- **mode**: one dimension of a tensor (not to be confused with the :class:`sequant::Index` that labels it).
- **batch** vs. **slice**: a *batch* is the loop that partitions a mode's *work*; a *slice* is one block of that mode's *data*.
- **role**: whether a mode is *contracted* (summed away at some node) or *external* (present in the final result). A contracted mode can
  be batched *locally*, right at the node that sums it. An external mode remains in the final result; batching it
  assembles disjoint result blocks. The optimizer chooses the contraction node that opens the external loop.
  External batching requires a separate opt-in (see ``batch_spectator_indices`` below).
- **external** includes both a free mode (present in one operand and the result) and a Hadamard mode (present in both
  operands and the result), following the index roles of `TAPP <https://arxiv.org/abs/2601.07827>`_ (Brandejs et al.).
  Here, *batch* means memory blocking; it does not mean a Hadamard mode, as in batched GEMM terminology.
- **combine**: partials over a contracted mode are summed; partials over an external mode occupy disjoint result blocks
  and are scattered into a destination.
- **persistent** vs. **volatile**: a persistent intermediate depends only on data that stays fixed across a calculation (e.g. two-electron
  integrals) and can be cached and reused; a volatile one depends on data that changes every iteration (e.g. coupled-cluster amplitudes)
  and must be recomputed each time.
- **peak memory**: the maximum, over the whole evaluation, of the combined size of every tensor simultaneously resident — the quantity
  the batched cost models minimize or bound.

``BatchPolicy``
------------------

:class:`sequant::BatchPolicy` is the single configuration object consulted both when :func:`sequant::optimize` decides *which* modes to
batch and when the runtime evaluator actually executes that decision. Its full field list is API-reference material, but seven fields
matter for everyday use:

- ``peak_threshold``: the master on/off switch — a peak-memory budget in bytes. The default, :math:`+\infty`, means *no batching at
  all*; setting a finite value turns batching on and sets the budget the optimizer searches under (see below for what
  happens when no schedule fits).
- ``is_batchable_contracted_index`` / ``is_batchable_external_index``: predicates selecting which index *spaces* may be sliced, split by
  role (see above). Both default to "nothing is batchable"; a caller opts specific spaces in explicitly.
- ``batch_target_size``: a per-index **upper bound** on block size, in elements; the runtime rounds down to whole backend tiles.
  A single tile larger than the bound forms its own batch, so the bound cannot override the backend's tile geometry.
- ``batch_spectator_indices``: opts *external* modes into batching as well (in addition to contracted ones); off by default.
- ``is_volatile_leaf`` / ``persistent_only``: identify amplitude-like leaves and, optionally, restrict batching over
  contracted modes to subtrees that contain none of them.

.. literalinclude:: /examples/user/batching.cpp
   :language: cpp
   :start-after: start-snippet-1
   :end-before: end-snippet-1
   :dedent: 2

Deciding what to batch
-------------------------

Passing a ``BatchPolicy`` to :func:`sequant::optimize` via ``OptimizeOptions::batch_policy``, together with the *batched*
objective function ``ObjectiveFunction::DenseTimeSpaceBatched`` (performance first, peak memory second), makes the
cost-based search consider slicing as one more way to reduce peak memory, in addition to choosing contraction order:

.. literalinclude:: /examples/user/batching.cpp
   :language: cpp
   :start-after: start-snippet-2
   :end-before: end-snippet-2
   :dedent: 2

The finite ``peak_threshold`` is a ceiling on the **modeled** peak. Among the schedules that fit it, the optimizer picks
the lowest performance cost, then the fewest sliced modes, then the lowest peak, so nothing is sliced unless slicing is
needed to fit the budget. If no schedule fits, it picks the lowest performance cost, then the lowest peak.

Exceeding the budget is not an error. The modeled peak is a per-term estimate from dense extents; backend workspace,
tile sizes, and intermediates shared across equations can change the realized peak, so a finite budget does not
guarantee that the process fits in that many bytes.

Executing a batched plan
---------------------------

Once ``optimize()`` has annotated an expression with a batching decision, actually *running* it numerically needs a real tensor backend
(TiledArray, BTAS, or TAPP) and is reached through the same :func:`sequant::evaluate` entry points used for any other evaluation, plus
:func:`sequant::make_evaluator` to adapt a ``BatchPolicy`` into the evaluator SeQuant's cache consults. Two execution strategies are
available, selected by ``BatchPolicy::scheduler``: the default ``BatchScheduler::forest_descent`` evaluates one equation at a time, while
``BatchScheduler::ordered`` fuses an entire forest of equations into a single schedule so compatible value forms can be reused across
equations within shared batch loops. Because these entry points are backend-dependent, they are outside the scope of a
backend-agnostic example here; the internal architecture behind both strategies is documented for contributors in
:doc:`/developer/batched_evaluation`, and the cost model driving the decisions above in :doc:`/developer/cost_model`.

To estimate the peak-memory impact of a ``BatchPolicy`` *before* running a real (and potentially expensive) calculation, SeQuant also
provides a dry-run meter (:func:`sequant::eval::dryrun::meter`) that replays the policy-selected scheduler against a
zero-data sizing backend, reporting the realized peak without allocating real tensors. Numerical accuracy still
requires a real evaluation. See :ref:`batched-evaluation-dry-runs`.
