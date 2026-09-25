Export Framework Architecture
=================================

:doc:`The user-facing export guide </user/guide/export>` covers the :class:`sequant::Generator` callback interface and lists the shipped
backends, noting that "the tree-walking, scalar-factor bookkeeping, and intermediate-reuse logic are handled once, centrally, for every
backend." This page documents that centralized machinery — the pipeline in ``SeQuant/core/export/export.hpp`` that turns a preprocessed
evaluation tree into a sequence of ``Generator`` callbacks — and what implementing a new backend concretely involves.

The export pipeline
-------------------------

``export_groups()`` (the driver behind the ``export_group()``/``export_expression()`` convenience wrappers) processes a range of
:class:`sequant::ExpressionGroup`\ s — named-or-unnamed collections of :class:`sequant::ExportNode` trees meant to be exported together,
e.g. as one function — in five stages:

1. **Validate.** If more than one group is given, the generator must report ``supports_named_sections()`` and every group must be named.
2. **Preprocess every tree** (``detail::preprocess_and_maybe_log`` / ``PreprocessVisitor``): scalar pruning, rebalancing, and
   intermediate renaming — see below — producing one ``PreprocessResult`` per tree (accumulated tensor/variable/index usage and reference
   counts).
3. **Global-scope declarations**: indices/tensors/variables that every tree agrees are known ahead of time are declared once, via
   ``generator.declare(...)``, for whichever kinds the generator reports at ``DeclarationScope::Global`` through
   ``index_declaration_scope()``/``variable_declaration_scope()``/``tensor_declaration_scope()``.
4. **Per group**: optional named-section bookkeeping, then ``Section``-scope declarations, then per tree: ``Expression``-scope
   declarations, ``generator.begin_expression()``, a tree walk, ``generator.end_expression()``.
5. ``generator.end_export()``.

The tree walk itself is driven by a ``GenerationVisitor`` in combined pre-/post-order: pre-order visits (leaves and, on the way down,
already-computed intermediates) call ``load_or_create()``, which decides ``create`` vs. ``load`` vs. ``set_to_zero`` from
``ExportContext``'s per-object ``LoadStrategy``/``ZeroStrategy`` and a ref-count map of what this walk has already loaded; post-order
visits call ``process_computation()``, which assembles the node's ``Product``/``Sum``/``Adjoint`` expression (multiplying back in any
pruned scalar factor) and calls ``generator.compute(...)``, then drops both children's reference counts, ``unload``-ing whichever drop to
zero. ``generator.persist()`` fires when a group's designated result node is reached.

Preprocessing: scalars, rebalancing, and renaming
------------------------------------------------------

``PreprocessVisitor`` does three separable jobs on each raw evaluation tree before it is walked (quoting its doc comment, which is the
best single explanation of *why*):

    "removing explicit appearances of scalar leafs. We don't want them to be represented in the tree. Instead, we keep track of them in
    a different way in order to be able to give scalar factors alongside the actual tensor contraction they are supposed to scale (this
    is necessary as there are backends which only support scaling in this context)

    rebalance the tree such that for any given non-leaf node, its left subtree is always larger (or equally large) than its right one.
    This ensures that we have to have the least amount of tensors loaded at the same time, when generating code for a backend which only
    supports stack-like memory allocations (e.g. when A is allocated before B, B must be deleted before A can be deleted).

    Rename intermediate tensors that have the same name and describe the same tensor block, which are required as two separate
    entities at the same time when evaluating the tree (thus a single tensor object is insufficient)."

How much scalar pruning actually happens is itself a capability the generator reports back via ``prunable_scalars()``
(``PrunableScalars::None``/``Constants``/``Variables``/``All``) — e.g. ITF only prunes ``Constants`` (it still needs ``Variable``
factors represented as tensors), while the Julia backends prune ``All``.

A related, narrower trick handles ``Sum`` nodes without creating an intermediate for every ``+``: rather than encode addition of a
non-leaf child explicitly, its computed result is "flushed downward" — the child ``ExportExpr`` takes on its parent's result
tensor/name directly (``ExportExpr::set_expr``) — so only leaf children still contribute an explicit ``compute()`` call, relying on
``+=`` semantics for the rest. A ``Sum`` node whose ``ComputeSelection`` selects neither child (``ComputeSelection::None``) exists purely
for tree connectivity and is skipped by ``process_computation()`` entirely.

Reuse vs. common-subexpression elimination
------------------------------------------------

Two different questions are easy to conflate here, and are handled in two different places:

- **"Has this same tensor/variable block already been computed in this run, and is a live copy still around?"** is export's job,
  answered during preprocessing. When ``detail::preprocess()`` sees a result it has produced before, it checks whether that earlier
  instance is still loaded: if so, the *new* occurrence is renamed (``detail::rename()``, appending a numeric suffix such as
  ``I`` → ``I2``) so the two can coexist, since both are needed live at once; if not, the same tensor object is simply reused and its
  ``ZeroStrategy`` is set to ``AlwaysZero`` — but note this only reuses the *name/slot*, not the *computation*: the value is still
  recomputed via a fresh ``compute()`` call. Within a single tree walk, ``GenerationVisitor``'s reference-count map additionally avoids
  redundant loads of the same already-live block across sibling subtrees.
- **"Do two different subtrees compute the same value at all, such that one of them can be elided entirely?"** — actual algebraic
  common-subexpression elimination — is *not* part of the export framework. It lives one layer down, in ``SeQuant-optimize``
  (``SeQuant/core/optimize/common_subexpression_elimination.hpp``, ``sequant::opt::cse``), operating on the same
  ``FullBinaryNode``/``EvalExpr`` tree shape export uses, and is meant to run *before* a tree ever reaches ``export_groups()``.
  ``utilities/external-interface/cse_step.cpp`` is the reference driver, running CSE as its own pipeline step ahead of the export step,
  filtered to subexpressions worth materializing (at least two tensors, so trivial subexpressions such as symmetrization terms are left
  alone).

This is the precise scope of the "intermediate-reuse logic ... handled once, centrally" statement in the user guide: bookkeeping around
already-decided intermediates, not the decision of which intermediates to introduce in the first place.

``GenerationOptimizer``: an optional peephole layer
----------------------------------------------------------

``GenerationOptimizer<MainGenerator, MainContext>`` wraps any ``Generator<C>`` and sits between the pipeline above and the real
backend. It buffers the sequence of lifecycle callbacks (``create``/``load``/``set_to_zero``/``unload``/``destroy``/``persist``/
``compute``) emitted for one expression, cancels immediately-adjacent ``load``+``unload`` pairs of the same object, and — within the
constraints of stack-like alloc/dealloc ordering (citing Saabas & Uustalu, ENTCS 190, 2007) — attempts to eliminate some further
redundant alloc/dealloc pairs by reordering. It changes only *when* objects are (re)loaded, never what gets computed, and is
deliberately conservative rather than cost-optimal about which reorderings it attempts. It is opt-in: wrap the generator
(``GenerationOptimizer<ItfGenerator<...>> generator(itf_generator);``) before passing it to ``export_groups()``.

Implementing a new backend
--------------------------------

Concretely, a new backend is:

- a ``Context`` type deriving from :class:`sequant::ExportContext` — or from :class:`sequant::ReorderingContext` to get
  cache-locality-aware index reordering (largest index space moved into the fastest-varying slot, for whichever ``MemoryLayout`` the
  backend uses) for free — carrying whatever naming/tagging state the target format needs. ``ItfContext`` (ITF/Molpro) is a heavier
  example, with per-space name and tag maps, explicit import maps for pre-existing tensors/variables, and its own index-batching
  storage; the Julia backends' contexts are much thinner, typically just a per-space tag and array-dimension-variable name.
- a :class:`sequant::Generator` subclass implementing its ~30 callbacks. Most are direct string rendering (``represent()`` for each of
  index/tensor/variable/constant/power, and the create/load/compute/unload/destroy/persist lifecycle per tensor and per variable); a
  handful are capability queries (``supports_named_sections()``, ``requires_named_sections()``, ``supports_index_batching()``,
  ``prunable_scalars()``, the three ``*_declaration_scope()`` getters) that steer what the centralized pipeline above does on the
  backend's behalf, rather than logic the backend implements itself.

:class:`sequant::TextGenerator` is the reference implementation to read first — its own doc comment describes it as "a dummy generator
producing a plain text representation of the code. Mostly intended for having a convenient backend for tests available but it is also
very useful for debugging or other cases in which a human-readable version of the code is required."

.. note::
   The adjoint/transpose case (``EvalOp::Adjoint``) currently exports only the index permutation, not elementwise conjugation: correct
   for a real-valued field, incomplete for a complex one, since the exported IR has no node to carry an explicit conjugation yet. A
   backend targeting complex arithmetic needs to be aware of this gap.

Debugging and tests
------------------------

``Logger::instance().export_equations`` optionally logs a TikZ rendering of each tree during preprocessing. The framework's test
coverage is data-driven: ``tests/unit/test_export.cpp`` runs every ``*.export_test`` fixture under ``tests/unit/export_tests/`` against
every registered generator, each fixture pairing a SeQuant expression with one expected-output block per backend format; a separate
``tests/unit/test_export_python.cpp`` covers the Python/NumPy/PyTorch backends. For a full real-world example of assembling
``ExpressionGroup``\ s, a ``Context``, and a (possibly ``GenerationOptimizer``-wrapped) ``Generator`` end to end — including running CSE
as a preceding pipeline step — see ``utilities/external-interface/export_step.cpp``.
