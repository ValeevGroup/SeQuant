Tensor Network Canonicalizer
===============================

:doc:`The canonicalization guide </user/guide/canonicalization>` motivates *why* canonical form matters and names the mechanism: a
colored-graph representation of the network, put into canonical form via the bundled `bliss <https://users.aalto.fi/~tjunttil/bliss/>`_
graph-automorphism library. This page picks up from there — how that graph is built and how its canonical form is translated back into a
concrete relabeling — for contributors working on :class:`sequant::TensorNetwork` (an alias for the current implementation,
``TensorNetworkV3``; ``V1``/``V2`` remain in the tree for comparison and are exercised by the same parametrized tests).

Graph construction
----------------------

``TensorNetworkV3::create_graph()`` turns the network into a ``bliss::Graph`` under one governing rule: **symmetries are encoded by
topology and color, and vertices that can be swapped must share a color** — so that a bliss automorphism of the graph corresponds exactly
to a permutation the network's declared symmetries allow. Concretely:

- Every index slot of every tensor gets a vertex; slots that may be freely permuted among themselves are bundled — a bra bundle, a ket
  bundle, a per-particle bra/ket ("braket") bundle for asymmetric tensors, a protoindex bundle for indices carrying protoindices — by
  introducing a vertex for the bundle and connecting it to each member's vertex. Each tensor additionally gets one "core" vertex, bundling
  its bra/ket (or braket) bundles.
- For a symmetric or antisymmetric tensor, all of its bra slot vertices share one color and all of its ket slot vertices share another —
  this coloring, not any special-cased logic, is what lets bliss's automorphism search permute bra (or ket) indices among themselves
  freely, and is therefore how permutational symmetry actually gets exploited. Column-symmetric (particle-symmetric) asymmetric tensors
  get matching colors across their braket-bundle vertices instead.
- ``Index`` vertices are added last and connected to whichever slot vertices reference them: an internal (contracted) index connects to
  exactly two slot vertices, an external one to a single slot vertex, and a shared auxiliary index (e.g. a Laplace-transform or
  density-fitting index) can connect to more than two. Indices are colored by their space plus their protoindices' colors, via
  ``VertexPainter`` (``SeQuant/core/tensor_network/vertex_painter.hpp``), which also deduplicates colors across a run. Vertex colors are a
  32-bit integer (``Graph::VertexColor``) because bliss maps them into RGB internally, which caps how many distinct colors a single graph
  can use.

Computing the canonical form
---------------------------------

The static ``TensorNetworkV3::canonicalize_graph(const Graph&, aut_hook)`` hands the constructed graph to bliss's ``canonical_form()``
(with the ``shs_fsm`` splitting heuristic), which returns a canonical labeling: a permutation of vertex ordinals that maps every network
isomorphic to this one (under the coloring/topology rules above) onto the same canonical graph. The labeling itself is determined only
up to the automorphisms of the graph, and which one bliss returns depends on the numbering of the input vertices. That is immaterial for
anonymous indices, which are renamed in canonical order below, but not for named ones: with ``ignore_named_index_labels`` all named
indices of a space share a color, an automorphism may exchange them, and the labeling would place their labels by input numbering
(`#666 <https://github.com/ValeevGroup/SeQuant/issues/666>`_). So, with labels ignored, the member ``canonicalize_graph`` selects among
the labelings the automorphisms allow: position by position in canonical order, each named-index position gets the index with the
smallest label (by the context's index comparer, as tensor canonicalizers order indices) that an automorphism fixing the positions
before it can bring there. The automorphisms fixing given vertices are found by rerunning bliss's automorphism search with those
vertices given colors of their own. The structure of the result does not depend on the labels; where the labels go does.

Translating that permutation back into an actual relabeling is the job of the (differently overloaded, same-named) *member* function
``canonicalize_graph(named_indices, ...)``. It walks the graph's vertices in canonical-rank order and, from each ``TensorBra``/
``TensorKet``/``TensorBraKet``/bundle vertex it encounters, records the canonical slot order (and, for antisymmetric tensors, whether that
reordering is an odd permutation) for the tensor that vertex belongs to. It then:

- relabels anonymous (dummy) indices via an ``IndexFactory``, in canonical-rank order, so isomorphic networks get identical dummy-index
  names;
- physically permutes each tensor's bra/ket/column slots per the recorded canonical order, accumulating a sign for antisymmetric
  permutations;
- sorts the tensors themselves by their core vertex's canonical rank, subject to a "do these tensors commute" relation, since a
  non-commuting product does not admit an arbitrary total order.

The by-product of all this sign bookkeeping is returned as ``nullptr`` (no sign flip) or ``ex<Constant>(-1)``. A lighter sibling,
``canonicalize_slots()``, builds and canonicalizes the same graph but stops short of physically reordering anything — it instead returns
``SlotCanonicalizationMetadata`` (a canonical named-index ordering plus the underlying ``bliss::Graph``, comparable via graph isomorphism)
for callers that only need to test two networks for equivalence, such as term matching or :doc:`Wick's theorem <wick>`.

The bliss call of the member ``canonicalize_graph`` also reports the generators of the automorphism group, which detect networks that
vanish by symmetry. An automorphism maps the network onto itself up to a phase, ``TensorNetworkV3::Graph::automorphism_phase()``: the
product, over the bra and ket bundles of every antisymmetric tensor (including fermionic normal operators), of the parity of the slot
permutation it induces. If that phase is -1 the network equals minus itself, so it is zero; e.g. in ``t{a1,a2;i1,i2}:S ã{a1,a2;i1,i2}``
the swap :math:`a_1 \leftrightarrow a_2` has phase :math:`(+1)(-1)`. Since the phase is a homomorphism of the group to
:math:`\{\pm 1\}`, a scored generator of phase -1 is enough to detect a zero. Generators that move a tensor, a named or external index, a
protoindex bundle, or an aux slot are not scored, so a zero can go undetected but is never invented. The phase is computed from the
graph alone: ``create_graph`` records, as it emits them, the slot vertices of each antisymmetric bundle, the index of each index vertex
and the vertices an automorphism must fix. When a generator of phase -1 is found the member ``canonicalize_graph`` returns
``ex<Constant>(0)`` and leaves the tensors as they were. :doc:`Wick's theorem <wick>` applies the same test to the input of its topology
analysis.

Topological vs. lexicographic canonicalization
----------------------------------------------------

``CanonicalizeOptions::method`` (``SeQuant/core/options.hpp``) selects between two independently useful strategies:

- ``Topological`` is the bliss-based procedure above — the only one guaranteed correct once the network contains indistinguishable
  tensors, at the cost of building and canonicalizing a graph.
- ``Lexicographic`` is a cheap label-based sort of tensors, slots, and indices. It produces a more human-readable ordering but is
  incomplete on its own: it can fail to recognize two networks as equal when they contain genuinely identical tensors.
- ``Complete`` (the default) runs ``Topological`` then ``Lexicographic``, combining correctness with a readable result; ``Rapid`` is
  ``Lexicographic`` alone, for when speed matters more than completeness.

Subtleties for contributors
--------------------------------

- ``BraKetSymmetry::Conjugate`` is not yet exploited by graph construction — indices of a conjugate-symmetric tensor are colored as
  ``Symm`` would be, since handling the associated complex conjugation correctly is, per the surrounding comment, "not entirely clear."
- Named vs. anonymous index handling is controlled by ``CreateGraphOptions::named_indices``/``distinct_named_indices``; a newer
  per-index "loop color" (``NamedIndexColorMap``, ``SeQuant/core/tensor_network/typedefs.hpp``) additionally lets same-space named
  indices belonging to different DAG-scope loops be told apart, which sliced/batched evaluation relies on.
- Auxiliary (``aux``) index slots have no permutational-symmetry support yet, both in graph construction and in the single-tensor
  ``DefaultTensorCanonicalizer::apply`` (marked with a ``TODO`` in both places).
- ``TensorNetworkV3::factorize()`` is unimplemented (aborts).
- Canonicalizing a *single* tensor's own bra/ket order — as opposed to a whole network — is a separate, deliberately pluggable concern:
  :class:`sequant::TensorCanonicalizer` is a base class that a contributor can implement against to customize how one tensor's slots get
  ordered, without touching the network-wide bliss machinery above. Instances are owned by the :class:`sequant::Context`, keyed by tensor
  label. Inside network canonicalization (``TensorNetworkV3::do_individual_canonicalization``) a tensor uses
  ``nondefault_tensor_canonicalizer_ptr(label)`` of a :func:`sequant::get_default_context_snapshot` taken once per network, i.e. the
  entry for exactly its own label, and otherwise the network's own canonicalizer (``DefaultTensorCanonicalizer`` or
  ``TensorBlockCanonicalizer``); the entry for the empty label is consulted only by ``Tensor::canonicalize()`` on a lone tensor. Since
  the lookup goes through the current context, a scoped context (see :doc:`/user/guide/context`) overrides it for the scope's duration,
  on the threads that see that scope.
  ``DefaultTensorCanonicalizer::apply`` is the reference implementation; it deliberately reimplements sort as a bubble sort (rather than
  using ``std::sort``) because it needs to count the transposition parity, and the standard sort algorithms make no guarantee about
  using swaps to get there.

Recorded canonical form
----------------------------

A fully canonicalized expression carries a mark that makes canonicalizing it again a no-op, so callers can request full
canonicalization without tracking whether it was already done. The free ``canonicalize()``, as well as ``Product::canonicalize()``,
``Sum::canonicalize()`` and ``Tensor::canonicalize()``, record their result with ``Expr::mark_canonical(opts)`` and return at once, with
no byproduct, for an expression that ``Expr::is_canonical(opts)``; rapid (``Lexicographic``-only) canonicalization never marks. Since the
summands of a ``Sum`` are marked by their own canonicalization, re-canonicalizing a ``Sum`` after some of its summands changed
canonicalizes only those.

A mark is valid only for the ``CanonicalizeOptions`` it was recorded under (named indices included) and only while the contexts in
effect canonicalize as they did when it was recorded, as reported by ``current_contexts_version()`` (``SeQuant/core/context.hpp``): it
combines the versions of the contexts in effect on the calling thread for all statistics. A context's version identifies its
canonicalization configuration (tensor canonicalizers, index comparers, cardinal tensor labels, canonicalization options, the index
space registry and the SP basis), so the value changes when a different configuration takes effect -- a default context is set or reset, a scoped one begins
or ends, or a setter of one of those runs -- but not for settings canonicalization does not read, and contexts with the same
configuration share it. It is a function of the contexts, not a counter, so marks recorded before a scoped context are valid again once
it ends, and marks recorded in the scope of one top-level ``WickTheorem`` are valid in that of the next. Canonicalization seals its result
with the value at which it started, so a change of the contexts while it runs leaves the result unmarked. In-place changes of a
canonicalizer or comparer object are not tracked (see ``Context::version()``).

The contract with ``Expr`` types is that every mutation of a node's own data clears its mark: ``Expr::reset_hash_value()`` does so, and
mutations that do not reset the hash, such as those of the ``Product`` scalar, call ``Expr::reset_canonical_mark()``; a new ``Expr``
type must do the same in each of its mutators. ``Operator`` and ``NormalOperatorSequence``, which are vectors of their elements, honor
it by resetting in their mutable element accessors (``operator[]``, ``at``, ``begin``, ``end``, ``push_back``, ``emplace_back``).
Mutation or replacement of a subexpression needs no such call: each node carries a stamp of its own data, and a mark digests the stamps
of the whole subtree, so merely iterating mutably, as ``visit()`` and ``expand()`` do, leaves marks intact. With assertions enabled,
checking a mark also revalidates the memoized hash of each leaf, which catches a leaf mutator that resets neither (the memoized hash
of a non-leaf is legitimately stale after in-place mutation of a subexpression). ``clone()`` keeps the mark.

Debugging and tests
------------------------

``Logger::instance().canonicalize``, ``.canonicalize_input_graph``, ``.canonicalize_dot``, and ``.tensor_network``
(``SeQuant/core/logger.hpp``) trace, respectively, canonicalization calls, the input graph before canonicalization, a dot-format dump of
the canonicalized graph, and general tensor-network construction. ``tests/unit/test_tensor_network.cpp`` runs a shared suite
parametrized over ``V1``/``V2``/``V3``, plus version-specific sections covering isomorphism via ``canonicalize_slots``, sign/phase
tracking, named-index ordering, and particle-symmetric/asymmetric/number-nonconserving cases.
