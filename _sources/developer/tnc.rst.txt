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

The static ``TensorNetworkV3::canonicalize_graph(const Graph&)`` hands the constructed graph to bliss's ``canonical_form()`` (with the
``shs_fsm`` splitting heuristic), which returns a permutation of vertex ordinals that is invariant under the graph's automorphism group —
this permutation *is* the canonical form; two isomorphic networks (under the coloring/topology rules above) always produce the same one.

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
  :class:`sequant::TensorCanonicalizer` is a registry base class (``register_instance``/``instance_ptr``, keyed by tensor label) that a
  contributor can implement against to customize how one tensor's slots get ordered, without touching the network-wide bliss machinery
  above. ``DefaultTensorCanonicalizer::apply`` is the reference implementation; it deliberately reimplements sort as a bubble sort
  (rather than using ``std::sort``) because it needs to count the transposition parity, and the standard sort algorithms make no guarantee
  about using swaps to get there.

Debugging and tests
------------------------

``Logger::instance().canonicalize``, ``.canonicalize_input_graph``, ``.canonicalize_dot``, and ``.tensor_network``
(``SeQuant/core/logger.hpp``) trace, respectively, canonicalization calls, the input graph before canonicalization, a dot-format dump of
the canonicalized graph, and general tensor-network construction. ``tests/unit/test_tensor_network.cpp`` runs a shared suite
parametrized over ``V1``/``V2``/``V3``, plus version-specific sections covering isomorphism via ``canonicalize_slots``, sign/phase
tracking, named-index ordering, and particle-symmetric/asymmetric/number-nonconserving cases.
