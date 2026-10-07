Wick's Theorem
=================

:doc:`The getting-started page </user/getting_started/wick>` covers *using* :class:`sequant::WickTheorem` — building operator products and
reading its normal-ordered/overlap output. This page covers its configuration surface and the contraction algorithm behind it, for anyone
extending or debugging that algorithm.

Configuration
----------------

:class:`sequant::WickTheorem` is built from a :class:`sequant::NormalOperatorSequence` or a general :class:`sequant::Expr`; either way
an index that appears once in the input is external and a repeated one is a dummy, unless the ``Context`` names the external indices
(see :ref:`context-canonicalization-options`). It is then tuned via a handful of fluent setters before calling ``compute()``:

- ``full_contractions(bool)`` (default ``true``): full contractions only, versus all contractions including partial ones.
- ``use_topology(bool)`` (default ``true``): when the input is an ``Expr``, treats the ``Op``\ s within a bra or ket of a
  ``NormalOperator`` that can be swapped as a symmetry of the input (see :ref:`below <wick-topology>`), and ``NormalOperator`` objects
  attached to the same tensor label, as topologically equivalent, so that contractions related by this equivalence are not separately
  enumerated.
- ``set_nop_connections()`` / ``set_nop_avoided_connections()``: force, or forbid, contraction between specific pairs of normal-operator
  ordinals. Under a ``Vacuum::MultiProduct`` vacuum a required connection is enforced after cumulant expansion rather than during
  the contraction recursion, since a cumulant connects operators that no pair does; an avoided pair is still rejected the moment a
  contraction between it is attempted (see :ref:`below <wick-extended>`).
- ``set_nop_partitions()`` / ``set_op_partitions()`` / ``make_default_op_partitions()``: declare explicit equivalence groups of normal
  operators, or of individual ``Op``\ s, so that contractions related by permuting within a group are counted once with a combinatorial
  degeneracy factor rather than enumerated redundantly — the general form of what ``use_topology()`` infers automatically.

``compute(count_only, skip_input_canonicalization)`` then applies the theorem and returns a ``Constant``, ``Product``, or ``Sum``. It is
**not reentrant, but is optionally threaded internally** (a ``Sum`` input has its summands processed concurrently, one ``WickTheorem``
instance per summand, merged into a shared accumulator under a mutex — see below); it throws if the input's vacuum does not match the
current :class:`sequant::Context`'s.

The contraction algorithm
------------------------------

``compute()`` (``SeQuant/core/wick.impl.hpp``) expands a general expression input into normal-operator-sequence form and fully
canonicalizes it, unless ``skip_input_canonicalization`` is set; for a ``Sum`` it then recurses per summand in parallel before merging
results. Canonicalization leaves an expression it already canonicalized alone (see :doc:`tnc`), so an input the caller has simplified
costs nothing to canonicalize again; a Product the caller has not simplified is canonicalized topologically, which is slower than the
rapid canonicalization it used to get. After canonicalizing the input, ``compute()`` scopes a context that registers null canonicalizers
for normal operators; a mark recorded outside that context is not valid within it, and is valid again once ``compute()`` returns. The
actual enumeration happens in ``compute_nontensor_wick``/``recursive_nontensor_wick``: it walks pairs of ``Op`` s across the flattened operator
sequence, left to right.
Under ``full_contractions_`` (the default) it only ever extends a contraction starting from the leftmost still-free operator — the
standard recursive formulation of full-contraction enumeration; with it disabled, all pairs are considered.

Because canonicalizing the produced normal operators would undo their normalization, a top-level ``compute()`` installs, after
canonicalizing its input, a scoped :class:`sequant::Context` in which the normal-operator labels map to ``NullTensorCanonicalizer``
for every statistics; a canonicalizer registered by the user for those labels is shadowed for the duration, not removed. The parallel
workers of the per-summand ``WickTheorem`` instances see that scope.

A candidate pair ``(left, right)`` contracts (``can_contract``) iff ``left`` is a quasiparticle annihilator, ``right`` is a quasiparticle
creator, and their quasiparticle spaces intersect (``is_qpannihilator``/``is_qpcreator``/``IndexSpaceRegistry::intersection``, from
``SeQuant/core/op.hpp``); under ``Vacuum::MultiProduct`` their actions must also differ (see :ref:`below <wick-extended>`). The
contraction *value* (``contract``) depends on whether those quasiparticle spaces are pure hole/particle subspaces or not:

- both pure: the contraction is a plain overlap :math:`s(L,R)`;
- neither pure: :math:`\delta(L,l)\, s(l,r)\, \delta(R,r)`;
- only the left space pure: :math:`s(L,r)\, \delta(R,r)`;
- only the right space pure: :math:`\delta(L,l)\, s(l,R)`,

where :math:`l = L \cap H` and :math:`r = R \cap P` (:math:`H`/:math:`P` the hole/particle subspaces), materialized via temporary indices
from ``Index::make_tmp_index`` when a space needs projecting. The middle factor comes from ``WickTheorem::contraction_value``, which
returns :math:`s` except under ``Vacuum::MultiProduct`` (see :ref:`below <wick-extended>`). For fermionic statistics, each contraction
additionally picks up a sign of :math:`-1` for every operator that sat strictly between ``left`` and ``right`` at the time of
contraction — the usual anticommutation rule.

.. _wick-topology:

Topological pruning and degeneracy
---------------------------------------

When ``use_topology_`` is set, ``is_topologically_unique()`` restricts contraction attempts within a partition group (from
``use_topology()``'s automatic grouping, or an explicit ``set_op_partitions()``/``set_nop_partitions()``) to the group's first free
member, skipping contractions that would only reproduce one already tried by symmetry. The automatic grouping puts two ``Op``\ s
in one partition only if swapping them alone — their index vertices together with the tensor slots these occupy, every other
vertex of the tensor network graph fixed — is an automorphism; an automorphism that also moves other indices, tensors or operators
does not make the ``Op``\ s independently interchangeable. The corresponding combinatorial weight —
equivalent to a multinomial coefficient over how many operators from each partition have already been paired off — is recovered by
``op_permutational_degeneracy()``, so the pruned enumeration and the exhaustive one agree on the final coefficient. They also agree for
an input that vanishes by symmetry, i.e. has an automorphism of phase -1 (see :doc:`tnc`), because such an input is returned as zero
before any contraction is attempted: by the up-front canonicalization of a ``Sum`` input, or by the topology analysis for a ``Product``.

For the spin-free case over a fermionic vacuum, an additional factor of :math:`2^{n}` is folded in per completed contraction, where
:math:`n` is the number of creation/annihilation "partner" cycles formed by that contraction — the generalized Wick's theorem for
spin-free operators (attributed, in the surrounding comment, to Kutzelnigg).

.. _wick-extended:

The extended Wick theorem (``Vacuum::MultiProduct``)
--------------------------------------------------------

Relative to a general (multiconfigurational) reference the quasiparticle picture is lost: an operator on a partially occupied
(*active*) orbital is both a quasiparticle creator and annihilator. SeQuant follows Kutzelnigg and Mukherjee's *generalized normal
order*: strings are normal-ordered so that their reference expectation value vanishes, and the theorem acquires, besides the pair
contractions (now valued :math:`\gamma` for ``cre·ann`` and :math:`\eta = \delta - \gamma` for ``ann·cre``), *multi-leg* contractions
of :math:`k \ge 2` creators and :math:`k` annihilators valued by the :math:`k`-body density cumulant :math:`\kappa_k`.

The reference must commute with the number operator (a state of definite particle number, or an ensemble of such states): only
then does every string with unequal numbers of creators and annihilators have a vanishing reference expectation value and cumulant,
so that the blocks above are the only multi-leg contractions (a reference that superposes particle numbers, such as a Bogoliubov
vacuum, would need anomalous contractions such as :math:`\langle a_p a_q \rangle` as well). Kutzelnigg and Mukherjee define generalized normal order for
number-conserving strings; SeQuant extends it to every string by the same rule, under which a string is the sum, over all sets of
disjoint internal contractions (pairs and cumulant blocks), of their values times the normal-ordered remainder. Hence
:math:`\{a_p\} = a_p`, the reference expectation value of every normal-ordered string vanishes, and the theorem holds for products
of such strings; the ``GNO products match the core-vacuum path`` test in ``tests/unit/test_mbpt.cpp`` checks it, built from that
rule, against the standard theorem under the core vacuum.

Rather than a second engine, this is layered on the standard one, which changes in three places:

- the classifiers in ``SeQuant/core/op.hpp`` gain a ``MultiProduct`` branch in which the "hole" space is the registry's
  *reference-occupied* space :math:`R` and the "particle" space its *vacuum-unoccupied* space :math:`U`; the active space is
  :math:`R \cap U`, and an active ``Op`` is classified as both. ``can_contract`` therefore also requires opposite actions, so that only
  ``cre·ann`` (over :math:`R`) and ``ann·cre`` (over :math:`U`) contract;
- ``WickTheorem::contraction_value`` returns ``γ`` (``left`` a creator) or ``η`` (``left`` an annihilator), from
  ``SeQuant/core/density.hpp``, as the middle factor of ``contract``. A ``γ`` over :math:`R` thus stands for :math:`\delta` on the core
  plus :math:`\gamma` on the active space, and likewise an ``η`` over :math:`U` for :math:`\delta` on the virtual space plus :math:`\eta`
  on the active space. An index with protoindices must not reach the active space (nonorthogonal active orbitals are not supported);
- the engine skips the required-connectivity filter (``WickTheorem::pairwise_connectivity()``), since a cumulant can connect
  operators that no pair does; it keeps rejecting contractions between avoided pairs early, because no cumulant can undo a direct
  contraction, and ``cumulant_expand`` filters the cumulant-mediated connections of both kinds.

Under this vacuum ``WickTheorem::compute`` hands its input and options to ``detail::extended_wick``
(``SeQuant/core/wick_extended.hpp``), which drives the rest; ``count_only`` and bosonic statistics are not supported there. For each
input term it applies the term's own :math:`\delta`\ s and overlaps over summed indices, since identifying two indices is not a
contraction (one may leave an index shared by two operators, and one left carries the indices of at most one operator),
canonicalizes, records which input ``NormalOperator`` every operator index comes from (its *provenance*; an index shared by two
operators is renamed apart, the :math:`\delta` binding the names multiplying the result), and runs the standard theorem
(the private ``WickTheorem::compute_contractions``) with partial contractions on the bare operator sequence, multiplying the c-number
factors back in afterwards. In that run every operator index is external, so no term is canonicalized and no surviving index is
renamed, which keeps the provenance valid. Topological equivalence needs dummy indices, though, so with ``use_topology`` the
partitions are computed on the canonicalized term before its shared indices are renamed (renaming leaves every operator at its
ordinal) and passed to that run with ``set_nop_partitions()``/``set_op_partitions()``; the statistics of the run accumulate in the
dispatching ``WickTheorem``'s ``stats()``. The terms of a Sum input are handled one after another, each with its own topology
analysis, where the standard path contracts them in parallel. Each resulting term then

- has its mixed-space ``γ``/``η`` and surviving operators split into pure pieces: a ``γ`` into a core :math:`\delta` plus an active
  ``γ``, an ``η`` into a virtual :math:`\delta` plus an active ``η``, a survivor onto its active part and, with partial contractions,
  its core and virtual parts. A part that is not a registered space is split over its base spaces. The pieces are reduced with the
  operator indices kept fixed, so each projected index stays :math:`\delta`-bound to an input index and inherits its provenance;
- is handed to ``detail::cumulant_expand``, which groups the surviving *active* operators into disjoint blocks of :math:`k`
  creators and :math:`k` annihilators, :math:`2 \le k \le` ``WickTheorem::max_cumulant_rank``, with legs from at least two
  input operators. Each block is valued :math:`\kappa_k` (``detail::block_value``), the term takes the parity of the permutation that
  pulls each block's legs, in order, to the front of the survivor string, and any remaining operators (partial contractions only) stay
  as a ``MultiProduct``-vacuum ``NormalOperator``. Terms that miss a required ``nop_connections`` pair or realize a
  ``nop_avoided_connections`` pair are dropped; a factor the theorem produces (a :math:`\gamma`, :math:`\eta`, :math:`\kappa`,
  :math:`\delta` or overlap) connects all the input operators its indices came from, any other tensor (e.g. a coefficient) none.

Finally every one-body ``η`` is optionally rewritten as :math:`\delta - \gamma` (``WickTheorem::eta_as_delta_minus_gamma``; a
multi-body ``η`` of the input is kept), the
:math:`\delta`\ s over summed indices are applied, and the result is simplified. In it ``γ{ann;cre}`` and ``η{ann;cre}`` are
one-body, Hermitian and column-symmetric, ``κ{ann…;cre…}`` is of rank :math:`\ge 2`, antisymmetric, Hermitian and
column-symmetric, and all their indices are active. These symmetries take part in the tensor hash, so a density spelled
differently would not merge with the engine's; hence ``γ``, ``η``, ``κ`` and the spin-free ``Γ`` are reserved labels
(``reserved::density_labels()``): a ``Tensor`` carrying one must have the symmetries ``density::symmetries()`` gives for its
label and rank, or its constructor (and ``set_label``) throws, and the parser supplies them; the permutational symmetry of a
one-body density is void, so any spelled-out one is accepted and normalized away. The one alternative is a perm-nonsymmetric
multi-body ``γ``, ``η`` or ``κ``, a spin component such as spin tracing produces. The factories in
``SeQuant/core/density.hpp`` build densities with them.

The no-double-counting rule is that a block is only ever built from operators that survive the standard theorem: each extended term
descends from exactly one partial-contraction term (the one carrying its pair contractions), so no :math:`1/k!` weights are needed.
Spin-free evaluation is not supported for this vacuum (``WickTheorem::compute`` throws); the hooks for it are
``WickTheorem::contraction_value`` and ``detail::block_value``/``detail::term_weight`` in ``SeQuant/core/wick_extended.hpp``.

``mbpt::ref_av`` has no branch of its own for this vacuum, only forcing full contractions. There the mbpt operators (``ã``) are
normal-ordered relative to the reference, so it evaluates a different quantity than the ``SingleProduct`` path, which
normal-orders them relative to the core; the two agree for products of elementary operators, which the tests verify for results
with cumulants up to :math:`\kappa_4`. :func:`sequant::mbpt::decompositions::cumulants_to_densities` converts cumulants of
every rank to densities. The tests are in ``tests/unit/test_wick_extended.cpp``, ``tests/unit/test_wick.cpp`` and the
``MRSO-MultiProduct`` section of ``tests/unit/test_mbpt.cpp``.

Reducing the result
------------------------

Raw contraction output can contain chains of Kronecker deltas and overlaps introduced by the space-projection cases above.
``WickTheorem::reduce()`` turns these into an index-replacement map and substitutes it through the expression, collapsing delta chains
and eliminating deltas wherever the internal/external status of the indices they bind allows it. The result is then, like any other
:class:`sequant::Product`/:class:`sequant::Sum`, put into canonical form by :doc:`the tensor-network canonicalizer <tnc>` so that like
terms collect correctly — Wick's-theorem correctness therefore also rests on canonicalization being correct.

Debugging and tests
------------------------

``Logger::instance().wick_harness``, ``.wick_contract``, and ``.wick_reduce`` (``SeQuant/core/logger.hpp``) trace, respectively, the
top-level expand/canonicalize/recurse flow, individual contraction attempts, and the delta/overlap reduction step. The algorithm's test
coverage lives in ``tests/unit/test_wick.cpp``.
