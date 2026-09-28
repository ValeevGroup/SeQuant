Wick's Theorem
=================

:doc:`The getting-started page </user/getting_started/wick>` covers *using* :class:`sequant::WickTheorem` — building operator products and
reading its normal-ordered/overlap output. This page covers its configuration surface and the contraction algorithm behind it, for anyone
extending or debugging that algorithm.

Configuration
----------------

:class:`sequant::WickTheorem` is built from a :class:`sequant::NormalOperatorSequence` (indices assumed external unless stated otherwise)
or a general :class:`sequant::Expr` (repeated indices assumed dummy), then tuned via a handful of fluent setters before calling
``compute()``:

- ``full_contractions(bool)`` (default ``true``): full contractions only, versus all contractions including partial ones.
- ``use_topology(bool)`` (default ``true``): treats ``Op``\ s of the same type within a ``NormalOperator``, and — when the input is an
  ``Expr`` — ``NormalOperator`` objects attached to the same tensor label, as topologically equivalent, so that contractions related by
  this equivalence are not separately enumerated when only the fully-contracted result (the vacuum average) is wanted.
- ``set_nop_connections()`` / ``set_nop_avoided_connections()``: force, or forbid, contraction between specific pairs of normal-operator
  ordinals.
- ``set_nop_partitions()`` / ``set_op_partitions()`` / ``make_default_op_partitions()``: declare explicit equivalence groups of normal
  operators, or of individual ``Op``\ s, so that contractions related by permuting within a group are counted once with a combinatorial
  degeneracy factor rather than enumerated redundantly — the general form of what ``use_topology()`` infers automatically.

``compute(count_only, skip_input_canonicalization)`` then applies the theorem and returns a ``Constant``, ``Product``, or ``Sum``. It is
**not reentrant, but is optionally threaded internally** (a ``Sum`` input has its summands processed concurrently, one ``WickTheorem``
instance per summand, merged into a shared accumulator under a mutex — see below); it throws if the input's vacuum does not match the
current :class:`sequant::Context`'s.

The contraction algorithm
------------------------------

``compute()`` (``SeQuant/core/wick.impl.hpp``) expands a general expression input into normal-operator-sequence form and, for a ``Sum``,
canonicalizes it and recurses per summand in parallel before merging results. The actual enumeration happens in
``compute_nontensor_wick``/``recursive_nontensor_wick``: it walks pairs of ``Op`` s across the flattened operator sequence, left to right.
Under ``full_contractions_`` (the default) it only ever extends a contraction starting from the leftmost still-free operator — the
standard recursive formulation of full-contraction enumeration; with it disabled, all pairs are considered.

A candidate pair ``(left, right)`` contracts (``can_contract``) iff ``left`` is a quasiparticle annihilator, ``right`` is a quasiparticle
creator, and their quasiparticle spaces intersect (``is_qpannihilator``/``is_qpcreator``/``IndexSpaceRegistry::intersection``, from
``SeQuant/core/op.hpp``). The contraction *value* (``contract``) depends on whether those quasiparticle spaces are pure hole/particle
subspaces or not:

- both pure: the contraction is a plain overlap :math:`s(L,R)`;
- neither pure: :math:`\delta(L,l)\, s(l,r)\, \delta(R,r)`;
- only the left space pure: :math:`s(L,r)\, \delta(R,r)`;
- only the right space pure: :math:`\delta(L,l)\, s(l,R)`,

where :math:`l = L \cap H` and :math:`r = R \cap P` (:math:`H`/:math:`P` the hole/particle subspaces), materialized via temporary indices
from ``Index::make_tmp_index`` when a space needs projecting. For fermionic statistics, each contraction additionally picks up a sign of
:math:`-1` for every operator that sat strictly between ``left`` and ``right`` at the time of contraction — the usual anticommutation
rule.

Topological pruning and degeneracy
---------------------------------------

When ``use_topology_`` is set, ``is_topologically_unique()`` restricts contraction attempts within a partition group (from
``use_topology()``'s automatic grouping, or an explicit ``set_op_partitions()``/``set_nop_partitions()``) to the group's first free
member, skipping contractions that would only reproduce one already tried by symmetry. The corresponding combinatorial weight —
equivalent to a multinomial coefficient over how many operators from each partition have already been paired off — is recovered by
``op_permutational_degeneracy()``, so the pruned enumeration and the exhaustive one agree on the final coefficient.

For the spin-free case over a fermionic vacuum, an additional factor of :math:`2^{n}` is folded in per completed contraction, where
:math:`n` is the number of creation/annihilation "partner" cycles formed by that contraction — the generalized Wick's theorem for
spin-free operators (attributed, in the surrounding comment, to Kutzelnigg).

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
