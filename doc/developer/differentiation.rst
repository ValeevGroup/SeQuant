Symbolic Differentiation
==========================

:doc:`The user guide page </user/guide/differentiation>` covers what :func:`sequant::deriv::diff` and
:func:`sequant::deriv::at_zero` compute, the convention of the result for a symmetric variable, and the scope rules;
this page assumes them. It describes how ``diff`` arrives at that result (``SeQuant/core/deriv/deriv.cpp``), and how
the tests under ``tests/unit/deriv`` decide that the result is right.

The engine is built from core calls only: it expands with :func:`sequant::expand`, writes the variable in core's
spelling of a single tensor with :class:`sequant::TensorBlockCanonicalizer`, counts indices with
``get_used_indices_with_counts``, renames with :func:`sequant::transform_expr`, and builds Kronecker deltas and
symmetrizers with :func:`sequant::make_kronecker`, ``make_antisymmetrizer`` and ``make_symmetrizer``. Its own parts are
the product rule, the scalar rule, the spellings of the variable, the pairing of a leaf's labels with a spelling's, the
factor for the like slots, and the checks.

The pipeline
--------------

``diff`` takes the variables one at a time. For each variable :math:`x` it runs three stages:

1. **Check and set up.** No calculus happens here. Products are distributed over sums (symmetrizers stay factors), every
   check runs and throws before the derivative by :math:`x` is taken, and in each term the summed indices that share a
   name with an index of :math:`x` (a slot index or a protoindex) are renamed to fresh temporary indices, so that the
   variable's names are free to be given to the leaf's summed indices. External indices are never renamed; a variable
   spelled with the name of one is rejected instead.
2. **Apply the differentiation rule.**

   .. math::

      \mathrm{diff}(E, x) = \sum_y s_y\, P_y\bigl[D_y(E)\bigr]

   where :math:`y` runs over the spellings of :math:`x` (next section) and :math:`s_y` is the sign of the operation that
   makes :math:`y` times the phase of core's spelling of :math:`x`, the sign that ``TensorBlockCanonicalizer`` returns.
   :math:`D_y` is the product rule over the factors of each term, with the elementary derivative of a leaf by the
   spelling :math:`y` (or the scalar rule for a ``Variable``); :math:`P_y` puts the symmetrizer over the like slots of
   :math:`y` on each term. The permutations of like slots are left to SeQuant's symmetrizers, which commute with
   differentiation: they are never expanded or enumerated.
3. **Clean up** with :func:`sequant::rapid_simplify`, which flattens, combines scalars and drops zero terms while
   leaving every index name as it is. :func:`sequant::simplify` is not used, for the reason the user page's note gives:
   it renames or reorders the indices a symmetrizer carries, and the result would no longer be an array over the
   variable's names.

The term of each spelling is the plain gradient by that spelling, not a transform of another term, so it is exact
whatever the other tensors and scalars of :math:`E` are. Each term carries its symmetrizer in the layout of its own
spelling (``Â{a;i}`` on one, ``Â{i;a}`` on another), so every contraction still pairs a bra slot with a ket slot.

The spellings of the variable
-------------------------------

The convention of the user page sums the derivative by each position over every *operation* of the variable: each
permutation of its slots that its declared symmetry allows, with its sign, alone and, where that is an operation,
followed by the exchange of bra and ket. The engine splits the operations in two. Those that permute like slots
(defined below) are summed by one symmetrizer factor. Every other operation writes the variable's indices in another
slot order; operations that differ only by a permutation of like slots give the same order, and each distinct order,
with the sign of an operation that makes it, is a *spelling* of the variable. The derivative is the sum of the rule
over the spellings.

The exchange of bra and ket is an operation where the variable has bra-ket symmetry, its adjoint spelling holds its
elements (``BraKetSymmetry::Symm``, or real index spaces) and its bra and ket have one nonzero rank. Otherwise, where
the adjoint spelling holds the elements, the two orders name one stored element, and the spellings of the adjoint are a
second family of spellings: so for a real tensor without bra-ket symmetry, whose adjoint carries another label
(``t⁺``), and for a bra-ket symmetric one with unequal bra and ket ranks, such as ``g{i;a,b}``.

The variable is written once in core's spelling of a single tensor, and its spellings are made from that; those of the
second family are made from the adjoint of that spelling as it is, since core's spelling of the adjoint of a bra-ket
symmetric tensor may be the tensor's own spelling again. The leaves of :math:`E` are compared with each spelling as they
are written; neither a leaf nor the whole expression is canonicalized. A leaf *meets* a spelling if it is an occurrence
of the variable in that spelling, as the user page defines it. Compared as written, a leaf over a larger index space
meets every spelling whose elements it holds and gives a term for each: ``h{p_1;p_2}`` holds the bra-ket symmetric
``h{i;a}`` both at ``(i, a)`` and at ``(a, i)``. Respelling a leaf first would also break the layout of the symmetrizer.
Core's spelling of a bra-ket symmetric leaf may exchange its bra and ket; the respelled leaf then matches the spelling
whose bra is the written leaf's ket, and the symmetrizer laid out by that spelling holds in its bra the indices that the
rest of the term holds in its bra, so a contraction would pair two bra slots.

``at_zero`` uses the same spellings, of each listed tensor and of its adjoint whatever the index spaces: it drops a
term with a leaf that meets some spelling and lies inside it, each slot within the index space of the spelling's slot,
and it rejects a leaf that meets some spelling and lies inside none, since setting the variable to zero would leave
that leaf neither zero nor unchanged.

Summed and external indices
-----------------------------

Whether a label of a leaf is renamed or keeps a Kronecker delta depends on whether it is summed in its term, and the
engine reads that from the term alone, with no list of external indices from the caller. Over the factors of the term
other than ``Â`` and ``Ŝ``, an index is summed if it stands in exactly two slots (bra, ket and auxiliary counted
together). An index in fewer slots is external: one in a single slot, and one that stands in no slot and appears only
as a protoindex, such as the pair label of a trace over that pair's virtual space. A protoindex of an external index is
external too, in however many slots it stands. An index in more than two slots that is not external makes the term no
valid expression, and the input is rejected; one that is external, a protoindex of an external index, is a label its
carrier fixes, as in a first derivative differentiated again. An auxiliary index counts like any other. The
symmetrizers are left out because their indices are the external indices of a residual: counted, they would make each
external index of ``Â{i_1;a_1} v{a_1,i_2;i_1,a_2}`` look summed.

The elementary derivative
---------------------------

All pairing of indices happens in one function, ``elementary``. The labels of a matching leaf are paired with those of
the spelling :math:`y`:

- the index of each slot with :math:`y`'s index in the same slot;
- a protoindex that is itself a slot index of the leaf: it is paired already, through that slot (:math:`y` must lay
  out its protoindices the same way, else the input is rejected);
- a protoindex that is not a slot, such as the pair label of a per-pair rotation, with :math:`y`'s protoindex in the
  same place of the list. A protoindex list is unordered, so every order of each distinct list of :math:`y` is taken
  and the terms are summed, over the Cartesian product of :math:`y`'s distinct lists. Where two labels of :math:`y`
  coincide, two of these terms coincide, which is the factor the convention gives an element that a symmetry maps to
  itself.

For each pair, a label that is summed in the term and whose basis includes :math:`y`'s is renamed to :math:`y`'s label
everywhere else in the term: this is the Kronecker delta summed away. Otherwise a slot index keeps the Kronecker delta
as a factor: when the label is external, since renaming it would change the free indices of the result; when the leaf
holds it twice (a trace), since it has taken one of :math:`y`'s names already; and when it runs over less than
:math:`y`'s label does, since summing the delta away would widen the remaining sum to all of :math:`y`'s index space.
A protoindex label that is not a slot must be renamed through at least one of its pairs, and only its other pairs, as
in a trace, keep their deltas; when it is renamed through none, because it is external or runs over less than
:math:`y`'s label, the input is rejected (see below). A Kronecker delta between two disjoint index spaces is zero, and
so is the term that would carry it: the trace ``h{p_1;p_1}`` by ``h{i;a}`` is zero.

A derivative is not a contraction
-----------------------------------

The deltas are built with ``make_kronecker``, never with ``make_identity``. A contraction pairs a vector with a covector
through the metric, so two orbitals of different bases contract to their overlap, which is what ``make_identity``
returns. A derivative pairs a coordinate with itself, which is a Kronecker delta in any coordinate system: orthonormal,
non-orthogonal or overcomplete. The pair natural orbitals (a virtual basis computed for each occupied pair) of all
pairs form an overcomplete set, but the amplitudes are not one object expanded in it: each pair has its own tensor in
its own basis, and elements of different pairs are independent. The redundancy shows up elsewhere: as the overlaps the
function itself holds where one pair's amplitudes enter another pair's equation, and as linearly dependent
stationarity equations. Neither is a rule of differentiation.

The symmetrizer of the variable
---------------------------------

SeQuant's ``Â`` and ``Ŝ`` are normalized: expanded, ``Â`` over two bra and two ket slots is a quarter of the signed sum
of its four permutations. Two slots of one side that the variable's symmetry relates are *like* if they have one basis
(index space), carry the same protoindices, and are both or neither in each protoindex list of the variable; two
columns of a column-symmetric variable are like if their bra slots are like and their ket slots are like. The sum over
the :math:`n` permutations of like slots is therefore :math:`n` times the symmetrizer: :math:`|B|!\,|K|!\,\hat{A}` for
a group :math:`B` of like bra slots and :math:`K` of like ket slots of an antisymmetric variable, and
:math:`|C|!\,\hat{S}` for a group :math:`C` of like columns of a column-symmetric one. A group may be part of a side,
such as the two occupied slots of a bra ``(i, i, a)``. The condition on the protoindex lists keeps a permutation of like
slots from changing which slot a protoindex refers to: the ket slots ``i_5`` and ``i_6`` of
``t{a_5<i_5>,a_6<i_5,i_6>;i_5,i_6}`` are not like, since only ``i_6`` is outside the list of ``a_5``, so their exchange
is a spelling of its own.

What a symmetrizer over like slots cannot express is rejected rather than approximated: ``Symmetry::Symm``, for which
SeQuant has no symmetrizer; two groups of like slots on one side; and a column-symmetric variable with unequal bra and
ket ranks, whose unpaired slots core also permutes.

Three more checks guard the matching itself. A leaf with the variable's label whose slots could hold the variable's
elements and that declares other symmetries would be read, position by position, as the same array with other
symmetries. A leaf must lay out its protoindices as each spelling it meets does (the same slots as protoindices, one
list where the spelling has one), since ``elementary`` pairs the lists slot by slot. And a protoindex label
of a leaf that is not a slot is not renamed when it is external or runs over less than the variable's label, and a
Kronecker delta would put it into one more slot: an external label that stands in a slot of the term would then look
summed (SeQuant reads an index in exactly two slots as summed), and a summed one would stand in three slots. The engine
rejects every such label, also an external one that stands in no slot, rather than tell these cases apart. A term whose
Kronecker delta joins disjoint index spaces is zero and is dropped before this check.

The judge
-----------

The tests do not trust hand-written expected values alone. Every result is also checked numerically by *the judge*,
an evaluator in ``tests/unit/deriv`` that knows nothing about derivatives: it fills each tensor with random data over a
small basis, with exactly its declared symmetry, and addresses protoindices by value. It computes the tangents with dual
numbers, numbers :math:`a + b\,\epsilon` with :math:`\epsilon^2 = 0`: raising an element by :math:`\epsilon` and
evaluating gives the derivative exactly, as the coefficient :math:`b`. Core's evaluation backends are not used for this:
they are optional builds, the dense ones have no protoindices, and none gives tangents.

The judge states the convention as one identity. With :math:`\mathrm{tangent}(J)` the derivative of :math:`E` when the
element :math:`J` of :math:`x` and its symmetry partners are raised together, :math:`\mathrm{orbit}(J)` the distinct
partners, and :math:`|G|` the number of operations of :math:`x`, counted by the evaluator's own enumeration of the
arrangements of :math:`x`:

.. math::

   \mathrm{diff}(E, x)(J) \cdot |\mathrm{orbit}(J)| = |G| \cdot \mathrm{tangent}(J), \qquad
   \mathrm{diff}(E, x)(gJ) = \mathrm{sign}(g) \cdot \mathrm{diff}(E, x)(J)

The tests come in two layers. The rules layer holds hand-written one-line inputs for each rule and each rejection. A
rule test has the judge's verdict and, where the result is short enough to write, its exact result as well; the others,
such as a three-slot cycle, are decided by the judge alone. A rejection test checks the message. The identities layer
takes no hand-written input: SeQuant generates coupled-cluster singles and doubles equations (with canonical and with
cluster-specific virtual orbitals; spin-orbital, closed-shell, and open-shell for canonical orbitals), and the expected
value is an identity: the judge's identity for every equation by every amplitude and multiplier and by a sample of the
other tensors that ``diff`` accepts, the model's own equations as derivatives of its Lagrangian, a Taylor expansion
about zero that checks ``diff`` and ``at_zero`` together, mixed partial derivatives that commute, and the validity of
each result's free indices. Each identity is shown to fail on one deliberately wrong result.
