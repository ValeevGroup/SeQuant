Symbolic Differentiation
==========================

:func:`sequant::deriv::diff` differentiates an expression built by SeQuant with respect to a scalar variable or a
tensor, and :func:`sequant::deriv::at_zero` sets variables to zero. Together they turn a Lagrangian into its
stationarity equations and a function of a perturbation strength into its response terms, without leaving SeQuant's
symbolic form. Both are declared in ``SeQuant/core/deriv/deriv.hpp`` and built as the CMake target ``SeQuant::deriv``,
which ``SeQuant::SeQuant`` links (``SeQuant::core`` does not).

The examples below assume this setup. Each ends with ``require_equivalent``, which throws unless the result and the
expected expression have the same canonical form:

.. literalinclude:: /examples/user/differentiation.cpp
   :language: cpp
   :start-after: start-snippet-0
   :end-before: end-snippet-0
   :dedent: 2

The two functions
-------------------

``diff(E, {x_1, ..., x_n})`` returns :math:`\partial^n E / \partial x_1 \cdots \partial x_n`. The variables are taken
from left to right, so ``diff(E, {x, y})`` is ``diff(diff(E, {x}), {y})``, and a list built at run time gives a
derivative of any order. A variable is a scalar :class:`sequant::Variable` or a single tensor, a *leaf* of the
expression tree: one :class:`sequant::Tensor`, not a product or a sum, and neither a symmetrizer nor a Kronecker delta.

``at_zero(E, {x_1, ..., x_n})`` returns ``E`` with every listed variable set to zero. It expands ``E``, drops each term
that has a factor that is one of the variables or its adjoint, or a positive power of either, and finishes with the same
light clean-up as ``diff`` (see `The form of the result`_): ``(x + y) (2 z)`` at ``x = 0`` comes back as ``2 y z``.
For a tensor variable, that factor is an occurrence of it (defined below), whatever its index names, that lies inside
the variable: in a spelling of the variable that it matches, each of its slots runs within the index space of the
spelling's slot. A leaf that holds elements of the variable and elements outside it, such as ``k{p_1;i_1}`` for the
variable ``k{a_2;i_2}``, would be neither zero nor unchanged, and is rejected. With several variables, a leaf must lie
inside one of them: ``k{p_1;i_1}`` is rejected for the variables ``k{m_2;i_2}`` and ``k{e_3;i_3}`` together, although
together they hold all of its elements.

A scalar variable is complex: ``x`` and its adjoint ``x^*`` (:func:`sequant::adjoint`) are independent variables, so
the derivative by ``x`` treats ``x^*`` as a constant.

.. literalinclude:: /examples/user/differentiation.cpp
   :language: cpp
   :start-after: start-snippet-1
   :end-before: end-snippet-1
   :dedent: 4

Tensor variables
------------------

In a term of ``E`` (a product of factors), an index is *summed* if it stands in exactly two slots of the term's
factors, not counting the symmetrizers ``Â`` and ``Ŝ``; every other index of the term is *external*, including one
that appears only as a protoindex. A *protoindex* is an index that another index carries, written ``a_1<i_1,i_2>``: it
says which basis ``a_1`` runs over, such as the virtual basis specific to the occupied pair ``(i_1, i_2)``. For
differentiation it is one more label of a tensor element, like the slot indices.

The derivative of one tensor element by another is 1 if all their labels agree, the index in every slot and the
protoindices each index carries, and 0 otherwise: a product of Kronecker deltas, one per label. Where the label of the
leaf in ``E`` is summed (over at least the variable's index space), its delta is summed away and the index takes the
variable's name in the rest of the term; where it is external, the Kronecker delta (:func:`sequant::make_kronecker`)
stays as a factor. The free indices of the result are therefore those of ``E`` plus those of the variable, under the
names the variable is written with. Tensor elements are independent coordinates even where their indices run over
different, non-orthogonal bases, so a derivative never produces an overlap; any overlap between two bases is a tensor
written in ``E``.

The *spelling* of a tensor is the way it is written: its label and the order of its slots. A tensor with a declared
symmetry can be spelled in several slot orders that hold the same elements, up to a sign: the antisymmetric
``g{i_3,a_3;i_4,a_4}`` is also ``-g{a_3,i_3;i_4,a_4}``. ``diff`` compares each leaf of ``E``, as it is written, slot by
slot with each spelling of the variable. The leaf is an *occurrence* of the variable in that spelling if it has the
variable's label and numbers of bra, ket and auxiliary slots, and in each slot an index whose index space overlaps that
of the spelling's index in that slot, that is not in a basis instance other than the variable's (see :doc:`context`)
and that carries as many protoindices. ``g{a_1,i_1;i_2,a_2}``, for example, is an occurrence of ``g{i_3,a_3;i_4,a_4}``
in the spelling ``-g{a_3,i_3;i_4,a_4}``, so its terms take the factor -1. A leaf with the same label that is an
occurrence in no spelling is another tensor, hence a constant.

The variable's symmetry
^^^^^^^^^^^^^^^^^^^^^^^^^

A tensor with a declared symmetry has fewer independent elements than positions. ``diff`` returns the derivative with
respect to the *independent* elements: the derivative by each position, summed over the symmetry operations of the
variable, each with its sign. At an element that no operation but the identity maps to itself, this is the plain
derivative by that element; in general it is the derivative along the element and its symmetry partners, times the
number of operations that leave the element unchanged.

The result writes the part of this sum that permutes the variable's *like slots* as SeQuant's normalized
antisymmetrizer ``Â`` over them, times the number of permutations it averages over. Like slots are slots of one side
(bra or ket) that run over the same index space, carry the same protoindices, and are both or neither in each
protoindex list of the variable. A variable that is symmetric only under exchanging its columns
(``ColumnSymmetry::Symm``, a column being a bra slot with the ket slot under it) gets the symmetrizer ``Ŝ`` over its
like columns instead. Every other operation, such as exchanging an occupied with a virtual slot, or bra with ket,
writes the variable in another slot order; each such spelling gives terms of its own, with the operation's sign. The
symmetrizer stays one factor of the result and is never expanded; symmetrizers already in ``E`` stay as they are. For
example, the antisymmetric ``t{a_1,a_2;i_1,i_2}`` by ``t{a_3,a_4;i_3,i_4}`` is
``4 Â{a_3,a_4;i_3,i_4}`` times four Kronecker deltas, the signed sum of the :math:`2! \cdot 2! = 4` arrangements of the
doubles amplitude; in ``1/4 g t`` the 4 cancels the :math:`1/4`:

.. literalinclude:: /examples/user/differentiation.cpp
   :language: cpp
   :start-after: start-snippet-2
   :end-before: end-snippet-2
   :dedent: 4

A leaf over a larger index space than the variable can hold the variable's elements in more than one slot order. The
bra-ket symmetric ``h{p_1;p_2}`` over all orbitals holds each element of ``h{i_3;a_4}`` twice, at ``h{i_3;a_4}`` and,
with bra and ket exchanged, at ``h{a_4;i_3}``. The leaf is an occurrence of the variable in both spellings, and the
result has one term for each:

.. literalinclude:: /examples/user/differentiation.cpp
   :language: cpp
   :start-after: start-snippet-7
   :end-before: end-snippet-7
   :dedent: 4

With this convention the coupled-cluster amplitude equations are the derivatives of the Lagrangian
:math:`\mathcal{L} = \langle 0 | (1 + \hat{\Lambda}) \bar{H} | 0 \rangle` with constant 1: the derivative by
:math:`\lambda_k` is ``CC::t()[k]`` and the derivative by :math:`t_k` is ``CC::λ()[k]`` in spin-orbital form. For
closed-shell spin-traced equations the same holds once the derivative is put through SeQuant's biorthogonal transform,
which :func:`sequant::mbpt::closed_shell_CC_spintrace_v2` applies to the residuals.

.. literalinclude:: /examples/user/differentiation.cpp
   :language: cpp
   :start-after: start-snippet-3
   :end-before: end-snippet-3
   :dedent: 4

Complex and real variables
----------------------------

The adjoint of a tensor ``x`` (:func:`sequant::adjoint`) is spelled with bra and ket exchanged, and its label gains
``⁺`` unless ``x`` is bra-ket symmetric. Whether the adjoint spelling ``x⁺`` holds other numbers than ``x`` decides
what the derivative includes:

- ``x⁺`` holds other numbers when ``x`` has no bra-ket symmetry, or ``BraKetSymmetry::Conjugate``, over complex index
  spaces (SeQuant's default). ``x`` and ``x⁺`` are then independent variables: leaves spelled ``x⁺`` are constants,
  and the derivative by ``x⁺`` is a separate call, ``diff(E, {adjoint(x)})``.
- ``x⁺`` holds the same elements as ``x`` when ``x`` is ``BraKetSymmetry::Symm`` (in any field), or when the index
  spaces of ``x`` are real. The result then includes the derivative by the leaves spelled ``x⁺``: a real ``t`` collects
  the terms of ``t⁺`` as well, and the derivative of :math:`G_{p_2 p_1} h_{p_1 p_2}` by a bra-ket symmetric ``h`` over
  one index space is :math:`G + G^T`.

Reality belongs to the index spaces, not to a tensor (see :doc:`context`): over real spaces every tensor is real. A
frequency-dependent response amplitude is complex even over real orbitals, so dynamic response equations are generated
over complex spaces.

.. literalinclude:: /examples/user/differentiation.cpp
   :language: cpp
   :start-after: start-snippet-4
   :end-before: end-snippet-4
   :dedent: 4

Responses
-----------

``diff`` takes the dependence on a variable from ``E`` as written. To differentiate with respect to the strength
:math:`\varepsilon` of a perturbation, write every quantity that depends on it as its series in :math:`\varepsilon`,
differentiate, and set :math:`\varepsilon = 0` with ``at_zero``. The :math:`n`-th derivative at zero is :math:`n!`
times the coefficient of :math:`\varepsilon^n`:

.. literalinclude:: /examples/user/differentiation.cpp
   :language: cpp
   :start-after: start-snippet-5
   :end-before: end-snippet-5
   :dedent: 4

Scope
-------

Three rules keep every input within what ``diff`` defines:

1. **Explicit dependence only.** A leaf with another label is a constant. If a quantity depends on the variable through
   a definition or a series, substitute it before calling ``diff``, as in the response example above.
2. **No variable under an inverse.** A :class:`sequant::TensorPower` of a leaf of the variable is rejected; a
   ``TensorPower`` of a tensor that does not depend on the variable, such as a denominator of orbital energies, is a
   constant factor. Where a variable enters through an inverse or an eigenvalue condition, keep the condition as a
   residual with a multiplier.
3. **A protoindex is a label of the element.** A pair-specific tensor, such as an amplitude, a per-pair rotation or a
   multiplier, of pair ``(i, j)`` and the same tensor of pair ``(k, l)`` are different variables. A tensor that is one
   object expanded in pair-specific bases (an integral, the Fock matrix, the overlap) is not; differentiate by it after
   :func:`sequant::mbpt::csv_transform`, where it is a plain tensor.

``E`` must also be a c-number expression: one without operators, such as the vacuum average of an operator product.

Rejected input throws :class:`sequant::Exception` before the derivative by the offending variable is taken. The
message starts with ``deriv::diff:`` (or ``deriv::at_zero:``) and names the exact condition. The conditions are of these
kinds:

- **What a variable is.** A variable is a ``Variable`` or a single tensor leaf, not a symmetrizer or a Kronecker delta,
  and a tensor variable puts each index in one slot only.
- **The variable's symmetries on its leaves.** A leaf with the variable's label whose slots could hold the variable's
  elements must declare the variable's symmetries; otherwise it would be read, slot by slot, as the same array with
  other symmetries.
- **What a symmetrizer over like slots can express.** The variable is not ``Symmetry::Symm`` (``Ŝ`` exchanges columns,
  and SeQuant has no symmetrizer that permutes bra slots alone without a sign), it has at most one group of like slots
  on each side, and a column-symmetric variable has as many bra as ket slots.
- **Protoindex labels as well-defined sets.** The protoindices of a slot of the variable are either all slots of the
  variable or none of them, and form an unordered list; the lists of two slots are the same or share no label. An
  occurrence of the variable lays out its protoindices as the spelling it matches does: the same slots as
  protoindices, and one list where the spelling has one. A protoindex label of an occurrence that is not one of its
  slots is summed in its term, over at least the index space of the variable's label; a term whose Kronecker delta
  joins disjoint index spaces is zero and is dropped before this check.
- **A valid function.** No index of a term of ``E`` stands in more than two slots, those of ``Â`` and ``Ŝ`` aside and
  auxiliary slots counted like bra and ket slots, since SeQuant reads an index in exactly two slots as summed and one in
  a single slot as external; the exception is a protoindex of an external index, a label its carrier fixes.
- **The scope above.** ``E`` holds operators, or the variable stands under a ``TensorPower``. Also rejected is a
  variable written with the name of an external index of ``E``, as happens when a result is differentiated again by
  the variable under the same index names; write the second variable with other index names.

``at_zero`` runs the checks that concern the variable alone (what a variable is, ``Symmetry::Symm``, protoindices that
mix the variable's slots with other labels), the check on the symmetries a leaf declares and the one on how an
occurrence lays out its protoindices, and adds two of its own: a variable under a zero or negative power, where setting
it to zero is undefined, and a leaf that holds elements of a variable and lies inside none of the variables.

.. literalinclude:: /examples/user/differentiation.cpp
   :language: cpp
   :start-after: start-snippet-6
   :end-before: end-snippet-6
   :dedent: 4

The form of the result
------------------------

The result is a flat sum of products with its scalars combined (:func:`sequant::rapid_simplify`); like terms are not
collected. Its free indices are those of ``E`` and of the variable, under the names they are written with. A summed
index of ``E`` that shares a name with an index of the variable is first renamed to one of SeQuant's temporary indices,
so a result may carry temporary index names as summed indices; a caller who resets the temporary-index counter
(``Index::reset_tmp_index()``) between two ``diff`` calls on a kept result can then get a name clash. The result is not
canonicalized; call :func:`sequant::simplify` on the final result for SeQuant's canonical form.

.. note::

   :func:`sequant::simplify` treats the indices that a symmetrizer carries as summed: it renames them or reorders the
   symmetrizer's slots with a sign. ``Â{a_3,a_4;i_3,i_4} g{i_3,i_4;a_3,a_4}``, for example, comes back as
   ``Â{a_1,a_2;i_1,i_2} g{i_1,i_2;a_1,a_2}``. After ``simplify`` the result no longer carries the variable's index
   names.

How ``diff`` computes the result is described in :doc:`the developer guide </developer/differentiation>`.
