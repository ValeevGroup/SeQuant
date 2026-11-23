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
number of operations that leave the element unchanged. A protoindex list is unordered, so the reorderings of each
protoindex list of the variable whose labels are not its slots, such as the pair ``(i, j)`` of a per-pair rotation,
count among the operations: where the two labels of such a pair coincide, the two orders give the same element, and
the derivative there is twice that by the element.

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
times the coefficient of :math:`\varepsilon^n`. The example writes the series by hand; ``mbpt::series`` writes it for
every tensor of a label (see `Series substitution`_):

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
- **A valid function.** For a tensor variable, no index of a term of ``E`` stands in more than two slots, those of ``Â``
  and ``Ŝ`` aside and auxiliary slots counted like bra and ket slots, since SeQuant reads an index in exactly two slots
  as summed and one in a single slot as external; the exception is a protoindex of an external index, a label its
  carrier fixes. The derivative by a scalar ``Variable`` reads no index, so it takes such a term.
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

Series substitution
---------------------

:func:`sequant::mbpt::series` (``SeQuant/domain/mbpt/series.hpp``) writes the dependence on a perturbation strength
into an expression: ``series(E, label, ε, n)`` replaces each leaf :math:`p` labelled ``label`` by its Taylor series in
the strength :math:`\varepsilon` through order :math:`n`,

.. math::

   p \to \sum_{k=0}^{n} \frac{1}{k!}\, \varepsilon^k p^{(k)}, \qquad p^{(0)} = p.

The coefficient :math:`p^{(k)}` is the :math:`k`-th derivative of :math:`p` at :math:`\varepsilon = 0`; the factor
:math:`1/k!` belongs to the series, not to the coefficient. It is an ordinary tensor: the leaf with the order added to
its label as SeQuant's perturbation operators write it (``t¹``, ``t²``), and with the same slots and symmetries. The
:math:`k`-th derivative of the series at zero is therefore the coefficient itself:

.. literalinclude:: /examples/user/differentiation.cpp
   :language: cpp
   :start-after: start-snippet-8
   :end-before: end-snippet-8
   :dedent: 4

The label is matched exactly: ``t`` matches ``t`` of every rank and every block, but not ``t¹``. A leaf with the adjoint
label, ``t⁺``, is matched too. The strength is real, so the series of ``t⁺`` is the adjoint of the series of ``t``, with
the coefficients ``t¹⁺``, ``t²⁺``, ... that :func:`sequant::adjoint` gives. Nothing is expanded or truncated: ``diff``
expands, and ``at_zero`` removes what is left of higher order. A series through order :math:`n` holds no coefficient
above :math:`n`, so a derivative of order above :math:`n` lacks the terms of the higher coefficients.

The strength is a scalar :class:`sequant::Variable` or a *strength tensor*: a tensor with one auxiliary index and no
other slot, such as ``ε{;;z_1}``, which holds several perturbations of one kind, such as the :math:`3N` coordinates of
:math:`N` nuclei. The index space of its index is the caller's choice, for example the batching space ``z`` that
:func:`sequant::mbpt::add_batching_spaces` registers. With a strength tensor the term of order :math:`k` is
:math:`\frac{1}{k!}\, p^{(k)}_{\ldots; z_1 \ldots z_k}\, \varepsilon_{z_1} \cdots \varepsilon_{z_k}`: the perturbation
indices are appended to the leaf's auxiliary slots, as in ``h¹{p_1;p_2;z_3} ε{;;z_3}``, and are new for every leaf,
named unlike every index of the expression. A derivative is then taken by an element such as ``ε{;;z_4}``, under index
names of the caller's choice, and ``at_zero`` by the strength tensor sets all of its elements to zero.

``series`` throws :class:`sequant::Exception`, with a message starting with ``mbpt::series:``, for an empty or adjoint
label, a strength that is neither an unconjugated ``Variable`` nor a tensor with exactly one auxiliary index and no
other slot, an order above 9, an expression that holds a ``Variable`` labelled like the label or an ``mbpt``
operator (expand its tensor form first), a leaf of the label under a ``TensorPower``, and an expression that already
holds a tensor labelled like one of the coefficients it would write.

.. note::

   Auxiliary slots have no symmetry in SeQuant, so a coefficient is not known to be symmetric in its perturbation
   indices. A second derivative by a strength tensor therefore comes back as
   ``1/2 h²{p_1;p_2;z_3,z_4} + 1/2 h²{p_1;p_2;z_4,z_3}``, not as ``h²{p_1;p_2;z_3,z_4}``: correct, but not merged.

Hartree-Fock orbital and nuclear derivatives
----------------------------------------------

The generators in ``sequant::mbpt::hf`` (``SeQuant/domain/mbpt/models/hf.hpp``) return the Hartree-Fock energy as a
function of an orbital rotation and of a perturbation strength; its derivatives are then a few calls of ``diff``,
``at_zero`` and ``series``. The generators belong to ``SeQuant::mbpt``, which does not link ``SeQuant::deriv``; code
that differentiates links both, as ``SeQuant::SeQuant`` does.

No equation is solved, and nothing is assumed about the orbitals. A term is removed only by the mathematics of the
differentiation (a derivative that is zero) or by ``at_zero``, never because of a physical condition: the Brillouin
condition (:math:`F_{ia} = 0` on the Hartree-Fock solution), stationarity and the :math:`2n+1` rule (on the solution,
the orbital response through order :math:`n` determines the energy derivatives through order :math:`2n+1`) are not
applied. A term that vanishes on the Hartree-Fock solution stays in the result; using such a condition is a separate
step for the caller.

Orbital rotation
^^^^^^^^^^^^^^^^^^

:func:`sequant::mbpt::hf::rotated_energy` returns the energy of the reference under an orbital rotation, through the
given order in the rotation,

.. math::

   E(\kappa) = \langle 0 | e^{-(\hat K - \hat K^\dagger)} \hat H e^{\hat K - \hat K^\dagger} | 0 \rangle
   - \langle 0 | \hat H | 0 \rangle, \qquad \hat K = \sum_{ai} \kappa_{ai}\, a^\dagger_a a_i,

generated with SeQuant's commutator expansion and Wick's theorem. :math:`\hat H` is ``op::H()``, normal ordered, so the
reference energy is not part of the result. :math:`\hat K` takes every orbital of the reference-occupied space into
every orbital of its complement, frozen and inactive orbitals included. The tensor of its amplitudes :math:`\kappa_{ai}`
carries the label ``κ̂``; the hat only keeps it apart from ``κ``, core's label of the cumulant. With two elements ``x``
and ``y`` of the rotation, the orbital gradient is ``at_zero(diff(E, {x}), {x})`` (``at_zero`` by ``κ̂`` sets ``κ̂⁺`` to
zero too), and the blocks of the orbital Hessian are :math:`A = \partial^2 E / \partial\kappa\, \partial\kappa^*`,
``diff(E, {x, adjoint(y)})``, and :math:`B = \partial^2 E / \partial\kappa\, \partial\kappa`, ``diff(E, {x, y})``. The
energy through second order needs no ``at_zero`` for these:

.. literalinclude:: /examples/user/differentiation.cpp
   :language: cpp
   :start-after: start-snippet-9
   :end-before: end-snippet-9
   :dedent: 4

Over complex index spaces ``κ̂`` and ``κ̂⁺`` are independent variables (see `Complex and real variables`_); the
derivatives by ``adjoint(x)`` are separate calls and give the complex conjugates. Over real index spaces the same calls
collect the terms of ``κ̂⁺`` as well: the gradient is ``2 f{i_1;a_1}`` and both second derivatives are
:math:`2(A + B)`. For a closed-shell reference over real spaces, the energy put through
:func:`sequant::mbpt::closed_shell_spintrace` gives the gradient ``4 f{i_1;a_1}`` and the closed-shell Hessian, with the
variable spelled as the spin-traced leaves are, ``κ̂{a_1;i_1}:N-N-S``.

Nuclear derivatives
^^^^^^^^^^^^^^^^^^^^^

At a displaced geometry, displaced by the strength :math:`\varepsilon`, the integrals change and the reference orbitals
are no longer orthonormal: their overlap is :math:`S(\varepsilon)`, with :math:`S(0) = 1`. In the symmetrically
orthonormalized orbitals, with :math:`X = S^{-1/2}`, the one-electron integrals are :math:`X h X` and the two-electron
integrals carry a factor :math:`X` on each of their four legs. The operators then refer to orthonormal orbitals at
every geometry, and the energy keeps its form, the orbitals orthonormalized first and rotated second:

.. math::

   E(\varepsilon, \kappa) = \langle 0 | e^{-(\hat K - \hat K^\dagger)} \tilde H(\varepsilon)
   e^{\hat K - \hat K^\dagger} | 0 \rangle.

Two generators write this function. :func:`sequant::mbpt::hf::energy` is the energy through the given order in the
rotation, written in the integrals ``h`` and ``g``: the reference energy ``h{i_1;i_1} + 1/2 g{i_1,i_2;i_1,i_2}`` plus
``rotated_energy`` with each Fock matrix element written as :math:`h_{pq} + \sum_i g_{pi,qi}`. ``h`` takes the
symmetries of the ``f`` it replaces and ``g`` those of the ``g`` of ``op::H()``, so all ``h`` and all ``g`` of the
result are one variable each. :func:`sequant::mbpt::hf::orthonormalized` replaces each ``h`` and ``g`` by its series
in the strength and puts :math:`X` on each leg, as the series :math:`X = \sum_k \binom{-1/2}{k} \Delta^k` in
:math:`\Delta = S(\varepsilon) - 1 = \varepsilon s^{(1)} + \frac{1}{2} \varepsilon^2 s^{(2)} + \ldots`. The overlap
coefficients ``s¹``, ``s²``, ... have one slot on each side and the bra-ket symmetry of ``h``; as a variable, ``s¹`` is
spelled ``s¹{p_1;p_2}:N-C-S``. Each leaf is cut at the given order in the strength, so a derivative of higher order
than that is not the derivative of the full function. The nuclear-repulsion energy, a function of the geometry alone,
is the caller's to add.

The orbitals follow the displacement: the rotation is a function :math:`\kappa(\varepsilon) = \varepsilon \kappa^{(1)} +
\frac{1}{2} \varepsilon^2 \kappa^{(2)} + \ldots` with :math:`\kappa(0) = 0`, and the total derivatives of the energy
are those of :math:`E(\varepsilon, \kappa(\varepsilon))`. They are taken by substituting the series of ``κ̂`` with
``series``, removing its term of order zero with ``at_zero``, and differentiating by the strength. For the
:math:`n`-th total derivative every order passed is at least :math:`n`: ``energy(n)``, ``orthonormalized`` through
:math:`n` and the series of ``κ̂`` through :math:`n`. The first-order
rotation :math:`\kappa^{(1)}` is the solution of the first-order response (coupled-perturbed Hartree-Fock) equation,
the total derivative of :math:`\partial E / \partial \kappa^*` at zero:

.. literalinclude:: /examples/user/differentiation.cpp
   :language: cpp
   :start-after: start-snippet-10
   :end-before: end-snippet-10
   :dedent: 4

The nuclear gradient is the derivative of the integrals with the density, minus the derivative of the overlap with the
Fock matrix, which on the Hartree-Fock solution is :math:`-\sum_{ij} s^{(1)}_{ij} F_{ji}`, the overlap derivative
contracted with the energy-weighted density (the occupied block of :math:`F`), plus
:math:`\sum_{ai} (F_{ia} \kappa^{(1)}_{ai} + \kappa^{(1)*}_{ai} F_{ai})`. That last term, the orbital gradient times
the first-order rotation, is zero on the solution, where :math:`F_{ia} = 0`, and is kept. The response equation is
:math:`\sum_{bj} (A_{ai,bj} \kappa^{(1)}_{bj} + B_{ai,bj} \kappa^{(1)*}_{bj}) + b_{ai} = 0` with the textbook
:math:`A_{ai,bj} = F_{ab} \delta_{ij} - F_{ji} \delta_{ab} + \langle aj \| ib \rangle` and
:math:`B_{ai,bj} = \langle ab \| ij \rangle`, the transpose of the ``A`` above and the conjugate of its ``B``, since
this derivative is by :math:`\kappa^*`; :math:`b` holds the derivatives of the integrals and of the overlap. The second
derivative, ``at_zero(diff(Et, {ε, ε}), {ε})``, likewise holds every term: the second derivative at fixed orbitals,
the mixed terms in :math:`\kappa^{(1)}`, the orbital Hessian between two :math:`\kappa^{(1)}`, and the orbital gradient
times :math:`\kappa^{(2)}`, which again vanishes on the solution. With a strength tensor in place of ``ε`` the same
calls give the gradient and the Hessian over all coordinates at once:

.. literalinclude:: /examples/user/differentiation.cpp
   :language: cpp
   :start-after: start-snippet-11
   :end-before: end-snippet-11
   :dedent: 4

The generators expect what follows:

- **A variable is spelled as the leaves are.** ``diff`` requires the leaves of a variable's label to declare the
  variable's symmetries (see `Scope`_), so a variable carries the symmetry flags the generators give:
  ``κ̂{a_1;i_1}:A-N-S``, ``h{p_1;p_2}:A-C-S``, ``g{p_1,p_2;p_3,p_4}:A-C-S``, ``s¹{p_1;p_2}:N-C-S``. It is also
  written over the index spaces of the context the function was generated in: the examples use the minimal spaces of
  :func:`sequant::mbpt::make_min_sr_spaces`, occupied ``i`` and unoccupied ``a``; with the frozen-core spaces of
  :func:`sequant::mbpt::make_sr_spaces` the rotation runs from ``m`` to ``e``.
- **The operator registry holds** ``κ̂``, for ``rotated_energy`` and ``energy``.
  :func:`sequant::mbpt::make_legacy_registry` and :func:`sequant::mbpt::make_minimal_registry` register it as an
  excitation operator; a custom registry must add it.
- **A single-determinant reference.** ``rotated_energy`` and ``energy`` throw :class:`sequant::Exception`, with a
  message starting with ``mbpt::hf::``, for a vacuum other than ``Vacuum::SingleProduct``, a reference that is not the
  vacuum, a registry without an unoccupied space, and pair-specific orbitals (``CSV::Yes``); ``energy`` also for a
  spin-free basis.
- **Integrals in h and g.** ``orthonormalized`` needs no rotation, only the complete space. It throws for an expression
  that still holds ``f`` (it takes the result of ``energy``, not of ``rotated_energy``), ``h⁺`` or ``g⁺``, or an
  ``mbpt`` operator, and for a leg of ``h`` or ``g`` that carries protoindices (pair-specific orbitals) or a spin.
