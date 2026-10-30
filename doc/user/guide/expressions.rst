The Expression Tree
====================

Every symbolic object that SeQuant manipulates — a tensor, a scalar, a sum, a product, an entire many-body equation — is represented
as a node in a tree of :class:`sequant::Expr` objects. Understanding this small, shared vocabulary makes the rest of the library predictable:
the same handful of node types recur everywhere, from the simple examples in :doc:`Getting started </user/getting_started/index>` to the
equations produced by :doc:`the coupled-cluster machinery <cc>`.

``Expr`` and ``ExprPtr``
------------------------

:class:`sequant::Expr` is the abstract base of every expression node; concrete node types (described below) derive from it. Expressions
are always managed through :class:`sequant::ExprPtr`, a `smart pointer <https://en.wikipedia.org/wiki/Smart_pointer>`_ (interchangeable with
``std::shared_ptr<Expr>``) that adds the
arithmetic operators (``+``, ``-``, ``*``) used to build up expressions programmatically. An ``ExprPtr`` is itself iterable over its
immediate subexpressions, so an expression tree can be walked, matched, or rewritten generically without knowing the concrete type of
every node.

Leaves: ``Tensor``, ``Constant``, ``Variable``
-----------------------------------------------

:class:`sequant::Tensor` is the most common leaf: a labeled tensor quantity with a bra and a ket index set (and, optionally, auxiliary
indices), carrying its own permutational and particle symmetry. Registering the :doc:`index spaces <context>` those indices belong to is
covered separately.

Scalar quantities are represented by two different leaf types, depending on whether their value is known while building the expression:

- :class:`sequant::Constant` holds a fixed complex-rational number, e.g. the :math:`\frac{1}{2}` prefactor of a commutator expansion.
- :class:`sequant::Variable` holds a *named* scalar, such as a perturbation strength :math:`\lambda`, whose value is not fixed by the
  expression itself.
- :class:`sequant::Power` represents ``base^exponent`` for a rational exponent, where the base is a ``Constant`` or a ``Variable``.
- :class:`sequant::TensorPower` represents the element-wise ``base^exponent`` of a tensor leaf for a rational exponent. Build either kind of
  power with :func:`sequant::pow`, which refuses a compound tensor-valued base (name it as a ``ResultExpr`` and raise the resulting leaf). An odd
  integer exponent keeps the base's permutational symmetry, an even one makes an antisymmetric base symmetric, and a non-integer exponent on an
  antisymmetric base is rejected. A ``TensorPower`` evaluates numerically on every backend; code export of one is supported by the text, Python and Julia TensorOperations generators (the ITF, Julia ITensor and Julia TensorKit generators refuse it; see :doc:`export`).

.. literalinclude:: /examples/user/expressions.cpp
   :language: cpp
   :start-after: start-snippet-1
   :end-before: end-snippet-1
   :dedent: 2

.. literalinclude:: /examples/user/expressions.cpp
   :language: cpp
   :start-after: start-snippet-4
   :end-before: end-snippet-4
   :dedent: 2

Composition: ``Product`` and ``Sum``
-------------------------------------

:class:`sequant::Product` (a scalar times zero or more factors) and :class:`sequant::Sum` (zero or more summands) combine other
expressions. Both are **associative and flatten automatically**: multiplying a ``Product`` into another ``Product``, or adding a
``Sum`` into another ``Sum``, does not create a deeper tree — the factors/summands are merged into the same node. A bare ``Constant``
factor is likewise folded straight into a ``Product``'s scalar prefactor rather than kept as a separate factor. This is why the ordinary
``+``/``*`` operators on ``ExprPtr`` are usually all that is needed to build up even large equations: the tree stays as flat as
possible without any manual bookkeeping.

.. literalinclude:: /examples/user/expressions.cpp
   :language: cpp
   :start-after: start-snippet-2
   :end-before: end-snippet-2
   :dedent: 2

.. note::
   Two expressions that are mathematically equal are not necessarily represented by *identical* trees until they have been put into
   canonical form; see :doc:`canonicalization`.

``ResultExpr``: naming an equation
------------------------------------

An :class:`sequant::Expr` on its own has no notion of "what it is equal to" — it is only ever the right-hand side of an implicit
equation. :class:`sequant::ResultExpr` closes that gap by pairing an expression with an explicit left-hand side, either a
:class:`sequant::Tensor` (for a tensorial result such as a residual :math:`R^{a_1}_{i_1}`) or a :class:`sequant::Variable` (for a scalar
result such as an energy). This is the representation handed to :doc:`optimize() <optimize>` and to the :doc:`code generators <export>`,
since both need to know the external index pairing and symmetry of the result, not just the expression that computes it.

.. literalinclude:: /examples/user/expressions.cpp
   :language: cpp
   :start-after: start-snippet-3
   :end-before: end-snippet-3
   :dedent: 2

A tensorial expression that a later equation needs as a factor, a Hadamard (element-wise) product such as
:math:`\ell^{a_1 a_2}_{i_1 i_2} = M^{a_1 a_2}_{i_1 i_2} / D^{a_1 a_2}_{i_1 i_2}` in a perturbative correction, is composed through heads: it
becomes its own ``ResultExpr``, whose head declares the external indices and symmetries once, and the later equation uses that head
as a leaf via :func:`ResultExpr::result_as_tensor() <sequant::ResultExpr::result_as_tensor>` rather than respelling it. This is also
how a compound expression is raised to a power: name it, then apply :func:`sequant::pow` to the leaf.

.. literalinclude:: /examples/user/expressions.cpp
   :language: cpp
   :start-after: start-snippet-5
   :end-before: end-snippet-5
   :dedent: 2

Under a head, the summands of a :class:`sequant::Sum` need not carry the same external indices: a summand may omit some of the
head's indices, meaning that it is constant along those modes and is broadcast over them, as in the orbital-energy denominator
``D{i_1,i_2;a_1,a_2} = ε{;a_1} + ε{;a_2} - ε{i_1;} - ε{i_2;}``. :func:`is_valid() <sequant::is_valid>` accepts such a ``ResultExpr``
if every summand's external indices are a subset of the head's and together the summands supply all of them; a bare ``Sum``, which
has no head to name the missing modes, must still be uniformly indexed. Aux modes cannot be broadcast over: a summand must carry
every aux index of the head. :func:`sequant::binarize` and :func:`sequant::optimize` of
a ``ResultExpr`` lower the broadcast first (:func:`sequant::lower_broadcast_sum`): each short summand becomes an outer product with
one reserved all-ones leaf ``𝟙`` (:func:`sequant::make_ones`) per missing index. Numerical evaluation requests a ``𝟙`` leaf from the
leaf evaluator like any other leaf, so the leaf evaluator must serve it (recognize it with :func:`sequant::is_ones`) as an all-ones
tensor of its mode's extent; for code export see :doc:`export`.

.. literalinclude:: /examples/user/expressions.cpp
   :language: cpp
   :start-after: start-snippet-6
   :end-before: end-snippet-6
   :dedent: 2

Systems of equations
~~~~~~~~~~~~~~~~~~~~~~~

Equations composed through heads, such as ``D``, ``ℓ`` and a consumer of ``ℓ``, can be binarized together:
:func:`sequant::binarize` of a ``std::span<const ResultExpr>`` returns an :class:`sequant::InlinedSystem` with one evaluation tree
per equation, in input order, a ``consumed`` flag for each equation whose head another equation used, and ``roots()``, the trees of
the equations no other equation used. Each leaf that names a head is replaced by that head's equation, binarized as it would be on
its own (so a broadcast sum is lowered and nested definitions are inlined in turn); the equations may be given in any order.

A :class:`sequant::Tensor` leaf names a head if its label, its bra, ket and aux ranks and the index space of every slot agree with
the head's; a leaf with the same label but a different shape remains an ordinary leaf, served by the leaf evaluator. A leaf that
names a head must declare the head's symmetries, otherwise binarization throws; spelling the use through
:func:`ResultExpr::result_as_tensor() <sequant::ResultExpr::result_as_tensor>` guarantees that. A :class:`sequant::Variable` leaf
names a scalar head by its label; a conjugated use throws, so define the conjugate as its own equation. Two equations with the same
label and shape, and a definition that reaches its own head, are errors as well; the API reference lists every exception.

An inlined equation enters its consumer as one factor with the layout and symmetries of the authored head: the consumer contracts
with ``ℓ`` as declared, not with the factors of its right-hand side, so the head, rather than the factors that define it, is what
the consumer's legality (for instance the context's strict bra-ket rule) is checked against. Repeated uses of a head yield
subtrees with equal hashes and are recognized as the same intermediate like any repeated subexpression. Batched evaluation
(see :doc:`batching`) slices through an inlined equation like through any other subtree, so a batch loop opened in the
consumer also bounds the intermediates of the definition: a denominator built from orbital energies need not be materialized
whole. Code export does not use the system overload; export each ``ResultExpr`` in order (see :doc:`export`), and a later
equation reads an earlier head like any other tensor.

.. literalinclude:: /examples/user/expressions.cpp
   :language: cpp
   :start-after: start-snippet-7
   :end-before: end-snippet-7
   :dedent: 2

The examples print expressions with ``<<`` and ``std::format``, which :file:`SeQuant/core/io/format.hpp` provides and which
default to SeQuant's text serialization format. Once built, an expression tree can be rendered to LaTeX or round-tripped through
that format; see :doc:`io` for both and for the format specifiers that select between them.
