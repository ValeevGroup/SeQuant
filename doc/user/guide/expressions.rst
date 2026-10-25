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

Since ``ExprPtr`` has pointer semantics, copying an ``ExprPtr`` does not copy the expression it points to: all copies refer to (and
modify) the same object. To obtain an independent copy, call ``clone()``.

``ExprHolder``: expressions as values
--------------------------------------

:class:`sequant::ExprHolder` stores a single expression with *value semantics*: copying an ``ExprHolder`` copies the held expression,
so no two holders ever share an expression and modifying one never affects another. It is a regular (``std::copyable``) type and can
thus be stored in standard containers and used with standard ranges and algorithms.

- An ``ExprHolder`` can be created from any expression object, e.g. ``ExprHolder h = Variable(L"x");``. Passing an rvalue moves the
  expression into the holder; creating one from an lvalue copies it and has to be spelled out explicitly (``ExprHolder h(expr);``).
- Creating an ``ExprHolder`` from an ``ExprPtr`` always copies the pointed-to expression (as other ``ExprPtr`` objects may still refer
  to it) and is therefore explicit as well. Conversely, an rvalue ``ExprHolder`` converts to an ``ExprPtr`` without copying.
- An ``ExprHolder`` converts to ``Expr &`` (``const Expr &`` if it is ``const``), so it can be passed to every function that takes an
  expression by reference; ``*h`` and ``h->`` give access to the held expression as well. Like ``ExprPtr``, it is iterable over the
  immediate subexpressions of the held expression.
- The arithmetic operators (``+``, ``-``, ``*``, and ``^`` for non-commutative products) work on ``ExprHolder`` and expression objects
  (as well as scalars) alike and return a new ``ExprHolder``; the operands remain untouched.
- After an ``ExprHolder`` has been moved from, the only valid operation on it is assigning a new expression to it.

``ExprHolder`` is a recent addition to SeQuant and most of its API still uses ``ExprPtr``.

Traversing an expression
------------------------

:func:`sequant::Expr::visit` calls a visitor on all subexpressions of an expression (with ``atoms_only = true``: only on its leaves).
The signature of the visitor determines what it may do:

- A visitor taking ``ExprPtr &`` can *replace* the visited subexpressions. Since the visited expression itself can't be replaced, it is
  not invoked on it.
- A visitor taking ``Expr &`` can *modify* the visited expressions in-place, and one taking ``const Expr &`` can only *inspect* them
  (it is the only kind of visitor accepted when visiting a ``const`` expression). Both are invoked on the visited expression itself
  as well.

Visitors taking ``const ExprPtr &`` are deprecated. They are still invoked on the visited expression itself, but receive a
*non-owning* pointer to it, which must not be retained beyond the visitor's invocation.

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

.. literalinclude:: /examples/user/expressions.cpp
   :language: cpp
   :start-after: start-snippet-1
   :end-before: end-snippet-1
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

The examples print expressions with ``<<`` and ``std::format``, which :file:`SeQuant/core/io/format.hpp` provides and which
default to SeQuant's text serialization format. Once built, an expression tree can be rendered to LaTeX or round-tripped through
that format; see :doc:`io` for both and for the format specifiers that select between them.
