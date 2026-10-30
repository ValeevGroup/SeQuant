Expression Types
=================

This page covers how SeQuant identifies and orders the concrete types of an expression tree, for contributors adding a
new type derived from :class:`sequant::Expr`, how a broadcast sum under a ``ResultExpr`` head is lowered, and how a system of
``ResultExpr`` equations is inlined into one evaluation tree per equation. It assumes the vocabulary of
:doc:`The Expression Tree </user/guide/expressions>`, including its account of which leaves name a head.

Type ids
-----------

A concrete ``Expr`` type implements ``Expr::type_id()``, usually by returning :func:`sequant::Expr::get_type_id` of
itself; ``CProduct`` and ``NCProduct`` instead report the id of their base ``Product``. Type ids serve two purposes:
equal ids mean "same type" to ``operator==`` and ``Expr::operator<``, and ``Expr::operator<`` orders expressions of
*unlike* types by their ids. The latter is visible in results wherever expressions of unlike types are sorted with
``operator<``. ``Expr::is<T>()`` also compares against ``get_type_id<T>()``, so the id of a type can be computed even if
no object reports it.

An id is computed once per type from two things the type may declare:

- ``static constexpr Expr::type_rank_type type_rank``: the coarse position of the type among unlike types (lower ranks
  sort first); a type that does not declare it gets ``Expr::default_type_rank``, and a value that does not fit is a
  compile-time error;
- ``static constexpr std::string static_type_name(std::type_identity<Self> = {})``, where ``Self`` is the declaring
  type: a name for the type; a class template composes it from its own name and the names of its template arguments,
  obtained via ``sequant::type_name_of<T>()``. SeQuant's types declare it ``constexpr``, which makes the name, and
  ``type_name_of``, usable in constant expressions; a non-``constexpr`` one is accepted too. A type that does not
  declare it gets a name derived from the compiler's spelling of the type.

The id is ``(rank << 56) | (fnv1a_64(name) >> 8)`` (``Expr::make_type_id``; ``Expr::type_rank_of`` recovers the rank),
so ranks order types first and the name only breaks ties within a rank. The ranks of SeQuant's own types are collected
in ``sequant::expr_type_rank`` (``SeQuant/core/expressions/expr.hpp``); in ascending order they are ``Tensor``,
``Product`` (which ``CProduct`` and ``NCProduct`` inherit), ``Constant``, ``Sum``, ``Variable``, ``Power``, ``TensorPower``, then every
type of default rank (e.g. ``NormalOperatorSequence`` and ``mbpt::Operator``), then ``BOperator``, ``FOperator``,
``BNOperator`` and ``FNOperator``. They are spaced so that a new type can be inserted between them.

Every ``Expr`` type in SeQuant declares its own name and a rank (``CProduct`` and ``NCProduct`` inherit ``Product``'s),
so its id, and hence the order of unlike types, is the same on every platform and in every program. A new type should
do the same: the compiler-derived name is deterministic for a given compiler, but differs between compilers, and so
does the relative order of types that rely on it. A derived type inherits its base's ``type_rank`` but not its name:
the ``std::type_identity<Self>`` parameter does not accept ``std::type_identity`` of a derived type, so a derived type
that declares no name of its own gets the compiler-derived one, and hence an id distinct from its base's.

Ids must be unique; ``get_type_id`` records each type's id and throws :class:`sequant::Exception` naming both types if a
second type arrives at an id already taken, for instance because two types declare the same rank and name. Types are
told apart by their mangled names (``typeid(T).name()``), so that a shared object built with hidden visibility, which
has its own ``std::type_info`` for a type, is not mistaken for a second type. Two types of the same name in unnamed
namespaces of different translation units have the same mangled name, so if both rely on the compiler-derived name they
silently share an id. A type in an unnamed namespace therefore declares ``static_type_name()`` if another translation
unit may define one of the same name.

Broadcast sums under a head
-------------------------------

Broadcast is a property of a ``ResultExpr``, not of a ``Sum``: only a head names the modes a summand lacks, and supplies their
slot (bra or ket) and extent. Hence :func:`is_valid(const ResultExpr&) <sequant::is_valid>` checks a Sum on the rhs per summand
(each summand's external indices are a subset of the head's, each carries every aux index of the head, since broadcast is over
bra and ket modes only, and together they cover it; a head index that a summand uses in any slot counts as supplied, since head
indices are never contracted), while ``is_valid`` of a bare ``Sum`` keeps requiring
uniform external indices. The symbolic layer keeps a broadcast Sum as written (``optimize(ResultExpr&)`` optimizes the summands
like any others); ``binarize(const ResultExpr&)`` wraps each summand that lacks head indices (by ``detail::supplied_result_indices`` on
the summand as written: unlike ``is_valid``, which strips symmetrizers first, it counts the slots of an Â or Ŝ factor,
because binarization keeps that factor as a leaf whose slots carry the head's indices) in a unary ``EvalOp::Broadcast`` node shaped like a ``Power`` node: the operand subtree on the left, a
``Constant`` sentinel on the right, and the head's indices, in slot order, as its canonical indices, so the Sum chain above it sees
uniformly indexed operands. The node's hash combines the operand's with the result's modes and the position of each operand mode
among them: one operand broadcast to different modes is a different value. At evaluation the ranges of the modes a node adds come
from ``BackendArrayOps::make_zeros`` on a rank-1 descriptor per mode, sliced to the current batch of a loop over that mode (by the
ordered executor, as the value's cell says it is sliced), and ``Result::broadcast`` fills the result; ``to_expr`` spells the node as
its operand. Export of the node is not supported yet.

Inlining a system of equations
-----------------------------------

``binarize`` of a system builds one inlining context over all equations (heads indexed by label) and binarizes each equation
through the same internal path as ``binarize(const ResultExpr&)``, with the context passed down. Inlining happens at the leaf:
where that path reaches a ``Tensor`` or ``Variable`` leaf that names a head, it binarizes the head's equation in place of the leaf
instead of making a leaf node, so lowering of a broadcast sum and inlining of further definitions apply to the subtree as to any
equation. A set of the equations under construction detects cycles. A consumed equation also keeps its own tree, and each use
binarizes its definition anew.

Before a definition is binarized at a use, its indices are renamed: the head's indices, slot by slot, to the leaf's (the head's
bra to the leaf's ket and vice versa for a swapped use of a bra-ket symmetric head), and every other index of its right-hand side
(a dummy, or every index under a scalar head) to a fresh index whose label occurs nowhere in the system, with proto indices
following the renaming. The subtree therefore cannot capture an index of its consumer or of another inlined subtree. (Renaming
goes through :func:`sequant::transform_expr`, which passes a scalar ``Power`` through like the other scalars.)

The root of an inlined subtree carries the renamed head and is flagged as an equation boundary
(``EvalExpr::is_equation_boundary``). Factor collection, which gathers the factors of a product chain to check and canonicalize the
contraction, stops at that node as it does at a ``Sum`` or unary node, rather than splicing the definition's factors into the
consumer's chain: the consumer contracts with the head's declared layout and symmetry, which is what the user authored, and its
legality checks see the head instead of, e.g., a Hadamard product with a dummy in a bra slot of two factors. For the same
reason :func:`sequant::to_expr` of a node above a boundary shows the head, not the expansion; the boundary node itself expands to
its right-hand side. The flag is not part of the node's hash, so an inlined subtree hashes like the same subtree written out,
and repeated uses of one head are recognized as one intermediate.

Export binarizes one ``ResultExpr`` at a time and never inlines; a tree from the system overload is not meant for export, since
its consumed heads would appear as intermediates without definitions.
