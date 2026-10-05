Expression Types
=================

This page covers how SeQuant identifies and orders the concrete types of an expression tree, for contributors adding a
new type derived from :class:`sequant::Expr`. It assumes the vocabulary of :doc:`The Expression Tree
</user/guide/expressions>`.

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
  type: a name for the type, usable in constant expressions; a class template composes it from its own name and the
  names of its template arguments, obtained via ``sequant::type_name_of<T>()``. A type that does not declare it gets a
  name derived from the compiler's spelling of the type.

The id is ``(rank << 56) | (fnv1a_64(name) >> 8)`` (``Expr::make_type_id``; ``Expr::type_rank_of`` recovers the rank),
so ranks order types first and the name only breaks ties within a rank. The ranks of SeQuant's own types are collected
in ``sequant::expr_type_rank`` (``SeQuant/core/expressions/expr.hpp``); in ascending order they are ``Tensor``,
``Product`` (which ``CProduct`` and ``NCProduct`` inherit), ``Constant``, ``Sum``, ``Variable``, ``Power``, then every
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
