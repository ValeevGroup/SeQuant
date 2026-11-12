Symbol Visibility
=================

All SeQuant libraries (``SeQuant-symb``, ``-eval``, ``-optimize``, ``-export``, ``-core``, ``-mbpt``, ``-eval-ta``) are compiled with
*hidden* symbol visibility (``CXX_VISIBILITY_PRESET hidden``, ``VISIBILITY_INLINES_HIDDEN``); this hides types (their RTTI and
vtables) as well as functions and variables. When they are built as shared libraries
(``BUILD_SHARED_LIBS=ON``), only explicitly annotated declarations are exported, on every platform. The embedded ``bliss`` is always a
static, hidden library.

Annotations
-----------

The macros are defined in the configured header ``SeQuant/api.hpp`` (see ``SeQuant/api.hpp.in``). Each library ``<MOD>`` has:

``SEQUANT_<MOD>_API``
   the public API of that library.

``SEQUANT_<MOD>_INTERNAL_API``
   not public API, but has to be exported for technical reasons: another SeQuant library, or inline/template code in a public header,
   uses it. It expands exactly like ``SEQUANT_<MOD>_API``; the separate name keeps the public API auditable by ``grep``.

``SEQUANT_<MOD>_TEMPLATE_DECL_API`` / ``SEQUANT_<MOD>_TEMPLATE_DEF_API``
   for explicit template instantiations that live in the library: the former goes on the ``extern template`` declaration in the
   header (which must precede any implicit instantiation), the latter on the ``template class``/``template`` definition in the
   ``.cpp`` file. (GCC/Clang want the attribute on the declaration, MSVC wants ``dllexport`` on the definition.)

Rules:

#. Functions and variables need an annotation only if they are *defined out of line* (in a ``.cpp`` file, or explicitly
   instantiated); header-only inline functions and function templates do not. Types are different: a type whose RTTI crosses
   libraries must be exported even if it is header-only, a class template through explicit instantiations (see `RTTI across
   libraries`_).
#. Annotate the first declaration, in the header. Out-of-line friend functions are not covered by a class-level annotation: annotate
   the ``friend`` declaration itself.
#. Export polymorphic classes and exception types at class level (``class SEQUANT_SYMB_API Foo``), so their vtable and RTTI are shared
   across libraries; the same goes for public classes whose out-of-line members are API. For internals of non-public classes, annotate
   the needed members only.
#. When one overload of a function is public, export the whole overload set.
#. Be conservative: being in namespace ``sequant`` does not make something public. Never export something just to make a unit test
   link (see below).

Process-wide state
------------------

Under hidden visibility every shared library gets its *own* copy of the function-local statics of inline functions, of inline static
data members, and of the statics of templates instantiated in it. Mutable state that must be process-wide (e.g. the thread's scoped
context overlays, the temporary-index counter, the number of threads, the ``Logger`` instance) therefore has to live in a single
library, behind an exported non-inline accessor. Likewise, keys/ids that are compared across libraries must not be addresses of
statics (use names, see ``Expr::get_type_id``, ``Result::id_for_type``, ``detail::implicit_context_key``).

RTTI across libraries
---------------------

A hidden type gets its own ``type_info`` in every library that uses it. libstdc++ and MSVC compare ``type_info`` by name, so this
goes unnoticed there, but libc++ (on Linux and on x86_64 macOS) compares by address: a ``dynamic_cast``, ``catch`` or ``typeid``
comparison then fails when the object was created, or the exception thrown, in another library. A type whose RTTI crosses library
boundaries must therefore be exported:

* a class is exported at class level with ``SEQUANT_<MOD>_API`` (or ``_INTERNAL_API``), which gives its typeinfo and vtable default
  visibility; this includes all exception types;
* a header-defined class template is exported through explicit instantiations in its library (see
  ``SEQUANT_<MOD>_TEMPLATE_DECL_API``), e.g. ``Operator<S>``, ``NormalOperator<S>``, ``NormalOperatorSequence<S>`` for both
  statistics;
* on Windows, ``dllexport`` of a class also exports, and hence fully instantiates, its template base classes. Not every class
  template supports that (e.g. ``container::svector``), so an exported class should hold a container as a member rather than derive
  from it (as ``Operator<S>`` and ``NormalOperatorSequence<S>`` do). Also avoid covariant return types in virtual members of exported
  explicit instantiations: MSVC rejects their thunks (C2908); and declare their copy/move constructors and assignment operators
  explicitly (``= default``): MSVC does not export implicitly declared move operations;
* the visibility of a template instantiation is also limited by that of its template arguments, so a user type used as such an
  argument (e.g. an irrep type) must not be hidden in the user's library either.

``Expr::is<T>``/``as<T>`` and ``Result::is<T>`` use name-derived type ids and do not depend on this. The unit test
``rtti-across-libraries`` checks the ``type_info`` identity of such types by address, which detects a missing export also with
libstdc++.

Third-party headers
-------------------

Hidden visibility also applies to the code of other libraries that a SeQuant library instantiates from their headers. Libraries
that keep process-wide state in statics of inline functions, or compare ``type_info`` objects, then get a private copy per
SeQuant library. TiledArray/MADNESS do both (e.g. the default ``World``), so SeQuant code must include TiledArray only via
``SeQuant/core/eval/backends/tiledarray/tiledarray.hpp``, which includes it with default visibility.

Unit tests
----------

In a shared build the unit tests link against the shared libraries and hence only see exported entities. Tests that exercise
non-exported (out-of-line) internals are compiled only in static builds: such test cases/sections are wrapped in
``#ifdef SEQUANT_UNIT_TESTS_INTERNALS`` (and test cases are tagged ``[internal]``); test files that exercise internals only are added
to their source list only ``if (NOT BUILD_SHARED_LIBS)``. CI covers both static and shared builds.
