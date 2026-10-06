Context and Configuration
==========================

Many SeQuant operations — normal ordering, Wick's theorem, the meaning of index labels, even how LaTeX output is typeset — depend on
settings that are impractical to pass explicitly to every function call. SeQuant instead keeps this configuration in an *implicit,
global context*: a global default that every function reads unless told otherwise. This page consolidates the two context
objects a user configures and the recommended way to set them up for many-body/quantum-chemistry work; the step-by-step manual setup of
index spaces is covered separately in :doc:`Getting started </user/getting_started/index_spaces>`.

The core ``Context``
---------------------

:class:`sequant::Context` bundles the settings that give meaning to an expression: the :class:`sequant::IndexSpaceRegistry` (the
vocabulary of index spaces in use, e.g. occupied/virtual), the ``Vacuum`` relative to which operators are normal-ordered
(``Vacuum::Physical`` — the true, particle-free vacuum — or ``Vacuum::SingleProduct`` — a single-determinant quasiparticle vacuum), the
``IndexSpaceMetric`` (whether the single-particle basis is orthonormal), and the ``SPBasis`` (spin-orbital vs. spin-free). It also
owns the :ref:`canonicalizer configuration <context-canonicalizer-configuration>`. It is accessed and replaced through
:func:`sequant::get_default_context`, :func:`sequant::set_default_context`, and :func:`sequant::reset_default_context`.

Constructing a ``Context`` from scratch and registering index spaces by hand, as shown in
:doc:`/user/getting_started/index_spaces`, is the right approach when a custom vocabulary of index spaces is needed. For standard
quantum-chemistry conventions, :func:`sequant::mbpt::load` is a one-line shortcut: it builds a ready-made
:class:`sequant::IndexSpaceRegistry` for one of a few common conventions (:class:`sequant::mbpt::Convention` — minimal,
single-reference, multi-reference, F12, ...; :class:`sequant::mbpt::SpinConvention` controls whether/how spin is tracked), sets it
and ``Vacuum::SingleProduct`` on a copy of the current default context, and installs that copy as the default context:

.. literalinclude:: /examples/user/context.cpp
   :language: cpp
   :start-after: start-snippet-1
   :end-before: end-snippet-1
   :dedent: 2

.. _context-canonicalizer-configuration:

Canonicalizer configuration
------------------------------

The ``Context`` owns the settings that govern how individual tensors are canonicalized (see :doc:`canonicalization`):

- the :class:`sequant::TensorCanonicalizer` objects, keyed by tensor label
  (:func:`sequant::Context::set_tensor_canonicalizer`). When tensors are canonicalized as part of a product (a tensor network), a
  tensor uses the canonicalizer registered for its own label, if any (:func:`sequant::Context::nondefault_tensor_canonicalizer_ptr`),
  and otherwise the built-in one that sorts bra/ket indices according to the tensor's symmetry. The canonicalizer keyed by the empty
  label (:func:`sequant::Context::tensor_canonicalizer_ptr`) is the one ``Tensor::canonicalize()`` applies to a lone tensor;
- the comparers that order indices and pairs of indices (:func:`sequant::Context::set_index_comparer`,
  :func:`sequant::Context::set_index_pair_comparer`);
- the *cardinal* tensor labels, which are given lexicographic preference during canonicalization
  (:func:`sequant::Context::set_cardinal_tensor_labels`).

Canonicalization reads this configuration from the default context for ``Statistics::Arbitrary`` only: a canonicalized product
mixes tensors that have no statistics with operators of either statistics, and must order all of them consistently. The configuration
of a context installed for a specific statistics is currently ignored. Keep it identical to that of the ``Statistics::Arbitrary``
context nonetheless (as :func:`sequant::set_scoped_modified_default_context` does), since a future version may consult the
statistics-specific context first.

These can be given up front through ``Context::Options`` or changed on an existing ``Context`` with the setters above; to change them
for the duration of a scope, use a scoped context (see below) with a modified canonicalizer configuration. The example
canonicalizes a product, since that is where a canonicalizer registered for a label is used (a sum canonicalizes each of its summands
separately):

.. literalinclude:: /examples/user/context.cpp
   :language: cpp
   :start-after: start-snippet-5
   :end-before: end-snippet-5
   :dedent: 2

The canonicalizer configuration is part of the ``Context`` value, so installing a new ``Context`` replaces it: one constructed from
``Context::Options`` or ``Context{}``, or the one :func:`sequant::reset_default_context` restores, carries the default configuration.
Set the canonicalizers and labels on the ``Context`` you install, or derive it from the current one
(``Context(get_default_context())``). :func:`sequant::mbpt::load` does the latter: it sets only the index space registry and the
vacuum on a copy of the current default context, so the rest of the configuration is kept. It installs the copy process-wide, hence
throws if the calling thread has a scoped context active, whose settings it would otherwise make permanent.
:func:`sequant::Context::set_cardinal_tensor_labels` takes the complete list of cardinal labels: the defaults (the reserved labels
for antisymmetrizers, symmetrizers and transpositions) are kept only if the list includes them, as the one returned by
:func:`sequant::mbpt::cardinal_tensor_labels` does.

The MBPT ``Context``
----------------------

Everything under :code:`SeQuant/domain/mbpt` (predefined operators, the :doc:`CC equation generator <cc>`, ...) additionally consults a
second, MBPT-specific context, :class:`sequant::mbpt::Context`. Its main content is an :class:`sequant::mbpt::OpRegistry`: a map from
operator labels (``"t"``, ``"f"``, ``"g"``, ...) to their :class:`sequant::mbpt::OpClass` (excitation, de-excitation, or general),
which is what lets SeQuant recognize e.g. that ``t`` is an excitation operator without inspecting its tensor form. It is configured
separately from the core ``Context``, via :func:`sequant::mbpt::set_default_mbpt_context`:

.. literalinclude:: /examples/user/context.cpp
   :language: cpp
   :start-after: start-snippet-2
   :end-before: end-snippet-2
   :dedent: 2

Most MBPT functions assert that this has been configured; forgetting this step is a common source of "OpRegistry is null"-type errors
when SeQuant code is first ported into a new program.

The other configurable field is ``CSV`` (default ``CSV::No``): when set to ``CSV::Yes``, operator tensors built by
:class:`sequant::mbpt::OpMaker` (e.g. cluster amplitudes ``t``) use cluster-specific virtuals — virtual-space indices that carry the
operator's occupied indices as proto-indices — instead of plain, independent virtual indices. :func:`sequant::mbpt::csv_transform`
expands such CSV-dependent tensors into an explicit basis (standard unoccupieds, PAOs, or AOs) when needed downstream.

.. literalinclude:: /examples/user/context.cpp
   :language: cpp
   :start-after: start-snippet-4
   :end-before: end-snippet-4
   :dedent: 2

Scoped context changes
------------------------

Both context types support `RAII <https://en.wikipedia.org/wiki/Resource_acquisition_is_initialization>`_-style, scoped overrides via :func:`sequant::set_scoped_default_context` (and its ``mbpt`` counterpart
:func:`sequant::mbpt::set_scoped_default_mbpt_context`): the returned resetter object restores the previous default context when it goes
out of scope, which is the safest way to temporarily change context for a single calculation without affecting surrounding code.
Scopes must end in the reverse order of their creation:

.. literalinclude:: /examples/user/context.cpp
   :language: cpp
   :start-after: start-snippet-3
   :end-before: end-snippet-3
   :dedent: 2

.. note::
   ``Context::Options`` fields are independent of one another and of the *previously active* context: constructing a new ``Context``
   (whether directly or through the scoped-override helpers) starts from the library defaults and applies only the fields explicitly
   given, rather than inheriting from whatever context happened to be active before.

Scoped contexts are per thread
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A scoped context is seen only by the thread that installed it and by the parallel work it launches through SeQuant: the parallel
primitives of ``SeQuant/core/runtime.hpp`` (``sequant::for_each`` and the others) make the calling thread's scoped contexts current in
their workers for the duration of the call. Threads created by the user see the process-wide default context instead.
:func:`sequant::set_default_context` changes that process-wide default; a thread that has scoped contexts active sees the change only
after those scopes end.

Because it replaces the process-wide default, a call of :func:`sequant::set_default_context` on any thread invalidates the reference
that :func:`sequant::get_default_context` returns on a thread without scoped contexts, and anything obtained from that reference by
reference. Code that keeps the context beyond a brief read, as the canonicalizers do for the duration of a canonicalization, holds a
copy made by :func:`sequant::get_default_context_snapshot` instead. The copy reflects every change of the process-wide default
completed before the call and is cheap: it shares the index space registry and the canonicalizer configuration with its source, and
takes the lock that guards the process-wide default only if that changed since the previous snapshot on the same thread.

Detecting changes
----------------------

Every ``Context`` has a version (:func:`sequant::Context::version`): a nonzero number, unique within the process, that changes whenever
the context is constructed, cloned or modified through one of its setters. :func:`sequant::current_context_version` returns the
version of the context in effect on the calling thread, which lets code that caches results derived from the context tell that the
cache is stale. The version tracks changes made through the ``Context`` interface; it does not track in-place mutation of an
:class:`sequant::IndexSpaceRegistry` shared with other contexts, nor of a canonicalizer or comparer object that the context refers to.
It identifies a context and its copies rather than their content: contexts that compare equal, such as two default-constructed ones,
may have different versions, so a cache keyed on the version can be invalidated without need, but never kept stale.
