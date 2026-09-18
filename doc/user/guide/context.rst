Context and Configuration
==========================

Many SeQuant operations — normal ordering, Wick's theorem, the meaning of index labels, even how LaTeX output is typeset — depend on
settings that are impractical to pass explicitly to every function call. SeQuant instead keeps this configuration in an *implicit,
thread-local context*: a global default that every function reads unless told otherwise. This page consolidates the two context
objects a user configures and the recommended way to set them up for many-body/quantum-chemistry work; the step-by-step manual setup of
index spaces is covered separately in :doc:`Getting started </user/getting_started/index_spaces>`.

The core ``Context``
---------------------

:class:`sequant::Context` bundles the settings that give meaning to an expression: the :class:`sequant::IndexSpaceRegistry` (the
vocabulary of index spaces in use, e.g. occupied/virtual), the ``Vacuum`` relative to which operators are normal-ordered
(``Vacuum::Physical`` — the true, particle-free vacuum — or ``Vacuum::SingleProduct`` — a single-determinant quasiparticle vacuum), the
``IndexSpaceMetric`` (whether the single-particle basis is orthonormal), and the ``SPBasis`` (spin-orbital vs. spin-free). It is accessed
and replaced through :func:`sequant::get_default_context`, :func:`sequant::set_default_context`, and
:func:`sequant::reset_default_context`.

Constructing a ``Context`` from scratch and registering index spaces by hand, as shown in
:doc:`/user/getting_started/index_spaces`, is the right approach when a custom vocabulary of index spaces is needed. For standard
quantum-chemistry conventions, :func:`sequant::mbpt::load` is a one-line shortcut: it builds a ready-made
:class:`sequant::IndexSpaceRegistry` for one of a few common conventions (:class:`sequant::mbpt::Convention` — minimal,
single-reference, multi-reference, F12, ...; :class:`sequant::mbpt::SpinConvention` controls whether/how spin is tracked) and installs it
as the default context with ``Vacuum::SingleProduct``:

.. literalinclude:: /examples/user/context.cpp
   :language: cpp
   :start-after: start-snippet-1
   :end-before: end-snippet-1
   :dedent: 2

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

Scoped context changes
------------------------

Both context types support `RAII <https://en.wikipedia.org/wiki/Resource_acquisition_is_initialization>`_-style, scoped overrides via :func:`sequant::set_scoped_default_context` (and its ``mbpt`` counterpart
:func:`sequant::mbpt::set_scoped_default_mbpt_context`): the returned resetter object restores the previous default context when it goes
out of scope, which is the safest way to temporarily change context for a single calculation without affecting surrounding code:

.. literalinclude:: /examples/user/context.cpp
   :language: cpp
   :start-after: start-snippet-3
   :end-before: end-snippet-3
   :dedent: 2

.. note::
   ``Context::Options`` fields are independent of one another and of the *previously active* context: constructing a new ``Context``
   (whether directly or through the scoped-override helpers) starts from the library defaults and applies only the fields explicitly
   given, rather than inheriting from whatever context happened to be active before.
