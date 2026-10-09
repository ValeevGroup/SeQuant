Context and Configuration
==========================

Many SeQuant operations — normal ordering, Wick's theorem, the meaning of index labels, even how LaTeX output is typeset — depend on
settings that are impractical to pass explicitly to every function call. SeQuant instead keeps this configuration in an *implicit,
global context*: a global default that every function reads unless told otherwise. This page consolidates the two context
objects a user configures and the recommended way to set them up for many-body/quantum-chemistry work; the step-by-step manual setup of
index spaces is covered separately in :doc:`Getting started </user/getting_started/index_spaces>`.

The core ``Context``
---------------------

:class:`sequant::Context` bundles the settings that give meaning to an expression: the :class:`sequant::IndexBasisRegistry` (the
vocabulary of index spaces in use, e.g. occupied/virtual), the ``Vacuum`` relative to which operators are normal-ordered
(``Vacuum::Physical`` — the true, particle-free vacuum —, ``Vacuum::SingleProduct`` — a single-determinant quasiparticle vacuum —, or
``Vacuum::MultiProduct`` — a general reference state, for which :class:`sequant::WickTheorem` applies the *extended* form of Wick's
theorem, with density cumulants), and the ``SPBasis`` (spin-orbital vs. spin-free); whether a basis is orthonormal is a property
of the basis (:func:`sequant::IndexBasis::metric`), not of the context. It also owns the :ref:`canonicalizer configuration <context-canonicalizer-configuration>`. It is
accessed and replaced through :func:`sequant::get_default_context`, :func:`sequant::set_default_context`, and
:func:`sequant::reset_default_context`.

A ``Context`` owns its registry, which its copies share and which cannot change while any context uses it: a registry given by value
is moved or copied in, and one given by ``std::shared_ptr`` is adopted if that is its only owner (e.g. a temporary, such as the result
of :func:`sequant::mbpt::make_sr_spaces`, or a moved-from pointer) and copied otherwise. Modifying a registry after giving it to a
context therefore does not affect the context, unless the registry was moved in or adopted and the caller kept another way to reach
it, such as a pointer to one of its spaces (see the warning of ``Context::set``). To change the registry of a context, copy the
registry, modify the copy and set it on a copy of the context.

Constructing a ``Context`` from scratch and registering index spaces by hand, as shown in
:doc:`/user/getting_started/index_spaces`, is the right approach when a custom vocabulary of index spaces is needed. For standard
quantum-chemistry conventions, :func:`sequant::mbpt::load` is a one-line shortcut: it builds a ready-made
:class:`sequant::IndexBasisRegistry` for one of a few common conventions (:class:`sequant::mbpt::Convention` — minimal,
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
  label (:func:`sequant::Context::tensor_canonicalizer_ptr`) is the one ``Tensor::canonicalize()`` applies to a lone tensor in which no
  index occurs more than once (one with a repeated index, protoindices included, is a tensor network and is canonicalized as one);
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
(``Context(get_default_context())``). :func:`sequant::mbpt::load` does the latter: it sets only the index basis registry and the
vacuum on a copy of the current default context, so the rest of the configuration is kept. It installs the copy process-wide, hence
throws if the calling thread has a scoped context active, whose settings it would otherwise make permanent.
:func:`sequant::Context::set_cardinal_tensor_labels` takes the complete list of cardinal labels: the defaults (the reserved labels
for antisymmetrizers, symmetrizers and transpositions) are kept only if the list includes them, as the one returned by
:func:`sequant::mbpt::cardinal_tensor_labels` does.

.. _context-canonicalization-options:

Canonicalization options
^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``Context`` also carries the :class:`sequant::CanonicalizeOptions` that :func:`sequant::canonicalize`,
:func:`sequant::simplify` and everything built on them use (``Context::Options::canonicalization_options``,
``Context::set(CanonicalizeOptions)``, :func:`sequant::Context::canonicalization_options`): the canonicalization method, which
indices are *named* (external), and whether the labels of the named indices affect the result. Like the rest of the configuration
they are read from the ``Statistics::Arbitrary`` context only, by canonicalization and by the theorem machinery alike, whatever the
statistics of the expression. By default the named indices of an
expression are deduced as those that occur once in it; naming them explicitly makes an index that occurs more than once external,
which is what the theorem machinery relies on as well: the named indices of the context are the external indices of a
:class:`sequant::WickTheorem`, the ones it does not sum over (an operator sequence given to it directly follows the same rule: an
index that appears in it twice is summed over unless the context names it). When the context names indices it names *every*
external index of the expressions canonicalized under it, since any other index is then a dummy. There are no per-call options;
to canonicalize under other options, scope a context that carries them (see :ref:`below <context-scoped>`).

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
expands such CSV-dependent tensors into an explicit basis (standard unoccupieds, PAOs, or AOs) when needed downstream; the
basis is given as an :class:`sequant::IndexSpace`, or as an :class:`sequant::IndexBasis` registered under a name, whose
registry entry says whether it is orthonormal (:func:`sequant::IndexBasis::metric`).

.. literalinclude:: /examples/user/context.cpp
   :language: cpp
   :start-after: start-snippet-4
   :end-before: end-snippet-4
   :dedent: 2

Several bases of one space can meet in an expression: the canonical and the localized orbitals of a perturbation theory,
the cluster-specific virtuals of the ground-state and of the perturbed amplitudes, and so on. An :class:`sequant::Index`
can therefore carry an optional *basis instance*, an opaque integer written after its proto indices, ``a_1<i_1,i_2;1>``
(``a_1<;1>`` without proto indices), unless the instance is registered under a name, which is then printed in place of
the space's label (:doc:`../getting_started/index_spaces`); an index without one is in its space's own basis, as before.
The name is part of the basis (:func:`sequant::IndexBasis::name`), as a space's label is of the space: an index parsed
from it or minted by SeQuant carries it, and so does one given the instance by number, which is resolved through the
registry. An index built from an :class:`sequant::IndexBasis` with a bare instance number is in an unnamed basis of its
own, which is not the named one (the two indices are different), unless that basis is first looked up with
:func:`sequant::default_registry_resolved`.
:func:`sequant::mbpt::add_pao_basis` registers the PAOs this way, as the named instance ``μ̃`` of the particle space
(an index in it prints as ``μ̃_1``), an alternative to the separate PAO space of :func:`sequant::mbpt::add_pao_spaces`.
Instances are granted per operator label and leg space with :func:`sequant::mbpt::OpRegistry::grant_basis`,
:class:`sequant::mbpt::OpMaker` mints a granted operator's legs with them, and the projectors of the :doc:`CC <cc>`
equations carry the grants of the amplitude being solved for. Integrals are never granted: in Wick's theorem their legs
take the instance of the leg they are contracted with.

.. literalinclude:: /examples/user/context.cpp
   :language: cpp
   :start-after: start-snippet-6
   :end-before: end-snippet-6
   :dedent: 2

.. _context-scoped:

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
completed before the call and is cheap: it shares the index basis registry and the canonicalizer configuration with its source, and
takes the lock that guards the process-wide default only if that changed since the previous snapshot on the same thread.
:func:`sequant::get_default_index_basis_registry` reads only the index basis registry in the same way, without copying the context.

Code that reads the contexts at several points of one operation, as :func:`sequant::canonicalize`, :func:`sequant::simplify` and
``WickTheorem::compute()`` do, pins them instead with :func:`sequant::pin_default_contexts`: on a thread without scoped contexts the
returned resetter scopes a copy of the contexts in effect, with their versions, so that a concurrent
:func:`sequant::set_default_context` is seen only once the pin ends; under scoped contexts, which cannot change under the caller, it
installs nothing.

Detecting changes
----------------------

Every ``Context`` has a version (:func:`sequant::Context::version`) that identifies its canonicalization configuration: the index
space registry, the tensor canonicalizers and index comparers (compared as objects, not by behavior), the cardinal tensor labels, the
canonicalization options and the single-particle basis (which determines the symmetry of normal operators). Two contexts share a
version if and only if canonicalization sees the same configuration in both; a version is never reused for another configuration. The
other settings (vacuum, first dummy index ordinal, typesetting, deserialization defaults) do not affect it.
:func:`sequant::current_context_version` returns the version of the context in effect on the calling thread, which lets code that
caches canonicalization results tell whether they are still valid. Such a cache is keyed on the versions for all statistics, since
canonicalization reads the single-particle basis from the context for the statistics of the normal operator at hand and the rest
from the context for arbitrary statistics. The version does not track in-place mutation of a canonicalizer or comparer object that
the context refers to.
