Coupled-Cluster Class
======================

Coupled-cluster (CC) theory is one of the most accurate and widely used quantum chemistry methods for describing electron correlation in molecular systems. It represents the many-electron wavefunction using an exponential ansatz:

.. math::

   |\Psi_{\text{CC}}\rangle = e^{\hat{T}}|\Phi_0\rangle

where :math:`|\Phi_0\rangle` is a reference determinant (typically Hartree-Fock), and :math:`\hat{T}` is a cluster operator that generates excited determinants. The cluster operator is typically expanded as:

.. math::

   \hat{T} = \hat{T}_1 + \hat{T}_2 + \hat{T}_3 + \ldots

where :math:`\hat{T}_n` generates :math:`n`-fold excited determinants. For computational tractability, the cluster operator is usually truncated. For example, CCSD includes only single and double excitations (:math:`\hat{T} = \hat{T}_1 + \hat{T}_2`).

The :class:`CC <sequant::mbpt::CC>` class provides a convenient interface for setting up and processing CC equations using SeQuant’s symbolic algebra engine. It supports various CC formulations, including traditional, unitary, and orbital-optimized ansätze.

Overview
--------

The :class:`CC <sequant::mbpt::CC>` class can be used to derive:

- Ground state amplitude equations
- λ (de-excitation) amplitude equations
- Equation-of-motion (EOM) CC equations for excited states
- Response equations for properties and perturbations

Expressions are generated in spin-orbital basis and can be post-processed using SeQuant's spin-tracing capabilities. See :ref:`cc-spin-tracing` for more details.


Ansatz Options
--------------

The :class:`CC <sequant::mbpt::CC>` class supports several CC ansätze through the :enum:`CC::Ansatz <sequant::mbpt::CC::Ansatz>` enum:

- ``Ansatz::T``: Traditional CC ansatz, where the wavefunction is represented as :math:`|\Psi_{\text{CC}}\rangle = e^{\hat{T}}|\Phi_0\rangle`. This is the standard approach used in most implementations.

- ``Ansatz::oT``: Orbital-optimized traditional ansatz. Singles amplitudes (:math:`\hat{T}_1`) are excluded from the cluster operator, with orbital optimization performed instead.

- ``Ansatz::U``: Unitary CC ansatz, where the wavefunction is represented as :math:`|\Psi_{\text{UCC}}\rangle = e^{\hat{T} - \hat{T}^\dagger}|\Phi_0\rangle`. Particularly useful for quantum computing applications.

- ``Ansatz::oU``: Orbital-optimized unitary ansatz, combines both unitary and orbital optimized ansätze.

Key Methods
-----------

Ground State Amplitudes
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   std::vector<ExprPtr> t(size_t commutator_rank = 4,
                          size_t pmax = std::numeric_limits<size_t>::max(),
                          size_t pmin = 0);

Derives the equations for the :math:`t` amplitudes (:math:`\langle \Phi_P|\bar{H}|\Phi_0 \rangle = 0`) up to specified excitation levels.

.. code-block:: cpp

   std::vector<ExprPtr> λ(size_t commutator_rank = 4);

Derives the equations for the :math:`\lambda` de-excitation amplitudes.

Coupled-Cluster Response
^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   std::vector<ExprPtr> tʼ(size_t rank = 1, size_t order = 1);
   std::vector<ExprPtr> λʼ(size_t rank = 1, size_t order = 1);

Derives perturbed amplitude equations for response theory calculations.

Equation-of-Motion Coupled-Cluster
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   std::vector<ExprPtr> eom_r(
       nₚ np, nₕ nh,
       const std::vector<std::size_t>& block_ranks = {},
       std::optional<UCCEOMAssembly> assembly = std::nullopt);
   std::vector<ExprPtr> eom_l(nₚ np, nₕ nh);

Derives equation-of-motion coupled-cluster (EOM-CC) equations for excited states. The ``eom_r`` method generates equations for the right eigenvectors, while ``eom_l`` generates equations for the left eigenvectors. Traditional CC always uses the connected H̄R product. For UCC, passing ``CC::UCCEOMAssembly::Commutator`` or ``CC::UCCEOMAssembly::ProjectedHbar`` explicitly selects the right-hand assembly. Without that argument, a Bernoulli expansion or non-empty ``block_ranks`` selects projected-H̄ assembly; otherwise, a BCH expansion selects commutator assembly.

For UCC, the optional ``block_ranks`` argument gives the per-block truncation orders as a row-major matrix over the EOM manifolds: nested-commutator order for BCH and :math:`\bar{H}^{k}` order for Bernoulli. Traditional CC does not support block ranks or projected-H̄ assembly; Bernoulli does not support commutator assembly.

Examples
--------

The following examples demonstrate how to use the :class:`CC <sequant::mbpt::CC>` class to derive CC equations for various ansätze and excitation levels.

From this point onward, assume the following namespaces are imported, and :class:`sequant::Context` and :class:`sequant::mbpt::Context` are set up as shown.

.. literalinclude:: /examples/user/cc.cpp
   :language: cpp
   :start-after: start-snippet-0
   :end-before: end-snippet-0
   :dedent: 2

Ground State CC Amplitude Equations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: /examples/user/cc.cpp
   :language: cpp
   :start-after: start-snippet-1
   :end-before: end-snippet-1
   :dedent: 2

EOM-CC Equations
^^^^^^^^^^^^^^^^

.. literalinclude:: /examples/user/cc.cpp
   :language: cpp
   :start-after: start-snippet-2
   :end-before: end-snippet-2
   :dedent: 2

Response and Perturbation Equations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: /examples/user/cc.cpp
   :language: cpp
   :start-after: start-snippet-3
   :end-before: end-snippet-3
   :dedent: 2

Advanced Usage
--------------

Truncating the Commutator Expansion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The similarity-transformed Hamiltonian is built by ``mbpt::lst``, see :ref:`mbpt-lst`. For traditional CC with two-body Hamiltonians, the commutator expansion is truncated at 4th order. However, for unitary CC or other Hamiltonians, you may need to explicitly set the commutator rank:

.. literalinclude:: /examples/user/cc.cpp
   :language: cpp
   :start-after: start-snippet-4
   :end-before: end-snippet-4
   :dedent: 2

For BCH expansions, ``Options::hbar_singles_comm_rank`` applies an additional singles-only similarity transform after the primary expansion. Unitary ansätze use :math:`\sigma_1 = T_1 - T_1^\dagger`, while non-unitary ansätze use :math:`T_1`. A zero rank disables this transform; Bernoulli expansions reject a positive rank.

.. _cc-hbar-connectivity:

Using :math:`\bar{H}` outside the CC class
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:func:`CC::hbar() <sequant::mbpt::CC::hbar>` is public, but for a non-unitary ansatz the expression it returns is **not** self-contained: each commutator is written as a connected product (see :ref:`mbpt-lst`), which only equals the commutator once you connect the operators when taking the expectation value. The class always does this for you — every non-unitary path passes ``default_op_connections()`` (or a superset) to ``ref_av``. A caller who does not will silently retain the disconnected terms:

.. code-block:: cpp

   using namespace sequant::mbpt;

   CC cc{2};

   // correct: op::ref_av defaults to default_op_connections()
   auto eqs = op::ref_av(op::P(nₚ(2)) * cc.hbar());

   // WRONG: empty connectivity keeps terms the commutator would have cancelled
   auto bad = op::ref_av(op::P(nₚ(2)) * cc.hbar(), {.connect = {}});

   // also correct when no additional singles transform is configured: ask lst
   // for the explicit commutator form, which needs no connectivity
   auto hbar = lst(op::H(), op::T(2), 4);

When ``hbar_singles_comm_rank`` is positive, reproducing a non-unitary :math:`\bar{H}` explicitly requires a second ``lst`` call for the configured singles transform.

Note that ``op::ref_av(expr)`` and ``op::ref_av(expr, {})`` are *not* the same call: the connectivity defaults to ``default_op_connections()`` only when the argument is omitted entirely, since ``EVOptions::connect`` is itself empty by default.

The two ``ref_av`` calls above swap roles for a unitary ansatz, so it is not simply "unaffected". There :func:`CC::hbar() <sequant::mbpt::CC::hbar>` already returns explicit commutators, and imposing connectivity on top of them would drop terms that must survive; a caller must therefore pass empty connections explicitly, as the class does internally:

.. code-block:: cpp

   CC ucc{2, {.ansatz = CC::Ansatz::U, .hbar_comm_rank = 3}};

   // correct for a unitary ansatz: hbar is already self-contained
   auto ueqs = op::ref_av(op::P(nₚ(2)) * ucc.hbar(), {.connect = {}});

.. _cc-spin-tracing:

Spin Tracing of Expressions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Equations generated by the :class:`CC <sequant::mbpt::CC>` class are in spin-orbital basis. SeQuant also provides capability to transform these equations into spin-traced forms.

Make sure to include the ``<SeQuant/domain/mbpt/spin.hpp>`` header to access spin-tracing functions.

.. literalinclude:: /examples/user/cc.cpp
   :language: cpp
   :start-after: start-snippet-5
   :end-before: end-snippet-5
   :dedent: 2
