Coupled-Cluster Class
======================

Coupled-cluster (CC) theory is one of the most accurate and widely used quantum chemistry methods for describing `electron correlation
<https://en.wikipedia.org/wiki/Electron_correlation>`_ in molecular systems. It represents the many-electron wavefunction using an exponential
`ansatz <https://en.wikipedia.org/wiki/Ansatz>`_:

.. math::

   |\Psi_{\text{CC}}\rangle = e^{\hat{T}}|\Phi_0\rangle

where :math:`|\Phi_0\rangle` is a reference determinant (typically `Hartree-Fock <https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method>`_), and :math:`\hat{T}` is a cluster operator that generates excited determinants. The cluster operator is typically expanded as:

.. math::

   \hat{T} = \hat{T}_1 + \hat{T}_2 + \hat{T}_3 + \ldots

where :math:`\hat{T}_n` generates :math:`n`-fold excited determinants. For computational tractability, the cluster operator is usually truncated. For example, CCSD includes only single and double excitations (:math:`\hat{T} = \hat{T}_1 + \hat{T}_2`).

The :class:`CC <sequant::mbpt::CC>` class provides a convenient interface for setting up and processing CC equations using SeQuant’s symbolic algebra engine. It supports various CC formulations, including traditional, unitary, and orbital-optimized ansätze.

Overview
--------

The :class:`CC <sequant::mbpt::CC>` class can be used to derive:

- Ground state amplitude equations
- λ (de-excitation) amplitude equations — Lagrange multipliers conjugate to the ground-state amplitudes, needed for properties, analytic gradients and
  perturbative corrections
- Equation-of-motion (EOM) CC equations for excited states
- Response equations for properties and `perturbations <https://en.wikipedia.org/wiki/Perturbation_theory_(quantum_mechanics)>`_

Expressions are generated in spin-orbital basis and can be post-processed using SeQuant's spin-tracing capabilities. See :ref:`cc-spin-tracing` for more details.

Multireference contexts are not fully supported yet: only ``hbar()``, ``energy()`` and ``t()`` are available, and only
with the BCH expansion.


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
       const std::vector<std::size_t>& block_ranks = {});
   std::vector<ExprPtr> eom_l(nₚ np, nₕ nh);

Derives equation-of-motion coupled-cluster (EOM-CC) equations for excited states. The ``eom_r`` method generates equations for the right eigenvectors, while ``eom_l`` generates equations for the left eigenvectors. Traditional CC always uses the connected H̄R product. UCC always uses Hamiltonian-matrix assembly, subtracting the scalar part of each diagonal block's Hamiltonian.

For UCC, the optional ``block_ranks`` argument gives the per-block truncation orders as a row-major matrix over the EOM manifolds: nested-commutator order for BCH and :math:`\bar{H}^{k}` order for Bernoulli. When omitted, every block uses the configured H̄ rank. A positive ``hbar_singles_comm_rank`` applies the additional singles transform even to rank-0 BCH blocks. Traditional CC does not support block ranks.

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

For unitary BCH expansions, ``Options::hbar_singles_comm_rank`` applies an additional singles-only similarity transform to H̄ generated by :math:`\sigma_1 = T_1 - T_1^\dagger` after the primary expansion. A zero rank disables this transform; a positive rank requires a unitary ansatz with singles amplitudes enabled and is not supported with Bernoulli expansions. Perturbed amplitude equations and RDMs are not supported with this option.

.. _cc-hbar-connectivity:

Using :math:`\bar{H}` outside the CC class
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:func:`CC::hbar() <sequant::mbpt::CC::hbar>` is public, but the form of the expression it returns depends on the ansatz
and on the reference. For a non-unitary ansatz with the reference equal to the Wick vacuum, each commutator is written as
a connected product (see :ref:`mbpt-lst`), which equals the commutator only once the operators are connected when taking
the expectation value; evaluated with empty connectivity it retains disconnected terms. For a unitary ansatz, or when the
reference differs from the Wick vacuum (where ``ref_av`` requires empty ``connect`` and ``do_not_connect`` lists), H̄ is
built from explicit commutators and is self-contained; imposing connectivity on it would drop terms that must survive.

:func:`CC::hbar_connections() <sequant::mbpt::CC::hbar_connections>` returns the connectivity that matches the form H̄
was built with, ``default_op_connections()`` or empty, and the class uses it (or a superset) for its own equations. Pass
it whenever you evaluate ``CC::hbar()`` yourself:

.. literalinclude:: /examples/user/cc.cpp
   :language: cpp
   :start-after: start-snippet-6
   :end-before: end-snippet-6
   :dedent: 2

The public operator-level and tensor-level ``ref_av`` and ``vac_av`` functions default to empty connectivity constraints;
omitting the options argument is equivalent to passing ``{}``. Alternatively, build H̄ with explicit commutators using
``mbpt::lst`` with its default options, which needs no connectivity.

.. _cc-spin-tracing:

Spin Tracing of Expressions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Equations generated by the :class:`CC <sequant::mbpt::CC>` class are in spin-orbital basis. SeQuant also provides capability to transform these equations into spin-traced forms; see :doc:`spin_tracing` for the underlying concepts (the ``Spin`` quantum number and the general ``spintrace()`` function).

Make sure to include the ``<SeQuant/domain/mbpt/spin.hpp>`` header to access spin-tracing functions.

.. literalinclude:: /examples/user/cc.cpp
   :language: cpp
   :start-after: start-snippet-5
   :end-before: end-snippet-5
   :dedent: 2
