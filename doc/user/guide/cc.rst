Coupled-Cluster Class
======================

Coupled-cluster (CC) theory is one of the most accurate and widely used quantum chemistry methods for describing electron correlation in molecular systems. It represents the many-electron wavefunction using an exponential ansatz:

.. math::

   |\Psi_{\text{CC}}\rangle = e^{\hat{T}}|\Phi_0\rangle

where :math:`|\Phi_0\rangle` is a reference determinant (typically Hartree-Fock), and :math:`\hat{T}` is a cluster operator that generates excited determinants. The cluster operator is typically expanded as:

.. math::

   \hat{T} = \hat{T}_1 + \hat{T}_2 + \hat{T}_3 + \ldots

where :math:`\hat{T}_n` generates :math:`n`-fold excited determinants. For computational tractability, the cluster operator is usually truncated. For example, CCSD includes only single and double excitations (:math:`\hat{T} = \hat{T}_1 + \hat{T}_2`).

The SeQuant ``CC`` class (see :class:`sequant::mbpt::CC`) provides a powerful symbolic algebra engine for deriving the complex equations that arise in CC theory for both ground and excited states. It supports various formulations of the CC method, including traditional, unitary, and orbital-optimized ansätze.

Overview
--------

The ``sequant::mbpt::CC`` class is a derivation engine for generating CC equations. It can derive:

- Ground state amplitude equations
- λ (de-excitation) amplitude equations
- Equation-of-motion (EOM) CC equations for excited states
- Response equations for properties and perturbations

Expressions are generated in spin-orbital basis and can be post-processed using SeQuant's spin-tracing capabilities. See :ref:`cc-spin-tracing` for more details.


Ansatz Options
--------------

The ``CC`` class supports several CC ansätze through the ``CC::Ansatz`` enum (see :enum:`sequant::mbpt::CC::Ansatz`):

- ``Ansatz::T``: Traditional CC ansatz, where the wavefunction is represented as :math:`|\Psi_{\text{CC}}\rangle = e^{\hat{T}}|\Phi_0\rangle`. This is the standard approach used in most implementations.

- ``Ansatz::oT``: Orbital-optimized traditional ansatz. Singles amplitudes (:math:`\hat{T}_1`) are excluded from the cluster operator, with orbital optimization performed instead.

- ``Ansatz::U``: Unitary CC ansatz, where the wavefunction is represented as :math:`|\Psi_{\text{UCC}}\rangle = e^{\hat{T} - \hat{T}^\dagger}|\Phi_0\rangle`. Particularly useful for quantum computing applications.

- ``Ansatz::oU``: Orbital-optimized unitary ansatz, combines both unitary and orbital optimized ansätze.

Key Methods
----------

Similarity Transformation
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   ExprPtr sim_tr(ExprPtr expr, size_t r);

Performs the similarity transformation of an operator. For traditional CC, this computes :math:`\bar{O} = e^{-\hat{T}} \hat{O} e^{\hat{T}}` expanded as a series of nested commutators truncated at order ``r``, where :math:`\hat{O}` is the operator to be transformed.

Ground State Amplitudes
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   std::vector<ExprPtr> t(size_t commutator_rank = 4,
                          size_t pmax = std::numeric_limits<size_t>::max(),
                          size_t pmin = 0);

Derives the equations for the :math:`t` amplitudes (:math:`\langle \Phi_P|\bar{H}|0 \rangle = 0`) up to specified excitation levels.

.. code-block:: cpp

   std::vector<ExprPtr> λ(size_t commutator_rank = 4);

Derives the equations for the :math:`\lambda` de-excitation amplitudes.

Coupled-Cluster Response
^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   std::vector<ExprPtr> t_pt(size_t order = 1, size_t rank = 1);
   std::vector<ExprPtr> λ_pt(size_t order = 1, size_t rank = 1);

Derives perturbed amplitude equations for response theory calculations.

Equation-of-Motion Coupled-Cluster
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   std::vector<ExprPtr> eom_r(nₚ np, nₕ nh);
   std::vector<ExprPtr> eom_l(nₚ np, nₕ nh);

Derives equation-of-motion coupled-cluster (EOM-CC) equations for excited states. The ``eom_r`` method generates equations for the right eigenvectors, while ``eom_l`` generates equations for the left eigenvectors.

Examples
--------

The following examples demonstrate how to use the ``CC`` class to derive CC equations for various ansätze and excitation levels.

From this point onward, assume the following namespaces are imported, and the ``sequant::Context`` is configured as shown.

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
-------------

Truncating the Commutator Expansion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For traditional CC with two-body Hamiltonians, the commutator expansion is truncated at the a 4th-order. However, for unitary CC or other Hamiltonians, you may need to explicitly set the commutator rank:

.. literalinclude:: /examples/user/cc.cpp
   :language: cpp
   :start-after: start-snippet-4
   :end-before: end-snippet-4
   :dedent: 2

.. _cc-spin-tracing:

Spin Tracing of Expressions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Equations generated by the ``CC`` class are in spin-orbital basis. SeQuant also provides capability to transform these equations into spin-traced forms.

Make sure to include the ``<SeQuant/domain/mbpt/spin.hpp>`` header to access spin-tracing functions.

.. literalinclude:: /examples/user/cc.cpp
   :language: cpp
   :start-after: start-snippet-5
   :end-before: end-snippet-5
   :dedent: 2
