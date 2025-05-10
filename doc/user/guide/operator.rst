Operator Interface for Symbolic Many-Body Algebra
=================================================

The ``mbpt::Operator`` class provides a convenient interface for symbolic many-body operators in SeQuant.
It is designed to represent, manipulate, and algebraically combine operators such as the Hamiltonian, excitation/deexcitation operators, and general tensor operators, all within a context-aware, quantum-number-tracking framework.

Context and Usage
-----------------
In SeQuant, the ``mbpt::Operator`` class is a symbolic object that encodes the quantum number transformation properties of many-body operators.
The ``mbpt::Operator`` class is used to construct fundamental building blocks of quantum chemistry Hamiltonians and other many-body operators.

For example, functions like ``H_()``, ``T_()``, ``Λ_()``, ``R_()``, and ``L_()``  in the ``mbpt::op`` namespace all return instances of ``Operator`` (or expressions built from them), each representing a specific type of many-body operator (Hamiltonian, excitation, deexcitation, etc).

When you call ``T_(K)`` or ``Λ_(K)`` in the ``mbpt::op`` namespace, you are constructing a symbolic excitation or deexcitation operator of rank ``K``. The returned object is an ``Operator`` that knows:

- Its label (e.g., "T" or "Λ")
- Corresponding :class:`sequant::Tensor` object (indices, rank, symmetry)
- How it changes the quantum numbers of a state it acts upon

This allows you to write code that manipulates operators at a high level, while the system automatically handles all the quantum number bookkeeping and algebraic rules.
This also enables efficient screening of terms using the quantum number tracking (e.g., number of creators and annihilators in each :class:`sequant::IndexSpace`), and context-aware operator algebra.


``QuantumNumberChange``: The Key Part of Operator Interface
-----------------------------------------------------------
The ``QuantumNumberChange`` class is a utility for representing how a many-body operator changes the quantum numbers of a state.
For an ``Operator`` object, this typically means tracking the number of particle and hole creators/annihilators in each relevant ``IndexSpace``.

How does this work?
^^^^^^^^^^^^^^^^^^^
- Each ``QuantumNumberChange`` instance stores, for each tracked subspace, an interval (e.g., "between 0 and 2 particles created").
- The number and meaning of these intervals is determined by the current ``sequant::Context`` (e.g., physical vacuum vs. single-product vacuum, number of base spaces etc).
- This allows the system to represent not just definite quantum number changes, but also operators with variable action.

How do we use it for operator algebra?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- ``QuantumNumberChange`` is the type used by the ``Operator`` class to track quantum number changes.
- When two operators are multiplied, their quantum number changes are combined using the rules encoded in ``QuantumNumberChange``.
- This enables SeQuant to **screen out expressions that would yield zero**, before applying Wick's Theorem without much computational overhead.
- ``QuantumNumberChange`` is also used **check if two operators commute with each other**, which is crucial for reducing the expression to its compact canonical form.

.. seealso::
  There are convenient helper functions available in ``sequant::mbpt`` namespace for constructing different types of ``QuantumNumberChange`` objects.
  For example, see ``sequant::mbpt::excitation_type_qns``.

Examples
--------

Let's look at some examples using the following SeQuant context:

.. literalinclude:: /examples/user/operator.cpp
   :language: cpp
   :start-after: start-snippet-0
   :end-before: end-snippet-0
   :dedent: 2

Constructing an ``Operator``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: /examples/user/operator.cpp
   :language: cpp
   :start-after: start-snippet-1
   :end-before: end-snippet-1
   :dedent: 2

The components of the ``Operator`` can be accessed as follows:

.. literalinclude:: /examples/user/operator.cpp
   :language: cpp
   :start-after: start-snippet-2
   :end-before: end-snippet-2
   :dedent: 2

Expression construction and screening
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: /examples/user/operator.cpp
   :language: cpp
   :start-after: start-snippet-3
   :end-before: end-snippet-3
   :dedent: 2

Vacuum averaging and final expression
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``sequant::mbpt::op::vac_av`` function can be used to compute the vacuum average of an operator level expression.

.. literalinclude:: /examples/user/operator.cpp
   :language: cpp
   :start-after: start-snippet-4
   :end-before: end-snippet-4
   :dedent: 2

.. note:: ``op::P`` is a projection operator which can be used to construct an excited bra or ket manifold based on the input.
