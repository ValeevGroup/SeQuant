Operator Interface for Symbolic Many-Body Algebra
=================================================

The ``mbpt::Operator`` class provides a convenient interface for symbolic many-body operators in SeQuant.
It is designed to represent, manipulate, and algebraically combine operators such as the Hamiltonian, excitation/deexcitation operators, and general tensor operators, all within a context-aware, quantum-number-tracking framework.

Context and Usage
-----------------
In SeQuant, the ``mbpt::Operator`` class is a symbolic object that encodes the quantum number transformation properties of many-body operators.
The ``mbpt::Operator`` class is used to construct fundamental building blocks of quantum chemistry Hamiltonians and other many-body operators.

For example, functions like ``H_()``, ``T_()``, ``Λ_()``, ``R_()``, and ``L_()``  in the ``mbpt::op`` namespace all return instances of ``Operator`` (or expressions built from them), each representing a specific type of many-body operator (Hamiltonian, excitation, deexcitation, etc.).

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

Illustrative Examples
---------------------
