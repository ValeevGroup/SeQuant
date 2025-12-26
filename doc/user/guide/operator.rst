Operator Interface for Symbolic Many-Body Algebra
=================================================

Motivation
----------

The :class:`mbpt::Operator <sequant::mbpt::Operator>` class represents an abstract Fock-space many-body operator, such as the Hamiltonian operator, that can be used within SeQuant expressions.
We already saw `earlier <../getting_started/operators.html>`_ that
the lower-level *tensor* representation of such operators involves specific index labels which makes tensor
representation  inconvenient when multiple instances of the same operator occur.
For example, consider a general 1-body operator

.. math::
    \hat{o}_1 \equiv o^{p_1}_{p_2} a_{p_2}^{p_1}.

The abstract representation :math:`\hat{o}_1` is convenient for efficient algebraic manipulations such as expansion and simplication:

.. math::
    \left(1 + \hat{o}_1\right)^2 \equiv = 1 + \hat{o}_1 + \hat{o}_1 + \hat{o}_1^2 = 1 + 2 \hat{o}_1 + \hat{o}_1^2

The tensor form is less convenient for such manipulations for 2 reasons:

- lowering multiple instances of :math:`\hat{o}_1` such as in its square :math:`\hat{o}_1^2 \equiv \hat{o}_1 \hat{o}_1` to the tensor form is context-sensitive since one must avoid index reuse:

  .. math::
      \hat{o}_1^2 = o^{p_1}_ {p_2} a_{p_2}^{p_1} o^{p_3}_{p_4} a_{p_4}^{p_3}
- operating on tensor forms of many-body operators involves frequent use of canonicalization, such as to be able to simplify :math:`\hat{o}_1 + \hat{o}_1` in the following tensor form

  .. math::
      \hat{o}_1 + \hat{o}_1 = o^{p_1}_ {p_2} a_{p_2}^{p_1} + o^{p_3}_{p_4} a_{p_4}^{p_3}

  to

  .. math::
      2 \hat{o}_1 = 2 o^{p_1}_ {p_2} a_{p_2}^{p_1}

  requires canonicalization of each term before attempting simplification. In the abstract form simplifying :math:`\hat{o}_1 + \hat{o}_1` to :math:`2 \hat{o}_1` is as simple as if :math:`\hat{o}_1` were a scalar.

To be useful for more than just algebraic simplification the abstract operator form
must include enough detail to determine what products of such operators do to quantum numbers (e.g., occupation numbers)
of a Fock-space state. Each ``mbpt::Operator`` thus denotes its effect on the quantum numbers. This allows efficient evaluation of
vacuum averages of products of ``mbpt::Operator`` objects by identifying
only the non-zero contributions.

``mbpt::Operator`` Anatomy
--------------------------

``mbpt::Operator`` is constructed from 3 callables (e.g., lambdas):

- ``void -> std::wstring_view``: returns the operator label such as "T" or "Λ"
- ``void -> ExptPtr``: returns the tensor form of the operator, and
- ``QuantumNumberChange<>& -> void``: encodes the action on the given input state by mutating its quantum numbers given as the argument to the lambda.

The quantum numbers of a state and their changes are represented by the same class, :class:`mbpt::QuantumNumberChange <sequant::mbpt::QuantumNumberChange>`.
For the Fock-space states the quantum numbers are the number of creators/annihilators. In general context the number of creators/annihilators
must be tracked separately for multiple ``IndexSpaces`` (e.g., for each hole and particle states in the particle-hole context). To support
general operators which can create a variable number of particles and holes the numbers of creators/annihilators must be encoded
by *intervals*. For example, if space :math:`\{p\}` is a union of hole (fully occupied in the vacuum) space :math:`\{i\}`
and particle (fully empty in the vacuum) space :math:`\{a\}` operator :math:`\hat{o}_1` includes 4 components,

.. math::
    o^{p_1}_{p_2} a_{p_2}^{p_1} = o^{i_1}_{i_2} a_{i_2}^{i_1} + o^{i_1}_{a_2} a_{a_2}^{i_1} + o^{a_1}_{i_2} a_{i_2}^{a_1} + o^{a_1}_{a_2} a_{a_2}^{a_1},

each of which is has 1 (particle or hole) creator and 1 (particle or hole) annihilator.
The net effect of this operator on a Fock-space state is to leave its number of hole/particle
creators/annihilators unchanged or increased by 1. E.g., its effect on a state with :math:`n`
particle creators is produce a state with :math:`[n,n+1]` particle creators; similarly,
action on a state with :math:`[n,m]` hole annihilators (:math:`n \leq m`) will produce a state with :math:`[n,m+1]` hole annihilators.

Finally, for a product of ``mbpt::Operator`` it is possible to determine its effect on the quantum numbers of the input state using simple rules based on Wick's theorem.
This enables SeQuant to **screen out expressions that would yield zero**, before applying Wick's theorem without much computational overhead.
This logic can also be used to **check if two operators commute with each other**, which is crucial for reducing the expression to their compact canonical form.

The use of interval arithmetic for quantum number changes allows to treat operators in their natural compact form as long as possible. E.g., :math:`\hat{o}_1` can be kept as a single term :math:`o^{p_1}_{p_2} a_{p_2}^{p_1}` rather than having to decompose
it into individual contributions, each with definite effect on quantum numbers. This is especially useful when the number of index spaces increases.

.. seealso::
    - Although the user can construct ``mbpt::Operator`` directly, SeQuant predefines factories for many commonly-used operators. For example, functions :class:`H_ <sequant::mbpt::op::H_>`, :class:`T_ <sequant::mbpt::op::T_>`, :class:`Λ_ <sequant::mbpt::op::Λ_>`, :class:`R_ <sequant::mbpt::op::R_>`, and :class:`L_ <sequant::mbpt::op::L_>`  in the ``mbpt::op`` namespace all return instances of ``Operator`` (or expressions built from them), each representing a specific type of many-body operator (Hamiltonian, excitation, deexcitation, etc).
    - There are convenient helper functions available in ``sequant::mbpt`` namespace for constructing different types of ``QuantumNumberChange`` objects. For example, see ``sequant::mbpt::excitation_type_qns``.

Operator Registry
-----------------

The :class:`mbpt::OpRegistry <sequant::mbpt::OpRegistry>` maintains a mapping from operator labels to their operator classes. In this section we describe how to use and extend the operator registry.

Each operator must be registered with one of three classes:

- :cpp:enumerator:`OpClass::ex <sequant::mbpt::OpClass::ex>`: Excitation operators (e.g., cluster excitation operator)
- :cpp:enumerator:`OpClass::deex <sequant::mbpt::OpClass::deex>`: De-excitation operators (e.g., cluster de-excitation operator)
- :cpp:enumerator:`OpClass::gen <sequant::mbpt::OpClass::gen>`: General operators (e.g., Hamiltonian-like operators)

.. note:: Some operator labels are reserved for SeQuant's internal MBPT and core functionalities and cannot be registered. See ``SeQuant/domain/mbpt/reserved.hpp`` for more details.

SeQuant provides two predefined registries:

- :func:`make_minimal_registry() <sequant::mbpt::make_minimal_registry>`: Contains only a minimal set of operators
- :func:`make_legacy_registry() <sequant::mbpt::make_legacy_registry>`: A backward-compatible registry that includes all operators used in legacy MBPT module.


Custom operators can be registered to extend SeQuant's vocabulary. Once registered, operators can be used with both :class:`OpMaker <sequant::mbpt::OpMaker>` and :class:`mbpt::Operator <sequant::mbpt::Operator>`:

.. literalinclude:: /examples/user/operator.cpp
   :language: cpp
   :start-after: start-snippet-5
   :end-before: end-snippet-5
   :dedent: 2

The registry is stored in the :class:`mbpt::Context <sequant::mbpt::Context>`. See the unit test cases for further manipulations with :class:`mbpt::Context <sequant::mbpt::Context>` and :class:`mbpt::OpRegistry <sequant::mbpt::OpRegistry>`.

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
If reference state differs from the Wick vacuum ``sequant::mbpt::op::ref_av`` function should be used instead to
compute the reference average.

.. literalinclude:: /examples/user/operator.cpp
   :language: cpp
   :start-after: start-snippet-4
   :end-before: end-snippet-4
   :dedent: 2

.. note:: ``op::P`` is a projection operator which can be used to construct an excited bra or ket manifold based on the input.
