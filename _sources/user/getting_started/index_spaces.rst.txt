Register Index Spaces
---------------------

Tensor expressions annotated by `abstract indices <https://en.wikipedia.org/wiki/Abstract_index_notation>`_ are common. In some contexts all tensor
modes refer to the same range or underlying vector space (as in all examples shown so far); then there is no need to distinguish modes of different
types. But in some contexts indices carry important semantic meaning. For example, the energy expression in the
`coupled-cluster method <https://doi.org/10.1017/CBO9780511596834>`_,

.. math::
   E_\mathrm{CC} =  F^{a_1}_ {i_1} t^{i_1}_ {a_1} + \frac{1}{4} \bar{g}^{a_1 a_2}_ {i_1 i_2} (t^{i_1 i_2}_ {a_1 a_2} + 2 t^{i_1}_ {a_1} t^{i_2}_ {a_2})

contains tensors with 2 types of modes, denoted by :math:`i` and :math:`a`, that represent single-particle (SP) states occupied and unoccupied in the
reference state, respectively. To simplify symbolic manipulation of such expressions SeQuant allows to define a custom vocabulary of index spaces and
to define their set-theoretic relationships. The following example illustrates the full space denoted by :math:`p` partitioned into occupied :math:`i`
and unoccupied :math:`a` base subspaces:

.. literalinclude:: /examples/user/getting_started/index_space_registry.cpp
   :language: cpp
   :start-after: start-snippet-1
   :end-before: end-snippet-1
   :dedent: 2

This and other vocabularies commonly used in quantum many-body context are supported out-of-the-box by SeQuant; their definitions are in :code:`SeQuant/domain/mbpt/convention.hpp`. The previous example is equivalent to the following:

.. literalinclude:: /examples/user/getting_started/index_space_registry.cpp
   :language: cpp
   :start-after: start-snippet-2
   :end-before: end-snippet-2
   :dedent: 2

Bitset representation of index spaces allows to define set-theoretic operations naturally. Bitset-based representation is used not only for index space *type* attribute (:code:`IndexSpace::Type`) but also for the *quantum numbers* attribute (:code:`IndexSpace::QuantumNumbers`). The latter can be used to represent spin quantum numbers, particle types, etc.
The main difference of the last example with the original example is that the :code:`make_min_sr_spaces()` factory changes the quantum numbers used by default (:code:`mbpt::Spin::any`) to make spin algebraic manipulations (like tracing out spin degrees of freedom) easier. Users can create their own definitions to suit their needs, but the vast majority of users will not need to venture outside of the predefined vocabularies.

Notice that the set-theoretic operations are only partially automated. It is the user's responsibility to define any and all unions and intersections of base spaces that they may encounter in their context. For this reason :class:`sequant::IndexSpaceRegistry` has its own :code:`unIon()` and :code:`intersection()` methods that perform error checking to ensure that only registered spaces are defined.

Quasiparticles
~~~~~~~~~~~~~~~

In most cases we are interested in using SeQuant to manipulate expressions involving operators in normal order relative to a vacuum state with a finite number of particles, rather than with respect to the genuine vacuum with zero particles. The choice of vacuum state as well as other related traits (whether the SP states are orthonormal, etc.) is defined by the implicit global context. The SeQuant programs until now used the genuine vacuum. The active context can be examined by calling :code:`get_default_context()`, changed via :code:`set_default_context()`, and reset to the default via :code:`reset_default_context()`:

.. literalinclude:: /examples/user/getting_started/index_spaces.cpp
   :language: cpp
   :start-after: start-snippet-1
   :end-before: end-snippet-1
   :dedent: 2

However, to deal with the single-product vacuum it is necessary to register at least one space and announce it as occupied in the vacuum state:

.. literalinclude:: /examples/user/getting_started/index_spaces.cpp
   :language: cpp
   :start-after: start-snippet-2
   :end-before: end-snippet-2
   :dedent: 2

It is also necessary to specify the *complete* space (union of all base spaces) so that the the space of unoccupied SP states can be determined:

.. literalinclude:: /examples/user/getting_started/index_spaces.cpp
   :language: cpp
   :start-after: start-snippet-3
   :end-before: end-snippet-3
   :dedent: 2

The Wick's theorem code itself is independent of the choice of vacuum:

.. literalinclude:: /examples/user/getting_started/index_spaces_wick.cpp
   :language: cpp

produces

.. math::
    {{\tilde{a}_ {{p_3}}}{\tilde{a}^{{p_1}}}{\tilde{a}_ {{p_4}}}{\tilde{a}^{{p_2}}}} = \bigl({\tilde{a}^{{p_1}{p_2}}_ {{p_3}{p_4}}} - {{s^{{p_1}}_ {{z_1}}}{s^{{z_1}}_ {{p_3}}}{\tilde{a}^{{p_2}}_ {{p_4}}}} - {{s^{{p_2}}_ {{z_1}}}{s^{{z_1}}_ {{p_4}}}{\tilde{a}^{{p_1}}_ {{p_3}}}} - {{s^{{p_1}}_ {{y_1}}}{s^{{y_1}}_ {{p_4}}}{\tilde{a}^{{p_2}}_ {{p_3}}}}
    + {{s^{{p_2}}_ {{z_1}}}{s^{{z_1}}_ {{p_3}}}{\tilde{a}^{{p_1}}_ {{p_4}}}} + {{s^{{p_1}}_ {{z_1}}}{s^{{p_2}}_ {{z_2}}}{s^{{z_1}}_ {{p_3}}}{s^{{z_2}}_ {{p_4}}}} + {{s^{{p_1}}_ {{y_1}}}{s^{{p_2}}_ {{z_1}}}{s^{{z_1}}_ {{p_3}}}{s^{{y_1}}_ {{p_4}}}}\bigr)

Note that:

* the tilde in :math:`\tilde{a}` denotes normal order with respect to single-product vacuum
* Einstein summation convention is implied, i.e., indices that appear twice in a given product (once in superscript, once in subscript) are summed over.
