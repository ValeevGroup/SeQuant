Operators
---------

Development of SeQuant is primarily motivated by the perturbative many-body methods, collectively referred to here as Many-Body Perturbation Theory (MBPT). Examples of such methods include the `coupled-cluster (CC) method <https://doi.org/10.1017/CBO9780511596834>`_ and `GW <https://doi.org/10.1103/PhysRev.139.A796>`_. The typical use case is to compute canonical forms of products of operators. For example, consider the coupled-cluster doubles (CCD) method.
*Amplitudes* :math:`t^{i_1 i_2}_ {a_1 a_2}` of the cluster operator,

.. math::
    \hat{t} \equiv \hat{t}_ 2 = \frac{1}{4} t^{i_1 i_2}_ {a_1 a_2} a_ {i_1 i_2}^{a_1 a_2},

are determined by solving the CCD equations:

.. math::
    \forall i_1, i_2, a_1, a_2: \quad 0 = \langle0\vert a^{i_1 i_2}_ {a_1 a_2} \exp(-\hat{t}_ 2) \hat{H} \exp(\hat{t}_ 2) \vert 0 \rangle = \langle0\vert a^{i_1 i_2}_ {a_1 a_2} \bigl( \hat{H} + [\hat{H}, \hat{t}_ 2] + \frac{1}{2} [[\hat{H}, \hat{t}_ 2], \hat{t}_ 2] \bigr) \vert 0 \rangle.

A pedestrian way to compose such expression is to define a cluster operator object using SeQuant tensors and normal-ordered operators:

.. literalinclude:: /examples/user/getting_started/ccd_pedestrian.cpp
   :language: cpp
   :start-after: start-snippet-1
   :end-before: end-snippet-1
   :dedent: 2

The normal-ordered Hamiltonian is defined similarly as a sum of 1- and 2-body contributions:

.. literalinclude:: /examples/user/getting_started/ccd_pedestrian.cpp
   :language: cpp
   :start-after: start-snippet-2
   :end-before: end-snippet-2
   :dedent: 2

Note that the compact definition of the Hamiltonian is due to the use of the union (:math:`p`) of base occupied (:math:`i`) and unoccupied (:math:`a`) spaces. Many other symbolic algebras only support use of nonoverlapping base spaces, in terms of which Hamiltonian and other tensor expressions would have a much more verbose form.
The commutator of the Hamiltonian and cluster operators is trivially composed:

.. literalinclude:: /examples/user/getting_started/ccd_pedestrian.cpp
   :language: cpp
   :start-after: start-snippet-3
   :end-before: end-snippet-3
   :dedent: 2

Note the use of :code:`simplify` to rewrite an expression in a simpler form. Its role will be emphasized later.
Unfortunately, we immediately run into the limitation of the "pedestrian" approach. Namely, the double commutator cannot be correctly obtained as

.. literalinclude:: /examples/user/getting_started/ccd_pedestrian.cpp
   :language: cpp
   :start-after: start-snippet-4
   :end-before: end-snippet-4
   :dedent: 2

due to the explicit use of specific dummy indices in the definition of :code:`t2`. Using it more than once in a given product will produce an expresson where each dummy index appears more than 2 times, breaking the Einstein summation convention.

The issue is actually not the reuse of the same of operator object, but the reuse of dummy indices. A straightforward, but brittle, solution is to ensure that each dummy index is only used once. E.g.,
to use :math:`\hat{t}_2` more than once in an expression we must make several versions of it, each with a separate set of dummy indices:

.. literalinclude:: /examples/user/getting_started/ccd_pedestrian.cpp
   :language: cpp
   :start-after: start-snippet-5
   :end-before: end-snippet-5
   :dedent: 2

This is too error-prone for making complex expressions. A better way is to represent :math:`\hat{t}_2` by an object that generates tensor form with unique dummy indices generated on the fly. in SeQuant such MBPT *operators* live in :code:`mbpt` namespace. The entire CCD amplitude equation is evaluated as follows:

.. literalinclude:: /examples/user/getting_started/ccd.cpp
   :language: cpp

The result is

.. math::
    \bigl({{{\frac{1}{4}}}{A^{{a_1}{a_2}}_ {{i_1}{i_2}}}{\bar{g}^{{i_1}{i_2}}_ {{a_1}{a_2}}}} - {{{\frac{1}{2}}}{A^{{a_1}{a_2}}_ {{i_1}{i_2}}}{f^{{a_3}}_ {{a_1}}}{\bar{t}^{{i_1}{i_2}}_ {{a_2}{a_3}}}} - {{A^{{a_1}{a_2}}_ {{i_1}{i_2}}}{\bar{g}^{{i_1}{a_3}}_ {{i_3}{a_1}}}{\bar{t}^{{i_2}{i_3}}_ {{a_2}{a_3}}}} + {{{\frac{1}{8}}}{A^{{a_1}{a_2}}_ {{i_1}{i_2}}}{\bar{g}^{{a_3}{a_4}}_ {{a_1}{a_2}}}{\bar{t}^{{i_1}{i_2}}_ {{a_3}{a_4}}}}
    + {{{\frac{1}{8}}}{A^{{a_1}{a_2}}_ {{i_1}{i_2}}}{\bar{g}^{{i_1}{i_2}}_ {{i_3}{i_4}}}{\bar{t}^{{i_3}{i_4}}_ {{a_1}{a_2}}}} + {{{\frac{1}{2}}}{A^{{a_1}{a_2}}_ {{i_1}{i_2}}}{f^{{i_1}}_ {{i_3}}}{\bar{t}^{{i_2}{i_3}}_ {{a_1}{a_2}}}} + {{{\frac{1}{16}}}{A^{{a_1}{a_2}}_ {{i_1}{i_2}}}{\bar{g}^{{a_3}{a_4}}_ {{i_3}{i_4}}}{\bar{t}^{{i_1}{i_2}}_ {{a_3}{a_4}}}{\bar{t}^{{i_3}{i_4}}_ {{a_1}{a_2}}}}
    - {{{\frac{1}{4}}}{A^{{a_1}{a_2}}_ {{i_1}{i_2}}}{\bar{g}^{{a_3}{a_4}}_ {{i_3}{i_4}}}{\bar{t}^{{i_1}{i_3}}_ {{a_3}{a_4}}}{\bar{t}^{{i_2}{i_4}}_ {{a_1}{a_2}}}} - {{{\frac{1}{4}}}{A^{{a_1}{a_2}}_ {{i_1}{i_2}}}{\bar{g}^{{a_3}{a_4}}_ {{i_3}{i_4}}}{\bar{t}^{{i_1}{i_2}}_ {{a_1}{a_3}}}{\bar{t}^{{i_3}{i_4}}_ {{a_2}{a_4}}}}
    - {{{\frac{1}{4}}}{A^{{a_1}{a_2}}_ {{i_1}{i_2}}}{\bar{g}^{{a_3}{a_4}}_ {{i_3}{i_4}}}{\bar{t}^{{i_1}{i_3}}_ {{a_3}{a_4}}}{\bar{t}^{{i_2}{i_4}}_ {{a_1}{a_2}}}} - {{{\frac{1}{4}}}{A^{{a_1}{a_2}}_ {{i_1}{i_2}}}{\bar{g}^{{a_3}{a_4}}_ {{i_3}{i_4}}}{\bar{t}^{{i_1}{i_2}}_ {{a_1}{a_3}}}{\bar{t}^{{i_3}{i_4}}_ {{a_2}{a_4}}}}
    + {{{\frac{1}{2}}}{A^{{a_1}{a_2}}_ {{i_1}{i_2}}}{\bar{g}^{{a_3}{a_4}}_ {{i_3}{i_4}}}{\bar{t}^{{i_1}{i_3}}_ {{a_1}{a_3}}}{\bar{t}^{{i_2}{i_4}}_ {{a_2}{a_4}}}}\bigr)

The use of MBPT operators rather than their tensor-level forms not only solves problems with the reuse of dummy indices, but also allows to implement additional optimizations such as algebraic simplifications of complex operator expressions and avoiding evaluation of operator products whose vacuum expectation values are guaranteed to vanish. This allows very efficient derivation of complex equations, e.g. CC equations through CCSDTQ are derived in a fraction of a second on a laptop:

.. code-block:: sh

    $ cmake --build build --target srcc
    $ time build/srcc  time ./srcc 4 t std so
    CC equations [type=t,rank=1,spinfree=false,screen=true,use_topology=true,use_connectivity=true,canonical_only=true] computed in 0.0027805 seconds
    R1(expS1) has 8 terms:
    CC equations [type=t,rank=2,spinfree=false,screen=true,use_topology=true,use_connectivity=true,canonical_only=true] computed in 0.012890749999999999 seconds
    R1(expS2) has 14 terms:
    R2(expS2) has 31 terms:
    CC equations [type=t,rank=3,spinfree=false,screen=true,use_topology=true,use_connectivity=true,canonical_only=true] computed in 0.039590500000000001 seconds
    R1(expS3) has 15 terms:
    R2(expS3) has 37 terms:
    R3(expS3) has 47 terms:
    CC equations [type=t,rank=4,spinfree=false,screen=true,use_topology=true,use_connectivity=true,canonical_only=true] computed in 0.107501417 seconds
    R1(expS4) has 15 terms:
    R2(expS4) has 38 terms:
    R3(expS4) has 53 terms:
    R4(expS4) has 74 terms:
    ./srcc 4 t std so  0.27s user 0.41s system 100% cpu 0.674 total
