Using SeQuant
=============

SeQuant is a general-purpose symbolic tensor algebra, but the primary use case is in quantum many-body physics. The following is a brief tutorial on
using SeQuant for this purpose.


Wick's Theorem
---------------

To get started let's use SeQuant to apply `Wick's theorem <https://en.wikipedia.org/wiki/Wick's_theorem>`_ to a simple product of elementary (creation
and annihilation) fermionic operators:

.. math::
   a_{p_3} a_{p_4} a^\dagger_{p_1} a^\dagger_{p_2}.

This is achieved by the following SeQuant program:

.. literalinclude:: /examples/user/getting_started/wick.cpp
   :language: cpp
   :start-after: start-snippet-1
   :end-before: end-snippet-1
   :dedent: 2

.. hint:: 
   All *core* user-facing SeQuant code lives in C++ namespace :code:`sequant`, and the shown code assumes this namespace has been imported via
   :code:`using namespace sequant`.

Running this program should produce a LaTeX expression for this formula:

.. math::
   {{a_ {{p_3}}}{a_ {{p_4}}}{a^{{p_1}}}{a^{{p_2}}}}
   = \bigl( - {{a^{{p_1}{p_2}}_ {{p_3}{p_4}}}} - {{s^{{p_1}}_ {{p_3}}}{s^{{p_2}}_ {{p_4}}}} + {{s^{{p_1}}_ {{p_3}}}{a^{{p_2}}_ {{p_4}}}}
     - {{s^{{p_1}}_ {{p_4}}}{a^{{p_2}}_ {{p_3}}}} + {{s^{{p_2}}_ {{p_3}}}{s^{{p_1}}_ {{p_4}}}} - {{s^{{p_2}}_ {{p_3}}}{a^{{p_1}}_ {{p_4}}}}
     + {{s^{{p_2}}_ {{p_4}}}{a^{{p_1}}_ {{p_3}}}}\bigr)

where the tensor notation is used to denote elementary and composite *normal-ordered* (or, shortly, *normal*) operators:
 
.. math::
   a^p \equiv a_p^\dagger

.. math::
   a^{p_1 p_2 \dots p_c}_ {q_1 q_2 \dots q_a} \equiv a_ {p_1}^\dagger  a_ {p_2}^\dagger \dots a_ {p_c}^\dagger a_ {q_a} \dots a_ {q_2} a_ {q_1}


:math:`s^p_q \equiv \langle q | p \rangle` denotes inner products ("overlaps") of 1-particle states. Wick's theorem can of course be applied directly
to products of normal composite operators, e.g,

.. literalinclude:: /examples/user/getting_started/wick.cpp
   :language: cpp
   :start-after: start-snippet-2
   :end-before: end-snippet-2
   :dedent: 2

produces

.. math::
   {{a^{{p_1}{p_2}}_ {{p_3}{p_4}}}{a^{␣\,{p_5}}_ {{p_6}{p_7}}}}
   = \bigl({a^{␣\,{p_1}{p_2}{p_5}}_ {{p_3}{p_4}{p_6}{p_7}}} - {{s^{{p_5}}_ {{p_4}}}{a^{␣\,{p_1}{p_2}}_ {{p_3}{p_6}{p_7}}}}
     + {{s^{{p_5}}_ {{p_3}}}{a^{␣\,{p_1}{p_2}}_ {{p_4}{p_6}{p_7}}}}\bigr)

where :math:`␣` is used in number-nonconserving operators to point out the empty "slots".

Same algebra can be performed for bosons:

.. literalinclude:: /examples/user/getting_started/wick.cpp
   :language: cpp
   :start-after: start-snippet-3
   :end-before: end-snippet-3
   :dedent: 2

.. math::
   {{b^{{p_1}{p_2}}_ {{p_3}{p_4}}}{b^{{p_5}{p_6}}_ {␣\,{p_7}}}}
   = \bigl( {b^{{p_5}{p_1}{p_2}{p_6}}_ {␣\,{p_3}{p_4}{p_7}}} + {{s^{{p_6}}_ {{p_3}}}{b^{{p_5}{p_2}{p_1}}_{␣\,{p_4}{p_7}}}}
     + {{s^{{p_5}}_{{p_3}}}{b^{{p_1}{p_2}{p_6}}_{␣\,{p_4}{p_7}}}}
     + {{s^{{p_6}}_{{p_4}}}{b^{{p_5}{p_1}{p_2}}_{␣\,{p_3}{p_7}}}} + {{s^{{p_5}}_{{p_4}}}{b^{{p_2}{p_1}{p_6}}_{␣\,{p_3}{p_7}}}}
     + {{s^{{p_5}}_{{p_3}}}{s^{{p_6}}_{{p_4}}}{b^{{p_1}{p_2}}_{␣\,{p_7}}}}
     + {{s^{{p_6}}_{{p_3}}}{s^{{p_5}}_{{p_4}}}{b^{{p_2}{p_1}}_{␣\,{p_7}}}} \bigr)

where :math:`b` denotes normal bosonic operators constructed analogously with the normal fermionic operators :math:`a`


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

..
    TODO: port to reST

    ```c++
    #include <SeQuant/core/context.hpp>
    #include <SeQuant/core/index.hpp>

    int main() {
      using namespace sequant;
      IndexSpaceRegistry isr;

      // base spaces
      isr.add(L"i", 0b01)
         .add(L"a", 0b10);
      // union of 2 base spaces
      // can create manually, as isr.add(L"p", 0b11) , or explicitly ...
      isr.add_union(L"p", {L"i", L"a"}); // union of i and a

      // can access unions and intersections of base and composite spaces
      assert(isr.unIon(L"i", L"a") == isr.retrieve(L"p"));
      assert(isr.intersection(L"p", L"i") == isr.retrieve(L"i"));
      
      // to use the vocabulary defined by isr use it to make a Context object and make it the default
      set_default_context(Context(std::move(isr)));
      
      // now can use space labels to construct Index objects representing said spaces
      Index i1(L"i_1");
      Index a1(L"a_1");
      Index p1(L"p_1");
      
      // set theoretic operations on spaces
      assert(i1.space().intersection(a1.space()) == IndexSpace::null);
    }
    ```

    This and other vocabularies commonly used in quantum many-body context are supported out-of-the-box by SeQuant; their definition is in `SeQuant/domain/mbpt/convention.hpp`. The previous example is equivalent to the following:

    ```c++
    #include <SeQuant/core/context.hpp>
    #include <SeQuant/domain/mbpt/convention.hpp>

    int main() {
      using namespace sequant;
      using namespace sequant::mbpt;

      // makes 2 base spaces, i and a, and their union
      auto isr = make_min_sr_spaces();
      set_default_context(Context(isr));

      // set theoretic operations on spaces
      auto i1 = Index(L"i_1");
      auto a1 = Index(L"a_1");
      assert(i1.space().attr().intersection(a1.space().attr()).type() ==
             IndexSpace::Type::null);
      assert(i1.space().attr().intersection(a1.space().attr()).qns() ==
             Spin::any);
    }
    ```

    Bitset representation of index spaces allows to define set-theoretic operations naturally. Bitset-based representation is used not only for index space _type_ attribute (`IndexSpace::Type`) but also for the _quantum numbers_ attribute (`IndexSpace::QuantumNumbers`). The latter can be used to represent spin quantum numbers, particle types, etc.
    The main difference of the last example with the original example is that the `make_min_sr_spaces()` factory changes the quantum numbers used by default (`mbpt::Spin::any`) to make spin algebraic manipulations (like tracing out spin degrees of freedom) easier. Users can create their own definitions to suit their needs, but the vast majority of users will not need to venture outside of the predefined vacbularies.

    Notice that the set-theoretic operations are only partially automated. It is the user's responsibility to define any and all unions and intersections of base spaces that they may encounter in their context. For this reason `IndexSpaceRegistry` has its own `unIon()` and `intersection()` methods that perform error checking to ensure that only registered spaces are defined.

    ### Quasiparticles

    In most cases we are interested in using SeQuant to manipulate expressions involving operators in normal order relative to a vacuum state with a finite number of particles, rather than with respect to the genuine vacuum with zero particles. The choice of vacuum state as well as other related traits (whether the SP states are orthonormal, etc.) is defined by the implicit global context. The SeQuant programs until now used the genuine vacuum. The active context can be examined by calling `get_default_context()`, changed via `set_default_context()`, and reset to the default via `reset_default_context()`:

    ```c++
    #include <SeQuant/core/context.hpp>

    int main() {
      using namespace sequant;

      // the default is to use genuine vacuum
      assert(get_default_context().vacuum() == Vacuum::Physical);
      // now set the context to a single product of SP states
      set_default_context(Context{Vacuum::SingleProduct, IndexSpaceMetric::Unit, BraKetSymmetry::symm});
      assert(get_default_context().vacuum() == Vacuum::SingleProduct);
      // reset the context back to the default
      reset_default_context();
      assert(get_default_context().vacuum() == Vacuum::Physical);
      
      return 0;
    }
    ```

    However, to deal with the single-product vacuum it is necessary to register at least one space and announce it as occupied in the vacuum state:

    ```c++
    isr.add(L"y", 0b01).vacuum_occupied_space(L"i");
    ```

    or, shorter,

    ```c++
    isr.add(L"y", 0b01, is_vacuum_occupied);
    ```

    It is also necessary to specify the _complete_ space (union of all base spaces) so that the the space of unoccupied SP states can be determined:

    ```c++
    isr.add(L"y", 0b01, is_vacuum_occupied)
       .add(L"z", 0b10)
       .add(L"p", 0b11, is_complete);
    ```


    The Wick's theorem code itself is independent of the choice of vacuum:

    ```c++
    #include <SeQuant/core/context.hpp>
    #include <SeQuant/core/wick.hpp>

    int main() {
      using namespace sequant;

      set_default_context(Context{IndexSpaceRegistry{}
                                      .add(L"y", 0b01, is_vacuum_occupied)
                                      .add(L"z", 0b10)
                                      .add(L"p", 0b11, is_complete),
                                  Vacuum::SingleProduct});

      auto cp1 = fcrex(L"p_1"), cp2 = fcrex(L"p_2");
      auto ap3 = fannx(L"p_3"), ap4 = fannx(L"p_4");

      std::wcout << to_latex(ap3 * cp1 * ap4 * cp2) << " = "
                 << to_latex(FWickTheorem{ap3 * cp1 * ap4 * cp2}
                                 .full_contractions(false)
                                 .compute())
                 << std::endl;

      return 0;
    }
    ```

    produces

    $$
    {{\tilde{a}_ {{p_3}}}{\tilde{a}^{{p_1}}}{\tilde{a}_ {{p_4}}}{\tilde{a}^{{p_2}}}} = \bigl({\tilde{a}^{{p_1}{p_2}}_ {{p_3}{p_4}}} - {{s^{{p_1}}_ {{z_1}}}{s^{{z_1}}_ {{p_3}}}{\tilde{a}^{{p_2}}_ {{p_4}}}} - {{s^{{p_2}}_ {{z_1}}}{s^{{z_1}}_ {{p_4}}}{\tilde{a}^{{p_1}}_ {{p_3}}}} - {{s^{{p_1}}_ {{y_1}}}{s^{{y_1}}_ {{p_4}}}{\tilde{a}^{{p_2}}_ {{p_3}}}}
    $$

    $$
    \qquad \qquad \qquad \quad + ~ {{s^{{p_2}}_ {{z_1}}}{s^{{z_1}}_ {{p_3}}}{\tilde{a}^{{p_1}}_ {{p_4}}}} + {{s^{{p_1}}_ {{z_1}}}{s^{{p_2}}_ {{z_2}}}{s^{{z_1}}_ {{p_3}}}{s^{{z_2}}_ {{p_4}}}} + {{s^{{p_1}}_ {{y_1}}}{s^{{p_2}}_ {{z_1}}}{s^{{z_1}}_ {{p_3}}}{s^{{y_1}}_ {{p_4}}}}\bigr) .
    $$

    Note that:

    - the tilde in $\tilde{a}$ denotes normal order with respect to single-product vacuum, and
    - Einstein summation convention is implied, i.e., indices that appear twice in a given product (once in superscript, once in a subscript) are summed over.

    ## Operators

    Development of SeQuant is primarily motivated by the perturbative many-body methods, collectively referred to here as Many-Body Perturbation Theory (MBPT). Examples of such methods include the [coupled-cluster (CC) method](https://doi.org/10.1017/CBO9780511596834) and [GW](https://doi.org/10.1103/PhysRev.139.A796). The typical use case is to compute canonical forms of products of operators. For example, consider the coupled-cluster doubles (CCD) method.
    _Amplitudes_ $t^{i_1 i_2}_ {a_1 a_2}$ of the cluster operator,

    $$
    \hat{t} \equiv \hat{t}_ 2 = \frac{1}{4} t^{i_1 i_2}_ {a_1 a_2} a_ {i_1 i_2}^{a_1 a_2},
    $$

    are determined by solving the CCD equations:

    $$
    \forall i_1, i_2, a_1, a_2: \quad 0 = \langle0\vert a^{i_1 i_2}_ {a_1 a_2} \exp(-\hat{t}_ 2) \hat{H} \exp(\hat{t}_ 2) \vert 0 \rangle = \langle0\vert a^{i_1 i_2}_ {a_1 a_2} \bigl( \hat{H} + [\hat{H}, \hat{t}_ 2] + \frac{1}{2} [[\hat{H}, \hat{t}_ 2], \hat{t}_ 2] \bigr) \vert 0 \rangle.
    $$

    A pedestrian way to compose such expression is to define a cluster operator object using SeQuant tensors and normal-ordered operators:

    ```c++
    auto t2 =
          ex<Constant>(rational(1, 4)) *
          ex<Tensor>(L"t", std::array{L"a_1", L"a_2"}, std::array{L"i_1", L"i_2"},
                     Symmetry::antisymm) *
          ex<FNOperator>(std::array{L"a_1", L"a_2"}, std::array{L"i_1", L"i_2"});
    ```

    The normal-ordered Hamiltonian is defined similarly as a sum of 1- and 2-body contributions:

    ```c++
    auto H = ex<Tensor>(L"f", std::array{L"p_1"}, std::array{L"p_2"},
                        Symmetry::nonsymm) *
             ex<FNOperator>(std::array{L"p_1"}, std::array{L"p_2"})
             +
             ex<Constant>(rational(1, 4)) *
             ex<Tensor>(L"g", std::array{L"p_1", L"p_2"},
                              std::array{L"p_3", L"p_4"}, Symmetry::antisymm) *
             ex<FNOperator>(std::array{L"p_1", L"p_2"},
                            std::array{L"p_3", L"p_4"});
    ```

    Note that the compact definition of the Hamiltonian is due to the use of the union ($p$) of base occupied ($i$) and unoccupied ($a$) spaces. Many other symbolic algebras only support use of nonoverlapping base spaces, in terms of which Hamiltonian and other tensor expressions would have a much more verbose form.

    Commutator of the Hamiltonian and cluster operator is trivially composed:

    ```c++
    inline auto commutator(auto op1, auto op2) {
      return simplify(op1 * op2 - op2 * op1);
    }

    auto c_ht = commutator(H, t2);
    ```

    Note the use of `simplify` to rewrite an expression in a simpler form. Its role will be emphasized later.

    Unfortunately, we immediately run into the limitation of the "pedestrian" approach. Namely,
    the double commutator cannot be correctly obtained as

    ```c++
    auto c_htt = ex<Constant>(rational(1, 2)) * commutator(commutator(H, t2), t2);
    ```

    due to the explicit use of specific dummy indices in the definition of `t2`. Using it more than once in a given product will produce an expresson where each dummy index appears more than 2 times, breaking the Einstein summation convention.

    The issue is actually not the reuse of the same of operator object, but the reuse of dummy indices. A straightforward, but brittle, solution is to ensure that each dummy index is only used once. E.g., 
    to use $\hat{t}_2$ more than once in an expression we must make several versions of it, each with a separate set of dummy indices:

    ```c++
    auto t2_0 =
          ex<Constant>(rational(1, 4)) *
          ex<Tensor>(L"t", std::array{L"a_1", L"a_2"}, std::array{L"i_1", L"i_2"},
                     Symmetry::antisymm) *
          ex<FNOperator>(std::array{L"a_1", L"a_2"}, std::array{L"i_1", L"i_2"});

    auto t2_1 =
          ex<Constant>(rational(1, 4)) *
          ex<Tensor>(L"t", std::array{L"a_3", L"a_4"}, std::array{L"i_3", L"i_4"},
                     Symmetry::antisymm) *
          ex<FNOperator>(std::array{L"a_3", L"a_4"}, std::array{L"i_3", L"i_4"});

    auto c_htt = ex<Constant>(rational(1, 4)) * commutator(commutator(H, t2_0), t2_1);
    ```

    This is too error-prone for making complex expressions. A better way is to represent $\hat{t}_2$ by an object that generates tensor form with unique dummy indices generated on the fly. in SeQuant such MBPT _operators_ live in `mbpt` namespace. The entire CCD amplitude equation is evaluated as follows:

    ```c++
    #include <SeQuant/core/context.hpp>
    #include <SeQuant/core/expr.hpp>
    #include <SeQuant/core/tensor.hpp>
    #include <SeQuant/core/op.hpp>
    #include <SeQuant/domain/mbpt/convention.hpp>
    #include <SeQuant/domain/mbpt/op.hpp>

    inline auto commutator(auto op1, auto op2) {
      return op1 * op2 - op2 * op1;
    }

    int main() {
      using namespace sequant;
      using namespace sequant::mbpt;
      set_default_context(Context(make_min_sr_spaces(), Vacuum::SingleProduct));

      auto hbar = H(2) + commutator(H(2), T_(2)) + ex<Constant>(rational(1,2)) * commutator(commutator(H(2), T_(2)), T_(2));
      auto ccd_eq = vac_av(P(2) * hbar);
      std::wcout << "<" << to_latex(P(2) * hbar) << "> = " << to_latex(ccd_eq) << std::endl;

      return 0;
    }
    ```

    The result is

    $$
    \bigl({{{\frac{1}{4}}}{A^{{a_1}{a_2}}_ {{i_1}{i_2}}}{\bar{g}^{{i_1}{i_2}}_ {{a_1}{a_2}}}} - {{{\frac{1}{2}}}{A^{{a_1}{a_2}}_ {{i_1}{i_2}}}{f^{{a_3}}_ {{a_1}}}{\bar{t}^{{i_1}{i_2}}_ {{a_2}{a_3}}}} - {{A^{{a_1}{a_2}}_ {{i_1}{i_2}}}{\bar{g}^{{i_1}{a_3}}_ {{i_3}{a_1}}}{\bar{t}^{{i_2}{i_3}}_ {{a_2}{a_3}}}} + {{{\frac{1}{8}}}{A^{{a_1}{a_2}}_ {{i_1}{i_2}}}{\bar{g}^{{a_3}{a_4}}_ {{a_1}{a_2}}}{\bar{t}^{{i_1}{i_2}}_ {{a_3}{a_4}}}}
    $$

    $$
    \qquad \quad + ~ {{{\frac{1}{8}}}{A^{{a_1}{a_2}}_ {{i_1}{i_2}}}{\bar{g}^{{i_1}{i_2}}_ {{i_3}{i_4}}}{\bar{t}^{{i_3}{i_4}}_ {{a_1}{a_2}}}} + {{{\frac{1}{2}}}{A^{{a_1}{a_2}}_ {{i_1}{i_2}}}{f^{{i_1}}_ {{i_3}}}{\bar{t}^{{i_2}{i_3}}_ {{a_1}{a_2}}}} + {{{\frac{1}{16}}}{A^{{a_1}{a_2}}_ {{i_1}{i_2}}}{\bar{g}^{{a_3}{a_4}}_ {{i_3}{i_4}}}{\bar{t}^{{i_1}{i_2}}_ {{a_3}{a_4}}}{\bar{t}^{{i_3}{i_4}}_ {{a_1}{a_2}}}}
    $$

    $$
    \quad - ~ {{{\frac{1}{4}}}{A^{{a_1}{a_2}}_ {{i_1}{i_2}}}{\bar{g}^{{a_3}{a_4}}_ {{i_3}{i_4}}}{\bar{t}^{{i_1}{i_3}}_ {{a_3}{a_4}}}{\bar{t}^{{i_2}{i_4}}_ {{a_1}{a_2}}}} - {{{\frac{1}{4}}}{A^{{a_1}{a_2}}_ {{i_1}{i_2}}}{\bar{g}^{{a_3}{a_4}}_ {{i_3}{i_4}}}{\bar{t}^{{i_1}{i_2}}_ {{a_1}{a_3}}}{\bar{t}^{{i_3}{i_4}}_ {{a_2}{a_4}}}}
    $$

    $$
     ~ + {{{\frac{1}{2}}}{A^{{a_1}{a_2}}_ {{i_1}{i_2}}}{\bar{g}^{{a_3}{a_4}}_ {{i_3}{i_4}}}{\bar{t}^{{i_1}{i_3}}_ {{a_1}{a_3}}}{\bar{t}^{{i_2}{i_4}}_ {{a_2}{a_4}}}}\bigr)
    $$

    The use of MBPT operators rather than their tensor-level forms not only solves problems with the reuse of dummy indices, but also allows to implement additional optimizations such as algebraic simplifications of complex operator expressions and avoiding evaluation of operator products whose vacuum expectation values are guaranteed to vanish. This allows very efficient derivation of complex equations, e.g. CC equations through CCSDTQ are derived in a fraction of a second on a  laptop:

    ```sh
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
    ```
