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
