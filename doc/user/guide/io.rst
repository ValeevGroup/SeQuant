***
I/O
***
SeQuant provides different forms of I/O for its objects. Generally, they can be found in the :code:`sequant::io` namespace. Different formats live in
different sub-namespaces. E.g. :func:`sequant::io::latex::to_string` yields a LaTeX representation of the respective object.

.. note::
   Frequently used I/O functions also have shorthands that can be brought in by including the :file:`sequant/core/io/shorthands.hpp` header. These
   shorthands live directly in the :code:`sequant` namespace. The example from above could be achieved via :func:`sequant::to_latex` once that
   header is included.

LaTeX
=====

The relevant function is :func:`sequant::io::latex::to_string`. SeQuant provides support for LaTeX conversion on most object types. If the shorthands
header is included, this functionality is also exposed as :func:`sequant::to_latex`.

The produced LaTeX code is not self-contained. It needs to be embedded in a suitable math environment. Furthermore, it assumes that the
`tensor <https://ctan.org/pkg/tensor>`_ has been loaded.

.. warning::
   Some objects also have a :func:`to_latex` member function. However, we don't provide any guarantees about which types do and they may also be
   removed in future SeQuant versions. Hence, you should always prefer using :func:`sequant::io::latex::to_string`.

.. note::
   The produced LaTeX code is typically not human-friendly (that is, hard to read) but should compile successfully once embedded in a suitable LaTeX
   document.



Serialization
=============

SeQuant provides `serialization <https://en.wikipedia.org/wiki/Serialization>`_ support. This means that SeQuant objects can be transformed into some
intermediary storable format (*serialization*), which can later be converted back (*deserialization*) to yield the original SeQuant object.

Functions associated with this live in the :code:`sequant::io::serialization` namespace. At the moment, there is only a text format implemented via
:func:`to_string` and :func:`from_string`. More formats may be added at a later point. It is important to note that in order to work with C++'s
type system, :func:`from_string` is a template that takes in the type of object it is supposed to produce. The most important variants are
:func:`from_string<ExprPtr>` and :func:`from_string<ResultExpr>` where the former also allows to deserialize any :class:`Expr` subclass.
Depending on the used template type, the expected format of the input will change accordingly.

As SeQuant's capabilities develop over time, it can be necessary to adapt the syntax of the serialization format.
By default, the abovementioned functions will always adhere to the latest parse syntax specification. However, they can be
instructed to work with a different version by explicitly specifying a :class:`SerializationSyntax` when calling them or using one of the versioned
function calls.

.. warning::
   All syntax versions except :class:`SerializationSyntax::Latest` are considered deprecated. Support for them will remain available for some time but
   might get removed in future versions of SeQuant.

If the shorthands header is included, text-based serialization and deserialization is available as :func:`sequant::serialize` and
:func:`sequant::deserialize<T>`.


Customizations
--------------

The exact serialization format can be customized to some extend by providing :class:`SerializationOptions` and :class:`DeserializationOptions`
instances when calling the respective functions. Among other things, they allow for specification of a specific :class:`SerializationSyntax`. For an
overview of all customization options, please refer to the documentation of those classes in the API reference.

.. note::
   Typically, you should leave these at their defaults to yield the most reliable behavior.


SerializationSyntax
-------------------

We will loosely follow `EBNF <https://en.wikipedia.org/wiki/Extended_Backus%E2%80%93Naur_form>`_ syntax for specifying the rules of the grammar
describing the parse syntax. However, we don't define a full formal grammar here as minor details and precedence issues may be left out/simplified in
order to increase readability.

V1
^^

.. note::
   Known limitations:

   - No support for second-quantized (normal-ordered) operators
   - No support for complex numbers
   - No support for representing operators (e.g. symmetrizers) explicitly

.. table::
   :widths: auto

   ==============  ===============================================================  ===========================================
   Component       Rule                                                             Note
   ==============  ===============================================================  ===========================================
   Result          (Tensor | Variable) ('=' | '<-') Expression                     
   Expression      Sum?                                                            
   Sum             Product ( ('+' | '-') Product)*                                 
   Product         Nullary ( '*'? Nullary )*                                        Explict '*' use is optional
   Nullary         '(' Sum ')' | Number | Tensor | Variable                        
   Number                                                                           Integer, Floating point or fraction
   Tensor          Name IndexGroup SymmetrySpec?                                   
   IndexGroup      | '{' IndexList? ( ';' IndexList? ( ';' IndexList? )? )? '}'     | Meaning is {<bra>;<ket>;<aux>}
                   | '^{' IndexList? '}_{' IndexList '}'                            | Meaning is ^{<ket>}_{<bra>} (no aux)
                   | '_{' IndexList? '}^{' IndexList '}'                            | Meaning is _{<bra>}^{<ket>} (no aux)
   IndexList       Index ( ',' Index )?
   Index           IndexSpaceName '_'? Integer
   IndexSpaceName                                                                    Name but no undersore allowed
   SymmetrySpec    ':' ( [ASN] ( '-' [SCN] ( '-' [SN] )? )? )                        :<Symmetry>-<BraKetSymmetr>-<ColumnSymmetry>
   Variable        Name
   Name                                                                              Single word (may include Unicode chars)
   ==============  ===============================================================  ===========================================

:func:`parse_expr` will start at rule :code:`Expression`, whereas :func:`parse_result_expr` will start at :code:`Result`.


Examples
""""""""

::

   R1{u1;i1} = f{u1;i1} - Ym1{u1;u2} f{u2;i1} - Ym1{u3;u2} * g{u1,u2;u3,i1} + 1/2 Ym2{u1,u4;u_2,u_3} g{u2,u3;u4,i1}:A-C-N

