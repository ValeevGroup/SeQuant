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

The relevant function is :func:`sequant::io::latex::to_string`. SeQuant provides support for `LaTeX <https://en.wikipedia.org/wiki/LaTeX>`_ conversion on most object types. If the shorthands
header is included, this functionality is also exposed as :func:`sequant::to_latex`.

The produced LaTeX code is not self-contained. It needs to be embedded in a suitable math environment. Furthermore, it assumes that the
`tensor <https://ctan.org/pkg/tensor>`_ has been loaded.

.. warning::
   Some objects also have a :func:`to_latex` member function. However, we don't provide any guarantees about which types do and they may also be
   removed in future SeQuant versions. Hence, you should always prefer using :func:`sequant::io::latex::to_string`.

.. note::
   The produced LaTeX code is typically not human-friendly (that is, hard to read) but should compile successfully once embedded in a suitable LaTeX
   document.



.. _io-Serialization:

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
   Nullary         '(' Sum ')' | Number | Tensor | Variable | RealImagPart     
   RealImagPart    ('Re' | 'Im') '[' Sum ']'                                        Real/imaginary part of the wrapped expression
   Number                                                                           Integer, Floating point or fraction
   Tensor          Name IndexGroup SymmetrySpec?                                   
   IndexGroup      | '{' IndexList? ( ';' IndexList? ( ';' IndexList? )? )? '}'     | Meaning is {<bra>;<ket>;<aux>}
                   | '^{' IndexList? '}_{' IndexList '}'                            | Meaning is ^{<ket>}_{<bra>} (no aux)
                   | '_{' IndexList? '}^{' IndexList '}'                            | Meaning is _{<bra>}^{<ket>} (no aux)
   IndexList       Index ( ',' Index )?
   Index           IndexSpaceName '_'? Integer
   IndexSpaceName                                                                    Name but no underscore allowed
   SymmetrySpec    ':' ( [ASN] ( '-' [SCNHA] ( '-' [SN] ( '-' [EON] )? )? )? )      :<Symmetry>-<BraKetSymmetry or Hermiticity>-<ColumnSymmetry>-<ConjugationParity>
   Variable        Name
   Name                                                                              Single word (may include Unicode chars)
   ==============  ===============================================================  ===========================================

The single-letter codes in ``SymmetrySpec`` abbreviate the corresponding enumerators, one letter per field: ``[ASN]`` is the tensor's
permutational :class:`Symmetry <sequant::Symmetry>` (``A`` = Antisymm, ``S`` = Symm, ``N`` = Nonsymm); ``[SCN]`` is its
:class:`BraKetSymmetry <sequant::BraKetSymmetry>` (``S`` = Symm, ``C`` = Conjugate, ``N`` = Nonsymm); and the final ``[SN]`` is its
:class:`ColumnSymmetry <sequant::ColumnSymmetry>` (``S`` = Symm, ``N`` = Nonsymm — this field has no Conjugate case). A tensor's trailing
``:A-C-S`` in the example below therefore reads "antisymmetric, conjugate bra-ket symmetric, symmetric column".

``RealImagPart`` spells the :class:`RealPart <sequant::RealPart>` and :class:`ImagPart <sequant::ImagPart>` expression nodes, which
:func:`simplify() <sequant::simplify>` emits when it folds a pair of complex-conjugate summands (``A + A*`` becomes ``2 Re[A]``,
``A - A*`` becomes ``2i Im[A]``). The square brackets follow the LaTeX these nodes render as (``\Re\left[...\right]``) and are what
keeps the spelling apart from a tensor or a variable named ``Re``/``Im``: neither of those admits a ``[``, so ``Re{i_1;a_1}`` is a
tensor, ``Re`` alone is a variable, and only ``Re[`` opens a real part. Deserialization goes through the same eager composition rules
as the :func:`real_part() <sequant::real_part>` / :func:`imaginary_part() <sequant::imaginary_part>` builders, so ``Re[Re[x]]`` yields
``Re[x]``, ``Im[Re[x]]`` yields ``0``, and a real scalar prefactor hoists out of the brackets (``Re[1/2 x]`` yields ``1/2 Re[x]``).

:func:`sequant::io::serialization::from_string<ExprPtr>` will start at rule :code:`Expression`, whereas
:func:`sequant::io::serialization::from_string<ResultExpr>` will start at :code:`Result`.


Examples
""""""""

The following parses a residual expression ``R1`` as a sum of four terms, the last one antisymmetric-bra-ket-symmetric-column (``:A-C-S``) and
scaled by the constant ``1/2``:

::

   R1{u1;i1} = f{u1;i1} - Ym1{u1;u2} f{u2;i1} - Ym1{u3;u2} * g{u1,u2;u3,i1} + 1/2 Ym2{u1,u4;u_2,u_3} g{u2,u3;u4,i1}:A-C-S

The following parses the real part of a fully contracted product, scaled by ``2`` -- the shape :func:`simplify() <sequant::simplify>`
leaves a folded conjugate pair in:

::

   2 Re[h{i_1;a_1} * t{a_1;i_1}]
