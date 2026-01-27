Parsing and Deparsing
=====================

SeQuant supports creating expression objects by parsing an appropriately formatted input string. This process is referred to as *parsing*.
Furthermore, SeQuant can also serialize expression objects into this string representation, which is referred to as *deparsing*.

The functions for performing these tasks are :func:`parse_expr(…)`, :func:`parse_result_expr(…)` and :func:`deparse(…)` respectively. The former two
differ in that one of them produces a :class:`ExprPtr` and the other a :class:`ResultExpr` object as their output. This implies that the parsed text
in the latter case is supposed to contain a result specification (see below).

As SeQuant's capabilities develop over time, it can be necessary to adapt the syntax of this text representation, referred to as the *parse syntax*.
By default, the abovementioned functions will always expect and produce text according to the latest parse syntax specification. However, they can be
instructed to work with a different version by explicitly specifying a :class:`SyntaxVersion` when calling them or using one of the versioned function
calls.

.. warning::
   All syntax versions except :class:`SyntaxVersion::Latest` are considered deprecated. Support for them will remain available for some time but might
   get removed in future versions of SeQuant.

Generally speaking, parse and deparse are inverse operations.


Customizations
--------------

The parse and deparse functions accept customization options of type :class:`ParseOptions` and :class:`DeparseOptions` respectively. They allow for
specification of a specific :class:`SyntaxVersion` and things like how to deal with symmetry annotations. For parsing, the latter refers to default
symmetries for objects that don't explicitly specify their symmetries and for deparsing this refers to whether or not those annotations are included
in the output.

For an overview of all customization options, please refer to the documentation of those classes in the API reference.


Syntax
------

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

:func:`parse_expr` will start at rule `Expression`, whereas :func:`parse_result_expr` will start at `Result`


Examples
""""""""

::

   R1{u1;i1} = f{u1;i1} - Ym1{u1;u2} f{u2;i1} - Ym1{u3;u2} * g{u1,u2;u3,i1} + 1/2 Ym2{u1,u4;u_2,u_3} g{u2,u3;u4,i1}:A-C-N

