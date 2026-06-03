External interface
==================

The external interface is supposed to allow interfacing with SeQuant from the outside, without having to write a C++ program that links to the SeQuant
library. The idea is to specify equations in text form and then submit them to SeQuant for processing. This works by using a JSON driver file that
contains instructions for what you want SeQuant to do.


.. _extint-input:

Input format
------------

See :ref:`io-Serialization`.

.. note::
   It is assumed that the input always specifies a result. That is, it is of the format :code:`lhs = rhs`. Furthermore, every input file may only
   contain a single result.


Driver file
-----------

The top-level entry in the JSON file specifies what action SeQuant should take. At the moment, only code-generation into the ITF format is supported.
Hence, every driver file currently has to start like this:

.. code-block:: json
   {
       "code_generation": {
           "output_format": "itf",
           ...
       }
   }

.. note::
   All paths specified in the driver file are understood to be relative to the JSON file's location (unless absolute paths are used, of course).

Beyond :code:`output_format`, the following top-level fields exist:

* :code:`output_path` (required): The path to the file the generated code is written to
* :code:`default_options` (optional): Allows specification of default processing options
* :code:`code_blocks` (required): Specifies a list of code blocks (groups of expressions)

Code block
^^^^^^^^^^
Every code block has to have a :code:`name` and a list of :code:`results`. The former is effectively the name of the to-be-generated function that
computes the individual results (expressions), whereas the latter is a list of expressions that shall be computed.

Every result has these mandatory fields:

* :code:`name`: The name of the tensor/scalar variable that shall hold the result of the computed expression
* :code:`equation_file`: Path to the file containing the input expression (cmp :ref:`extint-input`)

Additionally, the following *processing options* may be given. All of them may also be specified as part of the :code:`default_options` block in which
case those values are used, unless explicitly overwritten.

* :code:`density_fitting`: Whether to perform the density-fitting decomposition of the two-electron integral
* :code:`term_by_term`: Whether to split sums into individual summands for processing and code-generation. This yields to more readable but less
  performant code.
* :code:`optimize`: Whether to factorize the equations into a series of binary contractions
* :code:`subexpression_elimination`: Whether to eliminate common subexpressions (only possible when factorizing into binary contractions)
* :code:`expand_symmetrizer`: Whether to explicitly expand (write out) symmetrization operators
* :code:`spintracing`: What kind of spintracing to perform (if any). Possible options are

  * :code:`none`: Don't perform spintracing
  * :code:`closed_shell`: Apply spintracing using an algorithm suitable for closed-shell systems
  * :code:`rigorous`: Apply spintracing using an algorithm that should work for all cases, but is less efficient than :code:`closed_shell`

* :code:`projection`: What kind of projection/transformation to perform with the final result

  * :code:`primitive`: Don't do anything
  * :code:`biorthogonal`: Transform the result into a biorthogonal basis (only applicable to non-scalar results)


Index space specification
^^^^^^^^^^^^^^^^^^^^^^^^^

Every driver file has to contain the definition of index spaces that are used in expressions. This definition lives under the :code:`index_spaces`
element. Every specification has to provide the following attributes:

* :code:`name`: Name of the index space
* :code:`tag`: Tag for this index space (if any). Tags are used to encode the spaces of a tensor's indices in its name. May be empty.
* :code:`label`: Label used in expressions for indices in this space, e.g. in :code:`i1` the label is :code:`i`.
* :code:`size`: The size/dimension of indices in this index space. This affects factorization into binary contractions.

Example
^^^^^^^

.. code-block:: json

   {
       "code_generation": {
           "output_format": "itf",
           "default_options": {
               "density_fitting": false,
               "term_by_term": false,
               "spintracing": "closed_shell",
               "projection": "biorthogonal"
           },
           "output_path": "something.itfaa",
           "code_blocks": [
               {
                   "name": "First",
                   "results": [
                       {
                           "name": "One",
                           "equation_file": "first.inp",
                           "projection": "primitive"
                       }
                   ]
               },
               {
                   "name": "Second",
                   "results": [
                       {
                           "name": "Two",
                           "equation_file": "second.inp",
                           "projection": "primitive"
                       },
                       {
                           "name": "Three",
                           "equation_file": "third.inp",
                           "spintracing": "rigorous"
                       }
                   ]
               },
           ]
       },
       "index_spaces": [
           {
               "name": "virtual",
               "tag": "e",
               "label": "a",
               "size": 100
           },
           {
               "name": "active",
               "tag": "a",
               "label": "u",
               "size": 5
           },
           {
               "name": "occupied",
               "tag": "c",
               "label": "i",
               "size": 10
           }
       ]
   }
