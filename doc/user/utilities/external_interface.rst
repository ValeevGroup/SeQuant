External Interface
==================

The external interface is supposed to allow interfacing with SeQuant from the outside, without having to write a C++ program that links to the SeQuant
library. The idea is to specify equations in text form and then submit them to SeQuant for processing. This works by using a JSON driver file that
contains instructions for what you want SeQuant to do.


.. _extint-input:

Input Format
------------

See :ref:`io-Serialization` for the format in which the input equations are expected to be.

.. note::
   It is assumed that the input always specifies a result. That is, it is of the format :code:`lhs = rhs`. Furthermore, every input file may only
   contain a single result.



Driver File
-----------

There currently are two mostly independent implementations of the external interface available. Which version you want to use is determined by the
:code:`driver_format_version` key in the JSON driver. If you set it to a value of :code:`1`, you will get the old/legacy version of the interface. The
rest of the JSON file is expected to follow the syntax described in :ref:`external_interface_v1`. If you set it to :code:`2`, you will be using the
new, modular external interface. Its driver syntax is described in :ref:`external_interface_v2`. If the field is unset, the code defaults to :code:`1`
for reasons of backwards compatibility.

The main difference between the two versions is that the newer version is much more modular and flexible. Therefore, it is recommended that you use
that for all new tasks.

.. note::
   All paths specified in the driver file are understood to be relative to the JSON file's location (unless absolute paths are used, of course).


Common Syntax
^^^^^^^^^^^^^

.. _external_interface_idx_space:

Index Space Specification
"""""""""""""""""""""""""

Every driver file has to contain the definition of index spaces that are used in expressions. It lives under the top-level key :code:`index_spaces`
and is expected to be an array of objects. These objects can have the following properties

* :code:`label` (required, String): Label used in expressions for indices in this space, e.g. in :code:`i1` the label is :code:`i`.
* :code:`size` (required, Integer): The (approximate) size/dimension of indices in this index space. This affects things like factorization into
  binary contractions.
* :code:`real_valued` (Boolean): Whether the field used in this index space is real- rather than complex-valued. This affects for instance tensors
  with hermitian braket symmetry.
