.. _external_interface_v2:

Driver File v2
==============


The external interface implementation provides a series of individual processing steps that can be selected and chained by means of the JSON driver
file. Within this file, each processing step is represented as a JSON object with the following properties:

* :code:`kind` (required, String): This is a unique identifier for every different kind of steps and therefore determines what this step will do when
  exectuted.
* :code:`id` (String): Unique identifier for this particular step. This can be used to refer to (the output of) other steps.
* :code:`options` (Object): Set of options to tune the behavior of the step as needed. Some steps don't have options, whereas others require you to
  specify them. For most steps, options are optional. The properties of the options object depend on the step's kind - see
  :ref:`external_interface_step_kinds`.
* :code:`outputs` (Object): This can be used to give human-readable names to individual outputs of the current step in form of name-output pairs. See
  :ref:`external_interface_outputs` - the step ID in this case is implicitly the current step's ID and must not be included explicitly. These names
  are automatically propagated through the processing chain (for most kinds of steps).
* :code:`inputs` (String or Array of Strings): Specifies the inputs of this step, which must be outputs of other steps. See
  :ref:`external_interface_outputs`. 


The different steps are listed objects in the top-level :code:`steps` array. Steps are processed in order according to their order in this array.


.. _external_interface_outputs:

Referencing Outputs
-------------------

Outputs can be referenced by means of their IDs. An ID has the general format :code:`<step_id>.<output_id>`. :code:`<step_id>` must refer to the
:code:`id` of one of the steps that have been executed before. :code:`<output_id>` can be the human-readable names specified via a step's
:code:`outputs` property or an integer. In the latter case, the integer refers to the index of the output. Outputs are indexed starting from zero in
the order they are produced.

:code:`<output_id>` can also be an expression enclosed in square brackets in order to refer to multiple outputs at once. This can be a comma-delimited
(no spaces!) list of output names or indices, or a range of indices which are of the form :code:`<from>-<to>` like :code:`0-5`.

Examples of output IDs:

* :code:`my_step.0`
* :code:`my_step.some_name`
* :code:`my_step.[0-2]`
* :code:`my_step.[0-2,some_name,5]`

Finally, to refer to all outputs of a given step, just use the step's ID without anything appended to it. That is, `my_step` would automatically refer
to all outputs of the respective step.


.. external_interface_step_kinds:

Available Processing Step Kinds
-------------------------------


canonicalize
^^^^^^^^^^^^

Canonicalizes the input expressions.


cse
^^^

Performs common-subexpression elimination (CSE).


density_fitting
^^^^^^^^^^^^^^^

Inserts the density-fitting decomposition of the two-electron integrals.


export
^^^^^^

Exports the given expressions as code.


index_batching
^^^^^^^^^^^^^^

Configures the computation of results to happen in batches over certain result indices.


optimize
^^^^^^^^

Symbolically rewrites the expressions for an improved numerical evaluation. Most prominently, this factors tensor contractions into a series of binary
tensor contractions.


output
^^^^^^

Outputs expressions in the chosen markup style. Mainly intended for debugging purposes.


project
^^^^^^^

Performs the chosen projection with the inputs.


read_input
^^^^^^^^^^

Reads and parses expressions from files.


simplify
^^^^^^^^

Simplifies the given expressions.


spintracing
^^^^^^^^^^^

Spintraces the given expressions. That is, it performs spin-integration and potentially also spin-summation.


substitute
^^^^^^^^^^

Makes substitutions in the given expressions


to_export_tree
^^^^^^^^^^^^^^

Converts the given expressions into a tree data structure suitable for exports.


validate
^^^^^^^^

Validates the given expressions



Examples
--------

.. code-block:: json

   {
     "driver_format_version": 2,
     "index_spaces": [
       {
         "label": "a",
         "size": 1000,
         "real_valued": true
       },
       {
         "label": "u",
         "size": 5,
         "real_valued": true
       },
       {
         "label": "i",
         "size": 80,
         "real_valued": true
       },
       {
         "label": "F",
         "size": 1500,
         "real_valued": true
       }
     ],
     "steps": [
       {
         "id": "input",
         "kind": "read_input",
         "options": {
           "file_path": [
             "nevpt2/nevpt2_en0.inp",
             "nevpt2/nevpt2_en.inp",
   
             "nevpt2/nevpt2_res1_i1.inp",
             "nevpt2/nevpt2_res1_s0.inp",
             "nevpt2/nevpt2_res2_s1_singles.inp",
             "nevpt2/nevpt2_res1_s1.inp",
   
             "nevpt2/nevpt2_res2_p0.inp",
             "nevpt2/nevpt2_res2_p2.inp",
             "nevpt2/nevpt2_res2_i2.inp",
             "nevpt2/nevpt2_res2_p1.inp",
             "nevpt2/nevpt2_res2_s1.inp",
             "nevpt2/nevpt2_res2_s2.inp"
           ],
           "default_symmetry": "antisymmetric"
         },
         "outputs": {
           "ecc0": "0",
           "ecc": "1",
           "en": "0-1",
           "res1": "2-5",
           "res1_i1": "2",
           "res1_s0": "3",
           "res2_s1_singles": "4",
           "res1_s1": "5",
           "res2": "6-11",
           "res2_p0": "6",
           "res2_p2": "7",
           "res2_i2": "8",
           "res2_p1": "9",
           "res2_s1": "10",
           "res2_s2": "11",
           "res": "2-11"
         }
       },
       {
         "kind": "validate",
         "inputs": "input"
       },
       {
         "id": "DF",
         "kind": "density_fitting",
         "inputs": "input",
         "options": {
           "auxiliary_space": "F"
         }
       },
       {
         "id": "traced",
         "kind": "spintracing",
         "inputs": "DF",
         "options": {
           "algorithm": "closed_shell"
         }
       },
       {
         "id": "biorth",
         "kind": "project",
         "inputs": "traced.res",
         "options": {
           "method": "biorthogonal"
         }
       },
       {
         "id": "opt",
         "kind": "optimize",
         "inputs": [
           "traced.en",
           "biorth"
         ]
       },
       {
         "id": "treeify",
         "kind": "to_export_tree",
         "inputs": "opt"
       },
       {
         "kind": "export",
         "inputs": "treeify",
         "options": {
           "language": "itf",
           "optimize": true,
           "output": "nevpt2_v2.itfaa",
           "grouping": {
             "Energy0": "ecc0",
             "Energy": "ecc",
             "Residual": "res"
           },
           "relative_order": [
             "ecc0",
             "ecc",
             "res1_i1",
             "res1_s0",
             "res2_s1_singles",
             "res1_s1",
             "res2_p0",
             "res2_p2",
             "res2_i2",
             "res2_p1",
             "res2_s1",
             "res2_s2"
           ],
           "imports": {
             "R2{a1;i1}": "R2:ec"
           },
           "meta": {
             "index_spaces": {
               "a": {
                 "name": "External",
                 "tag": "e"
               },
               "u": {
                 "name": "Active",
                 "tag": "a"
               },
               "i": {
                 "name": "Closed",
                 "tag": "c"
               },
               "F": {
                 "name": "BasisMp2Fit",
                 "tag": "F"
               }
             },
             "min_index_id": 1
           }
         }
       }
   }

