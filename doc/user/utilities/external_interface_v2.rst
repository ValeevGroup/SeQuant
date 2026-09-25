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
* :code:`inherit_options_from` (String): The :code:`id` of another step of the same kind (anywhere in the :code:`steps` array) whose options
  shall be used as the base for this step's options. The local :code:`options` (if any) are applied on top as a
  `JSON merge patch <https://datatracker.ietf.org/doc/html/rfc7396>`_: nested objects are merged recursively, all other values replace the
  inherited ones and a value of :code:`null` removes the inherited option. The referenced step may itself inherit its options from yet another step.
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


How the Pipeline Executes
--------------------------

Steps are executed strictly in the order in which they appear in the :code:`steps` array - there is no separate dependency resolution based on
:code:`inputs`. Given IDs must be unique across the whole driver file.

Two kinds of data can flow between steps: *Expression* data (one or more named result expressions) and *ExportTree* data (the export-ready tree form
produced by the :code:`to_export_tree` step). Every step kind's documentation below states which of these it consumes and which it produces. Feeding a
step data of the wrong kind (e.g. passing plain expressions into a step that expects an :code:`ExportTree`) is an error.

Not every step produces output data - they exist purely for their side effects (e.g. aborting on invalid input, printing to the console or writing
something to a file) - and therefore cannot be referenced by another step's :code:`inputs`.

Errors encountered anywhere in the pipeline (be it while parsing a step's options or while executing it) abort the entire run; there is no
partial/best-effort processing.


.. _external_interface_step_kinds:

Available Processing Step Kinds
-------------------------------


batch_indices
^^^^^^^^^^^^^

Configures the computation of results to happen in batches over certain result indices. For every result, its indices are ranked by the size of the
index space they belong to. The smallest :code:`min_unbatched` indices are always left un-batched; among the remaining indices, up to
:code:`max_batched` become batch indices, chosen from the end of the size-ranked list selected by :code:`selection_strategy`. This step only annotates
which indices shall be batched - whether batching is actually realized (and how) is up to whichever code-generation backend is used by the subsequent
:code:`export` step, and not every backend supports it.

.. list-table::
   :header-rows: 1

   * - Input
     - Output
   * - ExportTree
     - ExportTree


.. list-table:: Options
   :header-rows: 1

   * - Option
     - Description
     - Default
     - Required
   * - :code:`min_unbatched`
     - Minimum number of (smallest-space) indices of a given result that are always left un-batched.
     - 2
     - No
   * - :code:`max_batched`
     - Upper limit on the number of indices that may be selected for batching.
     - unbounded
     - No
   * - :code:`selection_strategy`
     - Which end of the size-ranked index list batching candidates are taken from. :code:`largest` selects the largest-space indices (after excluding
       the :code:`min_unbatched` smallest ones); :code:`smallest` selects the smallest-space indices remaining after excluding the
       :code:`min_unbatched` largest ones.
     - :code:`largest`
     - No


canonicalize
^^^^^^^^^^^^

Canonicalizes the input expressions.

.. list-table::
   :header-rows: 1

   * - Input
     - Output
   * - Expression
     - Expression


This step does not accept any options.


cse
^^^

Performs common-subexpression elimination (CSE).

.. list-table::
   :header-rows: 1

   * - Input
     - Output
   * - ExportTree
     - ExportTree


.. list-table:: Options
   :header-rows: 1

   * - Option
     - Description
     - Default
     - Required
   * - :code:`min_usage`
     - Minimum number of times a given subexpression has to be used in order to be eligible for subexpression elimination.
     - 2
     - No
   * - :code:`merge_inputs`
     - Whether multiple inputs shall be merged together in order to undergo combined rather than individual CSE.
     - false
     - No



density_fitting
^^^^^^^^^^^^^^^

Inserts the density-fitting decomposition of the two-electron integrals. Every occurrence of the tensor named by :code:`integral_label` is rewritten
into a contraction of two three-index tensors (named by :code:`df_tensor_label`) carrying an additional auxiliary index from :code:`auxiliary_space`.
Expressions that don't contain the targeted tensor are passed through unmodified.

.. list-table::
   :header-rows: 1

   * - Input
     - Output
   * - Expression
     - Expression


.. list-table:: Options
   :header-rows: 1

   * - Option
     - Description
     - Default
     - Required
   * - :code:`auxiliary_space`
     - Label of the index space (as declared under the top-level :code:`index_spaces`) used for the auxiliary index introduced by the decomposition.
     - -
     - Yes
   * - :code:`integral_label`
     - Label of the tensor representing the two-electron integrals that shall be decomposed.
     - :code:`g`
     - No
   * - :code:`df_tensor_label`
     - Label given to the two three-index tensors produced by the decomposition.
     - :code:`DF`
     - No


export
^^^^^^

Exports the given expressions as code. Unlike most other steps, :code:`export` is not restricted to a single input - if several :code:`inputs` are
given, they are merged and jointly ordered before being exported together.

.. list-table::
   :header-rows: 1

   * - Input
     - Output
   * - ExportTree
     - None


.. list-table:: Options
   :header-rows: 1

   * - Option
     - Description
     - Default
     - Required
   * - :code:`language`
     - Target language/format to export to. Currently, the only supported value is :code:`itf`.
     - -
     - Yes
   * - :code:`output`
     - Path of the file the generated code is written to (relative to the driver file's location, unless given as an absolute path).
     - -
     - Yes
   * - :code:`optimize`
     - Whether to run a further optimization pass (eliminating redundant operations) over the intermediate representation before generating code.
     - true
     - No
   * - :code:`grouping`
     - Object mapping a group name to an input id/alias. All entries belonging to that id/alias are exported as part of the named group. Any entry not
       matched by any of these rules is put into an implicit group named :code:`Default`.
     - :code:`{}`
     - No
   * - :code:`relative_order`
     - Array of id/alias names giving a hint for the relative ordering of results within the generated code. Used together with the hard dependencies
       among results when determining the final order.
     - :code:`[]`
     - No
   * - :code:`imports`
     - Object mapping a serialized tensor expression to an import name. Tensor blocks matching the given expression are treated as externally supplied
       (imported) rather than being computed by the generated code.
     - :code:`{}`
     - No
   * - :code:`meta`
     - Additional, language-specific metadata. See below for the sub-schema accepted when :code:`language` is :code:`itf`.
     - :code:`{}`
     - No

For :code:`language: itf`, :code:`meta` accepts the following properties:

.. list-table::
   :header-rows: 1

   * - Option
     - Description
     - Default
     - Required
   * - :code:`index_spaces`
     - Object mapping an index space's label to an object with :code:`name` and :code:`tag` properties, overriding the display name and tag used for
       that space in the generated ITF code. Spaces not listed here (or fields left unset for a listed space) fall back to built-in defaults (which
       likely results in an error).
     - :code:`{}`
     - No
   * - :code:`min_index_id`
     - Non-negative integer offset that index IDs emitted into the generated ITF code start counting from.
     - 0
     - No


filter
^^^^^^

Filters the input expressions and assigns them into different named groups based on provided filter rules. For a summed expression, each summand is
tested against every group's rules independently; for a non-summed expression, the expression as a whole is tested. Matching summands/expressions are
accumulated into their group's output, so this step can (and typically does) produce more than one output per input - one for each configured group.

.. list-table::
   :header-rows: 1

   * - Input
     - Output
   * - Expression
     - Expression


.. list-table:: Options
   :header-rows: 1

   * - Option
     - Description
     - Default
     - Required
   * - :code:`groups`
     - Object mapping a group name to a filter specification (see below). At least one group must be given.
     - -
     - Yes
   * - :code:`keep_empty`
     - Whether a group that ends up empty for a given input is still emitted as an (empty) output, rather than being omitted entirely.
     - true
     - No

Each entry of :code:`groups` is a filter specification with the following properties:

.. list-table::
   :header-rows: 1

   * - Option
     - Description
     - Default
     - Required
   * - :code:`mode`
     - Whether a term has to satisfy :code:`all` of the group's rules (logical AND) or just :code:`any` of them (logical OR) in order to be considered
       a match.
     - :code:`all`
     - No
   * - :code:`rules`
     - Array of filter rules (see below) that define what this group matches.
     - -
     - Yes

Each entry of :code:`rules` is an object with the following properties:

.. list-table::
   :header-rows: 1

   * - Option
     - Description
     - Default
     - Required
   * - :code:`type`
     - Kind of rule. Currently, the only supported value is :code:`contains`, which matches if the term (or any of its subexpressions) contains a
       match for :code:`expr` or :code:`label`.
     - -
     - Yes
   * - :code:`negate`
     - If true, inverts the result of this rule.
     - false
     - No
   * - :code:`expr`
     - Serialized expression that is searched for (recursively, i.e. also within subexpressions) within the term being tested. Mutually exclusive with
       :code:`label`.
     - -
     - One of :code:`expr`/:code:`label`
   * - :code:`tensor_equality_mode`
     - Only used together with :code:`expr`. Controls how tensors are compared while matching: :code:`identity` requires exact index labels to match,
       whereas :code:`block`/:code:`shape` (synonyms) only compare index spaces, ignoring concrete index labels.
     - :code:`identity`
     - No
   * - :code:`label`
     - One (String) or several (Array of Strings) full-match regular expressions that are matched (recursively) against the labels of subexpressions
       within the term being tested. A match against any of the given patterns counts as a match for this rule. Mutually exclusive with :code:`expr`.
       Be aware that in order to enter a backslash into the regular expression, you have to escape it in a JSON string as :code:`\\`.
     - -
     - One of :code:`expr`/:code:`label`


optimize
^^^^^^^^

Symbolically rewrites the expressions for an improved numerical evaluation. Most prominently, this factors tensor contractions into a series of binary
tensor contractions, chosen so as to minimize the configured :code:`objective`.

.. list-table::
   :header-rows: 1

   * - Input
     - Output
   * - Expression
     - Expression


.. list-table:: Options
   :header-rows: 1

   * - Option
     - Description
     - Default
     - Required
   * - :code:`objective`
     - Cost metric to minimize when factorizing contractions. One of :code:`DenseFLOPs`, :code:`DenseSize`, :code:`DensePeakSize` or
       :code:`DensePeakSizeBatched`.
     - :code:`DenseFLOPs`
     - No
   * - :code:`reorder_sums`
     - Whether summands of a sum may be reordered so that terms sharing common intermediates end up next to each other.
     - true
     - No
   * - :code:`cse`
     - Whether to perform common-subexpression elimination while searching for a single term's evaluation order. :code:`subnet` recognizes equivalent
       subnetworks (more search effort, potentially fewer operations); :code:`none` disables this. Not to be confused with the standalone :code:`cse`
       step, which operates across a whole processing tree rather than within a single term.
     - :code:`none`
     - No
   * - :code:`intermediate_size_penalty`
     - Per-intermediate memory-footprint penalty added to the cost of a contraction (only consulted when :code:`objective` is :code:`DenseFLOPs`).
     - 0.0
     - No
   * - :code:`prune_outer_products`
     - Whether to prune disconnected (outer-product) subsets from the search space while looking for the best evaluation order.
     - true
     - No


output
^^^^^^

Outputs expressions in the chosen markup style. Mainly intended for debugging purposes.

.. list-table::
   :header-rows: 1

   * - Input
     - Output
   * - Expression
     - None


.. list-table:: Options
   :header-rows: 1

   * - Option
     - Description
     - Default
     - Required
   * - :code:`format`
     - Markup style to print in. :code:`latex` prints a LaTeX representation; :code:`serialize` prints SeQuant's own textual serialization format (see
       :ref:`io-Serialization`).
     - :code:`serialize`
     - No
   * - :code:`annotate_symmetry`
     - Only relevant for :code:`format: serialize`. Whether to include symmetry annotations in the printed representation.
     - true
     - No


project
^^^^^^^

Performs the chosen projection with the inputs.

.. list-table::
   :header-rows: 1

   * - Input
     - Output
   * - Expression
     - Expression


.. list-table:: Options
   :header-rows: 1

   * - Option
     - Description
     - Default
     - Required
   * - :code:`method`
     - Projection method to apply. Currently, the only supported value is :code:`biorthogonal`, which transforms the result into a biorthogonal basis.
     - -
     - Yes

.. note::
   Unlike in the :ref:`v1 driver <external_interface_v1>`, there is no no-op (:code:`primitive`) value for :code:`method` - simply omit this step if no
   projection is needed.


read_input
^^^^^^^^^^

Reads and parses expressions from files.

.. list-table::
   :header-rows: 1

   * - Input
     - Output
   * - None
     - Expression


This step does not accept any options.


simplify
^^^^^^^^

Simplifies the given expressions.

.. list-table::
   :header-rows: 1

   * - Input
     - Output
   * - Expression
     - Expression


This step does not accept any options.


spintracing
^^^^^^^^^^^

Spintraces the given expressions. That is, it performs spin-integration and potentially also spin-summation.

.. list-table::
   :header-rows: 1

   * - Input
     - Output
   * - Expression
     - Expression


.. list-table:: Options
   :header-rows: 1

   * - Option
     - Description
     - Default
     - Required
   * - :code:`algorithm`
     - Spintracing algorithm to use. :code:`rigorous` applies an algorithm that works for all cases; :code:`closed_shell` applies a cheaper algorithm
       that assumes a closed-shell reference.
     - :code:`rigorous`
     - No

.. note::
   There is no :code:`none` value for :code:`algorithm` - simply omit this step if spintracing isn't needed.


substitute
^^^^^^^^^^

Makes substitutions in the given expressions. Substitution targets and their replacements are given directly as serialized expressions in the options
(there is no mechanism to source them from another step's output); they are applied in the order in which they are declared, with each rule operating
on the result of the previous one.

.. list-table::
   :header-rows: 1

   * - Input
     - Output
   * - Expression
     - Expression


.. list-table:: Options
   :header-rows: 1

   * - Option
     - Description
     - Default
     - Required
   * - :code:`substitutions`
     - Object mapping a serialized expression identifying what to match (the key) to a serialized expression it shall be replaced with (the value).
       Multiple entries are applied in order.
     - -
     - Yes
   * - :code:`tensor_equality_mode`
     - Controls how tensors are compared while matching substitution targets: :code:`identity` requires exact index labels to match, whereas
       :code:`block`/:code:`shape` (synonyms) only compare index spaces, ignoring concrete index labels.
     - :code:`block`
     - No
   * - :code:`result_relabeling`
     - Object mapping an old result label to a new one. After substitutions have been applied to an expression, if its result label matches a key in
       this map, it is renamed to the associated value.
     - :code:`{}`
     - No

.. note::
   Unlike most other steps, :code:`substitute` does not propagate the named aliases of its inputs onto its outputs - downstream steps need to reference
   its outputs via the auto-generated :code:`<step_id>.<index>` scheme (or via this step's own :code:`outputs` property). Additionally, if none of the
   expressions belonging to a given input end up changed by any substitution, that input is dropped from the output entirely rather than being passed
   through unchanged.


to_export_tree
^^^^^^^^^^^^^^

Converts the given expressions into a tree data structure suitable for exports.

.. list-table::
   :header-rows: 1

   * - Input
     - Output
   * - Expression
     - ExportTree


This step does not accept any options.


validate
^^^^^^^^

Validates the given expressions, e.g. checking that indices occur at most twice within a product, that all summands of a sum share the same external
indices, and that a result's declared indices are consistent with those actually occurring in its expression. Processing is aborted as soon as the
first invalid expression is encountered.

.. list-table::
   :header-rows: 1

   * - Input
     - Output
   * - Expression
     - None


This step does not accept any options.

.. note::
   Since this step never produces any output, it is typically used without an explicit :code:`id`.



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
     ]
   }

This is an abridged excerpt. The full, runnable driver file can be found at
:file:`utilities/external-interface/examples/nevpt2_v2.json`; it additionally demonstrates the :code:`batch_indices` step and merging several inputs
into a single :code:`export` step.
