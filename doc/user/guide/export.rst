Code Generation
=================

Once an equation has been derived, :doc:`simplified <canonicalization>`, and :doc:`optimized <optimize>`, the final step is turning it
into code that a numerical tensor library or quantum-chemistry program can actually run. SeQuant's export framework
(:code:`SeQuant/core/export`) does this by walking a binarized evaluation tree and driving a pluggable code-generation backend, so
adding support for a new target language means implementing one small interface rather than writing a new tree-walker.

The ``Generator`` interface
------------------------------

A backend is a subclass of :class:`sequant::Generator`: a set of `callbacks <https://en.wikipedia.org/wiki/Callback_(computer_programming)>`_
(``create``, ``load``, ``compute``, ``unload``, ``declare``, ...) that get invoked, in the right order, while
:func:`sequant::export_expression` traverses a :doc:`ResultExpr <expressions>`'s
evaluation tree. Writing a new backend means implementing this callback interface; the tree-walking, scalar-factor bookkeeping, and
intermediate-reuse logic are handled once, centrally, for every backend — documented for contributors in :doc:`the developer guide
</developer/export>`. :class:`sequant::TextGenerator` is a minimal,
dependency-free backend that renders these callbacks as human-readable pseudocode — useful both as documentation of the callback
sequence and as a debugging aid for a real backend under development:

.. literalinclude:: /examples/user/export.cpp
   :language: cpp
   :start-after: start-snippet-1
   :end-before: end-snippet-1
   :dedent: 2

Available backends
---------------------

Beyond ``TextGenerator``, SeQuant ships several backends targeting real numerical tensor libraries and domain-specific formats:

- :class:`sequant::ItfGenerator` generates `ITF <https://doi.org/10.1002/wcms.82>`_ (Integrated Tensor Framework) code, a
  domain-specific tensor-contraction language used by parts of the Molpro quantum chemistry program.
- :class:`sequant::PythonEinsumGeneratorBase` and its NumPy/PyTorch specializations generate Python code using ``einsum``-style
  tensor contractions.
- :class:`sequant::JuliaTensorOperationsGenerator`, :class:`sequant::JuliaITensorGenerator`, and
  :class:`sequant::JuliaTensorKitGenerator` generate Julia code targeting, respectively, the
  `TensorOperations.jl <https://github.com/Jutho/TensorOperations.jl>`_,
  `ITensors.jl <https://itensor.github.io/ITensors.jl>`_, and `TensorKit.jl <https://jutho.github.io/TensorKit.jl>`_ packages.

Element-wise tensor powers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An evaluation tree may contain the element-wise power of a tensor (a :class:`sequant::TensorPower`, e.g. an orbital-energy
denominator raised to ``-1``). The text, Python (NumPy and PyTorch ``einsum``) and Julia TensorOperations backends export it,
spelling the exponent as follows; the operand is permuted to the result's index order where the two differ, and a scalar
prefactor, if any, is written in front of the power:

.. code-block:: text

   Compute R[i_1, a_1] += 2 D[a_1, i_1]^(-1)                                  (text)
   R_ov += 2 * np.einsum('ai->ia', D_vo, optimize=True)**(-1)                 (Python)
   R_ov .+= 2 * permutedims(D_vo, (2, 1)) .^ (-1)                             (Julia TensorOperations)

Julia spells a rational exponent as, e.g., ``(1//2)``. :class:`sequant::ItfGenerator` ("ITF has no element-wise tensor power"),
:class:`sequant::JuliaITensorGenerator` and :class:`sequant::JuliaTensorKitGenerator` ("<format> cannot express an element-wise
tensor power") have no way to express it and throw a :class:`sequant::Exception` instead. A custom backend opts in by overriding
:func:`sequant::Generator::represent` for :class:`sequant::TensorPower`; the default throws.

Broadcast sums
~~~~~~~~~~~~~~~~~

The evaluation tree of a broadcast sum under a ``ResultExpr`` head (see :doc:`expressions`) contains
:cpp:enumerator:`sequant::EvalOp::Broadcast` nodes, which no shipped backend exports yet: the export throws a
:class:`sequant::Exception`.

Each backend has its own ``Context`` type (e.g. ``ItfContext``) for backend-specific settings such as naming conventions; consult its
reference documentation for details. Numerically *evaluating* an expression directly in-process, rather than generating source code
for it, is a separate concern handled by :code:`SeQuant/core/eval` (see :doc:`evaluation` for its trace/profiling tooling).
