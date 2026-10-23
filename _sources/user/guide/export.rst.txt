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

Each backend has its own ``Context`` type (e.g. ``ItfContext``) for backend-specific settings such as naming conventions; consult its
reference documentation for details. Numerically *evaluating* an expression directly in-process, rather than generating source code
for it, is a separate concern handled by :code:`SeQuant/core/eval` (see :doc:`evaluation` for its trace/profiling tooling).
