Optimizing Evaluation Cost
============================

A symbolic tensor expression as built up from :doc:`Products and Sums <expressions>` says *what* to compute, but not *how*: a flat
:class:`sequant::Product` of several tensors leaves the pairwise contraction order unspecified, and a :class:`sequant::Sum` of many
terms leaves their evaluation order unspecified too. Both choices can change the number of floating-point operations (and the size of
intermediates) by orders of magnitude for the large tensor contractions typical of coupled-cluster and other many-body methods.
:func:`sequant::optimize` chooses good orderings for both, turning a symbolic expression into one that is also efficient to evaluate
or translate into code (see :doc:`export`).

The contraction search is per term. It does not jointly optimize the peak memory of a forest or the benefit of sharing
an intermediate across equations.

What it does
--------------

Given an expression (or a :doc:`ResultExpr <expressions>`), :func:`sequant::optimize`:

- picks a pairwise contraction order for every :class:`sequant::Product`, minimizing a cost metric (the total floating-point
  operation count, by default) using each index's basis extent (:func:`sequant::IndexBasis::extent`), and
- reorders the summands of every :class:`sequant::Sum` so that terms sharing common intermediates end up next to each other, which
  helps downstream common-subexpression elimination recognize them.

.. literalinclude:: /examples/user/optimize.cpp
   :language: cpp
   :start-after: start-snippet-1
   :end-before: end-snippet-1
   :dedent: 2

Here, contracting :code:`A` with :code:`B` first (rather than :code:`B` with :code:`C`) is cheaper given the relative sizes of the
occupied and virtual spaces involved, so :func:`sequant::optimize` groups them into an explicit sub-product before multiplying by
:code:`C` — visible above as the extra parenthesization in the optimized expression's structure.

Tuning
--------

:func:`sequant::optimize` takes an ``OptimizeOptions`` struct exposing further, more advanced controls. Its full field
list is API-reference material (see :class:`sequant::OptimizeOptions`), but a few matter for everyday use:

- ``objective_function`` selects the cost metric. ``DenseSize`` minimizes the total size of all intermediates, which
  can favor a different contraction order from minimizing the peak size of all simultaneously live tensors
  (``DenseSpaceTime``). ``DenseSpaceTime`` trades a small, bounded increase in peak (``peak_flops_tolerance``) for
  lower performance cost; ``DenseTimeSpace`` puts performance first and peak second. The batched variants of these
  two additionally slice a large mode into blocks to bound peak memory under a budget — substantial enough to have
  its own page: see :doc:`batching`.
- ``volatile_weight`` makes operations that must be rebuilt every iteration (those depending on a leaf marked by
  ``BatchPolicy::is_volatile_leaf``) count more than ones whose results can be cached. It is ignored by
  ``DenseSize``; under ``DenseSpaceTime`` it only affects the choice among schedules within ``peak_flops_tolerance``
  of the minimum peak.
- The optional ``roofline`` parameters make the performance cost account for data movement as well as arithmetic;
  see :ref:`cost-model-roofline`.
- Common-subexpression elimination during the contraction search (``CSEOptions::subnet``) is supported only by
  ``DenseFLOPs`` and ``DenseSize``; leave it disabled for the other objectives. Enabling it with another objective
  trips an assertion, and where assertions are compiled out (``SEQUANT_ASSERT_BEHAVIOR=IGNORE``) CSE is silently
  skipped.
