Optimizing Evaluation Cost
============================

A symbolic tensor expression as built up from :doc:`Products and Sums <expressions>` says *what* to compute, but not *how*: a flat
:class:`sequant::Product` of several tensors leaves the pairwise contraction order unspecified, and a :class:`sequant::Sum` of many
terms leaves their evaluation order unspecified too. Both choices can change the number of floating-point operations (and the size of
intermediates) by orders of magnitude for the large tensor contractions typical of coupled-cluster and other many-body methods.
:func:`sequant::optimize` chooses good orderings for both, turning a symbolic expression into one that is also efficient to evaluate
or translate into code (see :doc:`export`).

What it does
--------------

Given an expression (or a :doc:`ResultExpr <expressions>`), :func:`sequant::optimize`:

- picks a pairwise contraction order for every :class:`sequant::Product`, minimizing a cost metric (the total floating-point
  operation count, by default) using :class:`sequant::IndexSpace`'s approximate size for each index's extent, and
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

:func:`sequant::optimize` takes an ``OptimizeOptions`` struct exposing further, more advanced controls: alternative cost metrics
(e.g. minimizing intermediate storage or peak memory instead of raw flop count) and common-subexpression elimination across the whole
sum. One such alternative cost metric — minimizing *peak memory* by slicing a large mode into blocks — is substantial enough to have
its own page: see :doc:`batching`. For everything else, see the API reference for :class:`sequant::OptimizeOptions` for the full,
current set of options.
