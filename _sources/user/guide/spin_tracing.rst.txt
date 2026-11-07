Spin Tracing
==============

Equations derived by SeQuant's :doc:`Wick's theorem machinery </user/getting_started/wick>` and the :doc:`CC equation generator <cc>`
are, by default, expressed in a **spin-orbital** basis: every index implicitly ranges over spin-up (:math:`\alpha`) and spin-down
(:math:`\beta`) single-particle states alike. For a spin-unrestricted (open-shell) reference that is exactly what is needed. For a
spin-restricted, closed-shell reference, however, the spin degrees of freedom can be summed out analytically, yielding a smaller set
of equations over spatial orbitals alone — a transformation commonly called **spin tracing** or **spin adaptation**. This is a standard
step in deriving efficient closed-shell coupled-cluster and many-body perturbation theory equations for quantum chemistry.

The ``mbpt::Spin`` quantum number
------------------------------------

Spin is represented as an ordinary :doc:`index space quantum number <context>`, :class:`sequant::mbpt::Spin` (``alpha``, ``beta``, or
``any``/``none`` for "unspecified/spin-free"), attached to an :class:`sequant::Index` through its :class:`sequant::IndexSpace`. Helper
functions such as :func:`sequant::mbpt::make_spinalpha`, :func:`sequant::mbpt::make_spinbeta`, and :func:`sequant::mbpt::to_spin`
add, change, or query this quantum number on individual indices:

.. literalinclude:: /examples/user/spin_tracing.cpp
   :language: cpp
   :start-after: start-snippet-2
   :end-before: end-snippet-2
   :dedent: 2

Tracing spin out of an expression
------------------------------------

Rather than manipulating spin labels index-by-index, the usual entry point is :func:`sequant::mbpt::spintrace`: given a spin-orbital
expression, it sums over the spin cases of its internal indices and returns the spin-free (spatial-orbital) result. In the simplest case
— where none of the involved index spaces are already split into separate alpha/beta subspaces — this amounts to summing the (identical)
spin-conserving contributions:

.. literalinclude:: /examples/user/spin_tracing.cpp
   :language: cpp
   :start-after: start-snippet-1
   :end-before: end-snippet-1
   :dedent: 2

Specialized variants exist for the more demanding cases that arise in practice:

- :func:`sequant::mbpt::closed_shell_spintrace` is a more efficient alternative to ``spintrace`` specifically for closed-shell
  (spin-restricted) references, avoiding the exponential cost of the general algorithm.
- :func:`sequant::mbpt::closed_shell_CC_spintrace` additionally transforms the traced result into *biorthogonal* form, the
  representation typically wanted for closed-shell coupled-cluster equations (see :doc:`cc`); it also factors out and re-applies
  the necessary (anti)symmetrizers.
- :func:`sequant::mbpt::open_shell_spintrace` and :func:`sequant::mbpt::open_shell_CC_spintrace` produce, instead of a single spin-free
  result, one expression per distinct spin case, appropriate for an unrestricted (open-shell) reference.

All of these expect their input to already be in a specific normal form (a leading antisymmetrizer, produced by "complete"
canonicalization — see :doc:`canonicalization`, unless the expression doesn't have any external indices (e.g. energy expressions)); consult their
reference documentation for the exact preconditions and options before using them on a new class of equations. For a worked example applying
``closed_shell_CC_spintrace``/``open_shell_CC_spintrace`` to actual coupled-cluster amplitude equations, see :ref:`cc-spin-tracing` in the :doc:`cc`
page.
