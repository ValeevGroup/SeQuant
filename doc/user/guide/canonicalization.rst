Canonicalization
==================

Symbolic tensor expressions in many-body physics are riddled with *dummy* (summed-over) indices, and the same physical term can be
written down with many different, equally valid choices of dummy index labels. Two :class:`sequant::Expr` trees that are mathematically
equal are therefore not necessarily represented by *identical* trees — until they have been put into **canonical form**. Canonicalization
is what lets SeQuant recognize that two differently-labeled products are the same term (so they can be combined in a ``Sum``), and it
underlies essentially every simplification performed throughout the library, from :doc:`Wick's theorem </user/getting_started/wick>` to
:doc:`the coupled-cluster equation generator <cc>`.

Why it matters
----------------

Consider a term such as :math:`\bar{g}^{a_1 a_2}_{i_1 i_2} t^{i_1 i_2}_{a_1 a_2}`. Relabeling the summed indices
(:math:`i_1 \leftrightarrow i_2`, :math:`a_1 \leftrightarrow a_2`) yields an expression that is numerically identical, but whose
:class:`sequant::Product` tree differs in the index labels attached to its :class:`sequant::Tensor` factors. Without canonicalization,
adding this relabeled copy to the original would produce a two-term :class:`sequant::Sum` instead of a single term with a factor of 2 —
and, more importantly, larger derivations (e.g. the hundreds of terms arising in a CCSDT residual) would never collapse to their true,
much smaller size.

.. literalinclude:: /examples/user/canonicalization.cpp
   :language: cpp
   :start-after: start-snippet-1
   :end-before: end-snippet-1
   :dedent: 2

How it works, briefly
------------------------

Canonicalizing a single :class:`sequant::Tensor` — putting its own bra/ket indices in a fixed order consistent with its declared
permutational symmetry — is handled by a :class:`sequant::TensorCanonicalizer`. Canonicalizing a whole product of tensors (or of
normal-ordered operators) additionally requires choosing a consistent relabeling of the *dummy* indices shared between factors; SeQuant
does this by building a colored graph representation of the product (a *tensor network*) and computing its canonical form using the
bundled `bliss <https://users.aalto.fi/~tjunttil/bliss/>`_ `graph-automorphism <https://en.wikipedia.org/wiki/Graph_automorphism>`_ library. This machinery is what
``Expr::canonicalize()`` invokes internally, and what :func:`sequant::simplify` combines with cheap algebraic clean-up (flattening,
dropping zeros, ...) to fully reduce an expression:

.. literalinclude:: /examples/user/canonicalization.cpp
   :language: cpp
   :start-after: start-snippet-2
   :end-before: end-snippet-2
   :dedent: 2

None of this needs to be invoked explicitly in typical use: the :doc:`mbpt operator machinery <operator>` and
:doc:`CC equation generator <cc>` call ``simplify()``/canonicalization as needed while building up equations. Knowing that it happens —
and why two "different-looking" terms may in fact be the same one — is mainly useful for interpreting intermediate output and for
understanding why changing a tensor's declared symmetry can change how many terms an equation prints as. The implementation details of
the tensor-network canonicalizer are documented separately, for contributors, in :doc:`the developer guide </developer/tnc>`.
