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

A fully canonicalized expression remembers that it is in canonical form, so canonicalizing or simplifying it again before it is
modified does no work.

.. _canonicalization-internal-form:

One canonical form
--------------------

The library canonicalizes to one fixed form: the ``Complete`` method (:class:`sequant::CanonicalizationMethod`), with the library's
default tensor canonicalizers, index comparers and cardinal tensor labels (:func:`sequant::default_cardinal_tensor_labels`). No global
setting changes it; it depends on the :class:`sequant::Context` only through the index basis registry (the order of the index spaces)
and the single-particle basis (the symmetry of normal operators), see :doc:`context`. Code throughout the library relies on this
form, so it is the same whoever calls it and whatever else is configured.

Named indices
----------------

Canonicalization renames dummy indices, but must never rename the *named* (external) indices, whose labels carry meaning. Which
indices are named is therefore always specified explicitly; there is no default:

- a :class:`sequant::ResultExpr` carries them as the indices of its result, and ``canonicalize(ResultExpr&)`` /
  ``simplify(ResultExpr&)`` use them. These overloads also merge terms that differ only by a permutation of the result indices
  that the result's symmetry allows, provided an explicit (anti)symmetrizer of the result (``Ŝ``/``Â``) multiplies them: that
  the result has a symmetry does not make its terms have it (e.g. :math:`R = X + P(X)` is symmetric, :math:`2X` is not);
- for a plain ``ExprPtr`` the caller passes :data:`sequant::deduce_named_indices`, opting in to deducing the named indices as
  the indices that occur once in each product. This is wrong for theories in which an external index can occur several times
  (e.g. tensor hypercontraction), which is why it is never done implicitly. The library itself never deduces named indices.

The same holds for :class:`sequant::WickTheorem` and ``vac_av``/``ref_av``, and for the Python bindings
(``simplify(expr, deduce=True)``):

.. literalinclude:: /examples/user/canonicalization.cpp
   :language: cpp
   :start-after: start-snippet-3
   :end-before: end-snippet-3
   :dedent: 2

Nested products
-----------------

A product whose factors are tensors and (recursively) products of tensors, such as a product parenthesized for evaluation by
``optimize()``, is canonicalized as one tensor network of all its tensors: the dummy indices are named consistently across the
nested products, and the parenthesization is preserved; only the factors at each level are reordered.

.. _canonicalization-presentation:

Presentation
--------------

The canonical form is meant to identify equivalent expressions, not to be read. Other forms are requested explicitly:

- :func:`sequant::prettify` returns a copy of an expression in a form meant for reading: optionally canonicalized with the given
  presentation settings (:class:`sequant::CanonicalizationPresentation`: tensor canonicalizers per label, index comparers and
  cardinal tensor labels), then sorted lexicographically, which, unlike canonicalization, does not depend on the hashing algorithm.
  Named indices are never renamed, and the structure of the expression is preserved: products keep their parenthesization
  (e.g. the order of evaluation chosen by ``optimize()``), and sums are neither flattened nor are their summands combined.
  For a :class:`sequant::ResultExpr`, the terms of a sum that a symmetrizer of the result multiplies are presented alike if
  they are related by the symmetry of the result.
  :func:`sequant::mbpt::prettify_options` is a preset that presents mbpt expressions in their conventional order (the cardinal
  tensor labels of mbpt); ``prettify`` itself is part of the core library;
- the same presentation settings can be passed to ``canonicalize()`` and the ``Expr::canonicalize()`` members
  (``CanonicalizeOptions::presentation``). The result is not marked as canonical, and ``simplify()``, whose result must be the
  canonical form, rejects them;
- the LaTeX printers can print summands and factors in a lexicographic order (:class:`sequant::PrintOptions`, see :doc:`context`).

.. literalinclude:: /examples/user/canonicalization.cpp
   :language: cpp
   :start-after: start-snippet-4
   :end-before: end-snippet-4
   :dedent: 2

None of this needs to be invoked explicitly in typical use: the :doc:`mbpt operator machinery <operator>` and
:doc:`CC equation generator <cc>` call ``simplify()``/canonicalization as needed while building up equations. Knowing that it happens —
and why two "different-looking" terms may in fact be the same one — is mainly useful for interpreting intermediate output and for
understanding why changing a tensor's declared symmetry can change how many terms an equation prints as. The implementation details of
the tensor-network canonicalizer are documented separately, for contributors, in :doc:`the developer guide </developer/tnc>`.
