.. toctree::
   :caption: User
   :maxdepth: 2
   :hidden:

   Getting started <user/getting_started/index>
   Installation <user/getting_started/installing>
   User Guide <user/guide/index>

.. toctree::
   :caption: Developer
   :hidden:
   :maxdepth: 2

   Developer docs <developer/index>
   API reference <api>

.. toctree::
   :caption: See Also
   :hidden:

   TiledArray <https://valeevgroup.github.io/tiledarray>
   Massively Parallel Quantum Chemistry <https://valeevgroup.github.io/mpqc>
   Valeev Research Group <https://valeevgroup.github.io>


SeQuant: Symbolic Tensor Algebra in C++
=======================================

SeQuant is a framework for performing symbolic algebra of tensors over scalar fields (regular tensors) and over operator fields (tensor operators in,
e.g., quantum many-body physics). In addition to symbolic manipulation it can numerically evaluate (with an appropriate external tensor backend)
general tensor algebra expressions.

Computer algebra systems (CAS) like SeQuant are typically implemented within generic CAS like Mathematica or Maple, or using high-level languages like
Python. In fact, version 1 of SeQuant was written in Mathematica. However, the performance of high-level languages not sufficient for practical use
cases. SeQuant is written in C++ and is designed to be as efficient as possible without loss of generality.

.. grid:: 1 1 2 2
   :gutter: 2

   .. grid-item-card:: Getting started
      :link: user/getting_started/index
      :link-type: doc

      Quick primer on how to get started using SeQuant. This includes installation instructions as well as first examples.


   .. grid-item-card:: User guide
      :link: user/guide/index
      :link-type: doc

      A more in-depth view of the fundamental concepts of SeQuant and an overview over SeQuant's capabilities.


   .. grid-item-card:: Developer guide
      :link: developer/index
      :link-type: doc

      Explanation of some more low-level details and inner workings of SeQuant along with some guidance on the source code structure. Only needed when
      developing SeQuant itself.


   .. grid-item-card:: API reference
      :link: api
      :link-type: doc

      Complete API documentation generated from in-source doc comments

Developers
----------

The `Valeev Research Group (VRG) <https://valeevgroup.github.io>`_ in the Department of Chemistry at Virginia Tech kickstarted the initial design and
development of SeQuant. The ongoing development of SeQuant is driven by major
`contributions <https://github.com/ValeevGroup/SeQuant/graphs/contributors>`_ from VRG and the
`Köhn Group <https://www.itheoc.uni-stuttgart.de/research/koehn>`_ in the Department of Theoretical Chemistry at University of Stuttgart. 

Acknowledgement
---------------

Development of SeQuant has been possible thanks to the support of the US National Science Foundation (award 2217081) and the US Department of Energy
(awards DE-SC0022327 and DE-SC0022263)
