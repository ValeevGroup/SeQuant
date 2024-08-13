#######################################
SeQuant: Symbolic Tensor Algebra in C++
#######################################

.. image:: https://img.shields.io/badge/License-BSD_3--Clause-blue.svg
    :target: https://opensource.org/licenses/BSD-3-Clause

.. image:: https://github.com/ValeevGroup/SeQuant/actions/workflows/cmake.yml/badge.svg
    :target: https://github.com/ajay-mk/SeQuant/actions/workflows/cmake.yml

.. image:: https://github.com/ValeevGroup/SeQuant/actions/workflows/docs.yml/badge.svg
    :target: https://ajay-mk.github.io/SeQuant

.. note::

    This documentation is under development and may not be complete yet.


Synopsis
========
SeQuant is a framework for performing symbolic algebra of tensors over scalar fields (regular tensors) and over
operator fields (tensor operators in, e.g., quantum many-body physics).
In addition to symbolic manipulation it can numerically evaluate
(with an appropriate external tensor backend) general
tensor algebra expressions.

Computer algebra systems (CAS) like SeQuant are typically implemented within generic CAS like Mathematica or Maple, or
using high-level languages like Python. In fact, version 1 of SeQuant was written in Mathematica. However, the
performance of high-level languages not sufficient for practical use cases.
SeQuant is written in C++ and is designed to be as efficient as possible without loss of generality.

For getting started with SeQuant, refer to :doc:`source/using`. 

Installation
============
**The short version**:

- configure (from top SeQuant source dfirectory): :code:`cmake -DCMAKE_INSTALL_PREFIX=/path/to/install`

- build and test: :code:`cmake --build build --target install`

See :doc:`source/install` for detailed installation instructions.

Build Harness
=============
We will only consider how to use SeQuant from within an existing codebase that has a `CMake <https://cmake.org>`__ harness. 
If SeQuant has already been built and installed, it should be sufficient to add the SeQuant install prefix to :code:`CMAKE_INSTALL_PREFIX` CMake cache variable and adding the following lines to the :code:`CMakeLists.txt` file in your codebase:

.. code-block:: cmake

    find_package(SeQuant CONFIG)
    if (NOT TARGET SeQuant::SeQuant)
        cmake_minimum_required(VERSION 3.14.0)  # for FetchContent_MakeAvailable
        include(FetchContent)
        FetchContent_Declare(sequant
                GIT_REPOSITORY      https://github.com/ValeevGroup/SeQuant2.git
                )
        FetchContent_MakeAvailable(sequant)
        add_library(SeQuant::SeQuant ALIAS SeQuant)
    endif()
    target_link_libraries(your_executable_or_library_target PUBLIC SeQuant::SeQuant)

Developers
==========
SeQuant is developed by the `Valeev Research Group <https://valeevgroup.github.io>`__ in the Department of Chemistry at Virginia Tech.

Acknowledgement
===============
Development of SeQuant has been possible thanks to the support of the US National Science Foundation (award 2217081) and the US Department of Energy (awards DE-SC0022327 and DE-SC0022263)

.. toctree::
    :caption: Getting Started
    :hidden:
    :maxdepth: 2

    source/install
    source/using


.. toctree::
    :caption: Advanced Topics
    :maxdepth: 2
    :hidden:

    source/cc
    source/eval

.. toctree::
    :caption: API Reference
    :maxdepth: 1
    :hidden:

    genindex

.. toctree::
    :caption: See Also
    :hidden:

    TiledArray <https://valeevgroup.github.io/tiledarray>
    Massively Parallel Quantum Chemistry <https://valeevgroup.github.io/mpqc>
    Valeev Research Group <https://valeevgroup.github.io>
