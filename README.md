# SeQuant: Symbolic Tensor Algebra in C++

[![License: LGPL v3](https://img.shields.io/badge/License-LGPLv3-blue.svg)]()
[![Linux/MacOS Build](https://github.com/ValeevGroup/SeQuant/actions/workflows/cmake.yml/badge.svg)](https://github.com/ValeevGroup/SeQuant/actions/workflows/cmake.yml)
[![Docs](https://github.com/ValeevGroup/SeQuant/actions/workflows/docs.yml/badge.svg)](https://valeevgroup.github.io/SeQuant)


## Synopsis

SeQuant is a framework for performing symbolic algebra of tensors over scalar and
operator rings.  In addition to symbolic manipulation it can numerically evaluate
(with an appropriate external tensor backend) general  tensor algebra expressions.

Computer algebra systems (CAS) like SeQuant are typically implemented within generic CAS like Mathematica or Maple, or
using high-level languages like Python. In fact, version 1 of SeQuant was written in Mathematica. However, the
performance of high-level languages not sufficient for practical use cases.
SeQuant is written in C++ and is designed to be as efficient as possible without loss of generality.

See detailed documentation at [https://valeevgroup.github.io/SeQuant/](https://valeevgroup.github.io/SeQuant/).

## Installation

The short version:

- configure (from top SeQuant source dfirectory): `cmake -B build -S . -DCMAKE_INSTALL_PREFIX=/path/to/where/sequant/to/be/installed`
- build and install: `cmake --build build --target install`

For detailed instructions see [SeQuant: Installation Guide](https://valeevgroup.github.io/SeQuant/user/getting_started/installing.html).


### Build harness
We will only consider how to use SeQuant from within an existing codebase that has a [CMake](https://cmake.org) harness. If SeQuant has already been built and installed, it should be sufficient to add the SeQuant install prefix to `CMAKE_INSTALL_PREFIX` CMake cache variable and adding the following lines to the `CMakeLists.txt` file in your codebase:

```cmake
find_package(SeQuant CONFIG REQUIRED)
target_link_libraries(your_executable_or_library_target PUBLIC SeQuant::SeQuant)
```

It is often desirable to build SeQuant from source within a standalone codebase; this case be done using the [FetchContent CMake module](https://cmake.org/cmake/help/latest/module/FetchContent.html) as follows:

```cmake
find_package(SeQuant CONFIG)
if (NOT TARGET SeQuant::SeQuant)
    cmake_minimum_required(VERSION 3.14.0)  # for FetchContent_MakeAvailable
    include(FetchContent)
    FetchContent_Declare(sequant
            GIT_REPOSITORY      https://github.com/ValeevGroup/SeQuant.git
            )
    FetchContent_MakeAvailable(sequant)
    add_library(SeQuant::SeQuant ALIAS SeQuant)
endif()
target_link_libraries(your_executable_or_library_target PUBLIC SeQuant::SeQuant)
```

## Developers

The [Valeev Research Group (VRG)](https://valeevgroup.github.io) in the Department of Chemistry at Virginia Tech kickstarted the initial design and development of SeQuant. The ongoing development of SeQuant is driven by major [contributions](https://github.com/ValeevGroup/SeQuant/graphs/contributors) from VRG and the [Köhn Group](https://www.itheoc.uni-stuttgart.de/research/koehn) in the Department of Theoretical Chemistry at University of Stuttgart. 

## Acknowledgement

Development of SeQuant has been possible thanks to the support of the US National Science Foundation (award 2217081) and the US Department of Energy (awards DE-SC0022327 and DE-SC0022263)
