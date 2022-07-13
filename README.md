SeQuant: second quantization toolkit in C++
===========================================

# Synopsis

`SeQuant` is a framework for performing symbolic algebra designed specifically
for the algebra of second quantization in quantum many-body physics.
More abstractly it can symbolically transform and numerically
evaluate (with an appropriate external tensor backend) general
tensor algebra expressions.

# Installation

See file INSTALL.md .

# Getting started

## Build harness
To use SeQuant from within an existing codebase that has a CMake harness (the only case considered here) it should be sufficient to do this:
```cmake
find_package(SeQuant CONFIG REQUIRED)
target_link_libraries(your_executable_or_library_target PUBLIC SeQuant::SeQuant)
```
It is often desirable to build SeQuant from source within a standalone codebase; this case be done using the FetchContent CMake module as follows:
```cmake
find_package(SeQuant CONFIG)
if (NOT TARGET SeQuant::SeQuant)
    cmake_minimum_required(VERSION 3.14.0)  # for FetchContent_MakeAvailable
    include(FetchContent)
    FetchContent_Declare(sequant
            GIT_REPOSITORY      https://github.com/ValeevGroup/SeQuant2.git
            )
    FetchContent_MakeAvailable(sequant)
endif()
target_link_libraries(your_executable_or_library_target PUBLIC SeQuant::SeQuant)
```

## Using

To get started let's use SeQuant to apply Wick's theorem to a product of elementary (creation and annihilation)
operators:

```c++
#include <SeQuant/core/sequant.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/wick.hpp>

int main() {
  using namespace sequant;

  IndexSpace sp;
  Index p1(L"p_1", sp), p2(L"p_2", sp), p3(L"p_3", sp), p4(L"p_4", sp);

  auto cp1 = fcrex(p1), cp2 = fcrex(p2);
  auto ap3 = fannx(p3), ap4 = fannx(p4);

  std::wcout << to_latex(ap3 * ap4 * cp1 * cp2) << " = " << to_latex(FWickTheorem{ap3 * ap4 * cp1 * cp2}.set_external_indices(std::array{p1, p2, p3, p4}).full_contractions(false).compute()) << std::endl;
  
  return 0;
}
```
Running this program should produce ![this expression](doc/images/tut-img1.svg) .

# Developers

`SeQuant` is developed by the Valeev group at Virginia Tech.
