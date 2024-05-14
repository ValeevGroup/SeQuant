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

### Getting Started

To get started let's use SeQuant to apply Wick's theorem to a simple product of elementary (creation and annihilation)
 fermionic operators:

![a_{p_3} a_{p_4} a^\dagger_{p_1} a^\dagger_{p_2}](doc/images/tut-expr1.svg).

This is achieved by the following SeQuant program:

```c++
#include <SeQuant/core/wick.hpp>

int main() {
  using namespace sequant;

  IndexSpace sp;
  Index p1(L"p_1", sp), p2(L"p_2", sp), p3(L"p_3", sp), p4(L"p_4", sp);

  auto cp1 = fcrex(p1), cp2 = fcrex(p2);
  auto ap3 = fannx(p3), ap4 = fannx(p4);

  std::wcout << to_latex(ap3 * ap4 * cp1 * cp2) << " = "
             << to_latex(FWickTheorem{ap3 * ap4 * cp1 * cp2}
                             .set_external_indices(std::array{p1, p2, p3, p4})
                             .full_contractions(false)
                             .compute())
             << std::endl;
  
  return 0;
}
```

N.B. All _core_ user-facing SeQuant code lives in C++ namespace `sequant`, from now on we will assume this namespace has been imported via `using namespace sequant` and will omit the explicit namespace qualification for SeQuant components.

Running this program should produce a LaTeX expression for this formula:

![{{a^{}_{{p_3}}}{a^{}_{{p_4}}}{a^{{p_1}}_{}}{a^{{p_2}}_{}}} = { \bigl( - {{a^{{p_1}{p_2}}_{{p_3}{p_4}}}} - {{s^{{p_1}}_{{p_3}}}{s^{{p_2}}_{{p_4}}}} + {{s^{{p_1}}_{{p_3}}}{a^{{p_2}}_{{p_4}}}} - {{s^{{p_1}}_{{p_4}}}{a^{{p_2}}_{{p_3}}}} + {{s^{{p_2}}_{{p_3}}}{s^{{p_1}}_{{p_4}}}} - {{s^{{p_2}}_{{p_3}}}{a^{{p_1}}_{{p_4}}}} + {{s^{{p_2}}_{{p_4}}}{a^{{p_1}}_{{p_3}}}}\bigr) }](doc/images/tut-expr1-result1.svg)

, where the tensor notation is used to denote elementary and composite _normal_ operators,

![a^p \equiv \, & a_p^\dagger \\ a^{p_1 p_2 \dots p_c}_{q_1 q_2 \dots q_a} \equiv \, & a_{p_1}^\dagger  a_{p_2}^\dagger \dots a_{p_c}^\dagger a_{q_a} \dots a_{q_2} a_{q_1}](doc/images/tut-notation-eq1.svg)

, and

![s^p_q \equiv \langle q | p \rangle](doc/images/tut-notation-eq2.svg)

denotes 1-particle state inner products (overlaps). Wick's theorem can of course be applied to directly to products of normal composite operators, e.g,

```c++
  auto nop1 = ex<FNOperator>(std::array{p1, p2}, std::array{p3, p4});
  auto nop2 = ex<FNOperator>(std::array{p5}, std::array{p6, p7});

  std::wcout << to_latex(nop1 * nop2) << " = "
             << to_latex(FWickTheorem{nop1 * nop2}
                             .set_external_indices(
                                 std::array{p1, p2, p3, p4, p5, p6, p7})
                             .full_contractions(false)
                             .compute())
             << std::endl;
```

produces

![{{a^{{p_1}{p_2}}_{{p_3}{p_4}}}{a^{\textvisiblespace\,{p_5}}_{{p_6}{p_7}}}} = { \bigl({a^{\textvisiblespace\,{p_1}{p_2}{p_5}}_{{p_3}{p_4}{p_6}{p_7}}} - {{s^{{p_5}}_{{p_4}}}{a^{\textvisiblespace\,{p_1}{p_2}}_{{p_3}{p_6}{p_7}}}} + {{s^{{p_5}}_{{p_3}}}{a^{\textvisiblespace\,{p_1}{p_2}}_{{p_4}{p_6}{p_7}}}}\bigr) }](doc/images/tut-expr2-result1.svg)

, where `␣` is used in number-nonconserving operators to point out the empty "slots".

Of course, same manipulations can be used for bosons:

![{{b^{{p_1}{p_2}}_{{p_3}{p_4}}}{b^{{p_5}{p_6}}_{\textvisiblespace\,{p_7}}}} = { \bigl({b^{{p_5}{p_1}{p_2}{p_6}}_{\textvisiblespace\,{p_3}{p_4}{p_7}}} + {{s^{{p_6}}_{{p_3}}}{b^{{p_5}{p_2}{p_1}}_{\textvisiblespace\,{p_4}{p_7}}}} + {{s^{{p_5}}_{{p_3}}}{b^{{p_1}{p_2}{p_6}}_{\textvisiblespace\,{p_4}{p_7}}}} + {{s^{{p_6}}_{{p_4}}}{b^{{p_5}{p_1}{p_2}}_{\textvisiblespace\,{p_3}{p_7}}}} + {{s^{{p_5}}_{{p_4}}}{b^{{p_2}{p_1}{p_6}}_{\textvisiblespace\,{p_3}{p_7}}}} + {{s^{{p_5}}_{{p_3}}}{s^{{p_6}}_{{p_4}}}{b^{{p_1}{p_2}}_{\textvisiblespace\,{p_7}}}} + {{s^{{p_6}}_{{p_3}}}{s^{{p_5}}_{{p_4}}}{b^{{p_2}{p_1}}_{\textvisiblespace\,{p_7}}}}\bigr) }](doc/images/tut-expr3-result1.svg)

, where `b` denotes normal bosonic operators constructed analogously with the normal fermionic operators `a`, is obtained via

```c++
  auto nop3 = ex<BNOperator>(std::array{p1, p2}, std::array{p3, p4});
  auto nop4 = ex<BNOperator>(std::array{p5, p6}, std::array{p7});

  std::wcout << to_latex(nop3 * nop4) << " = "
             << to_latex(BWickTheorem{nop3 * nop4}
                             .set_external_indices(
                                 std::array{p1, p2, p3, p4, p5, p6, p7})
                             .full_contractions(false)
                             .compute())
             << std::endl;
```
### Register your relevant spaces
In scientific writing, it is common to index mathematical objects such as tensors and operators. These indices in a network of tensors encode a graph structure which can be used to deduce which modes of each tensor can be contracted. It is also common for these indices to indicate a block structure of a tensor in a subspace.
![](doc/images/fia.svg)
Most chemists would recognize the labeling of this tensor as the off-diagonal occupied, unoccupied block of the fock matrix. To understand this however, requires knowledge of what `i` and `a` mean. In other words, the labels `i` and `a` only have meaning in a context.
SeQuant now allows/requires the user to define the base labels and the set theoretics for the context they are in.

```c++
#include <SeQuant/core/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

int main() {
  using namespace sequant;
      IndexSpaceRegistry sample_ISR;

  IndexSpace x_space(L"x", 0b001);
  sample_ISR.add(x_space);

  IndexSpace y_space(L"y", 0b010);
  sample_ISR.add(y_space);

  IndexSpace z_space(L"z", 0b100);
  sample_ISR.add(z_space);

  IndexSpace xy_space(L"xy", 0b011);
  sample_ISR.add(xy_space);

  IndexSpace yz_space(L"yz", 0b110);
  sample_ISR.add(yz_space);

  IndexSpace xyz_space(L"xyz", 0b111);
  sample_ISR.add(xyz_space);

  assert(sample_ISR.unIon(z_space,y_space) == yz_space);
  assert(sample_ISR.intersection(xy_space,x_space) == x_space);
  
  Context new_cxt(Vacuum::SingleProduct,sample_ISR);
  set_default_context(new_cxt);
}
```
The above sample registry outlines the set theoreitcs for  indices x, y, and z using unsigned ints conveniently written as bitsets.
Notice that not all possible combinations need to be registered. It is the task of the user to predict any and all spaces that they may 
encounter in their context. Notice that `unIon()` and `intersection()` are part of the `IndexSpaceRegistry` class. This is because these operations should never produce an unregistered `IndexSpace`(except the `nullspace()`).
SeQuant of course provides a number of partitionings common for chemistry. These are located in `SeQuant/domain/mbpt/convention.hpp`.


### Quasiparticles

In most cases we are interested in using SeQuant to manipulate expressions involving operators in normal order relative to a vacuum state with a finite number of particles, rather than with respect to the genuine vacuum with zero particles. To make such composition easier SeQuant expressions depend on SeQuant _context_, which specifies things like the vacuum type, whether the single-particle (SP) basis is orthonormal, etc. The above SeQuant program used the default context, which assumes the genuine vacuum. The active context can be examined by calling `get_default_context()`, changed via `set_default_context()`, and reset to the default via `reset_default_context()`:

```c++
#include <SeQuant/core/context.hpp>

int main() {
  using namespace sequant;

  // the default is to use genuine vacuum
  assert(get_default_context().vacuum() == Vacuum::Physical);
  // now set the context to a single product of SP states
  set_default_context(Context{Vacuum::SingleProduct, IndexSpaceMetric::Unit, BraKetSymmetry::symm});
  assert(get_default_context().vacuum() == Vacuum::SingleProduct);
  // reset the context back to the default
  reset_default_context();
  assert(get_default_context().vacuum() == Vacuum::Physical);
  
  return 0;
}
```

However, to deal with the single-product vacuum it is necessary to specify which registered `IndexSpace` corresponds to vacuum occupied. This is because the definition of quasiparticle creators and annihilators depends on the operator's action on the single-product vacuum.
```c++
sample_ISR.vacuum_occupied_space(L"x");
```
Please keep in mind that at various points it is implied that spaces included in vacuum occupied have smaller `typeattr`(bitset/ unsigned int)
than those spaces with `typeattr` outside of this range.


we can evaluate Wick's theorem in single-product normal order:

```c++
  auto cp1 = fcrex(p1), cp2 = fcrex(p2);
  auto ap3 = fannx(p3), ap4 = fannx(p4);

  std::wcout << to_latex(ap3 * cp1 * ap4 * cp2) << " = "
             << to_latex(FWickTheorem{ap3 * cp1 * ap4 * cp2}
                             .set_external_indices(std::array{p1, p2, p3, p4})
                             .full_contractions(false)
                             .compute())
             << std::endl;
```

produces

![{{\tilde{a}_{{p_3}}}{\tilde{a}^{{p_1}}}{\tilde{a}_{{p_4}}}{\tilde{a}^{{p_2}}}} = { \bigl({\tilde{a}^{{p_1}{p_2}}_{{p_3}{p_4}}} - {{s^{{p_1}}_{{z_1}}}{s^{{z_1}}_{{p_3}}}{\tilde{a}^{{p_2}}_{{p_4}}}} - {{s^{{p_2}}_{{z_1}}}{s^{{z_1}}_{{p_4}}}{\tilde{a}^{{p_1}}_{{p_3}}}} - {{s^{{p_1}}_{{y_1}}}{s^{{y_1}}_{{p_4}}}{\tilde{a}^{{p_2}}_{{p_3}}}} + {{s^{{p_2}}_{{z_1}}}{s^{{z_1}}_{{p_3}}}{\tilde{a}^{{p_1}}_{{p_4}}}} + {{s^{{p_1}}_{{z_1}}}{s^{{p_2}}_{{z_2}}}{s^{{z_1}}_{{p_3}}}{s^{{z_2}}_{{p_4}}}} + {{s^{{p_1}}_{{y_1}}}{s^{{p_2}}_{{z_1}}}{s^{{z_1}}_{{p_3}}}{s^{{y_1}}_{{p_4}}}}\bigr) }](doc/images/tut-expr4-result1.svg)

. Note that:
- the tilde in `ã` denotes normal order with respect to single-product vacuum, and
- Einstein summation convention is implied, i.e., indices that appear twice in a given product are summed over.

## Operators
In the domain of Many Body Perturbation Theory (MBPT), there are a variety of common operators (`mbpt::Operator)`which 
are recognizable to chemists. For example: `H`, `f`, `T1`, `T2`, etc... More generally, operators are 
defined as some `Tensor` times a normal ordered operator sequence (`NormalOperator<Statistics::FermiDirac>`). These operators are meant 
provide an easy way to construct expressions common to chemistry.
Perhaps even more importantly, algebra at the higher level of operators is often 
more efficient. For example in the screening procedure for Wick's Theorem.

 Unfortunately, even within chemistry, there is no way to uniquely define many of the common operators. 
for example: `T1` in a multi-reference framework may look different from `T1` in a single-reference
picture. In fact, in order to define an excitation or deexcitation operator, we need to know where active(mobile) particles and holes come from.
Thus, in order for an `IndexSpaceRegistry` to be compatible with operator construction requires the user to choose which spaces contain active particles and
which contain active holes.
```c++
sample_ISR.active_particle_space(xy_space);
sample_ISR.active_hole_space(yz_space);
```
Notice that there is no requirement that the two spaces are orthogonal.
# Developers

`SeQuant` is developed by the Valeev group at Virginia Tech.
