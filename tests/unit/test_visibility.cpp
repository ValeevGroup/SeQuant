//
// Checks that types whose RTTI crosses library boundaries have a single
// type_info in the process (see doc/developer/visibility.rst).
//

#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/domain/mbpt/antisymmetrizer.hpp>
#include <SeQuant/domain/mbpt/op.hpp>

#include <typeinfo>

namespace {

/// @return the type_info of the first atom of @p expr that is a @p T
template <typename T>
const std::type_info* find_atom_type(const sequant::ExprPtr& expr) {
  const std::type_info* result = nullptr;
  expr->visit(
      [&](const sequant::ExprPtr& e) {
        // bind first: typeid(*e) warns about the side effects of operator*
        const sequant::Expr& atom = *e;
        if (!result && atom.is<T>()) result = &typeid(atom);
      },
      /* atoms_only = */ true);
  return result;
}

}  // namespace

// Comparing type_info objects by address, as libc++ does on Linux and x86_64
// macOS, detects a type that is hidden in one of the libraries also with
// libstdc++. On Windows every DLL has its own type_info objects (they are
// compared by name), so the identity check does not apply there.
TEST_CASE("rtti-across-libraries", "[visibility]") {
#if defined(_WIN32) || defined(__CYGWIN__)
  SUCCEED("type_info identity is per DLL on Windows");
#else
  using namespace sequant;

  SECTION("exception thrown by SeQuant-symb") {
    bool caught = false;
    try {
      [[maybe_unused]] auto e = io::serialization::from_string(L"t{i1;a1");
    } catch (const Exception& e) {
      caught = true;
      CHECK(&typeid(e) == &typeid(io::serialization::SerializationError));
    }
    REQUIRE(caught);
  }

  SECTION("NormalOperator created by SeQuant-mbpt") {
    // no contracted indices: only the operator's construction matters here
    const auto product =
        ex<Tensor>(L"t", bra{Index(L"i_1")}, ket{Index(L"a_1")}) *
        ex<FNOperator>(cre{Index(L"a_2"), Index(L"a_3")},
                       ann{Index(L"i_2"), Index(L"i_3")});
    const auto result = mbpt::antisymmetrize(product).result;
    const auto* type = find_atom_type<FNOperator>(result);
    REQUIRE(type != nullptr);
    CHECK(type == &typeid(FNOperator));
  }

  SECTION("mbpt::Operator created by SeQuant-mbpt") {
    const auto* type = find_atom_type<mbpt::op_t>(mbpt::op::T(1));
    REQUIRE(type != nullptr);
    CHECK(type == &typeid(mbpt::op_t));
  }
#endif
}
