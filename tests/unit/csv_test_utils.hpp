//
// Shared helpers for the basis-instance tests.
//

#ifndef SEQUANT_TESTS_UNIT_CSV_TEST_UTILS_HPP
#define SEQUANT_TESTS_UNIT_CSV_TEST_UTILS_HPP

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/index_basis.hpp>

#include <cstddef>
#include <map>

namespace sequant::tests::csv {

/// calls @p f on every bra/ket/aux slot of every tensor of @p expr
template <typename F>
void for_each_index(ExprPtr const& expr, F&& f) {
  expr->visit(
      [&f](ExprPtr const& x) {
        if (!x->is<AbstractTensor>()) return;
        for (auto const& idx : x->as<AbstractTensor>()._slots()) f(idx);
      },
      /* atoms_only = */ true);
}

/// @return the number of tensor slots of @p expr per basis instance
/// (`std::nullopt` counts the slots without one)
inline std::map<IndexBasis::optional_instance, std::size_t> instance_histogram(
    ExprPtr const& expr) {
  std::map<IndexBasis::optional_instance, std::size_t> result;
  for_each_index(expr, [&result](Index const& idx) {
    ++result[idx.basis().basis_instance()];
  });
  return result;
}

}  // namespace sequant::tests::csv

#endif  // SEQUANT_TESTS_UNIT_CSV_TEST_UTILS_HPP
