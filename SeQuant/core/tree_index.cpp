#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/tree_index.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <concepts>
#include <cstddef>
#include <initializer_list>
#include <ranges>
#include <string>
#include <type_traits>

namespace sequant {

TreeIndex::TreeIndex(container::svector<std::size_t> positions)
    : positions_(std::move(positions)) {}

TreeIndex::TreeIndex(std::initializer_list<std::size_t> positions)
    : TreeIndex(container::svector<std::size_t>(std::move(positions))) {}

template <typename ExprType>
ExprType &select_impl(ExprType &expr,
                      const container::svector<std::size_t> &positions) {
  constexpr bool is_expr_ptr =
      std::same_as<std::remove_cvref_t<ExprType>, ExprPtr>;

  ExprType *selected = &expr;

  for (std::size_t current : positions) {
    if (current >= selected->size()) {
      throw Exception("Position " + std::to_string(current) +
                      " is out of bounds for an expression of dimension " +
                      std::to_string(selected->size()));
    }

    using std::ranges::begin;

    auto it = begin(*selected) + current;

    if constexpr (is_expr_ptr) {
      selected = &(*it);
    } else {
      selected = &(*(*it));
    }
  }

  SEQUANT_ASSERT(selected);

  return *selected;
}

ExprPtr &TreeIndex::select_from(ExprPtr &expr) const {
  return select_impl(expr, positions_);
}

const ExprPtr &TreeIndex::select_from(const ExprPtr &expr) const {
  return select_impl(expr, positions_);
}

Expr &TreeIndex::select_from(Expr &expr) const {
  return select_impl(expr, positions_);
}

const Expr &TreeIndex::select_from(const Expr &expr) const {
  return select_impl(expr, positions_);
}

std::size_t TreeIndex::depth() const { return positions_.size(); }

}  // namespace sequant
