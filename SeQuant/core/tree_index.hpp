#ifndef SEQUANT_TREEINDEX_H
#define SEQUANT_TREEINDEX_H

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <compare>
#include <concepts>
#include <cstddef>
#include <initializer_list>
#include <iostream>
#include <ranges>

namespace sequant {

/// Represents an index (position) into a tree structure
/// @note Using this to select a subexpression from an expression tree takes
/// O(N) time where N = depth()
class TreeIndex {
 public:
  TreeIndex() = default;

  TreeIndex(container::svector<std::size_t> positions);

  TreeIndex(std::initializer_list<std::size_t> positions);

  template <std::ranges::range Positions>
    requires(std::integral<std::ranges::range_value_t<Positions>>)
  TreeIndex(Positions &&positions) {
    if constexpr (std::ranges::sized_range<Positions>) {
      positions_.reserve(std::ranges::size(positions));
    }

    for (auto val : positions) {
      SEQUANT_ASSERT(val >= 0);
      positions_.push_back(val);
    }
  }

  ExprPtr &select_from(ExprPtr &expr) const;
  const ExprPtr &select_from(const ExprPtr &expr) const;
  Expr &select_from(Expr &expr) const;
  const Expr &select_from(const Expr &expr) const;

  std::size_t depth() const;

  bool operator==(const TreeIndex &) const = default;
  std::strong_ordering operator<=>(const TreeIndex &) const = default;

  template <typename CharT>
  friend std::basic_ostream<CharT> &operator<<(
      std::basic_ostream<CharT> &stream, const TreeIndex &idx) {
    stream << "{";

    for (std::size_t i = 0; i < idx.positions_.size(); ++i) {
      stream << " " << idx.positions_[i] << " ";

      if (i + 1 < idx.positions_.size()) {
        stream << "->";
      }
    }

    stream << "}";

    return stream;
  }

 private:
  container::svector<std::size_t> positions_;
};

}  // namespace sequant

#endif  // SEQUANT_TREEINDEX_H
