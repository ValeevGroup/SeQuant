//
// Created by Robert Adam on 2023-09-27
//

#ifndef SEQUANT_CORE_UTILITY_INDICES_HPP
#define SEQUANT_CORE_UTILITY_INDICES_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/tensor.hpp>

#include <algorithm>
#include <tuple>
#include <type_traits>

namespace sequant {

namespace detail {
template <typename Range>
struct not_in {
  const Range &range;

  not_in(const Range &range) : range(range) {}

  template <typename T>
  bool operator()(const T &element) const {
    return std::find(range.begin(), range.end(), element) == range.end();
  }
};
}  // namespace detail

/// @returns Lists of non-contracted indices arising when contracting the two
/// given tensors in the order bra, ket, auxiliary
template <typename Container = container::svector<Index>>
std::tuple<Container, Container, Container> get_uncontracted_indices(
    const Tensor &t1, const Tensor &t2) {
  static_assert(std::is_same_v<typename Container::value_type, Index>);

  Container bra;
  Container ket;
  Container auxiliary;

  // Bra indices
  std::copy_if(t1.bra().begin(), t1.bra().end(), std::back_inserter(bra),
               detail::not_in{t2.ket()});
  std::copy_if(t2.bra().begin(), t2.bra().end(), std::back_inserter(bra),
               detail::not_in{t1.ket()});

  // Ket indices
  std::copy_if(t1.ket().begin(), t1.ket().end(), std::back_inserter(ket),
               detail::not_in{t2.bra()});
  std::copy_if(t2.ket().begin(), t2.ket().end(), std::back_inserter(ket),
               detail::not_in{t1.bra()});

  // Auxiliary indices
  std::copy_if(t1.auxiliary().begin(), t1.auxiliary().end(),
               std::back_inserter(auxiliary), detail::not_in{t2.auxiliary()});
  std::copy_if(t2.auxiliary().begin(), t2.auxiliary().end(),
               std::back_inserter(auxiliary), detail::not_in{t1.auxiliary()});

  return {std::move(bra), std::move(ket), std::move(auxiliary)};
}

}  // namespace sequant

#endif  // SEQUANT_CORE_UTILITY_INDICES_HPP
