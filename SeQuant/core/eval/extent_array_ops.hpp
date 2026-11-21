#ifndef SEQUANT_EVAL_EXTENT_ARRAY_OPS_HPP
#define SEQUANT_EVAL_EXTENT_ARRAY_OPS_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/backend_array_ops.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <algorithm>
#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <utility>

namespace sequant {

/// \brief Build a \c BackendArrayOps for a backend whose arrays are untiled,
/// from a per-basis-key extent map.
///
/// \param extent_of_base IndexBasis::base_key() -> the extent of that basis
/// \param make_zeros the backend's zero array of the given extents, callable as
///        `ResultPtr(container::svector<std::size_t> const& extents)`
///
/// \details \c make_zeros of the result sizes each mode of its descriptor by
/// its basis' extent; \c axis_batches chunks an axis into consecutive batches
/// of the target size (the last may be smaller). A nested (proto-bearing)
/// index is not supported.
/// \throw Exception from either closure if an index's basis is not in \p
///        extent_of_base, or from \c make_zeros if it is nested
template <typename MakeZeros>
[[nodiscard]] BackendArrayOps make_extent_array_ops(
    std::map<std::wstring, std::size_t> extent_of_base, MakeZeros make_zeros) {
  auto map = std::make_shared<std::map<std::wstring, std::size_t>>(
      std::move(extent_of_base));
  auto extent = [map](Index const& ix) {
    auto const it = map->find(std::wstring(ix.basis().base_key()));
    if (it == map->end())
      throw Exception("array ops: no extent for the basis of index " +
                      toUtf8(ix.full_label()));
    return it->second;
  };
  BackendArrayOps aops;
  aops.axis_batches = [extent](Index const& axis,
                               std::size_t target_batch_size) {
    container::svector<std::pair<std::size_t, std::size_t>> batches;
    auto const n = extent(axis);
    auto const step = std::max<std::size_t>(target_batch_size, 1);
    for (std::size_t lo = 0; lo < n; lo += step)
      batches.emplace_back(lo, std::min(lo + step, n));
    return batches;
  };
  aops.make_zeros = [extent, make_zeros = std::move(make_zeros)](
                        container::vector<Index> const& descriptor) {
    container::svector<std::size_t> extents;
    extents.reserve(descriptor.size());
    for (auto const& ix : descriptor) {
      if (ix.has_proto_indices())
        throw Exception("array ops: a nested zero array is not supported");
      extents.push_back(extent(ix));
    }
    return ResultPtr(make_zeros(extents));
  };
  return aops;
}

}  // namespace sequant

#endif  // SEQUANT_EVAL_EXTENT_ARRAY_OPS_HPP
