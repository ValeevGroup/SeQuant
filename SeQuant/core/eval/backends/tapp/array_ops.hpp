#ifndef SEQUANT_EVAL_BACKENDS_TAPP_ARRAY_OPS_HPP
#define SEQUANT_EVAL_BACKENDS_TAPP_ARRAY_OPS_HPP

#ifdef SEQUANT_HAS_TAPP

#include <SeQuant/core/eval/backends/tapp/result.hpp>
#include <SeQuant/core/eval/extent_array_ops.hpp>

#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <utility>

namespace sequant {

/// \brief Build a \c BackendArrayOps for the TAPP backend from a per-basis-key
/// extent map (see make_extent_array_ops).
///
/// \tparam T the \c TAPPTensor type of the zero arrays
/// \param extent_of_base IndexBasis::base_key() -> the extent of that basis
template <typename T>
[[nodiscard]] BackendArrayOps make_tapp_array_ops(
    std::map<std::wstring, std::size_t> extent_of_base) {
  return make_extent_array_ops(
      std::move(extent_of_base),
      [](container::svector<std::size_t> const& extents) -> ResultPtr {
        // a TAPPTensor is zero-initialized
        return eval_result<ResultTensorTAPP<T>>(T{
            container::svector<std::int64_t>(extents.begin(), extents.end())});
      });
}

}  // namespace sequant

#endif  // SEQUANT_HAS_TAPP

#endif  // SEQUANT_EVAL_BACKENDS_TAPP_ARRAY_OPS_HPP
