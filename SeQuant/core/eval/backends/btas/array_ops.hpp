#ifndef SEQUANT_EVAL_BACKENDS_BTAS_ARRAY_OPS_HPP
#define SEQUANT_EVAL_BACKENDS_BTAS_ARRAY_OPS_HPP

#ifdef SEQUANT_HAS_BTAS

#include <SeQuant/core/eval/backends/btas/result.hpp>
#include <SeQuant/core/eval/extent_array_ops.hpp>

#include <btas/btas.h>

#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace sequant {

/// \brief Build a \c BackendArrayOps for the BTAS backend from a per-basis-key
/// extent map (see make_extent_array_ops).
///
/// \tparam T the \c btas::Tensor type of the zero arrays
/// \param extent_of_base IndexBasis::base_key() -> the extent of that basis
template <typename T>
[[nodiscard]] BackendArrayOps make_btas_array_ops(
    std::map<std::wstring, std::size_t> extent_of_base) {
  return make_extent_array_ops(
      std::move(extent_of_base),
      [](container::svector<std::size_t> const& extents) -> ResultPtr {
        T zeros{btas::Range{
            std::vector<std::size_t>(extents.begin(), extents.end())}};
        zeros.fill(typename T::numeric_type{0});
        return eval_result<ResultTensorBTAS<T>>(std::move(zeros));
      });
}

}  // namespace sequant

#endif  // SEQUANT_HAS_BTAS

#endif  // SEQUANT_EVAL_BACKENDS_BTAS_ARRAY_OPS_HPP
