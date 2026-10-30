#ifndef SEQUANT_EVAL_BACKENDS_TILEDARRAY_ARRAY_OPS_HPP
#define SEQUANT_EVAL_BACKENDS_TILEDARRAY_ARRAY_OPS_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/backend_array_ops.hpp>
#include <SeQuant/core/eval/backends/tiledarray/result.hpp>
#include <SeQuant/core/index.hpp>

#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace sequant {

/// \brief Build a \c BackendArrayOps for the TiledArray backend from a
/// per-basis-key tiling map.
///
/// \tparam FlatArray the flat (Tensor-of-Scalars) \c TA::DistArray type;
/// \tparam ToTArray  the nested (Tensor-of-Tensor) \c TA::DistArray type.
///
/// \param tr1_of_base IndexBasis::base_key() -> the basis' full \c
/// TA::TiledRange1. \param world       the World the zero destinations are
/// built in; it must
///                    outlive the returned closures' use.
///
/// \details The two closures are the TA realization of external-axis batching's
/// backend needs (see \c BackendArrayOps): \c make_zeros builds a zero
/// destination -- flat or nested per the descriptor's proto structure, with a
/// nested result's inner tiles left empty and a sparse flat result's shape all
/// zero for the scatter writes to fill -- and \c axis_batches chunks an axis on
/// its basis' tile boundaries. Tiling is a property of the basis, so both are
/// sourced from \p tr1_of_base alone; no array in the DAG is consulted. Both
/// mpqc (from its orbital/basis registries) and unit tests (from the tranges
/// they build their leaves with) supply the map and call this.
template <typename FlatArray, typename ToTArray = FlatArray>
[[nodiscard]] BackendArrayOps make_ta_array_ops(
    std::map<std::wstring, TA::TiledRange1> tr1_of_base, TA::World& world) {
  auto map = std::make_shared<std::map<std::wstring, TA::TiledRange1>>(
      std::move(tr1_of_base));
  BackendArrayOps aops;
  aops.axis_batches = [map](Index const& axis, std::size_t target_batch_size) {
    return mode_batches_of_trange1(
        map->at(std::wstring(axis.basis().base_key())), target_batch_size);
  };
  aops.make_zeros =
      [map, &world](container::vector<Index> const& descriptor) -> ResultPtr {
    using numeric_type = typename FlatArray::numeric_type;
    std::vector<TA::TiledRange1> outer;
    bool nested = false;
    for (auto const& ix : descriptor) {
      if (ix.has_proto_indices()) {
        nested = true;  // an inner (nested) mode -- not an outer trange mode
        continue;
      }
      outer.push_back(map->at(std::wstring(ix.basis().base_key())));
    }
    TA::TiledRange otr(outer.begin(), outer.end());
    auto make_flat = [&]() -> ResultPtr {
      using shape_type = typename FlatArray::shape_type;
      if constexpr (!shape_type::is_dense()) {
        // An all-zero shape holds no tiles, so neither this nor a slice of it
        // allocates; the scatter's block assignment adds the tiles it writes
        // to the shape. A tile is zero only if its norm is below the threshold,
        // so under a zero threshold even a zero norm marks a tile nonzero: the
        // shape would claim tiles that are never set.
        if (shape_type::threshold() > 0)
          return eval_result<ResultTensorTA<FlatArray>>(FlatArray(
              world, otr, shape_type(typename shape_type::value_type(0), otr)));
      }
      FlatArray dest(world, otr);
      dest.fill_local(numeric_type(0));
      world.gop.fence();
      return eval_result<ResultTensorTA<FlatArray>>(std::move(dest));
    };
    if constexpr (std::is_same_v<FlatArray, ToTArray>) {
      // Flat-only backend (default ToTArray == FlatArray): descriptors are
      // flat.
      SEQUANT_ASSERT(!nested &&
                     "flat-only TA array-ops asked for a nested zero result");
      return make_flat();
    } else {
      if (!nested) return make_flat();
      // Nested: the outer trange from the non-proto indices, every outer tile
      // an empty-inner tensor (tot_inner_rank() == 0) -- a valid zero ToT whose
      // real inner tensors the scatter's write_into_slice fills in per batch.
      using OuterT = typename ToTArray::value_type;
      ToTArray dest(world, otr);
      for (auto it = dest.begin(); it != dest.end(); ++it)
        if (dest.is_local(it.index())) *it = OuterT{it.make_range()};
      world.gop.fence();
      return eval_result<ResultTensorOfTensorTA<ToTArray>>(std::move(dest));
    }
  };
  return aops;
}

}  // namespace sequant

#endif  // SEQUANT_EVAL_BACKENDS_TILEDARRAY_ARRAY_OPS_HPP
