#include "ops_count.hpp"

#include <SeQuant/core/tensor.hpp>

#include <cmath>
#include <numeric>
#include <set>

namespace sequant::factorize {

OpsCalcResult ops_count(const ExprPtr& prod, const rooted_tree& tree,  //
                        size_t nocc, size_t nvirt) {
  if (tree.children.empty()) {
    auto& tnsr = prod->at(tree.label)->as<Tensor>();
    auto result = OpsCalcResult{};
    result.flops = 0;
    result.indices = OpsCalcResult::idx_container_type(
        tnsr.const_braket().begin(), tnsr.const_braket().end());
    return result;
  }

  auto combine_results = [nocc, nvirt](const auto& lhs, const auto& rhs) {
    size_t ocount{0}, vcount{0};
    OpsCalcResult::idx_container_type unique_indices{};
    for (const auto& tt : {lhs.indices, rhs.indices})
      for (const auto& ii : tt) {
        auto success = unique_indices.insert(ii);
        if (!success.second) {
          unique_indices.erase(success.first);
          continue;
        }
        if (ii.space() == IndexSpace::active_occupied)
          ++ocount;
        else
          ++vcount;
      }
    auto thisOps = (ocount > 0 ? std::pow(nocc, ocount) : 1) *
                   (vcount > 0 ? std::pow(nvirt, vcount) : 1);

    auto ops = lhs.flops + rhs.flops;
    ops += thisOps > 1 ? (2 * thisOps) : 0;
    return OpsCalcResult{ops, unique_indices};
  };

  return std::accumulate(
      tree.children.begin(), tree.children.end(),

      ops_count(prod, rooted_tree{tree.label}, nocc, nvirt),

      [&combine_results, nocc, nvirt, &prod](const auto& x, const auto& y) {
        return combine_results(x, ops_count(prod, y, nocc, nvirt));
      }  //
  );
}

}  // namespace sequant::factorize
