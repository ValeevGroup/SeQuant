//
// Created by Bimal Gaudel on 3/27/25.
//

#ifndef SEQUANT_EVAL_FWD_HPP
#define SEQUANT_EVAL_FWD_HPP

#include <memory>

namespace sequant {

/// Backend-agnostic flag to control tensor de-nesting behavior during products.
/// When multiplying tensor-of-tensor types, this controls whether the result
/// should be "de-nested" (flattened) to a regular tensor or kept as nested.
enum class DeNest { True, False };

/// \brief Flavor of an \c EvalExpr::batch_axes entry: whether the index is
///        contracted away somewhere in the term (\c Contracted), or is an
///        external/spectator index that survives to the term's result
///        (\c External). See \c EvalExpr::batch_axes.
enum class AxisKind { Contracted, External };

template <typename TreeNode, bool force_hash_collisions = false>
class CacheManager;

class Result;

///
/// \brief Managed pointer to the result of an evaluation.
///
using ResultPtr = std::shared_ptr<Result>;

}  // namespace sequant

#endif  // SEQUANT_EVAL_FWD_HPP
