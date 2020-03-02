#ifndef SEQUANT_FACTORIZE_FACTORIZER_HPP
#define SEQUANT_FACTORIZE_FACTORIZER_HPP

///
/// Find optimal evaluation sequences of EvalTensor evaluation tree that result
/// into lower ops count.
///
/// @author Bimal Gaudel
/// @version Feb 2020
///

#include "eval_tensor_fwd.hpp"

#include <SeQuant/core/expr_fwd.hpp>
#include <tuple>

namespace sequant::factorize {

/// Generate the corresponding EvalTensor's of @c lexpr and @c rexpr. Such that
/// the resulting EvalTensor's have the maximum number of common sub-expressions
/// between them.
///
/// @param lexpr left sequant Expr to be fused.
/// @param rexpr right sequant Expr to be fused.
/// @return A tuple of shared pointers to EvalTensor's.
std::tuple<EvalTensorPtr, EvalTensorPtr> fuse_optimally(const ExprPtr& lexpr,
                                                        const ExprPtr& rexpr);

namespace detail {

/// Get the hash values of all nodes in a evaluation tensor.
///
/// @param tensor EvalTensor pointer.
/// @return A set of hash values of each node in the EvalTensor.
container::set<HashType> get_hash_values(const EvalTensorPtr& tensor);

///
/// Get the sum of operation counts for evaluations in a tensor whose hash
/// values do not exist in a set of hash_values passed as a parameter.
///
/// @param tensor A pointer to EvalTensor.
/// @param hash_values A set of hash_values to look at while dis-counting.
/// @return Sum of operations count of evaluations.
OpsCount get_unique_ops_count(const EvalTensorPtr& tensor,
                              const container::set<HashType>& hash_values);

}  // namespace detail

}  // namespace sequant::factorize

#endif /* ifndef SEQUANT_FACTORIZE_FACTORIZER_HPP */
