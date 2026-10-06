#ifndef SEQUANT_CORE_TENSOR_CANONICALIZER_FWD_HPP
#define SEQUANT_CORE_TENSOR_CANONICALIZER_FWD_HPP

#include <functional>
#include <utility>

namespace sequant {

class Index;
class TensorCanonicalizer;

/// compares Index objects during Tensor canonicalization
/// @sa TensorCanonicalizer::index_comparer_t
using tensor_index_comparer_t = std::function<bool(const Index&, const Index&)>;
/// @sa TensorCanonicalizer::index_pair_t
using tensor_index_pair_t = std::pair<const Index, const Index>;
/// compares pairs of Index objects during Tensor canonicalization
/// @sa TensorCanonicalizer::index_pair_comparer_t
using tensor_index_pair_comparer_t =
    std::function<bool(const tensor_index_pair_t&, const tensor_index_pair_t)>;

}  // namespace sequant

#endif  // SEQUANT_CORE_TENSOR_CANONICALIZER_FWD_HPP
