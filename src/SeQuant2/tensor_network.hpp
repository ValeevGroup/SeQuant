//
// Created by Eduard Valeyev on 2019-02-02.
//

#ifndef SEQUANT2_TENSOR_NETWORK_H
#define SEQUANT2_TENSOR_NETWORK_H

#include "../SeQuant2/tensor.hpp"

namespace sequant2 {

/// @brief A (non-directed) graph view of a set of Tensor objects

/// The main role of this is to canonize itself. Since Tensors can be connected
/// by multiple Index'es (thus edges are colored), what is canonized is actually
/// the graph of indices (i.e. the dual of the tensor graph), with Tensors
/// represented by one or more vertices.
class TensorNetwork {
 public:
  using IndexRef = std::reference_wrapper<Index>;

  template<typename ExprPtrRange>
  TensorNetwork(const ExprPtrRange &exprptr_range) {
    const bool contains_a_nontensor = ranges::any_of(
        exprptr_range,
        [](const ExprPtr &exprptr) {
          return !exprptr->is<Tensor>();
        });
    if (contains_a_nontensor)
      throw std::logic_error("TensorNetwork(exprptr_range): exprptr_range contains a non-Tensor");

    auto tsrptr_range = exprptr_range | ranges::view::transform([](const ExprPtr &ex) {
      return std::static_pointer_cast<Tensor>(ex);
    });
    assert(ranges::size(tsrptr_range)
               == ranges::size(exprptr_range)); // !!! FAILS !!! ranges::size(tsrptr_range) is always zero
    tensors_ = decltype(tensors_)(ranges::begin(tsrptr_range), ranges::end(tsrptr_range));
  }

  void canonicalize() { assert(false && "not yet implemented"); }

 private:
  // source tensors and indices
  container::svector<TensorPtr> tensors_;
  container::svector<IndexRef> indices_;

  // avalable after canonization
  container::svector<TensorPtr> canon_tensors_;
  container::svector<IndexRef> canon_indices_;
};

}  // namespace sequant2

#endif  // SEQUANT2_TENSOR_NETWORK_H
