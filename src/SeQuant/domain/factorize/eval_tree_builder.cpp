#include "eval_tree_builder.hpp"
#include "eval_tensor_fwd.hpp"
#include "intermediate_eval_tensor.hpp"
#include "leaf_eval_tensor.hpp"
#include "path_tree.hpp"

#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/tensor.hpp>
#include <boost/functional/hash.hpp>

#include <algorithm>

namespace sequant::factorize {

/// Getter for the built EvalTree.
const EvalTensorPtr& EvalTreeBuilder::get_eval_tree() const {
  return eval_tree_;
}

EvalTensorPtr EvalTreeBuilder::build_leaf(const ExprPtr& expr,
                                          bool real_valued) {
  auto& tnsr = expr->as<Tensor>();

  // if expr references to a real valued tensor check if swapping bra and ket
  // labels is required else set as not required
  bool swap_bk = real_valued ? swap_bra_ket_labels(expr) : false;

  // get the labels of indices in bra and ket
  IndexLabelContainer bra_index_labels, ket_index_labels;
  for (const auto& idx : tnsr.bra()) bra_index_labels.emplace_back(idx.label());
  for (const auto& idx : tnsr.ket()) ket_index_labels.emplace_back(idx.label());

  // swap labels as required
  if (swap_bk) std::swap(bra_index_labels, ket_index_labels);

  // collect bra and ket indices by appending ket indices to the bra index
  // container
  std::move(ket_index_labels.begin(), ket_index_labels.end(),
            std::back_inserter(bra_index_labels));

  // create an object to return
  auto leaf_tensor_ptr = std::make_shared<LeafEvalTensor>();
  // set indices of the object
  leaf_tensor_ptr->set_indices(bra_index_labels);
  // set hash value of the object
  leaf_tensor_ptr->set_hash_value(hash_leaf(expr, swap_bk));

  return leaf_tensor_ptr;
}

EvalTensorPtr EvalTreeBuilder::build_intermediate(const EvalTensorPtr& ltensor,
                                                  const EvalTensorPtr& rtensor,
                                                  Operation op) {
  if ((op == Operation::SUM) &&
      (ltensor->get_indices() != rtensor->get_indices())) {
    throw std::domain_error(
        "While summing two tensors, their indices must be identical!");
  }

  // create the intermediate
  auto imed = std::make_shared<IntermediateEvalTensor>();

  // set left and right tensors
  imed->set_left_ptr(ltensor);
  imed->set_right_ptr(rtensor);

  // set the operation type
  imed->set_operation(op);

  // fill the indices
  if (op == Operation::SUM) {
    imed->set_indices(ltensor->get_indices());
  } else {
    IndexLabelContainer lindices = ltensor->get_indices();
    IndexLabelContainer rindices = rtensor->get_indices();

    if (ltensor->is_leaf()) std::sort(lindices.begin(), lindices.end());

    if (rtensor->is_leaf()) std::sort(rindices.begin(), rindices.end());

    IndexLabelContainer imed_indices;
    std::set_symmetric_difference(lindices.begin(), lindices.end(),
                                  rindices.begin(), rindices.end(),
                                  std::back_inserter(imed_indices));
    imed->set_indices(imed_indices);
  }

  // set the hash value
  imed->set_hash_value(hash_intermediate(imed));

  // intermediate is built
  return imed;
}

HashType EvalTreeBuilder::hash_leaf(const ExprPtr& expr, bool swap_bk) {
  auto& tnsr = expr->as<Tensor>();
  HashType leaf_hash_value;

  boost::hash<std::wstring_view> label_hasher;
  leaf_hash_value = label_hasher(tnsr.label());

  using ispace_cast = size_t;
  boost::hash<ispace_cast> number_hasher;
  using index_container_type =
      decltype(tnsr.bra());  // const svector<Index, ..>

  auto hash_idx_space = [&leaf_hash_value,
                         &number_hasher](index_container_type& idx_vec) {
    for (const auto& idx : idx_vec)
      boost::hash_combine(
          leaf_hash_value,
          number_hasher(static_cast<ispace_cast>(idx.space().type())));
  };

  if (swap_bk) {
    hash_idx_space(tnsr.ket());
    hash_idx_space(tnsr.bra());
  } else {
    hash_idx_space(tnsr.bra());
    hash_idx_space(tnsr.ket());
  }

  return leaf_hash_value;
}

HashType EvalTreeBuilder::hash_intermediate(const EvalTensorPtr& evtensor) {
  // intermediate is hashed by the type of operation(sum, prod, etc.)
  // and the contracting/non-contracting index positions of the left and
  // the right tensor indices.
  // index positions are size_t by default
  // to combine the hashes cast Operation type to size_t as well

  auto imed_ptr = std::dynamic_pointer_cast<IntermediateEvalTensor>(evtensor);

  boost::hash<size_t> number_hasher;
  HashType imed_hash_value;

  // hashing the type of binary operation by casting operation into size_t
  imed_hash_value =
      number_hasher(static_cast<size_t>(imed_ptr->get_operation()));

  // hash the hash values of the left and right tensors
  // but do so in an order-agnostic manner
  HashType lhash_value = imed_ptr->get_left_tensor()->get_hash_value();
  HashType rhash_value = imed_ptr->get_right_tensor()->get_hash_value();

  if (lhash_value < rhash_value) {
    boost::hash_combine(imed_hash_value, lhash_value);
    boost::hash_combine(imed_hash_value, rhash_value);
  } else {
    boost::hash_combine(imed_hash_value, rhash_value);
    boost::hash_combine(imed_hash_value, lhash_value);
  }

  // the Operation is product ie. a tensor contraction, not only the left and
  // the right tensor's hash values are important but also the information
  // about which of their indices got  contracted is needed to uniquely
  // identify such an evaluation
  if (imed_ptr->get_operation() == Operation::PRODUCT) {
    // a lambda expression that combines the hashes of the uncontracting indices
    // of right or left tensor
    auto index_hash_combiner = [&](const EvalTensorPtr& lr_tensor) {
      // iterate through the imed_ptr indices
      for (const auto& idx : imed_ptr->get_indices()) {
        // iterate through the indices of the passed left/right tensor
        // by a counter
        for (auto i = 0; i < lr_tensor->get_indices().size(); ++i) {
          // if the index of the imed_ptr matches with the lr_tensor
          // then time to break as we found the contracted indices
          // no need to look further
          if (idx == lr_tensor->get_indices()[i]) {
            break;
          }
          // combine the hash of the counter of non-contracted index
          boost::hash_combine(imed_hash_value, number_hasher(i));
        }
      }
    };

    if (lhash_value < rhash_value) {
      index_hash_combiner(imed_ptr->get_left_tensor());
      index_hash_combiner(imed_ptr->get_right_tensor());
    } else {
      index_hash_combiner(imed_ptr->get_right_tensor());
      index_hash_combiner(imed_ptr->get_left_tensor());
    }
  }  // done combining hashes for product type evaluation
  return imed_hash_value;
}

bool EvalTreeBuilder::swap_bra_ket_labels(const ExprPtr& expr) {
  auto& tnsr = expr->as<Tensor>();

  if (tnsr.bra_rank() != tnsr.ket_rank())
    throw std::domain_error(
        "Only equal bra_rank() and ket_rank() are handled for now!");

  for (auto i = 0; i < tnsr.rank(); ++i) {
    if (tnsr.bra()[i].label()[0] == tnsr.ket()[i].label()[0]) continue;
    if (!(tnsr.bra()[i].label() > tnsr.ket()[i].label())) {
      return true;
    }  // if
  }    // for
  return false;
}

}  // namespace sequant::factorize
