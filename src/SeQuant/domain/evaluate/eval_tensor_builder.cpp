#include "eval_tensor_builder.hpp"
#include "eval_tensor.hpp"
#include "path_tree.hpp"

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/tensor.hpp>
#include <boost/functional/hash.hpp>

#include <algorithm>

namespace sequant::evaluate {

// ctor
EvalTensorBuilder::EvalTensorBuilder(bool complex_tensor_data)
    : complex_tensor_data{complex_tensor_data} {}

const EvalTensorPtr& EvalTensorBuilder::get_eval_tree() const {
  return eval_tree_;
}

void EvalTensorBuilder::build_eval_tree(const ExprPtr& expr) {
  // a lambda function to build product type evaltree of a single product
  auto prod_maker = [this](const ExprPtr& prod) {
    // @note using a trivial path tree
    // optimization can be done by finding an optimal path tree
    size_t i = 0;
    auto path_tree = std::make_shared<PathTree>(i);
    for (const auto& sum : *prod) {
      auto child = std::make_shared<PathTree>(i + 1);
      path_tree->add_child(child);
    }
    // because path tree has one label too extra
    path_tree->pop_last_child();
    return build_from_product(prod, path_tree);
  };

  // a lambda function to combine an eval tensor with a summand
  auto sum_accumulator = [&prod_maker, this](const EvalTensorPtr& left,
                                             const ExprPtr& summand) {
    // the @c left param is a binary evaluation tree such that it
    // is an intermediate type of evaluation tensor whose left child
    // is an antisymmetrization/symmetrization eval tensor and the right
    // child is any other kind of valid eval tensor
    auto left_imed = std::dynamic_pointer_cast<EvalTensorIntermediate>(left);
    auto& left_tensor = left_imed->get_right_tensor();
    //
    //
    auto right = prod_maker(summand);
    auto right_imed = std::dynamic_pointer_cast<EvalTensorIntermediate>(right);
    auto right_tensor = right_imed->get_right_tensor();
    auto result_right =
        build_intermediate(left_tensor, right_tensor, Operation::SUM);
    // result_left is the antisymmetrization or the symmetrization operation
    auto result_left = left_imed->get_left_tensor();
    return build_intermediate(left_tensor, right_tensor,
                              left_imed->get_operation());
  };
  // initialize the result with the first summand
  auto& sum = expr->as<Sum>();
  auto init = prod_maker(sum.summand(0));
  // build and store the tree
  eval_tree_ =
      std::accumulate(sum.begin() + 1, sum.end(), init, sum_accumulator);
}

EvalTensorPtr EvalTensorBuilder::build_from_product(
    const ExprPtr& expr, const PathTreePtr& path) const {
  const auto& prod = expr->as<Product>();

  // a lambda function to convert a linear container eg. as a vector, of
  // ExprPtr to sequant Tensors to an evaluation tensor
  auto product_accumulator = [&expr, &prod, this](const EvalTensorPtr& left,
                                                  const PathTreePtr& path) {
    std::shared_ptr<EvalTensor> right =
        path->is_leaf() ? build_leaf(prod.factor(path->get_label()))
                        : build_from_product(expr, path);

    return build_intermediate(left, right, Operation::PRODUCT);
  };

  auto left = build_leaf(prod.factor(path->get_label()));
  auto right =
      std::accumulate(path->get_children().begin(), path->get_children().end(),
                      left, product_accumulator);
  auto result = build_intermediate(left, right, Operation::PRODUCT);
  // @NOTE only real value is used for now
  result->set_scalar(prod.scalar().real());
  return result;
}

EvalTensorPtr EvalTensorBuilder::build_leaf(const ExprPtr& expr) const {
  auto& tnsr = expr->as<Tensor>();

  // if expr references to a real valued tensor check if swapping bra and ket
  // labels is required else set as not required
  bool swap_bk = need_bra_ket_swap(expr);

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
  auto leaf_tensor_ptr = std::make_shared<EvalTensorLeaf>();
  // set indices of the object
  leaf_tensor_ptr->set_indices(bra_index_labels);
  // set hash value of the object
  leaf_tensor_ptr->set_hash_value(hash_leaf(expr, swap_bk));

  return leaf_tensor_ptr;
}

HashType EvalTensorBuilder::hash_leaf(const ExprPtr& expr, bool swap_bk) const {
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

EvalTensorPtr EvalTensorBuilder::build_intermediate(
    const EvalTensorPtr& left_eval_expr, const EvalTensorPtr& right_eval_expr,
    Operation op) const {
  if ((op == Operation::SUM) &&
      (left_eval_expr->get_indices() != right_eval_expr->get_indices())) {
    throw std::domain_error(
        "While summing two tensors, their indices must be identical!");
  }

  // create the intermediate
  auto imed = std::make_shared<EvalTensorIntermediate>();

  // set left and right tensors
  imed->set_left_tensor(left_eval_expr);
  imed->set_right_tensor(right_eval_expr);

  // set the operation type
  imed->set_operation(op);

  // fill the indices
  if (op == Operation::SUM) {
    imed->set_indices(left_eval_expr->get_indices());
  } else if (op == Operation::PRODUCT) {
    IndexLabelContainer lindices = left_eval_expr->get_indices();
    IndexLabelContainer rindices = right_eval_expr->get_indices();

    // indices of intermediate tensors are canonicalized by sorting them
    if (left_eval_expr->is_leaf()) std::sort(lindices.begin(), lindices.end());

    if (right_eval_expr->is_leaf()) std::sort(rindices.begin(), rindices.end());

    IndexLabelContainer imed_indices;
    std::set_symmetric_difference(lindices.begin(), lindices.end(),
                                  rindices.begin(), rindices.end(),
                                  std::back_inserter(imed_indices));
    imed->set_indices(imed_indices);
  } else if ((op == Operation::ANTISYMMETRIZE) ||
             (op == Operation::SYMMETRIZE)) {
    // nothing special needs to be done
    // the evaluator should be smart enough
    // on handling (anti)symmetrization type operations
  }

  // set the hash value
  imed->set_hash_value(hash_intermediate(imed));

  // intermediate is built
  return imed;
}

HashType EvalTensorBuilder::hash_intermediate(
    const EvalTensorPtr& eval_expr) const {
  // intermediate is hashed by the type of operation(sum, prod, etc.)
  // and the contracting/non-contracting index positions of the left and
  // the right tensor indices.
  // index positions are size_t by default
  // to combine the hashes cast Operation type to size_t as well

  auto imed_ptr = std::dynamic_pointer_cast<EvalTensorIntermediate>(eval_expr);

  auto imed_operation = imed_ptr->get_operation();

  boost::hash<size_t> number_hasher;
  HashType imed_hash_value;

  // hashing the type of binary operation by casting operation into size_t
  imed_hash_value = number_hasher(static_cast<size_t>(imed_operation));

  // operation is of antisymmetrization type then its hash value is solely
  // determined by which indices are being permuted, this information is
  // carried by the left tensor -- the tensor with label 'A'
  if ((imed_operation == Operation::ANTISYMMETRIZE) ||
      (imed_operation == Operation::SYMMETRIZE)) {
    boost::hash_combine(imed_hash_value,
                        imed_ptr->get_left_tensor()->get_hash_value());
    return imed_hash_value;
  }

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

  // when Operation is product ie. a tensor contraction, not only the left and
  // the right tensor's hash values are important but also the information
  // about which of their indices got contracted is needed to uniquely
  // identify such an evaluation
  if (imed_operation == Operation::PRODUCT) {
    // a lambda expression that combines the hashes of the non-contracting
    // indices of right or left tensor
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
    };  // lambda function index_hash_combiner

    // when the same contraction occurs with left tensor switched to the right
    // we should obtain the same hash value. ie hashing should be agnostic to
    // the order of the two contracting tensors
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

bool EvalTensorBuilder::need_bra_ket_swap(const ExprPtr& expr) const {
  auto& tnsr = expr->as<Tensor>();

  // non-symmetric and bra-kets with invalid symmetries are handled
  if ((tnsr.braket_symmetry() == BraKetSymmetry::nonsymm) ||
      (tnsr.braket_symmetry() == BraKetSymmetry::invalid)) {
    return false;
  }

  // when working with complex tensors if it doesn't have symmetric bra-ket
  if (complex_tensor_data && (tnsr.braket_symmetry() != BraKetSymmetry::symm)) {
    return false;
  }

  if (tnsr.bra_rank() != tnsr.ket_rank())
    throw std::domain_error(
        "Only equal bra_rank() and ket_rank() are handled for now!");
  //

  for (auto i = 0; i < tnsr.rank(); ++i) {
    // if bra index and ket index at the same positions have the same
    // index space then we don't know what to do yet, so
    if (tnsr.bra()[i].space() == tnsr.ket()[i].space()) continue;

    // for the corresponding positions if the index space value of the ket
    // index is lower than that of the bra index, then swapping should be done
    if (tnsr.ket()[i].space() < tnsr.bra()[i].space()) return true;
  }

  // survived this far implies no swapping is needed
  return false;
}

}  // namespace sequant::evaluate
