#ifndef SEQUANT_EVALUATE_EVAL_TENSOR_BUILDER_HPP
#define SEQUANT_EVALUATE_EVAL_TENSOR_BUILDER_HPP

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/tensor.hpp>
#include <boost/functional/hash.hpp>

#include <algorithm>

#include "eval_tensor.hpp"

///
/// Evaluation tree builder from sequant Expr.
///
/// @author Bimal Gaudel
/// @version Feb 2020
///
namespace sequant::evaluate {
template <typename DataTensorType>
class EvalTensorBuilder {
 private:
  /// By default assume the tensor data is real -- not complex.
  const bool complex_tensor_data{false};

 public:
  /// Constructor.
  /// @param complex_tensor_data false by default. Set to true when working with
  /// complex-valued tensors.
  explicit EvalTensorBuilder(bool complex_tensor_data = false)
      : complex_tensor_data{complex_tensor_data} {};

  /// Build evaluation tree from a sequant expression.
  /// @note It is assumed that the expr is a sum and each summand is a product
  /// such that the first factor of the product is either an antisymmetrization
  /// tensor or symmetrization tensor.
  /// @param expr sequant ExprPtr.
  /// @return Pointer to the Evaluation Tensor
  /// @throws std::domain_error when the expr is neither Sum, Product or Tensor.
  EvalTensorPtr<DataTensorType> build_tree(const ExprPtr& expr) const {
    if (expr->is<Tensor>()) {
      return build_leaf(expr);
    }
    if (expr->is<Product>()) {
      return build_product(expr);
    }
    if (expr->is<Sum>()) {
      return build_sum(expr);
    }
    throw std::domain_error("Only sum, product or tensor is allowed!");
  }

 private:
  /// Build EvalTensor from a sequant Product.
  /// @param expr sequant ExprPtr to sequant Product.
  ///
  /// @return EvalTensor pointer.
  EvalTensorPtr<DataTensorType> build_product(const ExprPtr& expr) const {
    auto prod_accumulator = [this](const EvalTensorPtr<DataTensorType>& lexpr,
                                   const ExprPtr& factor) {
      return build_intermediate(lexpr, build_tree(factor), Operation::PRODUCT);
    };
    auto& prod = expr->as<Product>();
    auto& tnsr = prod.factor(0)->as<Tensor>();
    if ((tnsr.label() == L"A") || (tnsr.label() == L"P")) {
      auto right = std::make_shared<Product>(expr->begin() + 1, expr->end());
      right->scale(prod.scalar());
      return build_intermediate(build_tree(prod.factor(0)), build_tree(right),
                                tnsr.label() == L"A" ? Operation::ANTISYMMETRIZE
                                                     : Operation::SYMMETRIZE);
    }
    auto init = build_tree(prod.factor(0));
    init->set_scalar(prod.scalar().real());
    return std::accumulate(prod.begin() + 1, prod.end(), init,
                           prod_accumulator);
  }

  /// Build EvalTensor from a sequant Sum.
  /// @param expr sequant ExprPtr to sequant Sum.
  ///
  /// @return EvalTensor pointer.
  EvalTensorPtr<DataTensorType> build_sum(const ExprPtr& expr) const {
    auto sum_accumulator = [this](const EvalTensorPtr<DataTensorType>& lexpr,
                                  const ExprPtr& summand) {
      return build_intermediate(lexpr, build_tree(summand), Operation::SUM);
    };
    auto& sum = expr->as<Sum>();
    return std::accumulate(sum.begin() + 1, sum.end(),
                           build_tree(sum.summand(0)), sum_accumulator);
  }

  /// Build leaf EvalTensor from sequant tensor.
  ///
  /// @param expr sequant ExprPtr to a sequant Tensor.
  /// @return EvalTensorPtr to a EvalTensor.
  EvalTensorPtr<DataTensorType> build_leaf(const ExprPtr& expr) const {
    // create an object to return
    auto leaf_tensor_ptr =
        std::make_shared<EvalTensorLeaf<DataTensorType>>(expr);

    auto& tnsr = expr->as<Tensor>();

    // if expr references to a real valued tensor check if swapping bra and ket
    // labels is required else set as not required
    bool swap_bk = false;
    if((tnsr.label() ==  L"t") || (tnsr.label() ==  L"f"))
       swap_bk = need_bra_ket_swap(expr);

    // get the labels of indices in bra and ket
    IndexContainer bra_index_labels, ket_index_labels;
    for (const auto& idx : tnsr.bra()) bra_index_labels.emplace_back(idx);
    for (const auto& idx : tnsr.ket()) ket_index_labels.emplace_back(idx);

    // swap labels as required
    if (swap_bk) std::swap(bra_index_labels, ket_index_labels);

    // collect bra and ket indices by appending ket indices to the bra index
    // container
    std::move(ket_index_labels.begin(), ket_index_labels.end(),
              std::back_inserter(bra_index_labels));

    // set indices of the object
    leaf_tensor_ptr->set_indices(bra_index_labels);
    // set hash value of the object
    leaf_tensor_ptr->set_hash_value(hash_leaf(expr, swap_bk));

    return leaf_tensor_ptr;
  }

  /// Hash leaf tensor.
  ///
  /// @param expr sequant ExprPtr to a sequant expression.
  /// @param swap_bra_ket_labels Hash kets before bras if true. Should be set
  /// true while hashing real-valued tensors.
  /// @return Hash value of the tensor based on its kind and the index space of
  /// its bra and ket indices.
  HashType hash_leaf(const ExprPtr& expr, bool swap_bra_ket_labels) const {
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

    if (swap_bra_ket_labels) {
      hash_idx_space(tnsr.ket());
      hash_idx_space(tnsr.bra());
    } else {
      hash_idx_space(tnsr.bra());
      hash_idx_space(tnsr.ket());
    }

    return leaf_hash_value;
  }

  /// Build binary evaluation intermediate tensor from two evaluation tensors.
  ///
  /// @param left_eval_expr The left evaluation tensor.
  /// @param right_eval_expr The right evaluation tensor.
  /// @param op The binary evaluation type eg. sum, product, antisymmetrization
  /// @return EvalTensorPtr to a EvalTensor.
  /// @throw domain_error when the indices do not match for sum type evaluations
  /// of two tensors.
  EvalTensorPtr<DataTensorType> build_intermediate(
      const EvalTensorPtr<DataTensorType>& left_eval_expr,
      const EvalTensorPtr<DataTensorType>& right_eval_expr,
      Operation op) const {
    if ((op == Operation::SUM) && (left_eval_expr->get_indices().size() !=
                                   right_eval_expr->get_indices().size())) {
      throw std::domain_error(
          "While summing two tensors, their number of indices must match!");
    }

    // create the intermediate
    auto imed = std::make_shared<EvalTensorIntermediate<DataTensorType>>();

    // set left and right tensors
    imed->set_left_tensor(left_eval_expr);
    imed->set_right_tensor(right_eval_expr);

    // set the operation type
    imed->set_operation(op);

    // fill the indices
    if (op == Operation::SUM) {
      if (!((left_eval_expr->is_leaf() || right_eval_expr->is_leaf()))) {
        auto limed =
            std::dynamic_pointer_cast<EvalTensorIntermediate<DataTensorType>>(
                left_eval_expr);
        auto rimed =
            std::dynamic_pointer_cast<EvalTensorIntermediate<DataTensorType>>(
                right_eval_expr);
        if ((limed->get_operation() == Operation::ANTISYMMETRIZE) &&
            (rimed->get_operation() == limed->get_operation())) {
          return build_intermediate(
              limed->get_left_tensor(),
              build_intermediate(limed->get_right_tensor(),
                                 rimed->get_right_tensor(), Operation::SUM),
              Operation::ANTISYMMETRIZE);
        }
      }
      imed->set_indices(left_eval_expr->get_indices());
    } else if (op == Operation::PRODUCT) {
      IndexContainer lindices = left_eval_expr->get_indices();
      IndexContainer rindices = right_eval_expr->get_indices();

      // indices of intermediate tensors are canonicalized by sorting them
      if (left_eval_expr->is_leaf())
        std::sort(lindices.begin(), lindices.end());

      if (right_eval_expr->is_leaf())
        std::sort(rindices.begin(), rindices.end());

      IndexContainer imed_indices;
      std::set_symmetric_difference(lindices.begin(), lindices.end(),
                                    rindices.begin(), rindices.end(),
                                    std::back_inserter(imed_indices));
      imed->set_indices(imed_indices);
    } else if ((op == Operation::ANTISYMMETRIZE) ||
               (op == Operation::SYMMETRIZE)) {
      imed->set_indices(left_eval_expr->get_indices());
      // nothing special needs to be done
      // the evaluator should be smart enough
      // on handling (anti)symmetrization type operations
    }

    // set the hash value
    imed->set_hash_value(hash_intermediate(imed));

    // intermediate is built
    return imed;
  }

  /// Hash intermediate tensor
  ///
  /// @param eval_expr sequant ExprPtr to a sequant expression.
  /// @return Hash value of the tensor based on its kind and the index space of
  /// its bra and ket indices.
  HashType hash_intermediate(
      const EvalTensorPtr<DataTensorType>& eval_expr) const {
    // intermediate is hashed by the type of operation(sum, prod, etc.)
    // and the contracting/non-contracting index positions of the left and
    // the right tensor indices.
    // index positions are size_t by default
    // to combine the hashes cast Operation type to size_t as well

    auto imed_ptr =
        std::dynamic_pointer_cast<EvalTensorIntermediate<DataTensorType>>(
            eval_expr);

    auto imed_operation = imed_ptr->get_operation();

    boost::hash<size_t> number_hasher;
    HashType imed_hash_value;

    // hashing the type of binary operation by casting operation into size_t
    imed_hash_value = number_hasher(static_cast<size_t>(imed_operation));

    // operation is of (anti)symmetrization type then its hash value is solely
    // determined by which indices are being permuted, this information is
    // carried by the left tensor -- the tensor with label 'A'(or 'P')
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
      auto index_hash_combiner =
          [&](const EvalTensorPtr<DataTensorType>& lr_tensor) {
            // iterate through the imed_ptr indices
            for (const auto& idx : imed_ptr->get_indices()) {
              // iterate through the indices of the passed left/right tensor
              // by a counter
              for (auto i = 0; i < lr_tensor->get_indices().size(); ++i) {
                // if the index of the imed_ptr matches with the lr_tensor
                // then time to break as we found the contracted indices
                // no need to look further
                if (idx.label() == lr_tensor->get_indices()[i].label()) {
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

  /// Should we swap the whole bra indices with the whole ket indices?
  ///
  /// When a sequant Tensor represents a real valued data tensor extra
  /// canonicalization is done.
  /// @param expr sequant ExprPtr to a sequant expression.
  /// @return True if swapping is necessary.
  /// @throw domain_error if bra rank and ket rank do not match.
  bool need_bra_ket_swap(const ExprPtr& expr) const {
    auto& tnsr = expr->as<Tensor>();

    // non-symmetric and bra-kets with invalid symmetries are handled
    if ((tnsr.braket_symmetry() == BraKetSymmetry::nonsymm) ||
        (tnsr.braket_symmetry() == BraKetSymmetry::invalid)) {
      return false;
    }

    // when working with complex tensors if it doesn't have symmetric bra-ket
    if (complex_tensor_data &&
        (tnsr.braket_symmetry() != BraKetSymmetry::symm)) {
      return false;
    }

    if (tnsr.bra_rank() != tnsr.ket_rank())
      throw std::domain_error(
          "Only equal bra_rank() and ket_rank() are handled for now!");
    //

    for (auto i = 0; i < tnsr.rank(); ++i) {
      // if bra index and ket index at the same positions have the same
      // index space then we don't know what to do yet, so
      if (tnsr.bra().at(i).space() == tnsr.ket().at(i).space()) continue;

      // for the corresponding positions if the index space value of the ket
      // index is lower than that of the bra index, then swapping should be done
      else if (tnsr.ket().at(i).space() < tnsr.bra().at(i).space())
        return true;
      else
        return false;
    }

    // survived this far implies no swapping is needed
    return false;
  }

};  // class EvalTensorBuilder<DataTensorType>

}  // namespace sequant::evaluate

#endif /* ifndef SEQUANT_EVALUATE_EVAL_TENSOR_BUILDER_HPP */
