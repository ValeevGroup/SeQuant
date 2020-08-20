//
// Created by Bimal Gaudel on 8/5/20.
//

#ifndef SEQUANT_IMED_TENSOR_HPP
#define SEQUANT_IMED_TENSOR_HPP

#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/tensor_network.hpp>

#include "inferred_imed_data.hpp"

namespace sequant::factorize {

///
/// Represents the resulting tensor of a binary operation
/// between two Tensors, or two ImedTensors, or a Tensor
/// and an ImedTensor.
///
/// @author Bimal Gaudel
/// @version 05 August, 2020
///
class ImedTensor : public Tensor {
 private:
  bool is_sum_;

  ExprPtr left_;

  ExprPtr right_;

 public:
  ImedTensor() = default;

  virtual ~ImedTensor() = default;

  ImedTensor(ExprPtr left, ExprPtr right, const InferredImedData& data)
      : Tensor(L"I", data.bra_indices, data.ket_indices, data.symmetry,
               data.braket_symmetry, data.particle_symmetry),
        is_sum_{data.is_sum},
        left_{left},
        right_{right} {
    container::set<Expr::hash_type> operands_hash;

    //
    // ImedTensor::fill_operand_hash(this->left(), operands_hash);
    // ImedTensor::fill_operand_hash(this->right(), operands_hash);

    // auto hvalue = Expr::hash_type{};  // initialize hash value
    // boost::hash<Expr::hash_type> hash_value_hasher;
    // for (auto hh : operands_hash)
    //   boost::hash_combine(hvalue, hash_value_hasher(hh));

    // hashing this node
    // if (data.is_sum) {
    // } else {
    //   //
    // }
  }

  ExprPtr left() const { return left_; }

  ExprPtr right() const { return right_; }

  bool is_sum() const { return is_sum_; }

  static inline boost::hash<decltype(Tensor().label())> label_hasher{};

  static inline boost::hash<std::size_t> size_t_hasher{};

  static Expr::hash_type ispace_hasher(const Index& idx) {
    return size_t_hasher(static_cast<std::size_t>(idx.space().type()));
  }

  // hash and combine index spaces of Indexes from a container
  static Expr::hash_type ispace_container_hasher(
      const std::decay_t<decltype(Tensor().bra())>& cont) {
    auto result = Expr::hash_type{};
    for (const auto& idx : cont)
      boost::hash_combine(result, ImedTensor::ispace_hasher(idx));
    return result;
  }

  static Expr::hash_type ispace_based_tensor_hasher(const ExprPtr& expr) {
    auto& tnsr = expr->as<Tensor>();
    auto result = Expr::hash_type{};
    boost::hash_combine(result, label_hasher(tnsr.label()));
    boost::hash_combine(result,
                        ImedTensor::ispace_container_hasher(tnsr.bra()));
    boost::hash_combine(result,
                        ImedTensor::ispace_container_hasher(tnsr.ket()));
    return result;
  }

  ///
  /// @param expr ExprPtr to Tensor or ImedTensor
  ///
  static void fill_operand_hash(const ExprPtr& expr,
                                container::set<Expr::hash_type>& collect) {
    if (!expr->is<ImedTensor>()) {
      collect.insert(ImedTensor::ispace_based_tensor_hasher(expr));
      return;
    }
    // expr is ImedTensor
    auto& imed = expr->as<ImedTensor>();
    if (imed.is_sum()) {
      ImedTensor::fill_operand_hash(imed.left(), collect);
      ImedTensor::fill_operand_hash(imed.right(), collect);
    } else
      collect.insert(imed.hash_value());
  }
};

}  // namespace sequant::factorize

#endif  // SEQUANT_IMED_TENSOR_HPP
