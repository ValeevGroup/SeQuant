//
// Created by Bimal Gaudel on Jan 14, 2020
//

#ifndef SEQUANT_EVAL_EXPR_HPP
#define SEQUANT_EVAL_EXPR_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/tensor.hpp>

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>  // std::pair

namespace sequant::evaluate {

template <typename DataTensorType>  // DataTensorType eg: btas::Tensor<double>
using DTensor_Ptr = std::shared_ptr<DataTensorType>;

template <typename DataTensorType>
using Label_DTensorPtr_pair =
    std::pair<std::wstring, DTensor_Ptr<DataTensorType>>;

template <typename DataTensorType>
using Label_DTensorPtr_map =
    container::map<typename Label_DTensorPtr_pair<DataTensorType>::first_type,
                   typename Label_DTensorPtr_pair<DataTensorType>::second_type>;

namespace detail {

using ISpType_vec = container::svector<IndexSpace::Type>;

class CanonicalizedLabels {
  public:
  explicit CanonicalizedLabels(const ExprPtr& expr) : expr_{expr} {
    if (!expr->is<Tensor>())
      throw std::invalid_argument("ExprPtr should be to a Tensor");

    auto& tnsr = expr->as<Tensor>();
    bra_labels_ptr_ = std::make_shared<container::svector<std::wstring>>();
    bra_isptype_vec_ptr_ = std::make_shared<ISpType_vec>();
    ket_labels_ptr_ = std::make_shared<container::svector<std::wstring>>();
    ket_isptype_vec_ptr_ = std::make_shared<ISpType_vec>();
    for (const auto& b : tnsr.bra()) {
      bra_labels_ptr_->push_back(std::wstring(b.label()));
      bra_isptype_vec_ptr_->push_back(b.space().type());
    }
    for (const auto& k : tnsr.ket()) {
      ket_labels_ptr_->push_back(std::wstring(k.label()));
      ket_isptype_vec_ptr_->push_back(k.space().type());
    }
    canonicalize();
  }

  inline const auto& bra_labels() const {return bra_labels_ptr_;}
  inline const auto& ket_labels() const {return ket_labels_ptr_;}
  inline const auto& bra_space() const  {return bra_isptype_vec_ptr_;}
  inline const auto& ket_space() const  {return ket_isptype_vec_ptr_;}

  private:
  void canonicalize() {
    auto& tnsr = expr_->as<Tensor>();
    if ((tnsr.symmetry() != Symmetry::antisymm) &&
        (tnsr.braket_symmetry() != BraKetSymmetry::conjugate))
      throw std::domain_error(
          "Only Symmetry::antisymm and BraKetSymmetry::conjugate"
          "are handled for now!");
    if (tnsr.bra_rank() != tnsr.ket_rank())
      throw std::domain_error(
          "Only equal bra_rank() and ket_rank()"
          "are handled for now!");
    for (auto i = 0; i < tnsr.rank(); ++i) {
      if (tnsr.bra()[i].label()[0] == tnsr.ket()[i].label()[0])
        continue;  // found oo or vv
      else {
        if (!(tnsr.bra()[i] < tnsr.ket()[i])) {
          bra_isptype_vec_ptr_.swap(ket_isptype_vec_ptr_);
          bra_labels_ptr_.swap(ket_labels_ptr_);
        }
        break;
      }
    } // for
  } // canonicalize()

  std::shared_ptr<container::svector<std::wstring>> bra_labels_ptr_;
  std::shared_ptr<container::svector<std::wstring>> ket_labels_ptr_;
  std::shared_ptr<ISpType_vec> bra_isptype_vec_ptr_;
  std::shared_ptr<ISpType_vec> ket_isptype_vec_ptr_;
  const ExprPtr expr_;
};

template <typename DataTensorType>
class EvalExpr {
 public:
  EvalExpr<DataTensorType>()
      : dtensor_ptr_{nullptr}, isptype_vec_ptr_{nullptr} {}

  EvalExpr<DataTensorType>(const DTensor_Ptr<DataTensorType>& dtnsr_ptr,
                           const std::shared_ptr<ISpType_vec>& isptype_vec_ptr)
      : dtensor_ptr_{dtnsr_ptr}, isptype_vec_ptr_{isptype_vec_ptr} {}

 protected:
  DTensor_Ptr<DataTensorType> dtensor_ptr_;
  std::shared_ptr<ISpType_vec> isptype_vec_ptr_;
};

template <typename DataTensorType>
class EvalTensor : public EvalExpr<DataTensorType> {
 public:
  EvalTensor<DataTensorType>(
      const ExprPtr& tensor_ptr,
      const std::shared_ptr<Label_DTensorPtr_map<DataTensorType>>&
          label_tensor_map_ptr)
      : tensor_ptr_{tensor_ptr}, label_tensor_map_ptr_{label_tensor_map_ptr} {}

 private:
  ExprPtr tensor_ptr_;
  std::shared_ptr<Label_DTensorPtr_map<DataTensorType>> label_tensor_map_ptr_;
};

}  // namespace detail

}  // namespace sequant::evaluate

#endif  // SEQUANT_EVAL_EXPR_HPP
