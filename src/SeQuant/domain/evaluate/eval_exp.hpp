//
// Created by Bimal Gaudel on 1/14/20.
//

#ifndef SEQUANT_EVAL_FWD_HPP
#define SEQUANT_EVAL_FWD_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/tensor.hpp>

#include <memory>
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

using ISpaceType_vec = container::svector<IndexSpace::Type>;

template <typename DataTensorType>
class EvalTensor : public Expr {
 public:
  explicit EvalTensor<DataTensorType>(const ExprPtr& tnsr_ptr = nullptr)
      : dtensor_ptr_{tnsr_ptr} {
        auto& tnsr = tnsr_ptr->as<Tensor>();
      }

 private:
  ISpaceType_vec ispacetype_vec_;
  DTensor_Ptr<DataTensorType> dtensor_ptr_;
};

}  // namespace detail

}  // namespace sequant::evaluate

#endif  // SEQUANT_EVAL_FWD_HPP
