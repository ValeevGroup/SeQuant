//
// Created by Bimal Gaudel on 1/17/20.
//

#ifndef SEQUANT_EVAL_EXPR_HPP
#define SEQUANT_EVAL_EXPR_HPP

#include "eval_fwd.hpp"
#include "eval_tensor.hpp"

#include <SeQuant/core/expr_fwd.hpp>
#include <boost/functional/hash.hpp>

#ifdef SEQUANT_HAS_BTAS
#include <btas/btas.h>
#endif

#include <cstddef>
#include <iostream>
#include <memory>
#include <string>

namespace sequant::evaluate {

///
/// Specify the qualities of the data-tensors being calculated on.
///
/// Data-tensor examples are btas::Tensor<double>, TA::TArrayD.
///
/// @param label is a std::wstring and can be "g", "t", or "f"
///        as those are the labels used by sequant Tensor
///
/// @param spaces container of IndexSpace types. For example
///        if a order-two data-tensor is stored in such a way that its
///        first and second indices are occupied and virtual respectively
///        then {IndexSpace::active_occupied, IndexSpace::active_unoccupied}
///        could be a possible specification.
///
/// @method get_hash_value() returns hash_type by using the same hashing
///         technique as the EvalTensor constructor.
///
///        The only difference is, there the state of 'extra-canonicalization'
///        is deduced, i.e., for example, if bra-indices are all virtuals and
///        ket-indices are all occupied, then 'extra-canonicalization' swaps
///        bra<->ket index labels (as far as the sequant::Tensor has
///        Symmetry:antisymm and BraKetSymmetry:conjugate for this particular
///        example). However, no such deduction is made here. Party who is
///        calling this method is responsible for feeding the data-tensor with
///        expected IndexSpace ordering for a succesful evaluation. Reading
///        about how the 'extra-canonicalization' occurs while interpreting the
///        sequant::Tensor will be useful. Such a reading could be found
///        some place else under this project.
///
struct DataTensorSpecs {
  const std::wstring label;
  const container::svector<IndexSpace::Type> spaces;
  hash_type get_hash_value() const { return hash_value_; }

  DataTensorSpecs(const std::wstring& lbl,
                  const container::svector<IndexSpace::Type>& sp)
      : label{lbl}, spaces{sp} {
    boost::hash<std::wstring> hash_label;
    boost::hash<hash_type> hash_index_space;
    hash_value_ = hash_label(label);
    for (const auto& stype : spaces)
      boost::hash_combine(hash_value_, hash_index_space(stype));
  }

 private:
  hash_type hash_value_{0};
};

///
/// A priori information needed for evaluations.
///
template <typename DataTensorType>
struct EvalContext {
  ///
  /// Construct an evaluation context.
  /// @param hash_leaf_ptr_map A map of hash values of leaf tensors to the
  /// shared_ptr of their corresponding data tensors.
  /// @param hash_count_map A map of all tensor hashes to their number of
  /// appearance in the evaluation tree.
  ///
 public:
  EvalContext(const hash_to_dtensor_map<DataTensorType>& hash_leaf_ptr_map,
              const container::map<hash_type, std::size_t>& hash_count_map)
      : leaf_map_{hash_leaf_ptr_map} {
    for (const auto& item : hash_count_map) {
      if (item.second > 1) {
        imed_counts().insert(item);
        imed_map().insert(std::pair<hash_type, std::shared_ptr<DataTensorType>>(
            item.first, std::shared_ptr<DataTensorType>()));
      }
    }
  }

  explicit EvalContext(
      const hash_to_dtensor_map<DataTensorType>& hash_leaf_ptr_map)
      : leaf_map_{hash_leaf_ptr_map} {}

  /// Map of the hash value of leaf EvalTensor objects to
  /// their corresponding data-tensor pointers
  ///
  /// @return const hash_to_dtensor_map<DataTensorType>&
  const hash_to_dtensor_map<DataTensorType>& leaf_map() const {
    return leaf_map_;
  }

  /// Map of the hash value of intermediate EvalTensor objects to
  /// their corresponding data-tensor pointers
  ///
  /// @return const hash_to_dtensor_map<DataTensorType>&
  const hash_to_dtensor_map<DataTensorType>& imed_map() const {
    return imed_map_;
  }

  /// Mutable map of the hash value of intermediate EvalTensor objects to
  /// their corresponding data-tensor pointers
  ///
  /// @return hash_to_dtensor_map<DataTensorType>&
  hash_to_dtensor_map<DataTensorType>& imed_map() { return imed_map_; }

  /// Map of the hash value of intermediate EvalTensor objects to their
  /// corresponding counts.
  ///
  /// @note Only the intermediates that appear more than once are counted.
  const container::map<hash_type, std::size_t>& imed_counts() const {
    return imed_counts_;
  }

  /// Mutable map of the hash value of intermediate EvalTensor objects to their
  /// corresponding counts.
  container::map<hash_type, std::size_t>& imed_counts() { return imed_counts_; }

 private:
  using hash_to_count_map = container::map<hash_type, std::size_t>;

  hash_to_count_map imed_counts_;

  hash_to_dtensor_map<DataTensorType> leaf_map_;

  hash_to_dtensor_map<DataTensorType> imed_map_;
};

#ifdef SEQUANT_HAS_BTAS
template <typename DataTensorType>
DataTensorType eval_evtensor(const EvTensorPtr& evt_ptr,
                             EvalContext<DataTensorType>& context) {
  if (evt_ptr->is_leaf())
    return *context.leaf_map().find(evt_ptr->get_hash_value())->second;

  auto this_hash = evt_ptr->get_hash_value();
  auto this_hash_exists = context.imed_map().find(this_hash);

  if (this_hash_exists != context.imed_map().end()) {
    auto this_scalar = evt_ptr->get_scalar();
    assert(this_scalar.imag() == 0);

    // if the result already exists, return it
    // but, be wise, the pointer to the data tensor
    // could be a nullptr and scaling might need
    // to be done before return based on this
    // tensor's scalar
    //
    // pointer to data tensor
    auto dt_ptr = this_hash_exists->second;

    if (dt_ptr) {
      // data exists, return it by scaling if needed
      // also decrement the number of counts this data tensor is used
      //   - helpful to clean up when nobody needs the data
      context.imed_counts()[this_hash] -= 1;
      // if it is the last time the data is being called
      // move it to the local space and free it from the register
      if (context.imed_counts()[this_hash] == 0 || this_scalar.real() != 1) {
        DataTensorType result;

        // first case
        if (context.imed_counts()[this_hash] == 0) {
          result = std::move(*context.imed_map().find(this_hash)->second);
          context.imed_map()[this_hash] = std::shared_ptr<DataTensorType>();
        } else {  // part of the second case
          result = *(context.imed_map().find(this_hash)->second);
        }

        // second case
        if (this_scalar.real() != 1)
          btas::scal(this_scalar.real(), result);  // scaling

        return result;
      }
    } else {  // this is an useful intermediate but has no data yet
      auto result_ptr = std::make_shared<DataTensorType>();
      if (evt_ptr->get_op() == EvalTensor::Operation::Sum)
        *result_ptr = eval_evsum(evt_ptr, context);
      else
        *result_ptr = eval_evproduct(evt_ptr, context);
      context.imed_map()[this_hash] = result_ptr;

      if (this_scalar.real() == 1) return *result_ptr;
      auto result = *result_ptr;
      btas::scal(this_scalar.real(), result);
      return result;
    }
  }

  auto result = DataTensorType{};
  if (evt_ptr->get_op() == EvalTensor::Operation::Sum)
    result = eval_evsum(evt_ptr, context);
  else if (evt_ptr->get_op() == EvalTensor::Operation::Product)
    result = eval_evproduct(evt_ptr, context);
  else
    throw std::logic_error(
        "Only a Tensor, Sum or Product is handled from here.\n");
  return result;
}

template <typename DataTensorType>
DataTensorType eval_evsum(const EvTensorPtr& evt_ptr,
                          EvalContext<DataTensorType>& context) {
  auto lresult = eval_evtensor(evt_ptr->left_tensor(), context);
  auto rresult = eval_evtensor(evt_ptr->right_tensor(), context);

  btas::scal(evt_ptr->left_tensor()->get_scalar().real(), lresult);

  btas::scal(evt_ptr->right_tensor()->get_scalar().real(), rresult);

  return lresult + rresult;
}

template <typename DataTensorType>
DataTensorType eval_evproduct(const EvTensorPtr& evt_ptr,
                              EvalContext<DataTensorType>& context) {
  DataTensorType result{};
  auto lscalar = evt_ptr->left_tensor()->get_scalar();
  auto rscalar = evt_ptr->right_tensor()->get_scalar();

  assert((lscalar.imag() == 0) && (rscalar.imag() == 0));

  btas::contract(lscalar.real() * rscalar.real(),
                 eval_evtensor(evt_ptr->left_tensor(), context),
                 evt_ptr->left_tensor()->btas_indices(),
                 eval_evtensor(evt_ptr->right_tensor(), context),
                 evt_ptr->right_tensor()->btas_indices(), 0.0, result,
                 evt_ptr->btas_indices());

  return result;
}
#endif

void fill_hash_counts(const EvTensorPtr& evt_ptr,
                      container::map<hash_type, std::size_t>& counts_map) {
  if (evt_ptr->is_leaf()) return;

  counts_map[evt_ptr->get_hash_value()] += 1;

  fill_hash_counts(evt_ptr->left_tensor(), counts_map);
  fill_hash_counts(evt_ptr->right_tensor(), counts_map);
}

}  // namespace sequant::evaluate

#endif  // SEQUANT_EVAL_EXPR_HPP
