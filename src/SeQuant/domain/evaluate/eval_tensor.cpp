//
// Created by Bimal Gaudel on 1/17/20.
//

#include "eval_tensor.hpp"
#include "eval_fwd.hpp"

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <boost/functional/hash.hpp>

#include <utility>

namespace sequant::evaluate {

EvTensorPtr BinaryOpTypeBuilder::operator()(EvTensorPtr& evt_ptr,
                                            const ExprPtr& expr_ptr) const {
  auto right_evtptr = std::make_shared<EvalTensor>(EvalTensor(expr_ptr));
  if ((this->op == EvalTensor::Operation::Sum) &&
      evt_ptr->indices() != right_evtptr->indices()){
    std::wcout << __LINE__ << " " << to_latex(expr_ptr) << " " << evt_ptr->indices().size() << " " << right_evtptr->indices().size() << std::endl;
    throw std::domain_error(
        "While summing two tensors, their indices should be identical");
  }

  // std::wcout << __LINE__ << " " << to_latex(expr_ptr) << " " << evt_ptr->indices().size() << " " << right_evtptr->indices().size() << std::endl;
  auto result = std::make_shared<EvalTensor>();

  result->set_left_tensor(evt_ptr);
  result->set_right_tensor(right_evtptr);
  result->set_op(this->op);

  if (this->op == EvalTensor::Operation::Product) {
    auto compare_labels = [](const std::wstring& left,
                             const std::wstring& right) {
      // i_1 < i_2. a_5 < a_6
      if (left[0] == right[0]) return left < right;
      // but occupied are smaller than virtuals
      // i_1 < a_1, which will be true if
      // we do lexicographic comparison of the
      // two and flip the result
      return right < left;
    };

    auto lindices = label_container_type{result->left_tensor()->indices()};
    auto rindices = label_container_type{result->right_tensor()->indices()};

    std::sort(lindices.begin(), lindices.end(), compare_labels);
    std::sort(rindices.begin(), rindices.end(), compare_labels);

    // the set of indices the resulting EvalTensor gets
    // is the symmetric difference between the set of
    // indices of the left and the right
    std::set_symmetric_difference(
        lindices.begin(), lindices.end(), rindices.begin(), rindices.end(),
        std::back_inserter(result->indices_mutable()), compare_labels);

  } else if (this->op == EvalTensor::Operation::Sum) {
    std::copy(result->left_tensor()->indices().begin(),
              result->left_tensor()->indices().end(),
              std::back_inserter(result->indices_mutable()));
  } else {
    // no other kind of binary evaluation is handled from here
    throw std::domain_error("Operation can be Sum or Product only!");
  }

  // if at this point indices() of the result is empty
  // this is a scalar product of tensors which is not
  // handled by Operation::Product
  if (result->indices().empty())
    throw std::domain_error(
        "The number of non-contracting indices cannot be zero for a "
        "Operation::Product type evaluation");
  // hashing
  boost::hash<hash_type> hash_numbers;

  // hashing the type of binary operation
  hash_type this_hash = hash_numbers(static_cast<hash_type>(result->get_op()));
  // hash the hash values of the left and right tensors
  hash_type lhash_value = result->left_tensor()->get_hash_value();
  hash_type rhash_value = result->left_tensor()->get_hash_value();
  if (lhash_value < rhash_value) {
    boost::hash_combine(this_hash, lhash_value);
    boost::hash_combine(this_hash, rhash_value);
  } else {
    boost::hash_combine(this_hash, rhash_value);
    boost::hash_combine(this_hash, lhash_value);
  }

  if (this->op == EvalTensor::Operation::Product) {
    auto combine_idx_pos_hash = [&result, &hash_numbers](
                                    const EvTensorPtr& ev_ptr,
                                    hash_type& hash_holder) {
      for (const auto& idx : result->indices()) {
        for (hash_type i = 0; i < ev_ptr->indices().size(); ++i) {
          if (idx == ev_ptr->indices()[i]) {
            boost::hash_combine(hash_holder, hash_numbers(i));
            break;
          }  // if
        }    // for
      }      // for
    };       // lambda

    // the order of tensors while taking a product doesn't
    // matter so hash it further
    hash_type hash_holder1{0}, hash_holder2{0};
    combine_idx_pos_hash(result->left_tensor(), hash_holder1);
    combine_idx_pos_hash(result->right_tensor(), hash_holder1);

    combine_idx_pos_hash(result->right_tensor(), hash_holder2);
    combine_idx_pos_hash(result->left_tensor(), hash_holder2);

    if (hash_holder1 < hash_holder2) {
      boost::hash_combine(this_hash, hash_holder1);
      boost::hash_combine(this_hash, hash_holder2);
    } else {
      boost::hash_combine(this_hash, hash_holder2);
      boost::hash_combine(this_hash, hash_holder1);
    }

    // Note: The scalar of EvalTensor NOT hashed.
    // Consequently, hashes of two such tensors
    // match even if the tensors themselves differ
    // only by their scalars
  }

  result->set_hash_value(this_hash);
  return result;
}  // BinaryOpTypeBuilder::operator()

EvalTensor::EvalTensor(const Tensor& tnsr) {
  if (tnsr.bra_rank() != tnsr.ket_rank())
    throw std::domain_error(
        "Only equal bra_rank() and ket_rank()"
        "are handled for now!");

  // TRUE if index labels DON'T appear in a CANONICAL order
  // based on their IndexSpace
  bool swap_bra_ket_labels = false;
  if ((tnsr.symmetry() != Symmetry::antisymm) &&
      (tnsr.braket_symmetry() != BraKetSymmetry::conjugate))
    throw std::domain_error(
        "Only Symmetry::antisymm and BraKetSymmetry::conjugate"
        "are handled for now!");
  for (auto i = 0; i < tnsr.rank(); ++i) {
    if (tnsr.bra()[i].label()[0] == tnsr.ket()[i].label()[0]) continue;
    if (!(tnsr.bra()[i].label() > tnsr.ket()[i].label())) {
      swap_bra_ket_labels = true;
      break;
    }  // if
  }    // for

  // get the labels of indices in bra and ket
  label_container_type bra_index_labels, ket_index_labels;
  for (const auto& idx : tnsr.bra())
    bra_index_labels.emplace_back(std::wstring(idx.label()));
  for (const auto& idx : tnsr.ket())
    ket_index_labels.emplace_back(std::wstring(idx.label()));

  if (swap_bra_ket_labels) std::swap(bra_index_labels, ket_index_labels);
  std::move(ket_index_labels.begin(), ket_index_labels.end(),
            std::back_inserter(bra_index_labels));
  ket_index_labels.clear();
  // copy the labels to this->indices_
  std::move(bra_index_labels.begin(), bra_index_labels.end(),
            std::back_inserter(indices_mutable()));
  // hashing
  boost::hash<std::wstring_view> hash_idx_label;
  boost::hash<hash_type> hash_numbers;

  hash_type this_hash = hash_idx_label(tnsr.label());

  using index_container_type = container::svector<Index, 4>;
  auto hash_idx_space = [&this_hash,
                         &hash_numbers](const index_container_type& idx_vec) {
    for (const auto& idx : idx_vec)
      boost::hash_combine(
          this_hash, hash_numbers(static_cast<hash_type>(idx.space().type())));
  };

  if (swap_bra_ket_labels) {
    hash_idx_space(tnsr.ket());
    hash_idx_space(tnsr.bra());
  } else {
    hash_idx_space(tnsr.bra());
    hash_idx_space(tnsr.ket());
  }

  set_hash_value(this_hash);
}  // constructor from Tensor

EvalTensor::EvalTensor(const ExprPtr& expr) {
  if (expr->is<Tensor>()){
    auto result = EvalTensor(expr->as<Tensor>());
    *this = result;
  }
  else if (expr->is<Sum>()) {
    auto& sum = expr->as<Sum>();
    auto init =
        std::make_shared<EvalTensor>(EvalTensor(*sum.summands().begin()));
    auto builder = BinaryOpTypeBuilder{EvalTensor::Operation::Sum};
    auto result = std::accumulate(sum.summands().begin() + 1, sum.summands().end(),
                             init, builder);
    *this = *result;
  } else if (expr->is<Product>()) {
    auto& prod = expr->as<Product>();
    auto init =
        std::make_shared<EvalTensor>(EvalTensor(*prod.factors().begin()));
    auto builder = BinaryOpTypeBuilder{EvalTensor::Operation::Product};
    auto result = std::accumulate(prod.factors().begin() + 1, prod.factors().end(),
                             init, builder);
    result->set_scalar(prod.scalar());
    *this = *result;
  } else
    throw std::domain_error("Can handle Tensor, Sum, or Product only!");
}

constant_type EvalTensor::get_scalar() const { return scalar_; }

void EvalTensor::set_scalar(constant_type c) { scalar_ = c; }

EvalTensor::Operation EvalTensor::get_op() const { return operation_; }

void EvalTensor::set_op(EvalTensor::Operation op) { operation_ = op; }

const label_container_type& EvalTensor::indices() const { return indices_; }

label_container_type& EvalTensor::indices_mutable() { return indices_; }

hash_type EvalTensor::get_hash_value() const { return hash_value_; }

void EvalTensor::set_hash_value(hash_type h) { hash_value_ = h; }

const EvTensorPtr& EvalTensor::left_tensor() const { return left_tensor_; }

void EvalTensor::set_left_tensor(EvTensorPtr& lt) { left_tensor_ = lt; }

const EvTensorPtr& EvalTensor::right_tensor() const { return right_tensor_; }

void EvalTensor::set_right_tensor(EvTensorPtr& rt) { right_tensor_ = rt; }

bool EvalTensor::is_leaf() const { return left_tensor() == nullptr; }

#ifdef SEQUANT_HAS_BTAS
void EvalTensor::fill_btas_indices() {
  if (! btas_indices().empty()) return;
  /* BTAS threw assertion error while using hashes as indices*/
  /*   auto label_to_hash = [](const std::wstring& str) { */
  /*     boost::hash<std::wstring> wstr_hash; */
  /*     return wstr_hash(str); */
  /*   }; */

  /*   for (const auto& idx : indices()) */
  /*     btas_indices_.emplace_back(label_to_hash(idx)); */

  // convert the index labels to ordinals of size_t type
  auto wstring_to_number = [](const std::wstring& str) {
    size_t result = 0;
    for (const auto& ch : str) result += size_t(ch);
    return result;
  };

  for (const auto& idx : indices())
    btas_indices_.emplace_back(wstring_to_number(idx));

  if (is_leaf()) return;

  left_tensor()->fill_btas_indices();
  right_tensor()->fill_btas_indices();
}

const btas_index_container& EvalTensor::btas_indices() const {
  return btas_indices_;
}
#endif

}  // namespace sequant::evaluate


