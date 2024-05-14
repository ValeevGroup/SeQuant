//
// Created by Eduard Valeyev on 2019-03-26.
//

#ifndef SEQUANT_DOMAIN_MBPT_OP_IPP
#define SEQUANT_DOMAIN_MBPT_OP_IPP

#include <SeQuant/domain/mbpt/op.hpp>

namespace sequant {
namespace mbpt {

template <typename QuantumNumbers, Statistics S>
Operator<QuantumNumbers, S>::Operator() = default;

template <typename QuantumNumbers, Statistics S>
Operator<QuantumNumbers, S>::Operator(
    std::function<std::wstring_view()> label_generator,
    std::function<ExprPtr()> tensor_form_generator,
    std::function<void(QuantumNumbers&)> qn_action)
    : base_type(std::move(label_generator), std::move(tensor_form_generator)),
      qn_action_(std::move(qn_action)) {}

template <typename QuantumNumbers, Statistics S>
Operator<QuantumNumbers, S>::~Operator() = default;

template <typename QuantumNumbers, Statistics S>
QuantumNumbers Operator<QuantumNumbers, S>::operator()(
    const QuantumNumbers& qns) const {
  QuantumNumbers result(qns);
  this->apply_to(result);
  return result;
}

template <typename QuantumNumbers, Statistics S>
QuantumNumbers& Operator<QuantumNumbers, S>::apply_to(
    QuantumNumbers& qns) const {
  assert(qn_action_);
  if (is_vacuum(qns)) {  // action on vacuum is trivial ...
    qn_action_(qns);
  } else {  // action on a {operator. product of operators} = use Wick's theorem
    QuantumNumbers qns_this;
    qn_action_(qns_this);
    qns = combine(qns_this, qns);
  }
  return qns;
}

template <typename QuantumNumbers, Statistics S>
bool Operator<QuantumNumbers, S>::static_less_than(const Expr& that) const {
  assert(that.is<this_type>());
  auto& that_op = that.as<this_type>();

  // compare cardinal tensor labels first, then QN ranks
  auto& cardinal_tensor_labels = TensorCanonicalizer::cardinal_tensor_labels();
  const auto this_label = this->label();
  const auto that_label = that_op.label();
  if (this_label == that_label) return this->less_than_rank_of(that_op);
  const auto this_is_cardinal_it = ranges::find_if(
      cardinal_tensor_labels,
      [&this_label](const std::wstring& l) { return l == this_label; });
  const auto this_is_cardinal =
      this_is_cardinal_it != ranges::end(cardinal_tensor_labels);
  const auto that_is_cardinal_it = ranges::find_if(
      cardinal_tensor_labels,
      [&that_label](const std::wstring& l) { return l == that_label; });
  const auto that_is_cardinal =
      that_is_cardinal_it != ranges::end(cardinal_tensor_labels);
  if (this_is_cardinal && that_is_cardinal) {
    if (this_is_cardinal_it != that_is_cardinal_it)
      return this_is_cardinal_it < that_is_cardinal_it;
    else
      return this->less_than_rank_of(that_op);
  } else if (this_is_cardinal && !that_is_cardinal)
    return true;
  else if (!this_is_cardinal && that_is_cardinal)
    return false;
  else {  // !this_is_cardinal && !that_is_cardinal
    if (this_label == that_label)
      return this->less_than_rank_of(that_op);
    else
      return this_label < that_label;
  }
}

template <typename QuantumNumbers, Statistics S>
bool Operator<QuantumNumbers, S>::commutes_with_atom(const Expr& that) const {
  assert(that.is_cnumber() || that.is<this_type>());
  if (that.is_cnumber())
    return true;
  else {
    auto& that_op = that.as<this_type>();

    // if this has annihilators/creators in same space as that has
    // creator/annihilators return false

    auto delta_this = (*this)();
    auto delta_that = (that_op)();

    assert(this->size() % 2 == 0 && that.size() == this->size());

    return combine(delta_this, delta_that) == combine(delta_that, delta_this);
  }
}

template <typename QuantumNumbers, Statistics S>
void Operator<QuantumNumbers, S>::adjoint() {
  const auto dN = (*this)(QuantumNumbers{});
  using qnc_t = std::decay_t<decltype(dN)>;
  static_assert(std::is_same_v<QuantumNumbers, qnc_t>,
                "mbpt::Operator::adjoint: QuantumNumbers type mismatch");

  // grab label and update according to adjoint flag
  auto lbl = std::wstring(this->label());
  if (lbl.back() == sequant::adjoint_label) {
    assert(is_adjoint_);
    lbl.pop_back();
  } else {
    assert(!is_adjoint_);
    lbl.push_back(sequant::adjoint_label);
  }

  const auto tnsr = this->tensor_form();
  *this =
      Operator{[=]() -> std::wstring_view { return lbl; },  // label_generator
               [=]() -> ExprPtr {
                 return sequant::adjoint(tnsr);  // tensor_form_generator
               },
               [=](qnc_t& qn) {
                 qn += sequant::adjoint(dN);
                 return qn;  // qn_action
               }};
  this->is_adjoint_ = !this->is_adjoint_;  // toggle adjoint flag
}

template <typename QuantumNumbers, Statistics S>
bool Operator<QuantumNumbers, S>::less_than_rank_of(
    const this_type& that) const {
  return (*this)(QuantumNumbers{}) < that(QuantumNumbers{});
}

template <typename QuantumNumbers, Statistics S>
Expr::type_id_type Operator<QuantumNumbers, S>::type_id() const {
  return Expr::get_type_id<this_type>();
};

template <typename QuantumNumbers, Statistics S>
ExprPtr Operator<QuantumNumbers, S>::clone() const {
  return ex<this_type>(*this);
}

// Expresses general operators in human interpretable form. for example: \hat{T}_2 is a particle conserving 2-body excitation operator
// a non-particle conserving operator \hat{R}_2_1 implies that two particles are created followed by a single hole creation.
// conversely \hat{R}_1_2 implies the that only one particle is annihilated followed by two holes being created.
// The rule being, that for non-particle conserving operators, the first position indicates where the quasiparticle is going to and the second position indicates where it comes from.
// for the case of adjoint operators, the adjoint is represented by the symbol ⁺ and superscripting the quasi-particle numbers. for example: hat{R⁺}^{1,2}}
// For operators in which one or more quasi-particles has only partial coverage in the particle_space or hole_space, this notation is unsuitable, and we default to
// level printing of the operator.
template <typename QuantumNumbers, Statistics S>
std::wstring Operator<QuantumNumbers, S>::to_latex() const {
  return sequant::to_latex(*this);
}

template <typename QuantumNumbers, Statistics S>
Expr::hash_type Operator<QuantumNumbers, S>::memoizing_hash() const {
  using std::begin;
  using std::end;
  auto qns = (*this)(QuantumNumbers{});
  auto val = sequant::hash::value(qns);
  sequant::hash::combine(val, std::wstring(this->label()));
  this->hash_value_ = val;
  return *(this->hash_value_);
}

}  // namespace mbpt
}  // namespace sequant

#endif  // SEQUANT_DOMAIN_MBPT_OP_IPP
