//
// Created by Eduard Valeyev on 2019-03-26.
//

#ifndef SEQUANT_DOMAIN_MBPT_OP_HPP
#define SEQUANT_DOMAIN_MBPT_OP_HPP

#include <string>
#include <vector>

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/expr.hpp>

namespace sequant {
namespace mbpt {

/// enumerates the known Operator types
enum class OpType {
  h,    //!< 1-body Hamiltonian
  f,    //!< Fock operator
  g,    //!< 2-body Coulomb
  t,    //!< cluster amplitudes
  l,    //!< deexcitation cluster amplitudes
  A,    //!< antisymmetrizer
  L,    //!< left-hand eigenstate
  R,    //!< right-hand eigenstate
  R12,  //!< geminal kernel
  GR,   //!< GR kernel from f12 theory
  C     //!< cabs singles op
};

/// Operator character relative to Fermi vacuum
enum class OpClass { ex, deex, gen };

/// @return the tensor labels in the cardinal order
std::vector<std::wstring> cardinal_tensor_labels();

std::wstring to_wstring(OpType op);
OpClass to_class(OpType op);

//////////////////////////////

/// mbpt::Operator = Tensor * NormalOperator

/// It is an abstract (index-free) representation of a many-body operator
/// that tracks one or more quantum numbers (e.g., excitation level for
/// single-reference MBPT, particle numbers in each subspace for multi-reference
/// MBPT, etc.) and can be used for operator-level optimizations.
template <typename QuantumNumbers, Statistics S = Statistics::FermiDirac>
class Operator;

template <typename QuantumNumbers>
using FOperator = Operator<QuantumNumbers, Statistics::FermiDirac>;
template <typename QuantumNumbers>
using BOperator = Operator<QuantumNumbers, Statistics::BoseEinstein>;

/// Operator<void> does not track any quantum numbers, this is a base for
/// Operator
template <Statistics S>
class Operator<void, S> : public Expr {
 protected:
  Operator() = default;
  Operator(std::function<std::wstring()> label_generator,
           std::function<ExprPtr()> tensor_form_generator)
      : label_generator_(std::move(label_generator)),
        tensor_form_generator_(tensor_form_generator) {}

 public:
  virtual ~Operator() = default;

  virtual std::wstring label() const {
    assert(label_generator_);
    return label_generator_();
  }

  virtual ExprPtr tensor_form() const {
    assert(tensor_form_generator_);
    return tensor_form_generator_();
  }

 private:
  std::function<std::wstring()> label_generator_;
  std::function<ExprPtr()> tensor_form_generator_;
};

using FOperatorBase = FOperator<void>;
using BOperatorBase = BOperator<void>;

template <typename QuantumNumbers, Statistics S>
class Operator : public Operator<void> {
  using base_type = Operator<void>;

 protected:
  Operator() = default;

 public:
  Operator(std::function<std::wstring()> label_generator,
           std::function<ExprPtr()> tensor_form_generator,
           std::function<void(QuantumNumbers&)> qn_action)
      : base_type(std::move(label_generator), std::move(tensor_form_generator)),
        qn_action_(std::move(qn_action)) {}
  virtual ~Operator() = default;

  QuantumNumbers operator()(const QuantumNumbers& qns) const {
    QuantumNumbers result(qns);
    this->apply_to(result);
    return result;
  }

  virtual QuantumNumbers& apply_to(QuantumNumbers& qns) const {
    assert(qn_action_);
    qn_action_(qns);
    return qns;
  }

 private:
  std::function<void(QuantumNumbers&)> qn_action_;
};  // Operator

}  // namespace mbpt
}  // namespace sequant

#endif  // SEQUANT_DOMAIN_MBPT_OP_HPP
