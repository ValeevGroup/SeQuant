//
// Created by Eduard Valeyev on 2023-10-30.
//

#ifndef SEQUANT_DOMAIN_MBPT_VAC_AV_HPP
#define SEQUANT_DOMAIN_MBPT_VAC_AV_HPP

#include <SeQuant/domain/mbpt/op.hpp>

namespace sequant::mbpt {
inline namespace op {

/// defines the default op connections
inline OpConnections<std::wstring> default_op_connections() {
  static const OpConnections<std::wstring> defaults = {
      {L"h", L"t"},
      {L"f", L"t"},
      {L"f̃", L"t"},
      {L"g", L"t"},
      // NBs
      // - for unitary ansatz can also connect t^+ with Hamiltonian
      // - for exact (non-approximated) unitary ansatz will also need to connect
      // t^+ with t;
      //   for MR unitary ansatz can also connect t with t and t^+ with t^+ ,
      //   but since adjoint() does not change OpType all of these are expressed
      //   as same connection ... this points out the need to have a separate
      //   OpType for t^+ and t, and in general the contents of OpType must be
      //   customizable
      {L"t", L"h"},
      {L"t", L"f"},
      {L"t", L"f̃"},
      {L"t", L"g"}};
  return defaults;
}

/// concat 2 sets of op connections
template <typename T>
OpConnections<T> concat(const OpConnections<T>& a, const OpConnections<T>& b) {
  return ranges::concat_view(a, b) | ranges::to_vector;
}

/// @brief lowers an expression composed of Operators to tensor form
/// @param[in] expr input expression
/// @return expression with all Operators lowered to tensor form
/// @note mutates the input ExprPtr
inline ExprPtr lower_to_tensor_form(ExprPtr& expr) {
  auto op_lowerer = [](ExprPtr& leaf) {
    if (leaf.is<op_t>()) leaf = leaf.as<op_t>().tensor_form();
  };
  expr->visit(op_lowerer, /* atoms only = */ true);
  return expr;
}

///// @brief lowers an expression composed of Operators to tensor form
///// @param[in] expr_inp input expression
///// @return expression with all Operators lowered to tensor form
inline ExprPtr lower_to_tensor_form(const ExprPtr& expr_inp) {
  auto expr = expr_inp->clone();
  lower_to_tensor_form(expr);
  return expr;
}

// clang-format off
/// @brief computes the reference expectation value
/// @note equivalent to vac_av if the reference state is the Wick vacuum,
///       i.e. if `get_default_context().index_space_registry()->reference_occupied_space() == get_default_context().index_space_registry()->vacuum_occupied_space()`
/// @param[in] expr input expression
/// @param[in] opts controls the behavior, @see EVOptions
/// @note Uses `op::default_op_connections()` as default connectivity
// clang-format on
ExprPtr ref_av(ExprPtr expr, EVOptions<std::wstring> opts = {
                                 .connect = default_op_connections()});

/// @brief computes the vacuum expectation value
/// @internal evaluates only full contractions in  WickTheorem
/// @param[in] expr input expression
/// @param[in] opts controls the behavior, @see EVOptions
/// @note Uses `op::default_op_connections()` as default connectivity
/// @return the VEV
ExprPtr vac_av(ExprPtr expr, EVOptions<std::wstring> opts = {
                                 .connect = default_op_connections()});

}  // namespace op
}  // namespace sequant::mbpt
#endif  // SEQUANT_DOMAIN_MBPT_VAC_AV_HPP
