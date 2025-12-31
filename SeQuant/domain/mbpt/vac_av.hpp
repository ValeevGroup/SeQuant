//
// Created by Eduard Valeyev on 2023-10-30.
//

#ifndef SEQUANT_DOMAIN_MBPT_VAC_AV_HPP
#define SEQUANT_DOMAIN_MBPT_VAC_AV_HPP

#include <SeQuant/domain/mbpt/op.hpp>

namespace sequant::mbpt {
inline namespace op {

template <typename T>
using OpConnections = std::vector<std::pair<T, T>>;

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
inline OpConnections<std::wstring> concat(
    const OpConnections<std::wstring>& connections1,
    const OpConnections<std::wstring>& connections2) {
  return ranges::concat_view(connections1, connections2) | ranges::to_vector;
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

ExprPtr expectation_value_impl(
    ExprPtr expr, const OpConnections<std::wstring>& op_connections,
    bool use_topology, bool screen, bool skip_clone, bool full_contractions);

// clang-format off
/// @brief computes the reference expectation value
/// @note equivalent to vac_av if the reference state is the Wick vacuum,
///       i.e. if `get_default_context().index_space_registry()->reference_occupied_space() == get_default_context().index_space_registry()->vacuum_occupied_space()`
/// @param[in] expr input expression
/// @param[in] op_connections list of pairs of operator labels to be
/// connected; connections are defined left-to-right, i.e., pair
/// `{opL,opR}` declares that `opL` and `opR` are to be connected
/// when `opR` precedes `opL`, i.e. `opL` is to the left of `opR`
/// @param[in] use_topology if true, WickTheorem uses topological equivalence of
/// terms, default is true
/// @param[in] screen if true, expressions are screened before lowering to
/// Tensor level and calling WickTheorem, default is true
/// @param[in] skip_clone if true, will not clone the input expression
/// @return the reference expectation value
// clang_format on
ExprPtr ref_av(ExprPtr expr, const OpConnections<std::wstring>& op_connections = default_op_connections(),
               bool use_topology = true, bool screen = true,
               bool skip_clone = false);

/// @brief computes the vacuum expectation value
/// @internal evaluates only full contractions in  WickTheorem
/// @param[in] expr input expression
/// @param[in] op_connections list of pairs of operator labels to be
/// connected; connections are defined left-to-right, i.e., pair
/// `{opL,opR}` declares that `opL` and `opR` are to be connected
/// when `opR` precedes `opL`, i.e. `opL` is to the left of `opR`
/// @param[in] use_topology if true, WickTheorem uses topological equivalence of
/// terms, default is true
/// @param[in] screen if true, expressions are screened before lowering to
/// Tensor level and calling WickTheorem, default is true
/// @param[in] skip_clone if true, will not clone the input expression
/// @return the VEV
ExprPtr vac_av(ExprPtr expr, const OpConnections<std::wstring>& op_connections = default_op_connections(),
               bool use_topology = true, bool screen = true,
               bool skip_clone = false);

}  // namespace op
}  // namespace sequant::mbpt
#endif  // SEQUANT_DOMAIN_MBPT_VAC_AV_HPP
