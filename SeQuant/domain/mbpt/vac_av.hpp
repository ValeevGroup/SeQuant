//
// Created by Eduard Valeyev on 2023-10-30.
//
using namespace sequant::mbpt;

template <typename T>
using OpConnections = std::vector<std::pair<T, T>>;

/// defines the default op connections
inline OpConnections<mbpt::OpType> default_op_connections() {
  using mbpt::OpType;
  static const OpConnections<mbpt::OpType> defaults = {
      {OpType::h, OpType::t},
      {OpType::f, OpType::t},
      {OpType::f̃, OpType::t},
      {OpType::g, OpType::t},
      // NBs
      // - for unitary ansatz can also connect t^+ with Hamiltonian
      // - for exact (non-approximated) unitary ansatz will also need to connect
      // t^+ with t;
      //   for MR unitary ansatz can also connect t with t and t^+ with t^+ ,
      //   but since adjoint() does not change OpType all of these are expressed
      //   as same connection ... this points out the need to have a separate
      //   OpType for t^+ and t, and in general the contents of OpType must be
      //   customizable
      {OpType::t, OpType::h},
      {OpType::t, OpType::f},
      {OpType::t, OpType::f̃},
      {OpType::t, OpType::g}};
  return defaults;
}

/// concat 2 sets of op connections
template <typename T,
          typename = std::enable_if_t<std::is_same_v<T, mbpt::OpType> ||
                                      std::is_same_v<T, std::wstring>>>
inline OpConnections<T> concat(const OpConnections<T>& connections1,
                               const OpConnections<T>& connections2) {
  return ranges::concat_view(connections1, connections2) | ranges::to_vector;
}

/// lowers representation of op connections from mbpt::OpType to labels
/// @param[in] op_connections vector of pairs of operators to be connected
/// @return vector of pairs of operator labels to be connected
inline OpConnections<std::wstring> to_label_connections(
    const OpConnections<mbpt::OpType>& op_connections) {
  // convert mbpt::OpType to std::wstring
  using mbpt::optype2label;
  OpConnections<std::wstring> op_connect_wstr;
  for (const auto& [op1, op2] : op_connections) {
    op_connect_wstr.emplace_back(optype2label.at(op1), optype2label.at(op2));
  }
  return op_connect_wstr;
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

/// computes the vacuum expectation value (VEV)

/// @param[in] expr input expression
/// @param[in] op_connections list of pairs of operators to be
/// connected; connections are left-to-right, i.e., pair `{opL,opR}` declares
/// that `opL` and `opR` are to be connected when `opR` precedes `opL`, i.e.
/// `opL` is to the left of `opR`;
/// e.g., `{{OpType::h, OpType::t}}` will ensure that each
/// `OpType::h` will be connected to at least one `OpType::t` on its right
/// (preceding it). The default list of connections is given
/// by default_op_connections() .
/// @param[in] use_topology if true, WickTheorem uses topological equivalence of
/// terms, default is true
/// @param[in] screen if true, expressions are screened before lowering to
/// Tensor level and calling WickTheorem, default is true
/// @param[in] skip_clone if true, will not clone the input expression
/// @return the VEV
ExprPtr vac_av(ExprPtr expr,
               const OpConnections<mbpt::OpType>& op_connections =
                   default_op_connections(),
               bool use_topology = true, bool screen = true,
               bool skip_clone = false);

/// computes the vacuum expectation value (VEV)

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
ExprPtr vac_av(ExprPtr expr, const OpConnections<std::wstring>& op_connections,
               bool use_topology = true, bool screen = true,
               bool skip_clone = false);
