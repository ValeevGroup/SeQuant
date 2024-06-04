//
// Created by Eduard Valeyev on 2023-10-30.
//
using namespace sequant::mbpt;

/// defines the default op connections
inline std::vector<std::pair<mbpt::OpType, mbpt::OpType>>
default_op_connections() {
  using mbpt::OpType;
  static const std::vector<std::pair<mbpt::OpType, mbpt::OpType>> defaults = {
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
inline std::vector<std::pair<mbpt::OpType, mbpt::OpType>> concat(
    const std::vector<std::pair<mbpt::OpType, mbpt::OpType>> connections1,
    const std::vector<std::pair<mbpt::OpType, mbpt::OpType>> connections2) {
  return ranges::concat_view(connections1, connections2) | ranges::to_vector;
}

/// lowers representation of op connections from mbpt::OpType to labels
/// @param[in] op_connections vector of pairs of operators to be connected
/// @return vector of pairs of operator labels to be connected
inline std::vector<std::pair<std::wstring, std::wstring>> to_label_connections(
    const std::vector<std::pair<mbpt::OpType, mbpt::OpType>>& op_connections) {
  // convert mbpt::OpType to std::wstring
  using mbpt::optype2label;
  std::vector<std::pair<std::wstring, std::wstring>> op_connect_wstr;
  for (const auto& [op1, op2] : op_connections) {
    op_connect_wstr.emplace_back(optype2label.at(op1), optype2label.at(op2));
  }
  return op_connect_wstr;
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
/// @param[in] skip_clone if true, will not clone the input expression
/// @return the VEV
ExprPtr vac_av(ExprPtr expr,
               std::vector<std::pair<mbpt::OpType, mbpt::OpType>>
                   op_connections = default_op_connections(),
               bool skip_clone = false);

/// computes the vacuum expectation value (VEV)

/// @param[in] expr input expression
/// @param[in] op_connections list of pairs of operator labels to be
/// connected; connections are defined left-to-right, i.e., pair
/// `{opL,opR}` declares that `opL` and `opR` are to be connected
/// when `opR` precedes `opL`, i.e. `opL` is to the left of `opR`
/// @param[in] skip_clone if true, will not clone the input expression
/// @return the VEV
ExprPtr vac_av(
    ExprPtr expr,
    std::vector<std::pair<std::wstring, std::wstring>> op_connections,
    bool skip_clone = false);
