#ifndef SEQUANT_PYTHON_MBPT_H
#define SEQUANT_PYTHON_MBPT_H

#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/vac_av.hpp>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/wstring.hpp>

#include "python.h"

#include <cstdint>

namespace sequant::python::mbpt {

template <class F>
auto make_sr_op(F f) {
  auto op = [f](std::int64_t Rank) { return f(Rank); };
  return op;
}

// Overload for no op_connections
ExprPtr VacuumAverage(const ExprPtr& e) { return sequant::mbpt::op::vac_av(e); }

// overload  with string conversion
ExprPtr VacuumAverage(
    const ExprPtr& e,
    const std::vector<std::pair<std::string, std::string>>& op_connections) {
  sequant::mbpt::OpConnections<std::wstring> wop_connections;
  wop_connections.reserve(op_connections.size());
  for (const auto& [op1, op2] : op_connections) {
    wop_connections.emplace_back(sequant::to_wstring(op1),
                                 sequant::to_wstring(op2));
  }
  return sequant::mbpt::op::vac_av(e, wop_connections);
}

#define SR_OP(OP) \
  #OP, [](std::int64_t Rank) { return sequant::mbpt::OP(Rank); }, py::arg("Bra")

inline void __init__(py::module m) {
  sequant::mbpt::load(sequant::mbpt::Convention::Minimal);
  sequant::TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  // mbpt context setup
  sequant::mbpt::Context::Options opts;
  opts.op_registry_ptr = sequant::mbpt::make_legacy_registry();
  sequant::mbpt::set_default_mbpt_context(opts);

  m.def("F", &sequant::mbpt::F);
  m.def("H", &sequant::mbpt::H,
        "H(k = 2) returns a Hamiltonian operator with up to k-body terms",
        py::arg("k") = 2);

  m.def(SR_OP(A));
  m.def(SR_OP(T));
  m.def(SR_OP(T_));

  m.def("VacuumAverage", py::overload_cast<const ExprPtr&>(&VacuumAverage),
        py::arg("expr"));
  m.def("VacuumAverage",
        py::overload_cast<
            const ExprPtr&,
            const std::vector<std::pair<std::string, std::string>>&>(
            &VacuumAverage),
        py::arg("expr"), py::arg("op_connections"));
}

}  // namespace sequant::python::mbpt

#endif /* SEQUANT_PYTHON_MBPT_H */
