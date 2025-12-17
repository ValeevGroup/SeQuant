#ifndef SEQUANT_PYTHON_MBPT_H
#define SEQUANT_PYTHON_MBPT_H

#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/op.hpp>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>

#include "python.h"

#include <cstdint>

namespace sequant::python::mbpt {

template <class F>
auto make_sr_op(F f) {
  auto op = [f](std::int64_t Rank) { return f(Rank); };
  return op;
}

template <class... Args>
ExprPtr VacuumAverage(const ExprPtr& e, const Args&... args) {
  return sequant::mbpt::op::vac_av(e, args...);
}

#define SR_OP(OP) \
  #OP, [](std::int64_t Rank) { return sequant::mbpt::OP(Rank); }, py::arg("Bra")

inline void __init__(py::module m) {
  sequant::mbpt::load(sequant::mbpt::Convention::Minimal);
  sequant::TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  py::enum_<sequant::mbpt::OpType>(m, "OpType")
      .value("h", sequant::mbpt::OpType::h)
      .value("f", sequant::mbpt::OpType::f)
      .value("t", sequant::mbpt::OpType::t);

  m.def("F", &sequant::mbpt::F);
  m.def("H", &sequant::mbpt::H,
        "H(k = 2) returns a Hamiltonian operator with up to k-body terms",
        py::arg("k") = 2);

  m.def(SR_OP(A));
  m.def(SR_OP(T));
  m.def(SR_OP(t));

  m.def("VacuumAverage", &VacuumAverage<>);
  m.def("VacuumAverage",
        &VacuumAverage<std::vector<
            std::pair<sequant::mbpt::OpType, sequant::mbpt::OpType> > >);
}

}  // namespace sequant::python::mbpt

#endif /* SEQUANT_PYTHON_MBPT_H */
