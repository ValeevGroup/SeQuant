#ifndef SEQUANT_PYTHON_MBPT_H
#define SEQUANT_PYTHON_MBPT_H

#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/sr/sr.hpp>
#include "python.h"

namespace sequant::python::mbpt {

template <class F>
auto make_sr_op(F f) {
  auto op = [f](size_t Bra) { return f(Bra); };
  return op;
}

template <class... Args>
ExprPtr VacuumAverage(const ExprPtr& e, const Args&... args) {
  return sequant::mbpt::sr::so::vac_av(e, args...);
}

#define SR_OP(OP) \
#OP, [](size_t Bra) { return sequant::mbpt::sr::so::OP(Bra); }, py::arg("Bra")

inline void __init__(py::module m) {
  sequant::mbpt::set_default_convention();
  sequant::set_default_context(
      SeQuant(Vacuum::SingleProduct, IndexSpaceMetric::Unit,
              BraKetSymmetry::conjugate, SPBasis::spinorbital));
  sequant::TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  m.def("F", &sequant::mbpt::sr::so::F);
  m.def("H", &sequant::mbpt::sr::so::H, py::arg("antisymmetric") = false);

  m.def(SR_OP(A));
  m.def(SR_OP(T));

  m.def("VacuumAverage", &VacuumAverage<>);
  m.def("VacuumAverage", &VacuumAverage<std::vector<std::pair<int, int> > >);
}

}  // namespace sequant::python::mbpt

#endif /* SEQUANT_PYTHON_MBPT_H */
