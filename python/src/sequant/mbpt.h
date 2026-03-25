#ifndef SEQUANT_PYTHON_MBPT_H
#define SEQUANT_PYTHON_MBPT_H

#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/vac_av.hpp>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/string.hpp>

#include "python.h"

#include <cstdint>

namespace sequant::python::mbpt {

template <class F>
auto make_sr_op(F f) {
  auto op = [f](std::int64_t Rank) { return f(Rank); };
  return op;
}

using PyEVOptions = sequant::mbpt::EVOptions<std::string>;

ExprPtr VacuumAverage(const ExprPtr& e) { return sequant::mbpt::op::vac_av(e); }

ExprPtr VacuumAverage(const ExprPtr& e, const PyEVOptions& opts) {
  // helper for converting connections lists
  auto convert = [](const PyEVOptions::container_type& pairs) {
    sequant::mbpt::OpConnections<std::wstring> result;
    result.reserve(pairs.size());
    for (const auto& [a, b] : pairs) {
      result.emplace_back(sequant::toUtf16(a), sequant::toUtf16(b));
    }
    return result;
  };
  return sequant::mbpt::op::vac_av(e, {.connect = convert(opts.connect),
                                       .avoid = convert(opts.avoid),
                                       .screen = opts.screen,
                                       .use_topology = opts.use_topology,
                                       .skip_clone = opts.skip_clone});
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
  m.def(SR_OP(t));

  py::class_<PyEVOptions>(m, "EVOptions")
      .def(py::init<>())
      .def_readwrite("connect", &PyEVOptions::connect)
      .def_readwrite("avoid", &PyEVOptions::avoid)
      .def_readwrite("screen", &PyEVOptions::screen)
      .def_readwrite("use_topology", &PyEVOptions::use_topology)
      .def_readwrite("skip_clone", &PyEVOptions::skip_clone);

  m.def("VacuumAverage", py::overload_cast<const ExprPtr&>(&VacuumAverage),
        py::arg("expr"));
  m.def("VacuumAverage",
        py::overload_cast<const ExprPtr&, const PyEVOptions&>(&VacuumAverage),
        py::arg("expr"), py::arg("EVOptions"));
}

}  // namespace sequant::python::mbpt

#endif /* SEQUANT_PYTHON_MBPT_H */
