#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/tensor.hpp>

#include <iostream>
#include <string>
#include <vector>

#include "mbpt.h"

namespace py = pybind11;

namespace sequant::python {

inline std::vector<Index> make_index(std::vector<std::wstring> labels) {
  std::vector<Index> index;
  for (auto l : labels) {
    index.push_back(Index{l});
  }
  return index;
}

std::shared_ptr<Tensor> make_tensor(std::wstring label,
                                    std::vector<std::wstring> bra,
                                    std::vector<std::wstring> ket,
                                    std::vector<std::wstring> auxiliary) {
  return std::make_shared<Tensor>(label, make_index(bra), make_index(ket),
                                  make_index(auxiliary));
}

std::shared_ptr<Constant> make_constant(py::float_ number) {
  return std::make_shared<Constant>(to_rational(number.cast<double>()));
}

inline ExprPtr pow(const ExprPtr &b, int n) {
  ExprPtr e = b->clone();
  // auto p = *e;
  for (int i = 1; i < n; ++i) {
    e = e * b;
  }
  return e;
}

void index_space_register(py::kwargs kwargs) {
  for (auto kw : kwargs) {
    IndexSpace::register_instance(py::cast<std::wstring>(kw.first),
                                  IndexSpace::occupied);
  }
}

py::object summands(ExprPtr &expr) {
  if (auto s = std::dynamic_pointer_cast<Sum>(expr)) {
    return py::cast(s->summands());
  }
  return py::none();
}

py::object factors(ExprPtr &expr) {
  if (auto p = std::dynamic_pointer_cast<Product>(expr)) {
    return py::cast(p->factors());
  }
  return py::none();
}

// disambiguates sequant::simplify
ExprPtr &simplify(ExprPtr &expr) { return sequant::simplify(expr); }

py::object rational_to_fraction(const rational &r) {
  py::object Fraction = py::module::import("fractions").attr("Fraction");
  return Fraction(numerator(r), denominator(r));
}

std::string complex_to_string(const Complex<rational> &z) {
  return z.imag() != 0 ? (to_string(z.real()) + "1j * " + to_string(z.imag()))
                       : to_string(z.real());
}

}  // namespace sequant::python

NAMESPACE_BEGIN(PYBIND11_NAMESPACE)
NAMESPACE_BEGIN(detail)

template <>
struct is_holder_type<sequant::Expr, sequant::ExprPtr> : std::true_type {};

template <>
struct always_construct_holder<sequant::ExprPtr>
    : always_construct_holder<void> {};
template <>
class type_caster<sequant::ExprPtr>
    : public type_caster_holder<sequant::Expr, sequant::ExprPtr> {};

NAMESPACE_END(detail)
NAMESPACE_END(PYBIND11_NAMESPACE)

PYBIND11_MODULE(_sequant, m) {
  using namespace sequant;
  using namespace sequant::python;

#define SEQUANT_PYTHON_INDEXSPACE_TYPE_PROPERTY(TYPE)                \
  .def_property_static(                                              \
      #TYPE, [](py::object /* self */) { return IndexSpace::TYPE; }, \
      [](py::object /* self */, std::wstring s) {                    \
        IndexSpace::register_instance(s, IndexSpace::TYPE);          \
      })

  py::class_<IndexSpace>(m, "IndexSpace")
      SEQUANT_PYTHON_INDEXSPACE_TYPE_PROPERTY(active_occupied)
          SEQUANT_PYTHON_INDEXSPACE_TYPE_PROPERTY(occupied)
              SEQUANT_PYTHON_INDEXSPACE_TYPE_PROPERTY(unoccupied);

  py::class_<ExprPtr>(m, "ExprPtr")
      .def_property_readonly("latex", &ExprPtr::to_latex)
      .def("__add__", [](const ExprPtr &l, const ExprPtr &r) { return l + r; })
      .def("__sub__", [](const ExprPtr &l, const ExprPtr &r) { return l - r; })
      .def("__mul__", [](const ExprPtr &l, const ExprPtr &r) { return l * r; });

  py::class_<Expr, ExprPtr>(m, "Expr")
      .def_property_readonly("summands", &summands)
      .def_property_readonly("factors", &factors)
      .def_property_readonly("latex", &Expr::to_latex)
      .def("__add__", [](const ExprPtr &l, const ExprPtr &r) { return l + r; })
      .def("__sub__", [](const ExprPtr &l, const ExprPtr &r) { return l - r; })
      .def("__mul__", [](const ExprPtr &l, const ExprPtr &r) { return l * r; })
      .def("__pow__", [](const ExprPtr &b, int n) { return pow(b, n); });

  py::class_<Index, std::shared_ptr<Index>>(m, "Index")
      .def("__str__", &Index::label)
      .def("__repr__", &Index::label)
      .def_property_readonly("space", &Index::space);

  py::class_<Tensor, std::shared_ptr<Tensor>, Expr>(m, "Tensor")
      .def(py::init(&python::make_tensor))
      .def_property_readonly("label", &Tensor::label)
      .def_property_readonly("bra", &Tensor::bra)
      .def_property_readonly("ket", &Tensor::ket)
      .def_property_readonly("auxikiary", &Tensor::auxiliary)
      .def_property_readonly("braket",
                             [](const Tensor &t) {
                               auto braket = t.braket();
                               return std::vector<Index>(braket.begin(),
                                                         braket.end());
                             })
      .def_property_readonly("indices", [](const Tensor &t) {
        auto indices = t.indices();
        return std::vector<Index>(indices.begin(), indices.end());
      });

  py::class_<Complex<rational>>(m, "zRational")
      .def_property_readonly("real",
                             [](const Complex<rational> &r) {
                               return rational_to_fraction(r.real());
                             })
      .def_property_readonly("imag",
                             [](const Complex<rational> &r) {
                               return rational_to_fraction(r.imag());
                             })
      .def_property_readonly("latex", &Complex<rational>::to_latex)
      .def("__str__", &python::complex_to_string)
      .def("__repr__", &python::complex_to_string);

  py::class_<Constant, std::shared_ptr<Constant>, Expr>(m, "Constant")
      .def(py::init(&python::make_constant));

  py::class_<Product, std::shared_ptr<Product>, Expr>(m, "Product")
      .def_property_readonly("scalar", &Product::scalar);

  py::class_<Sum, std::shared_ptr<Sum>, Expr>(m, "Sum");

  m.def("simplify", &sequant::python::simplify);

  python::mbpt::__init__(m.def_submodule("mbpt"));
}
