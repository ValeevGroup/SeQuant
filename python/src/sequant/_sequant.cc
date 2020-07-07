#include <SeQuant/core/sequant.hpp>
#include <SeQuant/core/tensor.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/shared_ptr.h>
#include <pybind11/operators.h>
#include <pybind11/boost/container/small_vector.h>

#include <iostream>


namespace py = pybind11;

namespace sequant::python {

  inline std::vector<Index> make_index(std::vector<std::wstring> labels) {
    std::vector<Index> index;
    for (auto l : labels) {
      index.push_back(Index{l});
    }
    return index;
  }

  std::shared_ptr<Tensor> make_tensor(
    std::wstring label,
    std::vector<std::wstring> bra,
    std::vector<std::wstring> ket)
  {
    return std::make_shared<Tensor>(label, make_index(bra), make_index(ket));
  }

  void index_space_register(py::kwargs kwargs) {
    for (auto kw : kwargs) {
      IndexSpace::register_instance(
        py::cast<std::wstring>(kw.first),
        IndexSpace::occupied
      );
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

}

PYBIND11_MODULE(_sequant, m) {

  using namespace sequant;
  using namespace sequant::python;

#define SEQUANT_PYTHON_INDEXSPACE_TYPE_PROPERTY(TYPE)            \
  .def_property_static(                                          \
    #TYPE,                                                       \
    [](py::object) { return IndexSpace::occupied; },             \
    [](py::object, std::wstring s) {                             \
      IndexSpace::register_instance(s, IndexSpace::TYPE);        \
    }                                                            \
  )

  py::class_<IndexSpace>(m, "IndexSpace")
    SEQUANT_PYTHON_INDEXSPACE_TYPE_PROPERTY(occupied)
    ;

  py::class_<Expr, ExprPtr>(m, "Expr")
    .def_property_readonly("summands", &summands)
    .def_property_readonly("factors", &factors)
    .def_property_readonly("latex", &Expr::to_latex)
    .def("__add__", [](const ExprPtr &l, const ExprPtr &r) { return l+r; })
    .def("__sub__", [](const ExprPtr &l, const ExprPtr &r) { return l-r; })
    .def("__mul__", [](const ExprPtr &l, const ExprPtr &r) { return l*r; })
    ;

  py::class_<Tensor, std::shared_ptr<Tensor>, Expr>(m, "Tensor")
    .def(py::init(&python::make_tensor))
    ;

  py::class_<Product, std::shared_ptr<Product>, Expr>(m, "Product");

  py::class_<Sum, std::shared_ptr<Sum>, Expr>(m, "Sum");

}
