#ifndef SEQUANT_EXAMPLE_EVAL_READ_TENSOR_BTAS_HPP
#define SEQUANT_EXAMPLE_EVAL_READ_TENSOR_BTAS_HPP

#include "read_tensor.hpp"

#include <btas/btas.h>
#include <range/v3/view.hpp>

namespace sequant::eval {

btas::Tensor<double> read_tensor_btas(std::string_view fname);

}  // namespace sequant::example

#endif  // SEQUANT_EXAMPLE_EVAL_READ_TENSOR_BTAS_HPP
