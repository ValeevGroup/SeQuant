#ifndef SEQUANT_EXAMPLE_EVAL_READ_TENSOR_TA_HPP
#define SEQUANT_EXAMPLE_EVAL_READ_TENSOR_TA_HPP

#include "read_tensor.hpp"

#include <tiledarray.h>

namespace sequant::eval {

void read_tensor_ta(std::string_view fname, TA::TArrayD &tensor);

TA::TiledRange make_trange(size_t rank, size_t nobs);

} // namespace sequant::example

#endif // SEQUANT_EXAMPLE_EVAL_READ_TENSOR_TA_HPP
