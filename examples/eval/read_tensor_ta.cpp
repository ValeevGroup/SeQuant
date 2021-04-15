#include "read_tensor_ta.hpp"
#include "read_tensor.hpp"

#include <range/v3/view.hpp>
#include <tiledarray.h>

namespace sequant::example {

TA::TiledRange make_trange(size_t rank, size_t nobs) {
  using ranges::views::repeat;
  using ranges::views::take;
  auto tr1s = repeat(TA::TiledRange1{0, nobs}) | take(rank) | ranges::to_vector;
  return TA::TiledRange{tr1s.begin(), tr1s.end()};
}

void read_tensor_ta(std::string_view fname, TA::TArrayD &tensor) {
  // TODO assert tensor single tiled
  auto ta_tensor = TA::Tensor<double>{tensor.trange().make_tile_range(0)};
  read_tensor(fname, ta_tensor);
  *tensor.begin() = ta_tensor;
}

} // namespace sequant::example
