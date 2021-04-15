#include "read_tensor_btas.hpp"

namespace sequant::example {

btas::Tensor<double> read_tensor_btas(std::string_view fname) {
  using ranges::views::repeat;
  using ranges::views::take;

  auto const header = read_header(fname);
  auto const nobs = header.nocc + header.nvirt;
  auto range =
      btas::Range{repeat(nobs) | take(header.rank) | ranges::to_vector};
  auto tensor = btas::Tensor<double>{range};
  read_tensor(fname, tensor);
  return tensor;
}
}  // namespace sequant::example
