#include "read_tensor.hpp"
#include <sstream>

namespace sequant::eval {

input_header_t read_header(std::string_view fname) {
  auto ifs = std::ifstream{fname.data()};
  size_t rank = 0, nocc = 0, nvirt = 0;
  std::string header{};
  std::getline(ifs, header);
  auto hs = std::istringstream{header};
  hs >> rank;
  hs >> nocc;
  hs >> nvirt;
  return {rank, nocc, nvirt};
}

bool compatible_dims(std::string_view fname1, std::string_view fname2) {
  auto hf = read_header(fname1);
  auto hs = read_header(fname2);

  return hf.nocc == hs.nocc && hf.nvirt == hs.nvirt;
}

}  // namespace sequant::eval
