#ifndef SEQUANT_EVAL_EXAMPLE_READ_TENSOR_HPP
#define SEQUANT_EVAL_EXAMPLE_READ_TENSOR_HPP

#include <algorithm>
#include <array>
#include <fstream>
#include <string_view>

namespace sequant::eval {

struct input_header_t {
  size_t rank;
  size_t nocc;
  size_t nvirt;
};

input_header_t read_header(std::string_view fname);

bool compatible_dims(std::string_view fname1, std::string_view fname2);

template <typename Tensor_t>
void read_tensor(std::string_view fname, Tensor_t &tensor) {
  auto ifs = std::ifstream{fname.data()};
  // omit header line
  std::string header{};
  std::getline(ifs, header);
  header.clear();
  //
  // fill data to tensor
  double x;
  std::generate(std::begin(tensor), std::end(tensor), [&ifs, &x]() {
    ifs >> x;
    return x;
  });  // generate
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_EXAMPLE_READ_TENSOR_HPP
