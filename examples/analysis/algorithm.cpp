#include "algorithm.hpp"

#include <fmt/format.h>
#include <fstream>

namespace sequant {

container::set<Imed> intermediates(Sum const& sum) {
  auto ns =
      sum | ranges::views::transform(eval_node<EvalExpr>) | ranges::to_vector;
  return intermediates(ns);
}

container::vector<container::vector<size_t>> cse_graph(Sum const& sum) {
  using ranges::views::transform;
  auto const imeds = intermediates(sum);
  container::vector<container::vector<size_t>> result(
      sum.size(), std::vector<size_t>(sum.size(), 0));
  for (auto const& im : imeds) {
    auto const len = im.pos.size();
    if (len > 1) {
      auto x = 0;
      for (; x < len; ++x) {
        auto y = (x + 1) % len;
        result[im.pos[x]][im.pos[y]] = 1;
      }
    }
  }
  return result;
}

void write_cse_graph(
    std::filesystem::path const& file,
    container::vector<container::vector<size_t>> const& graph) {
  std::ofstream ofs{file};
  for (auto const& vec : graph)
    ofs << fmt::format(std::locale::classic(), "{}\n", fmt::join(vec, ","));
}

}  // namespace sequant
