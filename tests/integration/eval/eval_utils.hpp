//
// Created by Bimal Gaudel on 7/13/21.
//

#ifndef SEQUANT_EVAL_EVAL_UTILS_HPP
#define SEQUANT_EVAL_EVAL_UTILS_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

namespace sequant::eval {

template <typename R, typename F>
void cartesian_foreach(const std::vector<R>& rs, F f) {
  using container::svector;
  using It = decltype(std::begin(rs[0]));
  using T = typename R::value_type;
  svector<It> its, ends;
  for (const auto& r : rs) {
    its.push_back(std::begin(r));
    ends.push_back(std::end(r));
  }
  while (its.front() != ends.front()) {
    svector<T> s;
    s.reserve(its.size());
    for (auto& it : its) {
      s.push_back(*it);
    }
    f(s);
    size_t i = its.size();
    while (i > 0) {
      --i;
      ++its[i];
      if (i == 0) break;
      if (its[i] != ends[i]) break;
      its[i] = std::begin(rs[i]);
    }
  }
}

///
/// Maps the IndexSpace type of an index in the braket of a tensor
/// to nocc (for IndexSpace::active_occupied)
/// or nvirt (for IndexSpace::active_unoccupied)
///
/// \param tensor sequant::Tensor
/// \return View of an iterable with size_t-type elements.
///
auto range1_limits(sequant::Tensor const& tensor, size_t nocc, size_t nvirt) {
  auto isr = get_default_context().index_space_registry();
  static auto const ao = isr->retrieve(L"i");
  static auto const au = isr->retrieve(L"a");
  return ranges::views::transform(tensor.const_braket(),
                                  [nocc, nvirt](auto const& idx) {
                                    const auto& sp = idx.space();
                                    SEQUANT_ASSERT(sp == ao || sp == au);

                                    return sp == ao ? nocc : nvirt;
                                  });
}

template <typename Tensor_t>
void read_tensor(std::string_view fname, Tensor_t& tensor) {
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

#endif  // SEQUANT_EVAL_EVAL_UTILS_HPP
