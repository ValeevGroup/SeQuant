#ifndef SEQUANT_CORE_UTILITY_CS_CLEANUP_HPP
#define SEQUANT_CORE_UTILITY_CS_CLEANUP_HPP

#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/tensor.hpp>

#include <algorithm>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include <range/v3/range/conversion.hpp>
#include <range/v3/view.hpp>

namespace sequant {

template <typename TArray>
TArray cleanup_tensor(const TArray& array, const sequant::Tensor& tensor) {
  size_t total_rank = array.trange().rank();

  auto bra_rank = tensor.bra_rank();
  auto ket_rank = tensor.ket_rank();

  // there is no redundancy for r1 and r2, no need to apply hash filter
  if (total_rank <= 4) {
    return array;
  }

  auto bra_idx = ranges::views::iota(size_t{0}, bra_rank) | ranges::to_vector;
  auto ket_idx = ranges::views::iota(bra_rank, total_rank) | ranges::to_vector;

  auto ords_to_annot = [](auto const& ords) {
    using ranges::views::intersperse;
    using ranges::views::join;
    using ranges::views::transform;
    auto to_str = [](auto x) { return std::to_string(x); };
    return ords | transform(to_str) | intersperse(std::string{","}) | join |
           ranges::to<std::string>;
  };

  const auto l_annot =
      ords_to_annot(ranges::views::iota(size_t{0}, total_rank));

  TArray cleaned(array.world(), array.trange(), array.shape());
  cleaned.fill(0.0);

  TArray perm_sum(array.world(), array.trange(), array.shape());
  perm_sum.fill(0.0);

  rational norm_factor = rational(1, factorial(ket_rank));
  double norm_factor_d = static_cast<double>(norm_factor);

  std::vector<size_t> ket_perm = ket_idx;
  do {
    std::vector<size_t> full_perm;
    full_perm.reserve(total_rank);
    full_perm.insert(full_perm.end(), bra_idx.begin(), bra_idx.end());

    full_perm.insert(full_perm.end(), ket_perm.begin(), ket_perm.end());
    std::string perm_annot = ords_to_annot(full_perm);

    perm_sum(l_annot) += array(perm_annot);

  } while (std::next_permutation(ket_perm.begin(), ket_perm.end()));

  perm_sum(l_annot) = perm_sum(l_annot) * norm_factor_d;

  // cleaned = array - perm_sum
  cleaned(l_annot) = array(l_annot) - perm_sum(l_annot);

  return cleaned;
}

}  // namespace sequant

#endif  // SEQUANT_CORE_UTILITY_CS_CLEANUP_HPP
