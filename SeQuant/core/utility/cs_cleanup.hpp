#ifndef SEQUANT_CORE_UTILITY_CS_CLEANUP_HPP
#define SEQUANT_CORE_UTILITY_CS_CLEANUP_HPP

#include <SeQuant/core/rational.hpp>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace sequant {

template <typename TArray>
TArray cleanup_tensor(const TArray& tensor, const std::string& tensor_name) {
  size_t total_rank = tensor.trange().rank();
  if (total_rank <= 4) {
    std::cout << tensor_name << " has rank <= 2 or 4, no cleanup needed."
              << std::endl;
    return tensor;
  }

  size_t occ_rank = total_rank / 2;  // virtual indices = occ_rank

  std::cout << "\n=== Starting " << tensor_name << " Cleanup (permute last "
            << occ_rank << " indices) ===" << std::endl;
  auto norm_before = TA::norm2(tensor);
  std::cout << "Initial " << tensor_name << " norm: " << norm_before
            << std::endl;

  TArray cleaned(tensor.world(), tensor.trange(), tensor.shape());
  cleaned.fill(0.0);

  std::string virt_indices_str;
  std::string occ_indices_str;
  virt_indices_str.reserve(5 * occ_rank);
  occ_indices_str.reserve(5 * occ_rank);
  for (size_t i = 1; i <= occ_rank; ++i) {
    virt_indices_str += "a_" + std::to_string(i);
    occ_indices_str += "i_" + std::to_string(i);
    if (i < occ_rank) {
      virt_indices_str += ",";
      occ_indices_str += ",";
    }
  }
  std::string left_index_str =
      virt_indices_str + "," + occ_indices_str;  // the main tensor

  // double factorial = 1.0;
  // for (size_t i = 1; i <= occ_rank; ++i) {
  //     factorial *= static_cast<double>(i);
  // }
  // double inv_factor = 1.0 / factorial;

  rational inv_factor = rational(1, factorial(occ_rank));
  double inv_factor_d = static_cast<double>(inv_factor);

  // create a vector to generate permutations of [1, 2, ..., occ_rank]
  std::vector<size_t> perm_indices(occ_rank);  // occ_rank is the size
  std::iota(perm_indices.begin(), perm_indices.end(),
            1);  // fill with 1..occ_ra

  TArray perm_sum(tensor.world(), tensor.trange(), tensor.shape());
  perm_sum.fill(0.0);

  // loop over all possible permutations of occupied index labels
  do {
    std::string permuted_occ_str;
    permuted_occ_str.reserve(5 * occ_rank);
    for (size_t j = 0; j < occ_rank; ++j) {
      permuted_occ_str += "i_" + std::to_string(perm_indices[j]);
      if (j < occ_rank - 1) {
        permuted_occ_str += ",";
      }
    }
    std::string perm_indices_str = virt_indices_str + "," + permuted_occ_str;
    perm_sum(left_index_str) += tensor(perm_indices_str);

  } while (std::next_permutation(perm_indices.begin(), perm_indices.end()));

  perm_sum(left_index_str) = perm_sum(left_index_str) * inv_factor_d;

  // cleaned = tensor - perm_sum(including identity)
  cleaned(left_index_str) = tensor(left_index_str) - perm_sum(left_index_str);

  auto norm_after = TA::norm2(cleaned);
  std::cout << "=== " << tensor_name << " Cleanup Complete ===" << std::endl;
  std::cout << "Norm change: " << (norm_after - norm_before) << std::endl;
  std::cout << "Final " << tensor_name << " norm: " << norm_after << std::endl
            << std::endl;

  return cleaned;
}

}  // namespace sequant

#endif  // SEQUANT_CORE_UTILITY_CS_CLEANUP_HPP
