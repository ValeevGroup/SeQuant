//
// Created by Eduard Valeyev on 2019-02-27.
//

#ifndef SEQUANT_SPIN_HPP
#define SEQUANT_SPIN_HPP

#include <algorithm>

#include <codecvt>
#include <map>
#include <string>
#include <tuple>

#include "SeQuant/core/expr.hpp"
#include "SeQuant/core/tensor.hpp"

#define SPINTRACE_PRINT 0
#define PRINT_PERMUTATION_LIST 0

namespace sequant {

class spinIndex : public Index {
 private:
};

struct zero_result : public std::exception {};

inline bool tensor_symm(const Tensor& tensor) {
  bool result = false;
  assert(tensor.bra().size() == tensor.ket().size());
  // For each index check if QNS match.
  // return FALSE if one of the pairs does not match.
  auto iter_ket = tensor.ket().begin();
  for (auto&& bra : tensor.bra()) {
    if (IndexSpace::instance(bra.label()).qns() ==
        IndexSpace::instance(iter_ket->label()).qns()) {
      result = true;
    } else {
      return result;
    }
    iter_ket++;
  }
  return result;
}
//================================================
#if 1
// These functions are from:
// https://www.geeksforgeeks.org/number-swaps-sort-adjacent-swapping-allowed/
// They are required to calculate the number of inversions or adjacent swaps
// performed to get a (-1)^n sign for a permuted tensor
// TODO: Use STL containers and simplify for current application
int merge(int arr[], int temp[], int left, int mid, int right) {
  auto inv_count = 0;

  auto i = left; /* i is index for left subarray*/
  auto j = mid;  /* i is index for right subarray*/
  auto k = left; /* i is index for resultant merged subarray*/
  while ((i <= mid - 1) && (j <= right)) {
    if (arr[i] <= arr[j])
      temp[k++] = arr[i++];
    else {
      temp[k++] = arr[j++];

      /* this is tricky -- see above explanation/
        diagram for merge()*/
      inv_count = inv_count + (mid - i);
    }
  }

  /* Copy the remaining elements of left subarray
   (if there are any) to temp*/
  while (i <= mid - 1) temp[k++] = arr[i++];

  /* Copy the remaining elements of right subarray
   (if there are any) to temp*/
  while (j <= right) temp[k++] = arr[j++];

  /*Copy back the merged elements to original array*/
  for (i = left; i <= right; i++) arr[i] = temp[i];

  return inv_count;
}

int _mergeSort(int arr[], int temp[], int left, int right) {
  int mid, inv_count = 0;
  if (right > left) {
    /* Divide the array into two parts and call
      _mergeSortAndCountInv() for each of the parts */
    mid = (right + left) / 2;

    /* Inversion count will be sum of inversions in
       left-part, right-part and number of inversions
       in merging */
    inv_count = _mergeSort(arr, temp, left, mid);
    inv_count += _mergeSort(arr, temp, mid + 1, right);

    /*Merge the two parts*/
    inv_count += merge(arr, temp, left, mid + 1, right);
  }

  return inv_count;
}

int countSwaps(int arr[], int n) {
  int temp[n];
  return _mergeSort(arr, temp, 0, n - 1);
}
#endif
//================================================
// TODO: This MUST return an ExprPtr
ExprPtr expand_antisymm(const Tensor& tensor) {
  // bool result = false;
  assert(tensor.bra().size() == tensor.ket().size());

  // Generate a sum of asymmetric tensors if the input tensor is antisymmetric
  // AND more than one body otherwise, return the tensor
  if ((tensor.symmetry() == Symmetry::antisymm) && (tensor.bra().size() > 1)) {
    // auto n = tensor.ket().size();

    container::set<Index> bra_list;
    for (auto&& bra_idx : tensor.bra()) bra_list.insert(bra_idx);
    const auto const_bra_list = bra_list;

    container::set<Index> ket_list;
    for (auto&& ket_idx : tensor.ket()) ket_list.insert(ket_idx);

    Sum expr_sum{};
    auto p_count = 0;
    do {
      int permutation_int_array[bra_list.size()];
      auto counter_for_array = 0;

      // Store distance of each index in an array
      ranges::for_each(const_bra_list, [&](Index i) {
        auto pos = std::find(bra_list.begin(), bra_list.end(), i);
        int dist = std::distance(bra_list.begin(), pos);
        permutation_int_array[counter_for_array] = dist;
        counter_for_array++;
      });

      // Call function to count number of pair swaps
      auto permutation_count =
          countSwaps(permutation_int_array, sizeof(permutation_int_array) /
                                                sizeof(*permutation_int_array));

      // Generate tensor with new labels
      auto new_tensor = Tensor(tensor.label(), bra_list, ket_list);
      ExprPtr new_tensor_ptr = std::make_shared<Tensor>(new_tensor);

      // Tensor as product with (-1)^n, where n is number of adjacent swaps
      Product new_tensor_product{};
      new_tensor_product.append(std::pow(-1, permutation_count),
                                new_tensor_ptr);

      ExprPtr new_tensor_product_ptr =
          std::make_shared<Product>(new_tensor_product);
      expr_sum.append(new_tensor_product_ptr);
      p_count++;
    } while (std::next_permutation(bra_list.begin(), bra_list.end()));

    ExprPtr result = std::make_shared<Sum>(expr_sum);
    // std::wcout << tensor.to_latex() << ":\n" << result->to_latex() << "\n" << std::endl;
    return result;
  } else {
    ExprPtr result = std::make_shared<Tensor>(tensor);
    // std::wcout << tensor.to_latex() << ":\n" << result->to_latex() << "\n" << std::endl;
    return result;
  }
}

/// @param expr input expression
/// @param ext_index_groups groups of external indices
/// @return the expression with spin integrated out
ExprPtr spintrace(ExprPtr expr,
                  std::initializer_list<IndexList> ext_index_groups = {{}}) {
  // SPIN TRACE DOES NOT SUPPORT PROTO INDICES YET.
  auto check_proto_index = [&](const ExprPtr& expr) {
    if (expr->is<Tensor>()) {
      ranges::for_each(
          expr->as<Tensor>().const_braket(),
          [&](const Index& idx) { assert(!idx.has_proto_indices()); });
    }
  };
  expr->visit(check_proto_index);

  container::set<Index, Index::LabelCompare> grand_idxlist;
  auto collect_indices = [&](const ExprPtr& expr) {
    if (expr->is<Tensor>()) {
      ranges::for_each(expr->as<Tensor>().const_braket(),
                       [&](const Index& idx) { grand_idxlist.insert(idx); });
    }
  };
  expr->visit(collect_indices);

  //  std::cout << "grand_idxlist size: " << grand_idxlist.size() << std::endl;

  container::set<Index> ext_idxlist;
  for (auto&& idxgrp : ext_index_groups) {
    for (auto&& idx : idxgrp) {
      auto result = ext_idxlist.insert(idx);
      assert(result.second);
    }
  }

  container::set<Index> int_idxlist;
  for (auto&& gidx : grand_idxlist) {
    if (ext_idxlist.find(gidx) == ext_idxlist.end()) {
      int_idxlist.insert(gidx);
    }
  }

  // EFV: generate the grand list of index groups by concatenating list of
  // external index EFV: groups with the groups of internal indices (each
  // internal index = 1 group)
  using IndexGroup = container::svector<Index>;
  container::svector<IndexGroup> index_groups;

  for (auto&& i : int_idxlist) {
    index_groups.emplace_back(IndexGroup(1, i));
  }

  index_groups.insert(index_groups.end(), ext_index_groups.begin(),
                      ext_index_groups.end());

  //  for (auto&& i : index_groups)
  //    for (auto&& j : i) std::wcout << j.label() << " ";
  //  std::cout << std::endl;

  // EFV: for each spincase (loop over integer from 0 to 2^(n)-1, n=#of index
  // groups)
  // const uint64_t nspincases = std::pow(2, index_groups.size());
  const uint64_t nspincases = 5;  // TODO: remove this line before release

  Sum spin_expr_sum{};
  auto total_terms = 0;
  const auto tstart = std::chrono::high_resolution_clock::now();
  for (uint64_t spincase_bitstr = 0; spincase_bitstr != nspincases;
       ++spincase_bitstr) {
    // EFV:  assign spin to each index group => make a replacement list
    std::map<Index, Index> index_replacements;
    std::map<Index, Index> remove_spin_replacements;
    int64_t index_group_count = 0;
    for (auto&& index_group : index_groups) {
      auto spin_bit = (spincase_bitstr << (64 - index_group_count - 1)) >> 63;
      assert((spin_bit == 0) || (spin_bit == 1));

      for (auto&& index : index_group) {
        IndexSpace space;
        if (spin_bit == 0) {
          space = IndexSpace::instance(
              IndexSpace::instance(index.label()).type(), IndexSpace::alpha);
        } else {
          space = IndexSpace::instance(
              IndexSpace::instance(index.label()).type(), IndexSpace::beta);
        }
        auto subscript_label =
            index.label().substr(index.label().find(L'_') + 1);
        std::wstring subscript_label_ws(subscript_label.begin(),
                                        subscript_label.end());

        Index spin_index = Index::make_label_index(space, subscript_label_ws);

        index_replacements.emplace(std::make_pair(index, spin_index));
        remove_spin_replacements.emplace(std::make_pair(spin_index, index));
      }
      ++index_group_count;
    }

    // std::cout << "\nPermutation " << spincase_bitstr << ": " << std::endl;
    // for(auto&& i: index_replacements)
    //   std::wcout << i.first.label() << " -> " << i.second.label() << ",  " ;
    // std::cout << std::endl;

    // EFV:  for every term in the sum (else we only have 1 term)
    // EFV:    apply replacement
    // EFV:    screen out zeroes
    // EFV:    if nonzero: remove spins + add to the running total

    sequant::Tensor subs_expr;
    auto summand_count = 0;
    Sum temp_sum;

    // EFV:  for every term in the sum (else we only have 1 term)
    for (auto&& expr_summand : *expr) {
      const auto scaler_factor = expr_summand->as<Product>().scalar().real();
      Product temp_product{};
      for (auto&& expr_product : *expr_summand) {
        subs_expr = expr_product->as<Tensor>();

        // EFV:    apply replacement
        auto pass_mutated =
            subs_expr.transform_indices(index_replacements, false);

        ExprPtr subs_expr_ptr = std::make_shared<Tensor>(subs_expr);

        // EFV:    screen out zeroes
        if (!tensor_symm(subs_expr)) break;

        temp_product.append(1.0, subs_expr_ptr);

        auto antisymm_expansion = expand_antisymm(subs_expr);
        std::wcout << antisymm_expansion->to_latex() << "\n" << std::endl;
      }

      if (tensor_symm(subs_expr)) {
        // std::wcout << temp_product.to_latex() << std::endl;
        summand_count++;
      }

      // std::wcout << "temp_product: " << temp_product.to_latex() << std::endl;

      if ((*expr_summand).size() == temp_product.size()) {
        temp_product.scale(scaler_factor);
        ExprPtr subs_expr_product_ptr = std::make_shared<Product>(temp_product);
        temp_sum.append(subs_expr_product_ptr);
      }
    }
    ExprPtr subs_expr_sum_ptr = std::make_shared<Sum>(temp_sum);
    spin_expr_sum.append(subs_expr_sum_ptr);
    // std::wcout << "sum of terms for " << spincase_bitstr << " permutation: "
    // << to_latex_align(subs_expr_sum_ptr)
    //            << "\n";
    total_terms += summand_count;

  }  // Permutation FOR loop

  const auto tstop = std::chrono::high_resolution_clock::now();
  auto time_elapsed =
      std::chrono::duration_cast<std::chrono::milliseconds>(tstop - tstart);
  std::cout << "Time: " << time_elapsed.count() << " milli sec; " << total_terms
            << " terms after expansion." << std::endl;

  // std::wcout << spin_expr_sum.to_latex() << std::endl;
  return nullptr;
}

}  // namespace sequant

#endif  // SEQUANT_SPIN_HPP
