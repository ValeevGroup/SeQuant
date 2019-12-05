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
  const uint64_t nspincases = std::pow(2, index_groups.size());

  Sum spin_expr_sum{};
  auto total_terms = 0;
  const auto tstart = std::chrono::high_resolution_clock::now();
  for (uint64_t spincase_bitstr = 0; spincase_bitstr != nspincases;
       ++spincase_bitstr) {
    // EFV:  assign spin to each index group => make a replacement list
    std::map<Index, Index> index_replacements;
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

    for (auto&& expr_summand : *expr) {
      // std::wcout << "expr_summand_factor: " <<
      // expr_summand->as<Product>().scalar().real() << std::endl;
      const auto scaler_factor = expr_summand->as<Product>().scalar().real();
      Product temp_product{};
      for (auto&& expr_product : *expr_summand) {
        subs_expr = expr_product->as<Tensor>();
        auto pass_mutated =
            subs_expr.transform_indices(index_replacements, false);

        std::shared_ptr<Expr> subs_expr_ptr =
            std::make_shared<Tensor>(subs_expr);
        if (!tensor_symm(subs_expr)) break;
        temp_product.append(1.0, subs_expr_ptr);
      }
      if (tensor_symm(subs_expr)) {
        // std::wcout << temp_product.to_latex() << std::endl;
        summand_count++;
      }
      // std::wcout << "temp_product: " << temp_product.to_latex() << std::endl;
      if ((*expr_summand).size() == temp_product.size()) {
        temp_product.scale(scaler_factor);
        std::shared_ptr<Expr> subs_expr_product_ptr =
            std::make_shared<Product>(temp_product);
        temp_sum.append(subs_expr_product_ptr);
      }
    }
    std::shared_ptr<Expr> subs_expr_sum_ptr = std::make_shared<Sum>(temp_sum);
    spin_expr_sum.append(subs_expr_sum_ptr);
    // std::wcout << "sum for current permutation: " << to_latex_align(subs_expr_sum_ptr)
    //            << std::endl;
    total_terms += summand_count;
  }
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
