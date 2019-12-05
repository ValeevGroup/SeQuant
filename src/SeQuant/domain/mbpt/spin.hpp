//
// Created by Eduard Valeyev on 2019-02-27.
//

#ifndef SEQUANT_SPIN_HPP
#define SEQUANT_SPIN_HPP

#include <algorithm>

// #include <bitset>
#include <codecvt>
#include <map>
#include <string>
#include <tuple>

#include "SeQuant/core/expr.hpp"
#include "SeQuant/core/tensor.hpp"
//#include "SeQuant/core/wick.impl.hpp"

#define SPINTRACE_PRINT 0
#define PRINT_PERMUTATION_LIST 0

namespace sequant {

class spinIndex : public Index {
 private:
};

struct zero_result : public std::exception {};

/// @return true if bra and ket indices for a tensor in @c product have same
/// spin
inline bool check_spin_symm(const std::shared_ptr<Product>& product) {
  bool result = false;

  // Get each tensor in product
  for (auto&& n : *product) {
    assert(n->as<Tensor>().bra().size() == n->as<Tensor>().ket().size());

    // For each index check if QNS match.
    // return FALSE if one of the pairs does not match.
    auto iter_ket = n->as<Tensor>().ket().begin();
    for (auto&& bra : n->as<Tensor>().bra()) {
      if (IndexSpace::instance(bra.label()).qns() ==
          IndexSpace::instance(iter_ket->label()).qns()) {
        result = true;
      } else
        return result;
      iter_ket++;
    }
  }
  return result;
}

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
    } else
      return result;
    iter_ket++;
  }
  return result;
}

/// @param spin_label_set list of indices with alpha/beta spin qns
/// @return the same list of indices with null qns
inline const container::set<Index, Index::LabelCompare> remove_spin_label(
    const container::set<Index, Index::LabelCompare>& spin_label_set) {
  container::set<Index, Index::LabelCompare> result;
  for (auto&& i : spin_label_set) {
    auto subscript_label = i.label().substr(i.label().find(L'_') + 1);
    std::wstring subscript_label_ws(subscript_label.begin(),
                                    subscript_label.end());
    auto space = IndexSpace::instance(IndexSpace::instance(i.label()).type(),
                                      IndexSpace::nullqns);
    result.insert(result.end(),
                  Index::make_label_index(space, subscript_label_ws));
  }
  return result;
}
#if 0
/// @return true if made any changes
// template <typename T0, typename T1>
inline bool apply_index_replacement_rules(
    std::shared_ptr<Product>& product,
    const container::map<Index, Index>& const_replrules,
    const container::set<Index>& external_indices,
    container::set<Index, Index::LabelCompare>& all_indices) {
  // to be able to use map[]
  auto& replrules = const_cast<container::map<Index, Index>&>(const_replrules);

  expr_range exrng(product);

  /// this recursively applies replacement rules until result does not
  /// change removes the deltas that are no longer needed
  const bool tag_transformed_indices =
      true;  // to replace indices when maps image and domain overlap, tag
  // transformed indices
#ifndef NDEBUG
  // assert that tensors_ indices are not tagged if going to tag indices
  if (tag_transformed_indices) {
    for (auto it = ranges::begin(exrng); it != ranges::end(exrng); ++it) {
      const auto& factor = *it;
      if (factor->is<Tensor>()) {
        auto& tensor = factor->as<Tensor>();
        assert(ranges::none_of(tensor.const_braket(), [](const Index& idx) {
          return idx.tag().has_value();
        }));
      }
    }
  }
#endif
  bool mutated = false;
  bool pass_mutated = false;
  do {
    pass_mutated = false;

    for (auto it = ranges::begin(exrng); it != ranges::end(exrng);) {
      const auto& factor = *it;
      if (factor->is<Tensor>()) {
        bool erase_it = false;
        auto& tensor = factor->as<Tensor>();

        /// replace indices
        pass_mutated &=
            tensor.transform_indices(const_replrules, tag_transformed_indices);

        if (tensor.label() == L"S") {
          const auto& bra = tensor.bra().at(0);
          const auto& ket = tensor.ket().at(0);

          if (bra.proto_indices() == ket.proto_indices()) {
            const auto bra_is_ext = ranges::find(external_indices, bra) !=
                                    ranges::end(external_indices);
            const auto ket_is_ext = ranges::find(external_indices, ket) !=
                                    ranges::end(external_indices);

#ifndef NDEBUG
            const auto intersection_space =
                intersection(bra.space(), ket.space());
#endif

            if (!bra_is_ext && !ket_is_ext) {  // int + int
#ifndef NDEBUG
              if (replrules.find(bra) != replrules.end() &&
                  replrules.find(ket) != replrules.end())
                assert(replrules[bra].space() == replrules[ket].space());
#endif
              erase_it = true;
            } else if (bra_is_ext && !ket_is_ext) {  // ext + int
              if (includes(ket.space(), bra.space())) {
#ifndef NDEBUG
                if (replrules.find(ket) != replrules.end())
                  assert(replrules[ket].space() == bra.space());
#endif
                erase_it = true;
              } else {
#ifndef NDEBUG
                if (replrules.find(ket) != replrules.end())
                  assert(replrules[ket].space() == intersection_space);
#endif
              }
            } else if (!bra_is_ext && ket_is_ext) {  // int + ext
              if (includes(bra.space(), ket.space())) {
#ifndef NDEBUG
                if (replrules.find(bra) != replrules.end())
                  assert(replrules[bra].space() == ket.space());
#endif
                erase_it = true;
              } else {
#ifndef NDEBUG
                if (replrules.find(bra) != replrules.end())
                  assert(replrules[bra].space() == intersection_space);
#endif
              }
            } else {  // ext + ext
              if (bra == ket) erase_it = true;
            }

            if (erase_it) {
              pass_mutated = true;
              *it = ex<Constant>(1);
            }
          }  // matching proto indices
        }    // Kronecker delta
      }
      ++it;
    }
    mutated |= pass_mutated;
  } while (pass_mutated);  // keep replacing til fixed point

  // assert that tensors_ indices are not tagged if going to tag indices
  if (tag_transformed_indices) {
    for (auto it = ranges::begin(exrng); it != ranges::end(exrng); ++it) {
      const auto& factor = *it;
      if (factor->is<Tensor>()) {
        factor->as<Tensor>().reset_tags();
      }
    }
  }

  // update all_indices
  container::set<Index, Index::LabelCompare> all_indices_new;
  ranges::for_each(
      all_indices, [&const_replrules, &all_indices_new](const Index& idx) {
        auto dst_it = const_replrules.find(idx);
        auto insertion_result = all_indices_new.emplace(
            dst_it != const_replrules.end() ? dst_it->second : idx);
      });
  std::swap(all_indices_new, all_indices);

  return mutated;
}


template <typename T>
inline const T generate_alpha_idx(const T& all_indices) {
  T result;
  for (auto&& i : all_indices) {
    auto subscript_label = i.label().substr(i.label().find(L'_') + 1);
    std::wstring subscript_label_ws(subscript_label.begin(),
                                    subscript_label.end());

    auto space = IndexSpace::instance(IndexSpace::instance(i.label()).type(),
                                      IndexSpace::alpha);
    result.insert(result.end(),
                  Index::make_label_index(space, subscript_label_ws));
  }
  return result;
}

template <typename T>
inline const T generate_beta_idx(const T& all_indices) {
  T result;
  for (auto&& i : all_indices) {
    auto subscript_label = i.label().substr(i.label().find(L'_') + 1);
    std::wstring subscript_label_ws(subscript_label.begin(),
                                    subscript_label.end());

    auto space = IndexSpace::instance(IndexSpace::instance(i.label()).type(),
                                      IndexSpace::beta);
    result.insert(result.end(),
                  Index::make_label_index(space, subscript_label_ws));
  }
  return result;
}

/// @param source_index list all indices in tensor product
/// @param dest_index
/// @return boost map
template <typename T0, typename T1>
inline const container::map<Index, Index> map_constructor(
    const T0& source_index,
    const T1& dest_index) {
  container::map<Index, Index> result;
  auto dest_ptr = dest_index.begin();
  for (auto&& i : source_index) {
    result.emplace(i, *dest_ptr);
    dest_ptr++;
  }
  return result;
}


// This will return a list that should be used in the map function
// BUG:TODO: This function returns an incorrect map: the number sequence in
// label is not honored
/// @param list1 is the alpha index list
/// @param list2 is the beta index list
/// @return list of vectors containing permuted alpha and beta indices
template <typename T0, typename T1>
const std::vector<T1> generate_permutation_index_list(
    const T0& list1,
    const T1& list2) {

  std::vector<T1> result;
  result.push_back(list1);

//  ranges::for_each(list1, [&] (Index i) {std::wcout << i.to_latex() << " ";});
//  std::cout << std::endl;

  assert(list1.size() == list2.size());

  // Generate a string containing just 'A'
  std::string string_tuple(list1.size(), 'A');

  // Copy contents of list1 to a vector
  std::vector<Index> permuted_index_vector(list1.begin(), list1.end());

  // Replace the last 'A' in tuple with a 'B' and permute
  for (size_t i = 0; i < string_tuple.size(); i++) {
    string_tuple.replace(string_tuple.size() - i - 1, 1, "B");
    std::sort(string_tuple.begin(), string_tuple.end());

    do {
      // Iterator pointing to begining of l1 and l2
      auto alpha_iter = list1.begin();
      auto beta_iter = list2.begin();
      permuted_index_vector.clear();

      for (auto&& i : string_tuple) {
        if (i == 'A') {
          permuted_index_vector.insert(permuted_index_vector.end(), *alpha_iter);
        } else {
          permuted_index_vector.insert(permuted_index_vector.end(), *beta_iter);
        }
        alpha_iter++;
        beta_iter++;
      }

      // TODO: It is apparant from the following lines that the 'ORDERING' of indices in a set
      // must be edited such that the number in a label will have preference over the qns

      // std::cout << "V: "; // Vector
      // ranges::for_each(permuted_index_vector, [&] (Index i) {std::wcout << i.label() << " ";});
      // std::cout << std::endl;

      // Convert vector to container set and push back to result
      container::set<Index, Index::LabelCompare> permuted_index_list(
          permuted_index_vector.begin(), permuted_index_vector.end());

      // std::cout << "S: "; // Set
      // ranges::for_each(permuted_index_list, [&] (Index i) {std::wcout << i.label() << " ";});
      // std::cout << std::endl;

      result.push_back(permuted_index_list);
    } while (std::next_permutation(string_tuple.begin(), string_tuple.end()));
  }

  return result;
}
#endif

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

  for (auto&& i : index_groups)
    for (auto&& j : i) std::wcout << j.label() << " ";
  std::cout << std::endl;

  // EFV: for each spincase (loop over integer from 0 to 2^(n)-1, n=#of index
  // groups)
  const uint64_t nspincases = std::pow(2, index_groups.size());

  ExprPtr addition_list;
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

    std::cout << "\nPermutation " << spincase_bitstr << ": " << std::endl;
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
      // std::wcout << "expr_summand: " << expr_summand->as<Product>().scalar().real() << std::endl;
      const auto scaler_factor = expr_summand->as<Product>().scalar().real();
      Product temp_product{};
      for (auto&& expr_product : *expr_summand) {
        subs_expr = expr_product->as<Tensor>();
        auto pass_mutated =
            subs_expr.transform_indices(index_replacements, false);
        // std::wcout << " subs_expr: " << subs_expr.to_latex() << std::endl; // sub_expr is Tensor/Expr

        std::shared_ptr<Expr> subs_expr_ptr = std::make_shared<Tensor>(subs_expr);
        if (!tensor_symm(subs_expr))
          break;
        temp_product.append(1.0, subs_expr_ptr);
      }
      if (tensor_symm(subs_expr)) {
        // std::wcout << temp_product.to_latex() << std::endl;
        summand_count++;
      }
      // std::wcout << "temp_product: " << temp_product.to_latex() << std::endl;
      temp_product.scale(scaler_factor);
      std::shared_ptr<Expr> subs_expr_product_ptr = std::make_shared<Product>(temp_product);
      temp_sum.append(subs_expr_product_ptr);
    }
    std::shared_ptr<Expr> subs_expr_sum_ptr = std::make_shared<Sum>(temp_sum);
    std::wcout << "temp_sum: " << to_latex_align(subs_expr_sum_ptr) << std::endl;
    total_terms += summand_count;
  }
  const auto tstop = std::chrono::high_resolution_clock::now();
  auto time_elapsed =
      std::chrono::duration_cast<std::chrono::milliseconds>(tstop - tstart);
  std::cout << "Time: " << time_elapsed.count() << " milli sec; " << total_terms
            << " terms after expansion." << std::endl;

  return nullptr;
}

}  // namespace sequant

#endif  // SEQUANT_SPIN_HPP
