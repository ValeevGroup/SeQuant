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

/// @return true if made any changes
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

inline const container::set<Index, Index::LabelCompare> generate_alpha_idx(
    const container::set<Index, Index::LabelCompare>& all_indices) {
  container::set<Index, Index::LabelCompare> result;
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

inline const container::set<Index, Index::LabelCompare> generate_beta_idx(
    const container::set<Index, Index::LabelCompare>& all_indices) {
  container::set<Index, Index::LabelCompare> result;
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
inline const container::map<Index, Index> map_constructor(
    const container::set<Index, Index::LabelCompare>& source_index,
    container::set<Index, Index::LabelCompare>& dest_index) {
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
const std::vector<container::set<Index, Index::LabelCompare>>
generate_permutation_index_list(
    const container::set<Index, Index::LabelCompare>& list1,
    const container::set<Index, Index::LabelCompare>& list2) {

  std::vector<container::set<Index, Index::LabelCompare>> result;
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

      std::cout << "V: "; // Vector
      ranges::for_each(permuted_index_vector, [&] (Index i) {std::wcout << i.label() << " ";});
      std::cout << std::endl;

      // Convert vector to container set and push back to result
      container::set<Index, Index::LabelCompare> permuted_index_list(
          permuted_index_vector.begin(), permuted_index_vector.end());

      std::cout << "S: "; // Set
      ranges::for_each(permuted_index_list, [&] (Index i) {std::wcout << i.label() << " ";});
      std::cout << std::endl;

      result.push_back(permuted_index_list);
    } while (std::next_permutation(string_tuple.begin(), string_tuple.end()));
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

  container::set<Index> grand_idxlist;
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

  // Returns expr with a prefactor after spin tracing
  auto add_spin_labels = [&](ExprPtr& n) {
    auto result = n;
    if (n->is<Product>()) {
      std::wcout << (n->as<Product>()).to_latex() << std::endl;

      bool pass_mutated = false;
      // extract current indices
      container::set<Index, Index::LabelCompare> all_indices;
      ranges::for_each(*n, [&all_indices](const auto& factor) {
        if (factor->template is<Tensor>()) {
          ranges::for_each(factor->template as<const Tensor>().braket(),
                           [&all_indices](const Index& idx) {
                             auto result = all_indices.insert(idx);
                           });
        }
      });

      // Generate list with spins
      auto alpha_list = generate_alpha_idx(all_indices);
      auto beta_list = generate_beta_idx(all_indices);

      // generate maps TO the spin list
      auto i_to_alpha = map_constructor(all_indices, alpha_list);
      auto i_to_beta = map_constructor(all_indices, beta_list);

      for(auto &&i: i_to_alpha)
        std::wcout << "L" << __LINE__ << " " << i.first.label() << " -> " << i.second.label() << std::endl;

//      for(auto &&i: i_to_beta)
//        std::wcout << "L" << __LINE__ << " " << i.first.label() << " -> " << i.second.label() << std::endl;

      // Generates a vector of target index
      // BUG: The index replacement rules are incorrect
      auto tuple_vector_list =
          generate_permutation_index_list(alpha_list, beta_list);

      ranges::for_each(all_indices, [&] (Index i) {std::wcout << i.to_latex() << " ";});
      std::cout << std::endl;

      auto counter = 1;
      auto scaler_multiplier = 1;
      for (auto&& i : tuple_vector_list) {

        ranges::for_each(i, [&] (Index j) {std::wcout << j.to_latex() << " ";});
        std::cout << " " << counter << std::endl;

        // Generate replacement map
        auto index_map = map_constructor(all_indices, i);

        auto expr_cast = std::static_pointer_cast<Product>(n);

        // Ordering of labels is causing a problem
        pass_mutated = sequant::apply_index_replacement_rules(
            expr_cast, index_map, ext_idxlist, all_indices);

        // Check if spin qns on each pair of index is matching in a product
        auto spin_symm = check_spin_symm(expr_cast);

        // TODO: Use permut symm in expression to count -1.0
        auto spin_antisymm = false;

        if (spin_symm) {
          scaler_multiplier += 1;
          // std::wcout << "L" << __LINE__ << " " << n->to_latex() << std::endl;
          auto spin_removed_list = remove_spin_label(i);
          auto spin_free_map = map_constructor(i, spin_removed_list);
          pass_mutated = sequant::apply_index_replacement_rules(
              expr_cast, spin_free_map, ext_idxlist, all_indices);
          // std::wcout << "L" << __LINE__ << " " << n->to_latex() << std::endl;
        } else if (spin_antisymm) {

        }
        counter++;
      }

      auto multiplier = 1.0;
      multiplier *= scaler_multiplier;
      (n->as<Product>()).scale(multiplier);
      std::wcout << "L" << __LINE__ << " " << n->to_latex() << std::endl;
    }
  };
  expr->visit(add_spin_labels);

  std::wcout << "\n" << expr->to_latex() << std::endl;

#if SPINTRACE_PRINT
  //  List of all Index ( label + IndexSpace )
  for (auto&& idx : grand_idxlist) {
    std::wcout << idx.to_latex() << std::endl;
  }
  std::cout << std::endl;

//  for (auto i = grand_idxlist.begin(); i != grand_idxlist.end(); ++i)
//    std::wcout << i->to_latex() << std::endl;
#endif
  // TODO make list of groups, generate all spin tuples, apply and replace ...

  // Generate a list of tuples with A,B for spin
  auto tupleSize = 0;  // Counter for size of the list of tuples
  std::vector<std::string> tupleList;

  //  Create a list of A spins (the size of number of index)
  std::string strTuple(grand_idxlist.size(), 'A');

  //  Add list of A spins to tupleList
  tupleList.push_back(strTuple);
  tupleSize++;
#if PRINT_PERMUTATION_LIST
  std::cout << tupleSize << "\t" << strTuple << std::endl;
#endif

  //  Replace one spin in a tuple
  for (size_t i = 0; i < int_idxlist.size(); i++) {
    strTuple.replace(int_idxlist.size() - i - 1, 1, "B");
    std::sort(strTuple.begin(), strTuple.end());

    //  Permute the tuples
    do {
      tupleList.push_back(strTuple);
      tupleSize++;
#if PRINT_PERMUTATION_LIST
      std::cout << tupleSize << "\t" << strTuple << std::endl;
#endif
    } while (std::next_permutation(strTuple.begin(), strTuple.end()));
  }
  std::cout << "tupleList size: " << tupleList.size() << std::endl;

  // TODO: Merge tupleList elements and Index

  struct latex_visitor {
    void operator()(const std::shared_ptr<sequant::Expr>& expr) {
      result += expr->to_latex();
    }
    std::wstring result{};
  };
  //  latex_visitor v1{};
  //  expr->visit(v1);
  //  std::wcout << "v1.result = " << v1.result << std::endl;

  //  latex_visitor v2{};
  //  expr->visit(v2, true);
  //  std::wcout << "v2.result = " << v2.result << std::endl;

  // YOU DON't NEED A LIST OF BRA OR KET INDICES
  container::set<Index> bra_idxlist;
  auto collect_bra_indices = [&](const ExprPtr& expr) {
    if (expr->is<Tensor>()) {
      std::wcout << expr->as<Tensor>().to_latex() << std::endl;
      for (auto& i : expr->as<Tensor>().bra())
        std::wcout << sequant::to_latex(i) << std::endl;
      ranges::for_each(expr->as<Tensor>().bra(),
                       [&](const Index& idx) { bra_idxlist.insert(idx); });
    }
  };
  //  expr->visit(collect_bra_indices);

  //  Get sub expressions from an expression:
  auto print_subexpr = [](const auto& n) {
    std::wcout << "&" << n->to_latex() << L"\\\\" << std::endl;
  };
  //  ranges::for_each(expr->begin(), expr->end(), print_subexpr);

  // The bra and ket indices from each tensor are separated here
  auto print_subexpr_index = [&](ExprPtr& nn) {
    std::wcout << nn->to_latex() << " Expr: " << nn->is<Expr>() << std::endl;
    //    std::wcout << "Product: " << nn->is<Product>() << std::endl;
    //    std::wcout << "Sum: " << nn->is<Sum>() << std::endl;
    for (auto& n : *nn) {
      std::wcout << n->to_latex() << " tensor: " << n->is<Tensor>()
                 << std::endl;
      //    Get a list of indices for a tensor. BRA and KET needs to be separate
      //    lists. Add spin attribute to the index and regenetate the tensor
      container::set<Index> ketSpinIndexListA;
      container::set<Index> ketSpinIndexListB;
      for (auto& i : n->as<Tensor>().ket()) {
        auto subscript_label = i.label().substr(i.label().find(L'_') + 1);
        std::wstring subscript_label_ws(subscript_label.begin(),
                                        subscript_label.end());
        //        std::wcout << i.label() << " " << subscript_label <<
        //        std::endl;

        auto alpha_space = IndexSpace::instance(
            IndexSpace::instance(i.label()).type(), IndexSpace::alpha);
        ketSpinIndexListA.insert(
            ketSpinIndexListA.end(),
            Index::make_label_index(alpha_space, subscript_label_ws));
        auto beta_space = IndexSpace::instance(
            IndexSpace::instance(i.label()).type(), IndexSpace::beta);
        ketSpinIndexListB.insert(
            ketSpinIndexListB.end(),
            Index::make_label_index(beta_space, subscript_label_ws));
      }

      container::set<Index> braSpinIndexListA;
      container::set<Index> braSpinIndexListB;
      for (auto& i : n->as<Tensor>().bra()) {
        auto subscript_label = i.label().substr(i.label().find(L'_') + 1);
        std::wstring subscript_label_ws(subscript_label.begin(),
                                        subscript_label.end());
        //        std::wcout << i.label() << " " << subscript_label <<
        //        std::endl;

        auto alpha_space = IndexSpace::instance(
            IndexSpace::instance(i.label()).type(), IndexSpace::alpha);
        braSpinIndexListA.insert(
            braSpinIndexListA.end(),
            Index::make_label_index(alpha_space, subscript_label_ws));

        auto beta_space = IndexSpace::instance(
            IndexSpace::instance(i.label()).type(), IndexSpace::beta);
        braSpinIndexListB.insert(
            braSpinIndexListB.end(),
            Index::make_label_index(beta_space, subscript_label_ws));
      }

      auto tempA = Tensor((n->as<Tensor>()).label(), braSpinIndexListA,
                          ketSpinIndexListA);
      auto tempB = Tensor((n->as<Tensor>()).label(), braSpinIndexListB,
                          ketSpinIndexListB);

      std::wcout << to_latex(tempA) << " " << to_latex(tempB) << std::endl;

      // std::wcout << to_latex(tempA.as<Product>()) << std::endl;

      auto& temp_expr = tempA.as<Expr>();
      std::wcout << "tensor: " << tempA.is<Expr>()
                 << " Expr: " << temp_expr.is<Expr>() << std::endl;
    }
    // TODO: Form expression from individual tensors
    std::cout << std::endl;
  };
  //  expr->visit(print_subexpr_index);

  // TODO: Add same terms

  // TODO: Zero-out non-spinconserving terms

  // TODO: Reindex and return

  return nullptr;
}

}  // namespace sequant

#endif  // SEQUANT_SPIN_HPP
