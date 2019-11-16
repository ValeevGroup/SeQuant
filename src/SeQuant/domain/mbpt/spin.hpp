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

/// @return true if made any changes
inline bool apply_index_replacement_rules(
    std::shared_ptr<Product>& product,
    const container::map<Index, Index>& const_replrules,
    const container::set<Index>& external_indices,
    std::set<Index, Index::LabelCompare>& all_indices) {
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
  std::set<Index, Index::LabelCompare> all_indices_new;
  ranges::for_each(
      all_indices, [&const_replrules, &all_indices_new](const Index& idx) {
        auto dst_it = const_replrules.find(idx);
        auto insertion_result = all_indices_new.emplace(
            dst_it != const_replrules.end() ? dst_it->second : idx);
      });
  std::swap(all_indices_new, all_indices);

  return mutated;
}

const container::set<Index> generate_alpha_idx(
    const std::set<Index, Index::LabelCompare>& all_indices) {
  container::set<Index> result;
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

const container::set<Index> generate_beta_idx(
    const std::set<Index, Index::LabelCompare>& all_indices) {
  container::set<Index> result;
  for (auto&& i : all_indices) {
    auto subscript_label = i.label().substr(i.label().find(L'_') + 1);
    std::wstring subscript_label_ws(subscript_label.begin(),
                                    subscript_label.end());

    // TODO: The spin_type does not work
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
const container::map<Index, Index> map_constructor(
    const std::set<Index, Index::LabelCompare>& source_index,
    container::set<Index>& dest_index) {
  container::map<Index, Index> result;
  auto dest_ptr = dest_index.begin();
  for (auto&& i : source_index) {
    result.emplace(i, *dest_ptr);
    dest_ptr++;
  }
  return result;
}

const std::vector<container::map<Index, Index>> generate_map_list(
    container::set<Index>& source_index,
    container::set<Index>& destination_index) {
  assert(source_index.size() == destination_index.size());
  std::vector<container::map<Index, Index>> result;

  // TODO: generate tuple list
  std::cout << source_index.size() << std::endl;

  // TODO: call map_constructor

  // TODO: pushback to list

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

  auto add_spin_labels = [&](ExprPtr& n) {
    if (n->is<Product>()) {
      std::wcout << (n->as<Product>()).to_latex() << std::endl;

      bool pass_mutated = false;
      // extract current indices
      std::set<Index, Index::LabelCompare> all_indices;
      ranges::for_each(*n, [&all_indices](const auto& factor) {
        if (factor->template is<Tensor>()) {
          ranges::for_each(factor->template as<const Tensor>().braket(),
                           [&all_indices](const Index& idx) {
                             auto result = all_indices.insert(idx);
                           });
        }
      });

#if 0
      container::set<Index> alpha_list;
      container::set<Index> beta_list;
      {
        // For each index, copy .type() and add spin variable
        for (auto&& i : all_indices) {
          auto subscript_label = i.label().substr(i.label().find(L'_') + 1);
          std::wstring subscript_label_ws(subscript_label.begin(),
                                          subscript_label.end());

          auto alpha_space = IndexSpace::instance(
              IndexSpace::instance(i.label()).type(), IndexSpace::alpha);
          auto beta_space = IndexSpace::instance(
              IndexSpace::instance(i.label()).type(), IndexSpace::beta);
          alpha_list.insert(
              alpha_list.end(),
              Index::make_label_index(alpha_space, subscript_label_ws));
          beta_list.insert(
              beta_list.end(),
              Index::make_label_index(beta_space, subscript_label_ws));
        }
      }
#endif
      auto alpha_list = generate_alpha_idx(all_indices);
      auto beta_list = generate_beta_idx(all_indices);

      /*
          {
            std::vector<int> a1{1, 2, 3, 4, 5};
            std::vector<int> a2{2, 4, 6, 8, 10} ; //(a1.begin(), a1.end());

            container::map<int, int> map;
            assert(a1.size() == a2.size());
            for (size_t i = 0; i < a1.size(); ++i)
               map.emplace(a1[i], a2[i]); //  map[a1[i]] = a2[i];

            container::map<int, int>::iterator itr;
            for(itr = map.begin(); itr != map.end(); ++itr){
              std::cout << itr->first << " " << itr->second << std::endl;
            }
          }
      */

      /* {
         std::vector<Index> alpha_vector(alpha_list.begin(), alpha_list.end());
         std::vector<Index> beta_vector(beta_list.begin(), beta_list.end());
         std::vector<Index> allidx_vector(all_indices.begin(),
       all_indices.end());

         std::cout << __FILE__ << " " << __LINE__ << std::endl;

         container::map<Index, Index> spin_idx_map;
         assert(allidx_vector.size() == alpha_vector.size());
         for (size_t i = 0; i < allidx_vector.size(); ++i)
           spin_idx_map[allidx_vector[i]] = alpha_vector[i];

         container::map<Index, Index>::iterator vec_iter;
         for (vec_iter = spin_idx_map.begin(); vec_iter != spin_idx_map.end();
       ++vec_iter) { std::wcout << (vec_iter->first).label() << " " <<
       (vec_iter->second).label() << std::endl;
         }
       }*/

      //    std::cout << __FILE__ << " " <<  __LINE__ << std::endl;

      /*
      assert(all_indices.size() == alpha_list.size());
      auto ptr_t2 = alpha_list.begin();
      for(auto &&i : all_indices){
  //      spin_idx_map[*(ptr_t2)] = i;
        std::wcout << i.label() << " ";
        std::wcout << (*ptr_t2).label() << std::endl;
        ++ptr_t2;
      }
  */
      //    std::cout << __FILE__ << " " <<  __LINE__ << std::endl;

      //    container::map<Index,Index>::iterator itr;
      //    for(itr = spin_idx_map.begin(); itr != spin_idx_map.end(); ++itr)
      //      std::wcout << (itr->first).label() << " " << (itr->second).label()
      //      << std::endl;

#if 0
      container::map<Index, Index> i_to_alpha;
      container::map<Index, Index> i_to_beta;
      {
        auto alpha_list_ptr = alpha_list.begin();
        auto beta_list_ptr = beta_list.begin();
        for (auto&& i : all_indices) {
          i_to_alpha.emplace(i, *alpha_list_ptr);
          i_to_beta.emplace(i, *beta_list_ptr);
          std::wcout << i.label() << " " << (*alpha_list_ptr).label() << " "
                     << (*beta_list_ptr).label() << std::endl;
          alpha_list_ptr++;
          beta_list_ptr++;
        }
      }
#endif

      auto i_to_alpha = map_constructor(all_indices, alpha_list);
      auto i_to_beta = map_constructor(all_indices, beta_list);

      // std::shared_ptr<Product> expr_product_alpha (new
      // Product(n->as<Product>())); std::shared_ptr<Product> expr_product =
      // std::make_shared<Product>(Product(n->as<Product>()));
      auto expr_cast = std::static_pointer_cast<Product>(n);

      // TODO: Apply index replacement rules

      //      pass_mutated =
      //      sequant::apply_index_replacement_rules(expr_product,
      //      replacement_rules, ext_idxlist, all_indices);
      pass_mutated = sequant::apply_index_replacement_rules(
          expr_cast, i_to_alpha, ext_idxlist, all_indices);
      std::wcout << "mutated: " << pass_mutated << " -> "
                 << expr_cast->to_latex() << std::endl;

      //    std::shared_ptr<Product> expr_product_beta (new
      //    Product(n->as<Product>())); std::shared_ptr<Product>
      //    expr_product_beta =
      //    std::make_shared<Product>(Product(n->as<Product>())); auto
      //    expr_product_beta(expr_product_alpha);

      auto expr_cast2 = std::static_pointer_cast<Product>(n);
      pass_mutated = sequant::apply_index_replacement_rules(
          expr_cast2, i_to_beta, ext_idxlist, all_indices);
      std::wcout << "mutated: " << pass_mutated << " -> "
                 << expr_cast2->to_latex() << "\n"
                 << std::endl;
    }
  };
  expr->visit(add_spin_labels);

#if SPINTRACE_PRINT
  //  List of all Index ( label + IndexSpace )
  for (auto&& idx : grand_idxlist) {
    std::wcout << idx.to_latex() << std::endl;
  }
  std::cout << std::endl;

  for (auto i = grand_idxlist.begin(); i != grand_idxlist.end(); ++i)
    std::wcout << i->to_latex() << std::endl;

#endif

  // TODO make list of groups, generate all spin tuples, apply and replace ...

  // Generate list of spin quantum numbers from Grand Index List
  container::set<Index> alpha_list;
  container::set<Index> beta_list;
  {
    // For each index, copy .type() and add spin variable
    for (auto&& i : grand_idxlist) {
      auto alpha_space = IndexSpace::instance(
          IndexSpace::instance(i.label()).type(), IndexSpace::alpha);
      auto beta_space = IndexSpace::instance(
          IndexSpace::instance(i.label()).type(), IndexSpace::beta);
      alpha_list.insert(alpha_list.end(), Index::make_tmp_index(alpha_space));
      beta_list.insert(beta_list.end(), Index::make_tmp_index(beta_space));
    }

#if SPINTRACE_PRINT
    for (auto&& i : alpha_list) {
      std::wcout << i.label() << " {" << IndexSpace::instance(i.label()).type()
                 << "," << IndexSpace::instance(i.label()).qns() << "}"
                 << std::endl;
    }
    for (auto&& i : beta_list) {
      std::wcout << i.label() << " {" << IndexSpace::instance(i.label()).type()
                 << "," << IndexSpace::instance(i.label()).qns() << "}"
                 << std::endl;
    }
    std::cout << "alpha_list size: " << alpha_list.size() << std::endl;
    std::cout << "beta_list size: " << beta_list.size() << std::endl;
#endif
  }

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

  // TEST 1
  /*
  {
    for(auto&& i : tupleList){
      auto num = &i - &tupleList[0];
//      std::cout << num << "  " << std::endl;
      for(auto&& j: grand_idxlist){
        int num2 = &j - &grand_idxlist[0];
//      std::cout << tupleList.at(num) << std::endl;

      }
    }
  }


  std::vector<Index> grandIndex_vector;
  for(auto &&i : grand_idxlist){
    grandIndex_vector.push_back(i);
  }

  std::cout << "grandIndex_vector.size(): " << grandIndex_vector.size() <<
std::endl;

  for(auto i = grandIndex_vector.begin(); i<= grandIndex_vector.end(); ++i){
    std::cout << grandIndex_vector.at(i) << std::endl;
  }
*/

  //  TEST 2
  /*
  {
    auto temp_counter = 0;

    //  TODO: Create list of groups
    for (auto&& i : tupleList) {
      temp_counter++;
      auto num = &i - &tupleList[0];

      //      std::cout << num << " " << i << std::endl;

      // container::set<Index>::iterator it = alpha_list.begin();

      //  Loop over each character in the tuple
      for (auto&& j : i) {
        auto idx = &j - &i[0];

        if (j == 'A') {
        } else {
        }
      }
      //      std::cout << std::endl;
    }
    std::cout << "temp_counter: " << temp_counter << std::endl;
  }
*/
  //  auto expr2 = expr->clone();
  //  std::wcout << "expr2:\n" << expr2->to_latex() << std::endl;

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
