//
// Created by Eduard Valeyev on 2019-02-27.
//

#ifndef SEQUANT_SPIN_HPP
#define SEQUANT_SPIN_HPP

#include <algorithm>

#include <codecvt>
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

inline container::map<Index, Index> compute_index_replacement_rules(
    std::shared_ptr<Product> &product,
    const container::set<Index> &external_indices,
    const std::set<Index, Index::LabelCompare> &all_indices) {
  expr_range exrng(product);

  /// this ensures that all temporary indices have unique *labels* (not just
  /// unique *full labels*)
  auto index_validator = [&all_indices](const Index &idx) {
    return all_indices.find(idx) == all_indices.end();
  };
  IndexFactory idxfac(index_validator);
  container::map<Index /* src */, Index /* dst */> result;  // src->dst

  // computes an index in intersection of space1 and space2
  auto make_intersection_index = [&idxfac](const IndexSpace &space1,
                                           const IndexSpace &space2) {
    const auto intersection_space = intersection(space1, space2);
    if (intersection_space == IndexSpace::null_instance()) throw zero_result{};
    return idxfac.make(intersection_space);
  };

  // transfers proto indices from idx (if any) to img
  auto proto = [](const Index &img, const Index &idx) {
    if (idx.has_proto_indices()) {
      if (img.has_proto_indices()) {
        assert(img.proto_indices() == idx.proto_indices());
        return img;
      } else
        return Index(img, idx.proto_indices());
    } else {
      assert(!img.has_proto_indices());
      return img;
    }
  };

  // adds src->dst or src->intersection(dst,current_dst)
  auto add_rule = [&result, &proto, &make_intersection_index](const Index &src,
                                                              const Index &dst) {
    auto src_it = result.find(src);
    if (src_it == result.end()) {  // if brand new, add the rule
      auto insertion_result = result.emplace(src, proto(dst, src));
      assert(insertion_result.second);
    }
    else {  // else modify the destination of the existing rule to the
      // intersection
      const auto &old_dst = src_it->second;
      assert(old_dst.proto_indices() == src.proto_indices());
      if (dst.space() != old_dst.space()) {
        src_it->second =
            proto(make_intersection_index(old_dst.space(), dst.space()), src);
      }
    }
  };

  // adds src1->dst and src2->dst; if src1->dst1 and/or src2->dst2 already
  // exist the existing rules are updated to map to the intersection of dst1,
  // dst2 and dst
  auto add_rules = [&result, &idxfac, &proto, &make_intersection_index](
      const Index &src1, const Index &src2, const Index &dst) {
    // are there replacement rules already for src{1,2}?
    auto src1_it = result.find(src1);
    auto src2_it = result.find(src2);
    const auto has_src1_rule = src1_it != result.end();
    const auto has_src2_rule = src2_it != result.end();

    // which proto-indices should dst1 and dst2 inherit? a source index without
    // proto indices will inherit its source counterpart's indices, unless it
    // already has its own protoindices: <a_ij|p> = <a_ij|a_ij> (hence replace p
    // with a_ij), but <a_ij|p_kl> = <a_ij|a_kl> != <a_ij|a_ij> (hence replace
    // p_kl with a_kl)
    const auto &dst1_proto =
        !src1.has_proto_indices() && src2.has_proto_indices() ? src2 : src1;
    const auto &dst2_proto =
        !src2.has_proto_indices() && src1.has_proto_indices() ? src1 : src2;

    if (!has_src1_rule && !has_src2_rule) {  // if brand new, add the rules
      auto insertion_result1 = result.emplace(src1, proto(dst, dst1_proto));
      assert(insertion_result1.second);
      auto insertion_result2 = result.emplace(src2, proto(dst, dst2_proto));
      assert(insertion_result2.second);
    } else if (has_src1_rule &&
        !has_src2_rule) {  // update the existing rule for src1
      const auto &old_dst1 = src1_it->second;
      assert(old_dst1.proto_indices() == dst1_proto.proto_indices());
      if (dst.space() != old_dst1.space()) {
        src1_it->second = proto(
            make_intersection_index(old_dst1.space(), dst.space()), dst1_proto);
      }
      result[src2] = src1_it->second;
    } else if (!has_src1_rule &&
        has_src2_rule) {  // update the existing rule for src2
      const auto &old_dst2 = src2_it->second;
      assert(old_dst2.proto_indices() == dst2_proto.proto_indices());
      if (dst.space() != old_dst2.space()) {
        src2_it->second = proto(
            make_intersection_index(old_dst2.space(), dst.space()), dst2_proto);
      }
      result[src1] = src2_it->second;
    } else {  // update both of the existing rules
      const auto &old_dst1 = src1_it->second;
      const auto &old_dst2 = src2_it->second;
      const auto new_dst_space =
          (dst.space() != old_dst1.space() || dst.space() != old_dst2.space())
          ? intersection(old_dst1.space(), old_dst2.space(), dst.space())
          : dst.space();
      if (new_dst_space == IndexSpace::null_instance()) throw zero_result{};
      Index new_dst;
      if (new_dst_space == old_dst1.space()) {
        new_dst = old_dst1;
        if (new_dst_space == old_dst2.space() && old_dst2 < new_dst) {
          new_dst = old_dst2;
        }
        if (new_dst_space == dst.space() && dst < new_dst) {
          new_dst = dst;
        }
      } else if (new_dst_space == old_dst2.space()) {
        new_dst = old_dst2;
        if (new_dst_space == dst.space() && dst < new_dst) {
          new_dst = dst;
        }
      } else if (new_dst_space == dst.space()) {
        new_dst = dst;
      } else
        new_dst = idxfac.make(new_dst_space);
      result[src1] = proto(new_dst,dst1_proto);
      result[src2] = proto(new_dst,dst2_proto);
    }
  };

  /// this makes the list of replacements ... we do not mutate the expressions
  /// to keep the information about which indices are related
  for (auto it = ranges::begin(exrng); it != ranges::end(exrng);
       ++it) {
    const auto &factor = *it;
    if (factor->type_id() == Expr::get_type_id<Tensor>()) {
      const auto &tensor = static_cast<const Tensor &>(*factor);
      if (tensor.label() == L"S") {
        assert(tensor.bra().size() == 1);
        assert(tensor.ket().size() == 1);
        const auto &bra = tensor.bra().at(0);
        const auto &ket = tensor.ket().at(0);
        assert(bra != ket);

        const auto bra_is_ext = ranges::find(external_indices, bra) !=
            ranges::end(external_indices);
        const auto ket_is_ext = ranges::find(external_indices, ket) !=
            ranges::end(external_indices);

        const auto intersection_space = intersection(bra.space(), ket.space());
        assert(intersection_space != IndexSpace::null_instance());

        if (!bra_is_ext && !ket_is_ext) {  // int + int
          const auto new_dummy = idxfac.make(intersection_space);
          add_rules(bra, ket, new_dummy);
        } else if (bra_is_ext && !ket_is_ext) {  // ext + int
          if (includes(ket.space(), bra.space())) {
            add_rule(ket, bra);
          } else {
            add_rule(ket, idxfac.make(intersection_space));
          }
        } else if (!bra_is_ext && ket_is_ext) {  // int + ext
          if (includes(bra.space(), ket.space())) {
            add_rule(bra, ket);
          } else {
            add_rule(bra, idxfac.make(intersection_space));
          }
        }
      }
    }
  }

  return result;
}

/// @return true if made any changes
inline bool apply_index_replacement_rules(
    std::shared_ptr<Product> &product,
    const container::map<Index, Index> &const_replrules,
    const container::set<Index> &external_indices,
    std::set<Index, Index::LabelCompare> &all_indices) {
  // to be able to use map[]
  auto &replrules = const_cast<container::map<Index, Index> &>(const_replrules);

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
      const auto &factor = *it;
      if (factor->is<Tensor>()) {
        auto &tensor = factor->as<Tensor>();
        assert(ranges::none_of(tensor.const_braket(), [](const Index &idx) {
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
      const auto &factor = *it;
      if (factor->is<Tensor>()) {
        bool erase_it = false;
        auto &tensor = factor->as<Tensor>();

        /// replace indices
        pass_mutated &=
            tensor.transform_indices(const_replrules, tag_transformed_indices);

        if (tensor.label() == L"S") {
          const auto &bra = tensor.bra().at(0);
          const auto &ket = tensor.ket().at(0);

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
              if (bra == ket)
                erase_it = true;
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
      const auto &factor = *it;
      if (factor->is<Tensor>()) {
        factor->as<Tensor>().reset_tags();
      }
    }
  }

  // update all_indices
  std::set<Index, Index::LabelCompare> all_indices_new;
  ranges::for_each(
      all_indices, [&const_replrules, &all_indices_new](const Index &idx) {
        auto dst_it = const_replrules.find(idx);
        auto insertion_result = all_indices_new.emplace(dst_it != const_replrules.end() ? dst_it->second
                                                                                        : idx);
      });
  std::swap(all_indices_new, all_indices);

  return mutated;
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

  /*
  std::wcout << "G: " << grand_idxlist.size() << " E: " << ext_idxlist.size() << " I: " << int_idxlist.size() << std::endl;
  for(auto&&i : grand_idxlist)
    std::wcout << " G: " << i.label();
  std::cout << std::endl;
  for(auto&&i : ext_idxlist)
    std::wcout << " E: " << i.label();
  std::cout << std::endl;
  for(auto&&i : int_idxlist)
    std::wcout << " I: " << i.label();
  std::cout << std::endl;
*/

  auto add_spin_labels = [&] (ExprPtr& n){
  if(n->is<Product>()){
    std::wcout << (n->as<Product>()).to_latex() << std::endl;

    bool pass_mutated = false;
    // extract current indices
    std::set<Index, Index::LabelCompare> all_indices;
    ranges::for_each(*n, [&all_indices](const auto &factor) {
      if (factor->template is<Tensor>()) {
        ranges::for_each(factor->template as<const Tensor>().braket(),
                         [&all_indices](const Index &idx) {
                           auto result = all_indices.insert(idx);
                         });
      }
    });

    for(auto&& i: all_indices)
      std::wcout << i.label() << " ";

    std::cout << std::endl;
    container::set<Index> alpha_list;
    container::set<Index> beta_list;
    {
      // For each index, copy .type() and add spin variable
      for (auto &&i : all_indices) {
        auto subscript_label = i.label().substr(i.label().find(L'_') + 1);
        std::wstring subscript_label_ws(subscript_label.begin(), subscript_label.end());

        auto alpha_space = IndexSpace::instance(
            IndexSpace::instance(i.label()).type(), IndexSpace::alpha);
        auto beta_space = IndexSpace::instance(
            IndexSpace::instance(i.label()).type(), IndexSpace::beta);
        alpha_list.insert(alpha_list.end(), Index::make_label_index(alpha_space, subscript_label_ws));
        beta_list.insert(beta_list.end(), Index::make_label_index(beta_space, subscript_label_ws));
      }
    }
      for(auto&& i: alpha_list)
        std::wcout << i.label() << " ";

      for(auto&& i: beta_list)
        std::wcout << i.label() << " ";
      std::cout << std::endl;

/*
    container::map<Index, Index> int_to_alpha;
    container::map<Index, Index> int_to_beta;

    std::transform(beta_list.begin(), beta_list.end(),alpha_list.begin(), alpha_list.end(), std::inserter(int_to_alpha, int_to_alpha.end()),[] (Index a, Index b)
                   { return std::make_pair(a,b);});
    */

    // TODO: Index replacement rules
    std::shared_ptr<Product> expr_product (new Product(n->as<Product>()));
    const auto replacement_rules = sequant::compute_index_replacement_rules(expr_product, ext_idxlist, all_indices);

    // TODO: Apply index replacement rules
//    if(!replacement_rules.empty()){

      pass_mutated = sequant::apply_index_replacement_rules(expr_product, replacement_rules, ext_idxlist, all_indices);
//    }

    std::wcout << "mutated: " << pass_mutated << " -> " << expr_product->to_latex() << std::endl;


  }
  };
  expr->visit(add_spin_labels);

 std::cout << __FILE__ << " " <<  __LINE__ << std::endl;



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
