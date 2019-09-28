//
// Created by Eduard Valeyev on 2019-02-27.
//

#ifndef SEQUANT_SPIN_HPP
#define SEQUANT_SPIN_HPP

#include <algorithm>
#include <string>
#include <tuple>
#include "../../core/expr.hpp"
#include "../../core/tensor.hpp"

#define SPINTRACE_PRINT 0
#define PRINT_PERMUTATION_LIST 0

namespace sequant {

class spinIndex : public Index {
 private:
};

/// @param expr input expression
/// @param ext_index_groups groups of external indices
/// @return the expression with spin integrated out
ExprPtr spintrace(ExprPtr expr,
                  std::initializer_list<IndexList> ext_index_groups = {{}}) {
  std::cout << "In sequant::spintrace: " << std::endl;

  container::set<Index> grand_idxlist;
  auto collect_indices = [&](const ExprPtr& expr) {
    if (expr->is<Tensor>()) {
      ranges::for_each(expr->as<Tensor>().const_braket(),
                       [&](const Index& idx) { grand_idxlist.insert(idx); });
    }
  };
  expr->visit(collect_indices);

  std::cout << "grand_idxlist size: " << grand_idxlist.size() << std::endl;

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

#if SPINTRACE_PRINT
  //  List of all Index ( label + IndexSpace )
  for (auto&& idx : grand_idxlist) {
    std::wcout << idx.to_latex() << std::endl;
  }
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

  // TODO: Add same terms

  // TODO: Zero-out non-spinconserving terms

  // TODO: Reindex and return

  return nullptr;
}

}  // namespace sequant

#endif  // SEQUANT_SPIN_HPP
