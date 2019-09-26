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

  for (auto&& idx : grand_idxlist) {

    // Checking IndexSpace
    if(IndexSpace::instance(idx.label()) == IndexSpace::nullqns){
//      std::cout << "nullqns" << std::endl;
      std::cout << IndexSpace::instance(idx.label()).qns() << std::endl;
    }
    std::wcout << idx.to_latex() << std::endl;
  }

  container::set<Index> ext_idxlist;
  for (auto&& idxgrp : ext_index_groups) {
    for (auto&& idx : idxgrp) {
      auto result = ext_idxlist.insert(idx);
      assert(result.second);
    }
  }

  std::cout << std::endl;

  container::set<Index> int_idxlist;
  for (auto&& gidx : grand_idxlist) {
    if (ext_idxlist.find(gidx) == ext_idxlist.end()) {
      int_idxlist.insert(gidx);
    }
  }

  for (auto&& idx : int_idxlist) {
    std::wcout << idx.to_latex() << std::endl;
  }

  // TODO make list of groups, generate all spin tuples, apply and replace ...

  // TODO: Add spin quantum numbers and generate tuples

  //   Using tuples
  //  auto spinTuple = std::make_tuple('A', 'B');
  //  std::cout << "spinTuple: " << std::get<0>(spinTuple) << " "
  //            << std::get<1>(spinTuple) << std::endl;

  //  Using vector
  std::vector<std::string> spinVector = {"A", "B"};
  std::cout << "spinVector: " << spinVector[0] << " " << spinVector[1]
            << std::endl;
  std::cout << "Length of index list: " << int_idxlist.size() << std::endl;

  //  std::vector<Index> tupleList(int_idxlist.size());
  //  std::cout << "tupleList size: " << tupleList.size() << std::endl;

  auto tupleSize = 0;  // Counter for size of the list of tuples
  std::vector<std::string> tupleList;

  //  Create a list of alpha spins the size of number of index
  std::string strTuple;
  for (size_t i = 0; i < int_idxlist.size(); i++) {
    strTuple.append(spinVector[0]);
  }
  tupleSize++;
  tupleList.push_back(strTuple);

  //  do {
  //    tupleSize++;
  //    tupleList.push_back(strTuple);
  //        std::cout << strTuple << '\n';
  //  } while (std::next_permutation(strTuple.begin(), strTuple.end()));

  // Replace one spin at a time and permute
  for (size_t i = 0; i < int_idxlist.size(); i++) {
    strTuple.replace(int_idxlist.size() - i - 1, 1, spinVector[1]);
    std::sort(strTuple.begin(), strTuple.end());
    do {
      tupleSize++;
      tupleList.push_back(strTuple);
      //          std::cout << strTuple << '\n';
    } while (std::next_permutation(strTuple.begin(), strTuple.end()));
    //        std::cout << std::endl;
    //    std::cout << i << " " << strTuple << std::endl;
  }

  std::cout << "tupleList.size(): " << tupleList.size() << std::endl;

  // TODO: Merge tupleList elements and Index

  for (auto&& i : tupleList) {
    //    std::cout << i << std::endl;
    std::vector<char> spinIndexList(i.begin(), i.end());
    for(auto&&j : spinIndexList){

    }
  }

  // TODO: Add same terms

  // TODO: Zero-out non-spinconserving terms

  // TODO: Reindex and return

  return nullptr;
}

}  // namespace sequant

#endif  // SEQUANT_SPIN_HPP
