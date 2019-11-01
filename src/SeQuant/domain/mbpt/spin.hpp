//
// Created by Eduard Valeyev on 2019-02-27.
//

#ifndef SEQUANT_SPIN_HPP
#define SEQUANT_SPIN_HPP

#include <algorithm>

#include <string>
#include <tuple>

// #include <range/v3/algorithm/for_each.hpp>

#include "SeQuant/core/expr.hpp"
#include "SeQuant/core/tensor.hpp"

#define SPINTRACE_PRINT 1
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
  std::cout << "In sequant::spintrace." << std::endl;

  container::set<Index> grand_idxlist;
  auto collect_indices = [&](const ExprPtr& expr) {
    if (expr->is<Tensor>()) {
      //  std::wcout << expr->to_latex() << std::endl;
      //  std::wcout << expr->as<Tensor>().bra_rank() << " " <<
      //  expr->as<Tensor>().ket_rank() << std::endl;
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
  std::cout << std::endl;

  for (auto i = grand_idxlist.begin(); i != grand_idxlist.end(); ++i)
    std::wcout << i->to_latex() << std::endl;

#endif

  // TODO make list of groups, generate all spin tuples, apply and replace ...

  // Generate list of spin quantum numbers from Grand Index List
  container::set<Index> alpha_list;
  container::set<Index> beta_list;
  {
    //    Index::reset_tmp_index();

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

  std::wcout << "Before: " << expr->to_latex() << std::endl;
  expand(expr);
  //  std::wcout << "After: " << expr->to_latex() << std::endl;
  //  for (auto &&i : expand(expr)){
  //    std::wcout << i << std::endl;
  //  }

  //  TEST 2
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

  //  auto expr2 = expr->clone();
  //  std::wcout << "expr2:\n" << expr2->to_latex() << std::endl;

  struct latex_visitor {
    void operator()(const std::shared_ptr<sequant::Expr>& expr) {
      result += expr->to_latex();
    }
    std::wstring result{};
  };

  latex_visitor v1{};
  expr->visit(v1);
  //  std::wcout << "v1.result = " << v1.result << std::endl;

  latex_visitor v2{};
  expr->visit(v2, true);
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
  expr->visit(collect_bra_indices);

  //  Get sub expressions from an expression:
  auto print_subexpr = [](const auto& n) {
    std::wcout << "&" << n->to_latex() << L"\\\\" << std::endl;
  };
  ranges::for_each(expr->begin(), expr->end(), print_subexpr);

  // The bra and ket indices from each tensor are separated here
  auto print_subexpr_index = [&](ExprPtr& n) {
    if (n->is<Tensor>()) {
      std::wcout << "&" << n->to_latex() << L"\\\\" << std::endl;

      //    Get a list of indices for a tensor. Don't separate into bra/ket. Add
      //    spin attribute to the index and regenetate the tensor

      for (auto& i : n->as<Tensor>().const_braket()) {
        std::wcout << sequant::to_latex(i) << " "
                   << IndexSpace::instance(i.label()).type();
        //        std::wcout << temp.label() << " {" <<
        //        IndexSpace::instance(temp.label()).type()
        //                   << "," << IndexSpace::instance(temp.label()).qns()
        //                   << "}"
        //                   << std::endl;
      }
      //      for (auto& i : n->as<Tensor>().ket()) {
      //        std::wcout << sequant::to_latex(i) << " ";
      //      }
      //      std::cout << std::endl;
      //      for (auto& i : n->as<Tensor>().bra())
      //        std::wcout << sequant::to_latex(i) << " ";
      std::cout << std::endl;
    }
  };
  expr->visit(print_subexpr_index);
  //  ranges::for_each(expr->begin(), expr->end(), print_subexpr_index);

  // TODO: Add same terms

  // TODO: Zero-out non-spinconserving terms

  // TODO: Reindex and return

  return nullptr;
}

}  // namespace sequant

#endif  // SEQUANT_SPIN_HPP
