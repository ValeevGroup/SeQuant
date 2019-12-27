//
// Created by Bimal Gaudel on 12/20/19.
//

#include "factorizer.hpp"
#include "path_tree.hpp"

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/tensor.hpp>

#include <algorithm>
#include <iostream>
#include <numeric>

namespace sequant {
namespace factorize {
namespace detail {

using vec_path_ptr = sequant::container::svector<std::shared_ptr<PathTree>>;

CostResult::CostResult(
    const ExprPtr& ptr) {  // ptr is a sequant::Tensor pointer
  const auto& tensor = ptr->as<sequant::Tensor>();
  std::copy(tensor.bra().begin(), tensor.bra().end(),
            std::back_inserter(indices));
  std::copy(tensor.ket().begin(), tensor.ket().end(),
            std::back_inserter(indices));
}

CostResult CostCounter::operator()(const CostResult& res1,
                                   const CostResult& res2) {
  using index_vec = sequant::container::svector<sequant::Index>;

  auto result = CostResult{};

  std::set_symmetric_difference(res1.indices.begin(), res1.indices.end(),
                                res2.indices.begin(), res2.indices.end(),
                                std::back_inserter(result.indices));

  index_vec indx_union;
  std::set_union(res1.indices.begin(), res1.indices.end(), res2.indices.begin(),
                 res2.indices.end(), std::back_inserter(indx_union));

  unsigned long flops = 1;
  for (const auto& i : indx_union) {
    if (i.space() == sequant::IndexSpace::active_occupied)
      flops *= nocc;
    else if (i.space() == sequant::IndexSpace::active_unoccupied)
      flops *= nvirt;
    else
      throw "Only know about occupied and unoccupied Index spaces.";
  }
  result.flops += res1.flops + flops + res2.flops;
  return result;
}

CostResult compute_cost(const std::shared_ptr<PathTree>& tree,
                        const sequant::Product& product,
                        const CostCounter& counter) {
  if (tree->is_leaf())
    return CostResult{product.factors().at(tree->get_label())};

  auto children = sequant::container::svector<CostResult>{};
  for (const auto& i : tree->get_children()) {
    children.push_back(compute_cost(i, product, counter));
  }
  return std::accumulate(children.begin(), children.end(),
                         CostResult{product.factors().at(tree->get_label())},
                         counter);
}

void path_scanner(vec_path_ptr& paths, const sequant::Product& product,
                  const CostCounter& counter) {
  if (paths.size() == 1) {
    std::wcout << paths[0]->print_tree() << L"  "
               << compute_cost(paths[0], product, counter).flops << "\n";
  }

  for (auto i = 0; i < paths.size(); ++i) {
    for (auto j = i + 1; j < paths.size(); ++j) {
      // form a product of ith and jth tensors
      paths[i]->add_child(paths[j]);

      // set up new arguments for the recursive call
      vec_path_ptr new_args;
      // newly formed prod. is the first element
      new_args.push_back(paths[i]);

      // all tensors from paths argument should be
      // in the new_args except the ith and the jth
      // since they are absorbed in the product formed above
      bool skip_recursive_call =
          false;  // unlike in Python there's no 'for else' loop in C++
      for (auto k = 0; k < paths.size(); ++k) {
        if ((k != i) && (k != j)) {
          // remove redundancy by lexicographic comparison
          if ((!paths[k]->is_leaf()) &&
              (paths[k]->get_label() < paths[i]->get_label())) {
            skip_recursive_call = true;
            new_args.clear();
            break;
          }
          new_args.push_back(paths[k]);
        }
      }
      if (!skip_recursive_call) path_scanner(new_args, product, counter);

      paths[i]->pop_last_child();
    }  // for j
  }    // for i
}

}  // namespace detail
}  // namespace factorize
}  // namespace sequant
