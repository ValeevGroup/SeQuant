//
// Created by Bimal Gaudel on 12/20/19.
//

#include "factorizer.hpp"
#include "path_tree.hpp"

#include <SeQuant/core/container.hpp>
#include <memory>
#include <SeQuant/core/tensor.hpp>

#include <algorithm>
#include <iostream>
#include <numeric>

namespace sequant::factorize {
namespace detail {

using vec_path_ptr = container::svector<std::shared_ptr<PathTree>>;

ContractionCostResult::ContractionCostResult(
    const ExprPtr& ptr) {  // ptr is a Tensor pointer
  const auto& tensor = ptr->as<Tensor>();
  std::copy(tensor.bra().begin(), tensor.bra().end(),
            std::back_inserter(indices));
  std::copy(tensor.ket().begin(), tensor.ket().end(),
            std::back_inserter(indices));
  std::sort(indices.begin(), indices.end(), sequant::Index::LabelCompare());
}

ContractionCostResult ContractionCostCounter::operator()(const ContractionCostResult& res1,
                                   const ContractionCostResult& res2) const {
  // TODO
  // rewrite the ContractionCostCounter such that the operator() only
  // deals with the reduced ratio of the different sizes of the IndexSpaces
  // eg. imagine calculating 20*20*20*20*100*100*100*100 vs 1*1*1*1*5*5*5*5
  // we care more about relative flops counts than absolute ones
  using index_vec = container::svector<Index>;
  auto label_comp = sequant::Index::LabelCompare();

  auto result = ContractionCostResult{};
  std::set_symmetric_difference(res1.indices.begin(), res1.indices.end(),
                                res2.indices.begin(), res2.indices.end(),
                                std::back_inserter(result.indices),
                                label_comp);

  index_vec indx_union;
  std::set_union(res1.indices.begin(), res1.indices.end(), res2.indices.begin(),
                 res2.indices.end(), std::back_inserter(indx_union),
                 label_comp);

  FlopsType result_flops = 1;
  for (const auto& i : indx_union)
    result_flops *= (map_ptr->find(i.space().type()))->second;

  result.flops = res1.flops + result_flops + res2.flops;
  return result;
}

ContractionCostResult compute_path_cost(const std::shared_ptr<PathTree>& tree,
                        const Product& product,
                        const ContractionCostCounter& counter) {
  if (tree->is_leaf()) {
      return ContractionCostResult{product.factors().at(tree->get_label())};
  }

  auto children = container::svector<ContractionCostResult>{};
  for (const auto& i : tree->get_children())
    children.push_back(compute_path_cost(i, product, counter));

  return std::accumulate(children.begin(), children.end(),
      ContractionCostResult{product.factors().at(tree->get_label())}, counter);
}

void optimal_path(vec_path_ptr& paths, const Product& product,
                  const ContractionCostCounter& counter,
                  std::shared_ptr<PathCostResult>& running_cost) {
  if (paths.size() == 1) {
    auto current_cost = compute_path_cost(paths[0], product, counter);
    if (current_cost.flops < running_cost->flops) {
      running_cost->flops = current_cost.flops;
      running_cost->path = std::make_shared<PathTree>(*paths[0]);
    }
    return;
  }

  for (auto i = 0; i < paths.size(); ++i) {
    for (auto j = i + 1; j < paths.size(); ++j) {
      // form a product of ith and jth tensors
      paths[i]->add_child(paths[j]);

      // set up new arguments for the recursive call
      vec_path_ptr new_args;
      // newly formed prod. is the first element
      new_args.push_back(paths[i]);

      // unlike in Python there's no 'for else' loop in C++
      bool skip_recursive_call = false;

      // all tensors from paths argument should be
      // in the new_args except the ith and the jth
      // since they are absorbed in the product formed above
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
      if (!skip_recursive_call) optimal_path(new_args, product, counter, running_cost);

      paths[i]->pop_last_child();
    }  // for j
  }    // for i
}

ExprPtr path_to_product(const std::shared_ptr<PathTree>& path,
                                 const Product& product) {
  ProductPtr result(new Product{});
  if (path->is_leaf()) {
    result->append(product.factors()[path->get_label()]);
    return result;
  }

  result->append(product.factors()[path->get_label()]);
  for (const auto& i : path->get_children()) {
    auto res = path_to_product(i, product);
    auto& res_product = res->as<Product>();
    // product of single factors is flattened
    if (res_product.factors().size() == 1) result->append(1.0, std::move(res_product.factors()[0]));
    // else append as it is
    else result->append(std::move(res));
  }
  return result;
}

}  // namespace detail

ExprPtr factorize_product(const Product& product,
    const std::shared_ptr<container::map<IndexSpace::Type, size_t>>& ispace_size) {

  using detail::vec_path_ptr;
  using detail::PathTree;

  auto initial_path = [&](Product prod){
    auto fsize = prod.factors().size();
    auto result = vec_path_ptr{};
    for (size_t i=0; i < fsize; ++i) {
      result.push_back(std::make_shared<PathTree>(i));
    }
    return result;
  };
  auto paths = initial_path(product);
  auto running_path_result = std::make_shared<detail::PathCostResult>();

  detail::optimal_path(paths, product,
      detail::ContractionCostCounter{ispace_size},
      running_path_result);
  auto factorized_expr = detail::path_to_product(running_path_result->path, product);
  auto& factorized_prod = factorized_expr->as<Product>();
  factorized_prod.scale(product.scalar());
  return factorized_expr;
}
/// factorize an Expr
/// Note: antisymmetrization Exprs will be gone from the returned result
/// @return ExprPtr
ExprPtr factorize_expr(const ExprPtr& expr_ptr,
                                const std::shared_ptr<container::map<IndexSpace::Type, size_t>>& ispace_map,
                                bool factorize){
  if (expr_ptr->is<Product>()){
    // form a new product by omitting antisymmetrization tensors
    // and return the result of the sequant::factorize::product_factorize on it

    ProductPtr new_prod(new Product{});
    auto& prod = expr_ptr->as<Product>();
    for (auto &&fac: prod) {
      // CATCH: no further checking if the factors are non-Tensor type
      auto& tnsr = fac->as<Tensor>();
      // skip antisym. tensors
      if (tnsr.label() != L"A") new_prod->append(1, fac);
    }
    new_prod->scale(prod.scalar());

    if (factorize)
      return factorize_product(new_prod->as<Product>(), ispace_map);
    else return new_prod;

  } else if (expr_ptr->is<Sum>()){

    SumPtr new_sum(new Sum{});
    auto& sum = expr_ptr->as<Sum>();
    for (auto &&sumand: sum)
      new_sum->append(factorize_expr(sumand, ispace_map, factorize));
    return new_sum;
  } else {
    return expr_ptr;
  }
}
}  // namespace sequant::factorize
