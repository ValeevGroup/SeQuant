#include "single_term_opt.hpp"
#include <SeQuant/core/tensor.hpp>

namespace sequant::factorize {

ExprPtr sto_exhaustive_scan(const ExprPtr& expr, size_t nocc, size_t nvirt) {
  return repack_prod(expr, OptimalRootedTree{expr, nocc, nvirt}.tree());
}

ExprPtr repack_prod(const ExprPtr& expr, const factorize::rooted_tree& tree,
                    bool scale) {
  if (tree.children.empty()) {
    auto prod = ex<Product>(Product{expr->at(tree.label)});
    if (scale) {
      auto& prod_deref = prod->as<Product>();
      prod_deref.scale(expr->as<Product>().scalar());
    }
    return prod;
  }
  auto prod = std::make_shared<Product>();

  if (scale) prod->scale(expr->as<Product>().scalar());

  prod->append(1., expr->at(tree.label));
  for (const auto& tt : tree.children) {
    auto child_prod = repack_prod(expr, tt, false);
    if (child_prod->size() > 1)
      prod->append(child_prod);
    else
      prod->append(1., child_prod);
  }
  return prod;
}

OptimalRootedTree::OptimalRootedTree(const ExprPtr& prod, size_t nocc,
                                     size_t nvirt)
    : prod{prod},
      flops{std::numeric_limits<decltype(flops)>::max()},
      nocc{nocc},
      nvirt{nvirt} {}

void OptimalRootedTree::operator()(const rooted_tree& tree) {
  if (auto&& res = ops_count(prod, tree, nocc, nvirt); res.flops < flops) {
    tree_ = tree;
    flops = res.flops;
  }
}

rooted_tree OptimalRootedTree::tree() {
  if (tree_.has_value()) return tree_.value();
  if (prod->is<Tensor>()) return rooted_tree{0};

  auto init = std::vector<rooted_tree>(prod->size());
  for (size_t ii = 0; ii < prod->size(); ++ii) init[ii] = rooted_tree{ii};

  enumerate_eval_sequence(init, std::ref(*this));

  return tree_.value();
}

}  // namespace sequant::factorize
