#include "single_term_opt.hpp"
#include <SeQuant/core/tensor.hpp>

namespace sequant::factorize {

ExprPtr single_term_opt(
    const ExprPtr& expr, size_t nocc, size_t nvirt,
    const std::function<ExprPtr(const ExprPtr&, size_t, size_t)>& backend) {
  if (expr->is<Product>()) {
    container::svector<ExprPtr> operators, nonOperators;
    ranges::partition_copy(*expr, ranges::back_inserter(operators),
                           ranges::back_inserter(nonOperators),
                           [](const auto& x) {  //
                             if (!x->template is<Tensor>()) return true;
                             auto l = x->template as<Tensor>().label();
                             return l == L"A" || l == L"P";
                           });

    auto factors =
        ranges::views::concat(
            operators,
            *backend(ex<Product>(nonOperators.begin(), nonOperators.end()),
                     nocc, nvirt)) |
        ranges::to<container::svector<ExprPtr>>;

    return ex<Product>(
        Product{expr->as<Product>().scalar(), factors.begin(), factors.end()});
  }

  auto sum = Sum{};
  for (const auto& x : *expr)
    sum.append(single_term_opt(x, nocc, nvirt, backend));

  return sum.clone();
}

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
