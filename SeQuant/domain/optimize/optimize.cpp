#include "optimize.hpp"

namespace sequant::optimize {

bin_eval_expr_flops_counter::bin_eval_expr_flops_counter(
    size_t no, size_t nv, const container::set<size_t>& imeds)
    : counter{no, nv}, imed_hashes{imeds} {}

size_t bin_eval_expr_flops_counter::operator()(
    utils::binary_node<utils::eval_expr> const& node) const {
  return 0;
}

size_t bin_eval_expr_flops_counter::operator()(
    utils::binary_node<utils::eval_expr> const& node, size_t lflops,
    size_t rflops) const {
  if (imed_hashes.contains(node->hash())) return 0;

  return counter(node, lflops, rflops);
}

sto_result single_term_opt(Product const& prod, size_t nocc, size_t nvirt,
                           container::set<size_t> const& imeds_hash,
                           bool canon) {
  using seq_t = utils::eval_seq<size_t>;

  struct {
    Product const& facs;

    ExprPtr operator()(size_t pos) { return facs.factor(pos)->clone(); }

    ExprPtr operator()(ExprPtr lf, ExprPtr rf) {
      auto p = Product{};
      p.append(lf);
      p.append(rf);

      return ex<Product>(std::move(p));
    }
  } fold_prod{prod};  // struct

  const auto counter = bin_eval_expr_flops_counter{nocc, nvirt, imeds_hash};

  auto result = sto_result{std::numeric_limits<size_t>::max(), {}};

  auto finder = [&result, &fold_prod, &counter, &prod, canon](const auto& seq) {
    auto expr = seq.evaluate(fold_prod);
    if (canon) {
      expr->canonicalize();
      pull_scalar(expr);

      if (prod.scalar() != 1.) {
        if (!expr->template is<Product>())  // in case expr is non-product
          expr = ex<Product>(Product{expr});
        *expr *= Constant{prod.scalar()};
      }
    }  // if canon
    auto node = utils::binarize_expr(expr);

    auto flops = node.evaluate(counter);

    if (flops == result.ops) {
      result.optimal_seqs.emplace_back(std::move(node));
    } else if (flops < result.ops) {
      result.optimal_seqs.clear();
      result.optimal_seqs.emplace_back(std::move(node));
      result.ops = flops;
    } else {
      // flops > optimal flops. do nothing.
    }
  };  // finder

  auto init_seq = ranges::views::iota(size_t{0}) |
                  ranges::views::take(prod.size()) |
                  ranges::views::transform([](auto x) { return seq_t{x}; }) |
                  ranges::to_vector;

  utils::enumerate_eval_seq<size_t>(init_seq, finder);

  return result;
}

void pull_scalar(sequant::ExprPtr expr) noexcept {
  using sequant::Product;
  if (!expr->is<Product>()) return;
  auto& prod = expr->as<Product>();

  auto scal = prod.scalar();
  for (auto&& x : *expr)
    if (x->is<Product>()) {
      auto& p = x->as<Product>();
      scal *= p.scalar();
      p.scale(1.0 / p.scalar());
    }

  prod.scale(1.0 / prod.scalar());
  prod.scale(scal);
}

}  // namespace sequant::optimize
