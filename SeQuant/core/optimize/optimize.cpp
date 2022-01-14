#include <SeQuant/core/optimize.hpp>
#include <SeQuant/core/clone_packed.hpp>

namespace sequant {

namespace opt {

namespace detail {
std::vector<size_t> on_bits_pos(size_t n, size_t num_bits) {
  auto result = std::vector<size_t>{};
  result.reserve(num_bits);
  for (auto i = 0; i < num_bits; ++i)
    if (n & (1 << i)) result.push_back(i);
  return result;
}
}  // namespace detail

ExprPtr tail_factor(ExprPtr const& expr) noexcept {
  if (expr->is<Tensor>())
    return expr->clone();

  else if (expr->is<Product>()) {
    auto scalar = expr->as<Product>().scalar();
    auto facs = ranges::views::tail(*expr);
    return ex<Product>(Product{scalar, ranges::begin(facs), ranges::end(facs)});
  } else {
    // sum
    auto summands = *expr | ranges::views::transform(
                                [](auto const& x) { return tail_factor(x); });
    return ex<Sum>(Sum{ranges::begin(summands), ranges::end(summands)});
  }
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

double log_flops(container::vector<Index> const& commons,
                 container::vector<Index> const& diffs, double log_nocc,
                 double log_nvirt) {
  double res = 0;
  for (auto&& idx : ranges::views::concat(commons, diffs))
    if (idx.space() == IndexSpace::active_occupied)
      res += log_nocc;
    else if (idx.space() == IndexSpace::active_unoccupied)
      res += log_nvirt;
    else
      throw std::runtime_error(
          "Unexpected IndexSpace encountered while computing flops.");
  return res;
}

eval_seq_t single_term_opt_v2(TensorNetwork const& network, size_t nocc,
                              size_t nvirt) {
  // number of terms
  auto const nt = network.tensors().size();
  if (nt == 1) return eval_seq_t{0};
  if (nt == 2) return eval_seq_t{0,1,-1};
  auto nth_tensor_indices = container::vector<container::vector<Index>>{};
  nth_tensor_indices.reserve(nt);

  for (auto i = 0; i < nt; ++i) {
    auto bk = container::vector<Index>{};
    for (auto idx : braket(*network.tensors().at(i))) bk.push_back(idx);

    ranges::sort(bk, Index::LabelCompare{});
    nth_tensor_indices.emplace_back(std::move(bk));
  }
  double const log_nocc = std::log10(nocc);
  double const log_nvirt = std::log10(nvirt);

  container::vector<OptRes> result((1<<nt), OptRes{{},0,{}});

  // power_pos is used, and incremented, only when the
  // result[1<<0]
  // result[1<<1]
  // result[1<<2]
  // and so on are set
  size_t power_pos = 0;
  for (auto n = 1; n < (1 << nt); ++n) {
    double cost = std::numeric_limits<double>::max();
    auto const on_bits = detail::on_bits_pos(n, nt);
    size_t p1 = 0, p2 = 0;
    container::vector<Index> tindices{};
    detail::scan_biparts_some_bits(
        on_bits, [&result = std::as_const(result), &tindices, log_nocc,
                  log_nvirt, &cost, &p1, &p2](auto p1_, auto p2_) {
          auto [commons, diffs] =
              common_indices(result[p1_].indices, result[p2_].indices);
          auto new_cost = log_flops(commons, diffs, log_nocc, log_nvirt) +
                          result[p1_].flops + result[p2_].flops;
          if (new_cost < cost) {
            cost = new_cost;
            tindices = std::move(diffs);
            p1 = p1_;
            p2 = p2_;
          }
        });  //

    auto seq = eval_seq_t{};
    if (tindices.empty()) {
      cost = 0;
      tindices = std::move(nth_tensor_indices[power_pos]);
      // seq = make_sequence(++power_pos);
      seq = eval_seq_t{static_cast<int>(power_pos++)};
    } else {
      // cost set
      // tindices set
      seq = ranges::views::concat(result[p1].sequence, result[p2].sequence) | ranges::to<eval_seq_t>;
      seq.push_back(-1);
    }

    result[n].flops = cost;
    result[n].indices = std::move(tindices);
    result[n].sequence = std::move(seq);
  }

  return result[(1 << nt) - 1].sequence;
}

ExprPtr single_term_opt_v2(Product const& prod, size_t nocc, size_t nvirt) {
  if (prod.factors().size() < 3) return clone_packed(prod);
  auto seq = single_term_opt_v2(TensorNetwork{prod}, nocc, nvirt);
  auto result = container::vector<ExprPtr>{};
  for (auto i: seq)
    if (i==-1){
      auto rexpr = *result.rbegin();
      result.pop_back();
      auto lexpr = *result.rbegin();
      result.pop_back();
      auto p = Product{};
      p.append(lexpr);
      p.append(rexpr);
      result.push_back(clone_packed(p));
    } else {
      result.push_back(prod.at(i));
    }

  (*result.rbegin())->as<Product>().scale(prod.scalar());
  return *result.rbegin();
}

}  // namespace opt

EvalNode optimize(const ExprPtr& expr) {
  static const size_t NOCC = 10;
  static const size_t NVIRT = 100;

  using ranges::views::transform;
  if (expr->is<Tensor>())
    return to_eval_node(expr);
  else if (expr->is<Product>()) {
    // return *(opt::single_term_opt(expr->as<Product>()).optimal_seqs.begin());
    auto opt_expr = opt::single_term_opt_v2(expr->as<Product>(), NOCC, NVIRT);
    return to_eval_node(opt_expr);
  } else if (expr->is<Sum>()) {
    auto smands = *expr | transform([](auto const& s) {
      return to_expr(optimize(s));
    }) | ranges::to_vector;

    return to_eval_node(ex<Sum>(Sum{smands.begin(), smands.end()}));
  } else
    throw std::runtime_error{"optimization attempted on unsupported Expr type"};
}

}  // namespace sequant
