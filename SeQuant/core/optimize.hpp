#ifndef SEQUANT_OPTIMIZE_OPTIMIZE_HPP
#define SEQUANT_OPTIMIZE_OPTIMIZE_HPP

#include <ios>
#include <limits>
#include <utility>

#include <SeQuant/core/clone_packed.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/eval_seq.hpp>
#include <SeQuant/core/tensor_network.hpp>

namespace sequant {
/// Optimize an expression assuming the number of virtual orbitals
/// greater than the number of occupied orbitals.

/// \param expr Expression to be optimized.
/// \return EvalNode object.
// EvalNode optimize(ExprPtr const& expr);

namespace opt {

namespace detail {
///
/// For n number of bits, assuming all of them are on,
/// bipartition them into all possibilities, except for
/// the trivial (all zero bits, all one bits) partition
/// eg. for n = 3
///       (001, 110)
///       (010, 101)
///       (011, 100)
///
template <typename F,
          std::enable_if_t<std::is_invocable_v<F, size_t, size_t>, bool> = true>
void scan_biparts_all_bits(size_t n, F&& scanner) {
  auto ulim = (1 << n);
  for (auto i = 1; i < ulim / 2; ++i)
    std::invoke(std::forward<F>(scanner), i, (ulim - 1 - i));
}

///
/// given positions of bits that are on,
/// P = {i, j,...,m}
/// generate binary partitions of the positions
/// such as {m} and {i, j, ...} (= P - {m}) and so on
/// except for the trivial {i, j,...,m} {} partition
///
template <typename F,
          std::enable_if_t<std::is_invocable_v<F, size_t, size_t>, bool> = true>
void scan_biparts_some_bits(std::vector<size_t> const& bs, F&& scanner) {
  scan_biparts_all_bits(bs.size(), [&scanner, &bs](size_t a1, size_t _) {
    size_t p1{0}, p2{0};
    for (auto i = 0; i < bs.size(); ++i) {
      // if ith bit is on, ith elem in bs is included in p1
      //                       else it is included in p2
      if ((1 << i) & a1)
        p1 |= (1 << bs[i]);
      else
        p2 |= (1 << bs[i]);
    }
    std::invoke(std::forward<F>(scanner), p1, p2);
  });
}

/// given a number @c n, return a vector of ON bit positions
/// only first num_bits bits will be checked from right to left
/// in the bit representation of @c n
/// By default the value of @c num_bits is the total number of bits in
/// the representation of @c n.
std::vector<size_t> on_bits_pos(size_t n, size_t num_bits = sizeof(size_t) * 8);

}  // namespace detail

///
/// Omit the first factor from the top level product from given expression.
/// Intended to drop "A" and "S" tensors from CC amplitudes as a preparatory
/// step for evaluation of the amplitudes.
///
ExprPtr tail_factor(ExprPtr const& expr) noexcept;

///
///
/// Pulls out scalar to the top level from a nested product.
/// If @c expr is not Product, does nothing.
void pull_scalar(sequant::ExprPtr expr) noexcept;

///
/// Result of the single term optimization of a term.
/// Holds operations count.
///
/// Iterable of one or more binary_expr<EvalExpr> root node pointers
/// that lead to the same operations count.
///
/// ie. degenerate evaluations leading to the minimal operations count are
/// stored as binary tree nodes.
///
struct STOResult {
  AsyCost cost;

  container::vector<EvalNode> optimal_seqs;
};

///
/// \tparam IdxToSz map-like {IndexSpace : size_t}
/// \param idxsz see @c IdxToSz
/// \param commons Index objects
/// \param diffs   Index objects
/// \return flops count
/// @note @c commons and @c diffs have unique indices individually as well as
///       combined
template <typename IdxToSz,
          std::enable_if_t<std::is_invocable_r_v<size_t, IdxToSz, Index>,
                           bool> = true>
double ops_count(IdxToSz const& idxsz, container::vector<Index> const& commons,
                 container::vector<Index> const& diffs) {
  double ops = 1.0;
  for (auto&& idx : ranges::views::concat(commons, diffs))
    ops *= std::invoke(idxsz, idx);
  // ops == 1.0 implies both commons and diffs empty
  return ops == 1.0 ? 0 : ops;
}

/// Perform single term optimization on a product. Deprecated.
///
/// @return STOResult
template <typename F = std::function<bool(EvalNode const&)>,
          std::enable_if_t<std::is_invocable_r_v<bool, F, EvalNode const&>,
                           bool> = true>
[[deprecated]] STOResult single_term_opt(
    Product const& prod, F&& pred = [](auto const&) { return true; }) {
  using ranges::to_vector;
  using ranges::views::iota;
  using ranges::views::take;
  using ranges::views::transform;
  using seq_t = EvalSeq<size_t>;

  struct {
    Product const& facs;

    ExprPtr operator()(size_t pos) { return facs.factor(pos)->clone(); }

    ExprPtr operator()(ExprPtr lf, ExprPtr rf) {
      auto p = Product{};
      p.append(std::move(lf));
      p.append(std::move(rf));

      return ex<Product>(std::move(p));
    }
  } fold_prod{prod};  // struct

  auto result = STOResult{AsyCost::max(), {}};

  auto finder = [&result, &fold_prod, &pred, prod](auto const& seq) {
    auto expr = seq.evaluate(fold_prod);

    if (prod.scalar() != 1.) {
      if (!expr->template is<Product>())  // in case expr is non-product
        expr = ex<Product>(Product{std::move(expr)});
      *expr *= Constant{prod.scalar()};
    }

    auto node = to_eval_node(expr);

    auto cost = asy_cost(node, std::forward<F>(pred));

    if (cost == result.cost) {
      result.optimal_seqs.emplace_back(std::move(node));
    } else if (cost < result.cost) {
      result.optimal_seqs.clear();
      result.optimal_seqs.emplace_back(std::move(node));
      result.cost = std::move(cost);
    } else {
      // cost > optimal cost. do nothing.
    }
  };  // finder

  auto init_seq = iota(size_t{0}) | take(prod.size()) |
                  transform([](auto x) { return seq_t{x}; }) |
                  ranges::to_vector;

  enumerate_eval_seq<size_t>(init_seq, finder);
  return result;
}

///
/// any element in the vector belongs to the integral range [-1,N]
/// where N is the length of the [Expr] (ie. the iterable of expressions)
///   * only applicable for binary evaluations
///   * the integer -1 can appear in certain places: it implies the binary
///     operation between the last two expressions
///   * eg.
///         * {0,1,-1,2,-1} => ( (e[0], e[1]), e[2])
///         * {0,1,-1,2,3,-1,-1} => ((e[0], e[1]), (e[2],e[3]))
///
using eval_seq_t = container::vector<int>;

///
/// Represents a result of optimization on a range of expressions
/// for a binary evaluation
///
struct OptRes {
  /// Free indices remaining upon evaluation
  container::vector<sequant::Index> indices;

  /// The flops count of evaluation
  double flops;

  /// The evaluation sequence
  eval_seq_t sequence;
};

/// returns a pair of index vectors
/// first element of the pair is the vector of common indices compared by labels
/// second element of the pair is the set symmetric difference of the input
/// index vectors if either of the input index container is empty, the result is
/// a pair of empty vectors
/// @note I1 and I2 containers are assumed to be sorted by using
/// Index::LabelCompare{};
template <typename I1, typename I2>
std::pair<container::vector<Index>, container::vector<Index>> common_indices(
    I1 const& idxs1, I2 const& idxs2) {
  container::vector<Index> i1vec(ranges::begin(idxs1), ranges::end(idxs1)),
      i2vec(ranges::begin(idxs2), ranges::end(idxs2));
  if (i1vec.empty() || i2vec.empty()) return {{}, {}};

  container::vector<Index> commons, diffs;
  std::set_intersection(std::begin(i1vec), std::end(i1vec), std::begin(i2vec),
                        std::end(i2vec), std::back_inserter(commons),
                        Index::LabelCompare{});
  std::set_symmetric_difference(
      std::begin(i1vec), std::end(i1vec), std::begin(i2vec), std::end(i2vec),
      std::back_inserter(diffs), Index::LabelCompare{});
  return {commons, diffs};
}

template <typename IdxToSz,
          std::enable_if_t<std::is_invocable_r_v<size_t, IdxToSz, Index>,
                           bool> = true>
eval_seq_t single_term_opt_v2(TensorNetwork const& network,
                              IdxToSz const& idxsz) {
  // number of terms
  auto const nt = network.tensors().size();
  if (nt == 1) return eval_seq_t{0};
  if (nt == 2) return eval_seq_t{0, 1, -1};
  auto nth_tensor_indices = container::vector<container::vector<Index>>{};
  nth_tensor_indices.reserve(nt);

  for (auto i = 0; i < nt; ++i) {
    auto bk = container::vector<Index>{};
    for (auto idx : braket(*network.tensors().at(i))) bk.push_back(idx);

    ranges::sort(bk, Index::LabelCompare{});
    nth_tensor_indices.emplace_back(std::move(bk));
  }

  container::vector<OptRes> result((1 << nt), OptRes{{}, 0, {}});

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
        on_bits, [&result = std::as_const(result), &tindices, &idxsz, &cost,
                  &p1, &p2](auto p1_, auto p2_) {
          auto [commons, diffs] =
              common_indices(result[p1_].indices, result[p2_].indices);
          auto new_cost = ops_count(idxsz, commons, diffs) + result[p1_].flops +
                          result[p2_].flops;
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
      seq = eval_seq_t{static_cast<int>(power_pos++)};
    } else {
      // cost set
      // tindices set
      seq = ranges::views::concat(result[p1].sequence, result[p2].sequence) |
            ranges::to<eval_seq_t>;
      seq.push_back(-1);
    }

    result[n].flops = cost;
    result[n].indices = std::move(tindices);
    result[n].sequence = std::move(seq);
  }

  return result[(1 << nt) - 1].sequence;
}

///
/// \param prod  Product to be optimized.
/// \param idxsz An invocable object that maps an Index object to size.
/// \return Parenthesized product expression.
///
/// @note @c prod is assumed to consist of only Tensor expressions
///
template <typename IdxToSz,
          std::enable_if_t<std::is_invocable_v<IdxToSz, Index>, bool> = true>
ExprPtr single_term_opt_v2(Product const& prod, IdxToSz const& idxsz) {
  if (prod.factors().size() < 3) return clone_packed(prod);
  auto seq = single_term_opt_v2(TensorNetwork{prod}, idxsz);
  auto result = container::vector<ExprPtr>{};
  for (auto i : seq)
    if (i == -1) {
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

///
/// \param expr  Expression to be optimized.
/// \param idxsz An invocable object that maps an Index object to size.
/// \return Optimized expression converted to EvalNode.
template <typename IdxToSz,
          std::enable_if_t<std::is_invocable_v<IdxToSz, Index>, bool> = true>
EvalNode optimize(const ExprPtr& expr, IdxToSz const& idxsz) {
  using ranges::views::transform;
  if (expr->is<Tensor>())
    return to_eval_node(expr);
  else if (expr->is<Product>()) {
    auto opt_expr = opt::single_term_opt_v2(expr->as<Product>(), idxsz);
    return to_eval_node(opt_expr);
  } else if (expr->is<Sum>()) {
    auto smands = *expr | transform([&idxsz](auto const& s) {
      return to_expr(optimize(s, idxsz));
    }) | ranges::to_vector;

    return to_eval_node(ex<Sum>(Sum{smands.begin(), smands.end()}));
  } else
    throw std::runtime_error{"optimization attempted on unsupported Expr type"};
}

}  // namespace sequant

#endif  // SEQUANT_OPTIMIZE_OPTIMIZE_HPP
