#ifndef SEQUANT_OPTIMIZE_OPTIMIZE_HPP
#define SEQUANT_OPTIMIZE_OPTIMIZE_HPP

#include <ios>
#include <iostream>
#include <limits>
#include <utility>

#include <SeQuant/core/clone_packed.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/tensor_network.hpp>

namespace sequant {
/// Optimize an expression assuming the number of virtual orbitals
/// greater than the number of occupied orbitals.

/// \param expr Expression to be optimized.
/// \return EvalNode object.
// EvalNode optimize(ExprPtr const& expr);

namespace opt {

namespace {

///
/// Non-trivial, unique bipartitions of the bits in an unsigned integral.
///  eg.
///  decimal (binary) => [{decimal (binary)}...]
///  -------------------------------------------
///         3 (11) => [{1 (01), 2 (10)}]
///      11 (1011) => [{1  (0001), 10 (1010)},
///                     {2 (0010), 9 (1001)},
///                     {3 (0011), 8 (1000)}]
///      0 (0)     => [] (empty: no partitions possible)
///      2 (10)    => [] (empty)
///      4 (100)   => [] (empty)
///
/// \tparam I Unsigned integral type.
/// \param n Represents a bit set.
/// \param func func is function that takes two arguments of type I.
///
template <
    typename I, typename F,
    typename = std::enable_if_t<std::is_integral_v<I> && std::is_unsigned_v<I>>,
    typename = std::enable_if_t<std::is_invocable_v<F, I, I>>>
void biparts(I n, F&& func) {
  if (n == 0) return;
  I const h = static_cast<I>(std::floor(n / 2.0));
  for (auto n_ = 1; n_ <= h; ++n_) {
    auto const l = n & n_;
    auto const r = (n - n_) & n;
    if ((l | r) == n) std::invoke(std::forward<F>(func), l, r);
  }
}

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
double ops_count(IdxToSz const& idxsz, container::svector<Index> const& commons,
                 container::svector<Index> const& diffs) {
  double ops = 1.0;
  for (auto&& idx : ranges::views::concat(commons, diffs))
    ops *= std::invoke(idxsz, idx);
  // ops == 1.0 implies both commons and diffs empty
  return ops == 1.0 ? 0 : ops;
}

///
/// any element in the vector belongs to the integral range [-1,N)
/// where N is the length of the [Expr] (ie. the iterable of expressions)
///   * only applicable for binary evaluations
///   * the integer -1 can appear in certain places: it implies the binary
///     operation between the last two expressions
///   * eg.
///         * {0,1,-1,2,-1} => ( (e[0], e[1]), e[2])
///         * {0,1,-1,2,3,-1,-1} => ((e[0], e[1]), (e[2],e[3]))
///
using eval_seq_t = container::svector<int>;

///
/// Represents a result of optimization on a range of expressions
/// for a binary evaluation
///
struct OptRes {
  /// Free indices remaining upon evaluation
  container::svector<sequant::Index> indices;

  /// The flops count of evaluation
  double flops;

  /// The evaluation sequence
  eval_seq_t sequence;
};

///
/// returns a vector of Index objects that are common in @c idxs1 and @c idxs2.
///
/// @note I1 and I2 containers are assumed to be sorted by using
/// Index::LabelCompare{};
template <typename I1, typename I2>
container::svector<Index> common_indices(I1 const& idxs1, I2 const& idxs2) {
  using std::back_inserter;
  using std::begin;
  using std::end;
  using std::set_intersection;

  auto result = container::svector<Index>{};

  set_intersection(begin(idxs1), end(idxs1), begin(idxs2), end(idxs2),
                   back_inserter(result), Index::LabelCompare{});
  return result;
}

///
/// returns a vector of Index objects that are common in @c idxs1 and @c idxs2.
///
/// @note I1 and I2 containers are assumed to be sorted by using
/// Index::LabelCompare{};
template <typename I1, typename I2>
container::svector<Index> diff_indices(I1 const& idxs1, I2 const& idxs2) {
  using std::back_inserter;
  using std::begin;
  using std::end;
  using std::set_symmetric_difference;

  auto result = container::svector<Index>{};

  set_symmetric_difference(begin(idxs1), end(idxs1), begin(idxs2), end(idxs2),
                           back_inserter(result), Index::LabelCompare{});
  return result;
}

/// T is integral type
/// TODO: Use C++20 <bit> header when possible
template <typename T>
bool has_single_bit(T x) noexcept {
  return x != 0 && (x & (x - 1)) == 0;
}

///
/// \tparam IdxToSz
/// \param network A TensorNetwork object.
/// \param idxsz An invocable on Index, that maps Index to its dimension.
/// \return Optimal evaluation sequence that minimizes flops. If there are
///         equivalent optimal sequences then the result is the one that keeps
///         the order of tensors in the network as original as possible.
template <typename IdxToSz,
          std::enable_if_t<std::is_invocable_r_v<size_t, IdxToSz, Index>,
                           bool> = true>
eval_seq_t single_term_opt(TensorNetwork const& network, IdxToSz const& idxsz) {
  // number of terms
  auto const nt = network.tensors().size();
  if (nt == 1) return eval_seq_t{0};
  if (nt == 2) return eval_seq_t{0, 1, -1};
  auto nth_tensor_indices = container::svector<container::svector<Index>>{};
  nth_tensor_indices.reserve(nt);

  for (auto i = 0; i < nt; ++i) {
    auto const& tnsr = *network.tensors().at(i);
    auto bk = container::svector<Index>{};
    bk.reserve(bra_rank(tnsr) + ket_rank(tnsr));
    for (auto&& idx : braket(tnsr)) bk.push_back(idx);

    ranges::sort(bk, Index::LabelCompare{});
    nth_tensor_indices.emplace_back(std::move(bk));
  }

  container::svector<OptRes> results((1 << nt), OptRes{{}, 0, {}});

  // power_pos is used, and incremented, only when the
  // result[1<<0]
  // result[1<<1]
  // result[1<<2]
  // and so on are set
  size_t power_pos = 0;
  for (size_t n = 1; n < (1 << nt); ++n) {
    double curr_cost = std::numeric_limits<double>::max();
    std::pair<size_t, size_t> curr_parts{0, 0};
    container::svector<Index> curr_indices{};

    // function to find the optimal partition
    auto scan_parts = [&curr_cost,                              //
                       &curr_parts,                             //
                       &curr_indices,                           //
                           & results = std::as_const(results),  //
                       &idxsz](                                 //
                          size_t lpart, size_t rpart) {
      auto commons =
          common_indices(results[lpart].indices, results[rpart].indices);
      auto diffs = diff_indices(results[lpart].indices, results[rpart].indices);
      auto new_cost = ops_count(idxsz,           //
                                commons, diffs)  //
                      + results[lpart].flops     //
                      + results[rpart].flops;
      if (new_cost <= curr_cost) {
        curr_cost = new_cost;
        curr_parts = decltype(curr_parts){lpart, rpart};
        curr_indices = std::move(diffs);
      }
    };

    biparts(n, scan_parts);

    auto& curr_result = results[n];
    if (has_single_bit(n)) {
      assert(curr_indices.empty());
      // evaluation of a single atomic tensor
      curr_result.flops = 0;
      curr_result.indices = std::move(nth_tensor_indices[power_pos]);
      curr_result.sequence = eval_seq_t{static_cast<int>(power_pos++)};
    } else {
      curr_result.flops = curr_cost;
      curr_result.indices = std::move(curr_indices);
      auto const& first = results[curr_parts.first].sequence;
      auto const& second = results[curr_parts.second].sequence;

      curr_result.sequence =
          (first[0] < second[0] ? ranges::views::concat(first, second)
                                : ranges::views::concat(second, first)) |
          ranges::to<eval_seq_t>;
      curr_result.sequence.push_back(-1);
    }
  }

  return results[(1 << nt) - 1].sequence;
}

}  // namespace

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
/// \param prod  Product to be optimized.
/// \param idxsz An invocable object that maps an Index object to size.
/// \return Parenthesized product expression.
///
/// @note @c prod is assumed to consist of only Tensor expressions
///
template <typename IdxToSz,
          std::enable_if_t<std::is_invocable_v<IdxToSz, Index>, bool> = true>
ExprPtr single_term_opt(Product const& prod, IdxToSz const& idxsz) {
  if (prod.factors().size() < 3) return clone_packed(prod);
  auto seq = single_term_opt(TensorNetwork{prod}, idxsz);
  auto result = container::svector<ExprPtr>{};
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

/////
///// \param expr  Expression to be optimized.
///// \param idxsz An invocable object that maps an Index object to size.
///// \return Optimized expression for lower evaluation cost.
template <typename IdxToSize,
          typename =
              std::enable_if_t<std::is_invocable_r_v<size_t, IdxToSize, Index>>>
ExprPtr optimize(ExprPtr const& expr, IdxToSize&& idx2size) {
  using ranges::views::transform;
  if (expr->is<Tensor>())
    return expr->clone();
  else if (expr->is<Product>())
    return opt::single_term_opt(expr->as<Product>(),
                                std::forward<IdxToSize>(idx2size));
  else if (expr->is<Sum>()) {
    auto smands = *expr | transform([&idx2size](auto&& s) {
      return optimize(s, std::forward<IdxToSize>(idx2size));
    }) | ranges::to_vector;
    return ex<Sum>(Sum{smands.begin(), smands.end()});
  } else
    throw std::runtime_error{"Optimization attempted on unsupported Expr type"};
}

}  // namespace sequant

#endif  // SEQUANT_OPTIMIZE_OPTIMIZE_HPP
