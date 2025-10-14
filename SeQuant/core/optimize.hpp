#ifndef SEQUANT_OPTIMIZE_OPTIMIZE_HPP
#define SEQUANT_OPTIMIZE_OPTIMIZE_HPP

#include <cmath>
#include <cstddef>
#include <functional>
#include <iterator>
#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>

#if __cplusplus >= 202002L
#include <bit>
#endif

namespace {

///
/// \tparam T integral type
/// \return true if @c x has a single bit on in its bit representation.
///
template <typename T>
bool has_single_bit(T x) noexcept {
#if __cplusplus < 202002L
  return x != 0 && (x & (x - 1)) == 0;
#else
  return std::has_single_bit(x);
#endif
}

}  // namespace

namespace sequant {
/// Optimize an expression assuming the number of virtual orbitals
/// greater than the number of occupied orbitals.

class Tensor;

/// \param expr Expression to be optimized.
/// \return EvalNode object.
// EvalNode optimize(ExprPtr const& expr);

namespace opt {

///
/// \param idxsz An invocable that returns size_t for Index argument.
/// \param idxs Index objects.
/// \return flops count
///
template <typename IdxToSz, typename Idxs,
          std::enable_if_t<std::is_invocable_r_v<size_t, IdxToSz, const Index&>,
                           bool> = true>
double ops_count(IdxToSz const& idxsz, Idxs const& idxs) {
  auto oixs = tot_indices(idxs);
  double ops = 1.0;
  for (auto&& idx : ranges::views::concat(oixs.outer, oixs.inner))
    ops *= std::invoke(idxsz, idx);
  // ops == 1.0 implies zero flops.
  return ops == 1.0 ? 0 : ops;
}

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
void biparts(I n, F const& func) {
  if (n == 0) return;
  I const h = static_cast<I>(std::floor(n / 2.0));
  for (I n_ = 1; n_ <= h; ++n_) {
    auto const l = n & n_;
    auto const r = (n - n_) & n;
    if ((l | r) == n) func(l, r);
  }
}

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
  EvalSequence sequence;
};

///
/// Returns a vector of Index objects that are common in @c idxs1 and @c idxs2
/// that are sorted using Index:LabelCompare{}.
///
/// @note I1 and I2 containers are assumed to be sorted by using
/// Index::LabelCompare{};
///
template <typename I1, typename I2, typename Comp = std::less<Index>>
container::svector<Index> common_indices(I1 const& idxs1, I2 const& idxs2) {
  using std::back_inserter;
  using std::begin;
  using std::end;
  using std::set_intersection;

  SEQUANT_ASSERT(std::is_sorted(begin(idxs1), end(idxs1), Comp{}));
  SEQUANT_ASSERT(std::is_sorted(begin(idxs2), end(idxs2), Comp{}));

  container::svector<Index> result;

  set_intersection(begin(idxs1), end(idxs1), begin(idxs2), end(idxs2),
                   back_inserter(result), Comp{});
  return result;
}

///
/// Returns a vector of Index objects that are common in @c idxs1 and @c idxs2
/// that are sorted using Index::LabelCompare{}.
///
/// @note I1 and I2 containers are assumed to be sorted by using
/// Index::LabelCompare{};
///
template <typename I1, typename I2, typename Comp = std::less<Index>>
container::svector<Index> diff_indices(I1 const& idxs1, I2 const& idxs2) {
  using std::back_inserter;
  using std::begin;
  using std::end;
  using std::set_symmetric_difference;

  SEQUANT_ASSERT(std::is_sorted(begin(idxs1), end(idxs1), Comp{}));
  SEQUANT_ASSERT(std::is_sorted(begin(idxs2), end(idxs2), Comp{}));

  container::svector<Index> result;

  set_symmetric_difference(begin(idxs1), end(idxs1), begin(idxs2), end(idxs2),
                           back_inserter(result), Comp{});
  return result;
}

///
/// \tparam IdxToSz
/// \param network A TensorNetwork object.
/// \param idxsz An invocable on Index, that maps Index to its dimension.
/// \return Optimal evaluation sequence that minimizes flops. If there are
///         equivalent optimal sequences then the result is the one that keeps
///         the order of tensors in the network as original as possible.
template <typename IdxToSz,
          std::enable_if_t<std::is_invocable_r_v<size_t, IdxToSz, const Index&>,
                           bool> = true>
EvalSequence single_term_opt(TensorNetwork const& network,
                             IdxToSz const& idxsz) {
  using ranges::views::concat;
  using IndexContainer = container::svector<Index>;
  // number of terms
  auto const nt = network.tensors().size();
  if (nt == 1) return EvalSequence{0};
  if (nt == 2) return EvalSequence{0, 1, -1};
  auto nth_tensor_indices = container::svector<IndexContainer>{};
  nth_tensor_indices.reserve(nt);

  for (std::size_t i = 0; i < nt; ++i) {
    auto const& tnsr = *network.tensors().at(i);

    nth_tensor_indices.emplace_back();
    auto& ixs = nth_tensor_indices.back();
    for (auto&& j : slots(tnsr)) ixs.emplace_back(j);

    ranges::sort(ixs, std::less<Index>{});
  }

  container::svector<OptRes> results((1 << nt), OptRes{{}, 0, {}});

  // power_pos is used, and incremented, only when the
  // result[1<<0]
  // result[1<<1]
  // result[1<<2]
  // and so on are set
  size_t power_pos = 0;
  for (size_t n = 1; n < (1ul << nt); ++n) {
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
      auto new_cost = ops_count(idxsz,                   //
                                concat(commons, diffs))  //
                      + results[lpart].flops             //
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
      SEQUANT_ASSERT(curr_indices.empty());
      // evaluation of a single atomic tensor
      curr_result.flops = 0;
      curr_result.indices = std::move(nth_tensor_indices[power_pos]);
      curr_result.sequence = EvalSequence{static_cast<int>(power_pos++)};
    } else {
      curr_result.flops = curr_cost;
      curr_result.indices = std::move(curr_indices);
      auto const& first = results[curr_parts.first].sequence;
      auto const& second = results[curr_parts.second].sequence;

      curr_result.sequence = (first[0] < second[0] ? concat(first, second)
                                                   : concat(second, first)) |
                             ranges::to<EvalSequence>;
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
  using ranges::views::filter;
  using ranges::views::reverse;

  if (prod.factors().size() < 3)
    return ex<Product>(Product{prod.scalar(), prod.factors().begin(),
                               prod.factors().end(), Product::Flatten::No});
  auto const tensors =
      prod | filter(&ExprPtr::template is<Tensor>) | ranges::to_vector;
  auto seq = single_term_opt(TensorNetwork{tensors}, idxsz);
  auto result = container::svector<ExprPtr>{};
  for (auto i : seq)
    if (i == -1) {
      auto rexpr = *result.rbegin();
      result.pop_back();
      auto lexpr = *result.rbegin();
      result.pop_back();
      auto p = Product{1, ExprPtrList{lexpr, rexpr}, Product::Flatten::No};
      result.push_back(ex<Product>(Product{
          1, p.factors().begin(), p.factors().end(), Product::Flatten::No}));
    } else {
      result.push_back(tensors.at(i));
    }

  auto& p_ = (*result.rbegin()).as<Product>();
  for (auto&& v : prod | reverse | filter(&Expr::template is<Variable>))
    p_.prepend(1, v, Product::Flatten::No);

  p_.scale(prod.scalar());
  return *result.rbegin();
}

///
/// \brief Create clusters out of positions of terms in a sum that share common
///        intermediates.
///
/// \param expr A Sum to find clusters in.
/// \return A vector of clusters (vectors of position index of terms
///         in \c expr).
///
container::vector<container::vector<size_t>> clusters(Sum const& expr);

///
/// \brief Reorder summands so that terms having common intermediates appear
///        closer.
///
/// \param sum Expression to reorder.
/// \return Expression with summands re-ordered.
///
Sum reorder(Sum const& sum);

///
/// \param expr  Expression to be optimized.
/// \param idxsz An invocable object that maps an Index object to size.
/// \param reorder_sum If true, the summands are reordered so that terms with
///                    common sub-expressions appear closer to each other.
/// \return Optimized expression for lower evaluation cost.
template <typename IdxToSize, typename = std::enable_if_t<std::is_invocable_r_v<
                                  size_t, IdxToSize, const Index&>>>
ExprPtr optimize(ExprPtr const& expr, IdxToSize const& idx2size,
                 bool reorder_sum) {
  using ranges::views::transform;
  if (expr->is<Product>()) {
    if (ranges::all_of(*expr, [](auto&& x) {
          return x->template is<Tensor>() || x->template is<Variable>();
        }))
      return opt::single_term_opt(expr->as<Product>(), idx2size);
    else {
      auto const& prod = expr->as<Product>();

      container::svector<ExprPtr> non_tensors(prod.size());
      container::svector<ExprPtr> new_factors;

      for (auto i = 0; i < prod.size(); ++i) {
        auto&& f = prod.factor(i);
        if (f.is<Tensor>() || f.is<Variable>())
          new_factors.emplace_back(f);
        else {
          non_tensors[i] = f;
          auto target_idxs = get_unique_indices(f);
          new_factors.emplace_back(
              ex<Tensor>(L"I_" + std::to_wstring(i), bra(target_idxs.bra),
                         ket(target_idxs.ket), aux(target_idxs.aux)));
        }
      }

      auto result = opt::single_term_opt(
          Product(prod.scalar(), new_factors, Product::Flatten::No), idx2size);

      auto replacer = [&non_tensors](ExprPtr& out) {
        if (!out->is<Tensor>()) return;
        auto const& tnsr = out->as<Tensor>();
        auto&& label = tnsr.label();
        if (label.at(0) == L'I' && label.at(1) == L'_') {
          size_t suffix = std::stoi(std::wstring(label.data() + 2));
          out = non_tensors[suffix].clone();
        }
      };

      result->visit(replacer, /* atoms_only = */ true);
      return result;
    }
  } else if (expr->is<Sum>()) {
    auto smands = *expr | transform([&idx2size](auto&& s) {
      return optimize(s, idx2size, /* reorder_sum= */ false);
    }) | ranges::to_vector;
    auto sum = Sum{smands.begin(), smands.end()};
    return reorder_sum ? ex<Sum>(opt::reorder(sum)) : ex<Sum>(std::move(sum));
  } else
    return expr->clone();
}

}  // namespace opt

///
/// Optimize the expression using IndexSpace::aproximate_size() for reference
/// index extent.
///
/// \param expr  Expression to be optimized.
/// \param reorder_sum If true, the summands are reordered so that terms with
///                    common sub-expressions appear closer to each other.
///                    True by default.
/// \return Optimized expression for lower evaluation cost.
ExprPtr optimize(ExprPtr const& expr, bool reorder_sum = true);

/// Optimize the expression using IndexSpace::aproximate_size() for reference
/// index extent.
///
/// \param expr  Expression to be optimized.
/// \param reorder_sum If true, the summands are reordered so that terms with
///                    common sub-expressions appear closer to each other.
///                    True by default.
/// \return Optimized expression for lower evaluation cost.
ResultExpr& optimize(ResultExpr& expr, bool reorder_sum = true);

/// Optimize the expression using IndexSpace::aproximate_size() for reference
/// index extent.
///
/// \param expr  Expression to be optimized.
/// \param reorder_sum If true, the summands are reordered so that terms with
///                    common sub-expressions appear closer to each other.
///                    True by default.
/// \return Optimized expression for lower evaluation cost.
[[nodiscard]] ResultExpr& optimize(ResultExpr&& expr, bool reorder_sum = true);

}  // namespace sequant

#endif  // SEQUANT_OPTIMIZE_OPTIMIZE_HPP
