#ifndef SEQUANT_CORE_OPTIMIZE_SINGLE_TERM_DETAIL_HPP
#define SEQUANT_CORE_OPTIMIZE_SINGLE_TERM_DETAIL_HPP

// Shared single-term-optimization detail: cost counters, OptRes/EvalSequence,
// init_results, and the subset/bipartition DP scaffolding. Split out of
// single_term.hpp so that cost_model.hpp can include just these helpers
// without the public single_term_opt<Metric> definitions (which themselves
// depend on cost_model.hpp). This breaks the single_term.hpp <-> cost_model.hpp
// include cycle and lets cost_model.hpp be included standalone.

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/optimize/options.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/external/bliss/graph.hh>

#include <range/v3/algorithm/contains.hpp>
#include <range/v3/algorithm/find.hpp>
#include <range/v3/view/concat.hpp>

#include <range/v3/view/filter.hpp>
#include <range/v3/view/indirect.hpp>
#include <range/v3/view/move.hpp>
#include <range/v3/view/reverse.hpp>
#include <range/v3/view/transform.hpp>

#include <algorithm>
#include <bit>
#include <cstdlib>
#include <functional>
#include <limits>
#include <type_traits>

namespace sequant::opt {

template <typename F>
concept has_index_extent = std::is_invocable_r_v<size_t, F, Index const&>;

namespace detail {

/// \brief Whether an "option-like" callable is engaged.
///
/// \c std::function and function pointers report their non-empty state via
/// contextual conversion to bool; any other callable (e.g. a plain lambda) has
/// no empty state and is always engaged. This lets cost helpers accept either a
/// (possibly empty) \c std::function or a bare lambda for optional callables
/// such as \c inner_pow, instead of requiring a bool-convertible type.
template <typename F>
constexpr bool option_engaged(F const& f) {
  if constexpr (std::is_constructible_v<bool, F const&>)
    return static_cast<bool>(f);
  else
    return true;
}

/// \brief Composite-aware extent product of a \ref tot_indices split.
///
/// Multiplies the extents of the outer indices (via \p ixex) by the extents of
/// the inner (CSV/PNO tensor-of-tensor) composites. Shared by all cost counters
/// (\ref flops_counter, \ref memsize_counter, \ref footprint_counter) so every
/// objective sizes CSV/PNO composites identically.
///
/// \param tot_idxs a \ref tot_indices split (\c outer plus \c inner
/// composites). \param ixex     invocable mapping an Index to its extent (sizes
/// \c outer and,
///                 when \p inner_pow is empty, \c inner too).
/// \param inner_pow optional invocable <tt>(composite, k) -> double</tt> giving
///                 the size of a composite that shares its proto-index set with
///                 \c k composites in the same tensor (the k-th power mean of
///                 the per-pair domain). When non-empty, composites are grouped
///                 by proto-index set and each member of a k-group is sized by
///                 \c inner_pow(composite, k); the outer \c nocc^N times the
///                 group product then equals the true block-sparse volume
///                 \c Sum_pairs d^k.
/// \return the extent product (element count).
/// \note When \p inner_pow is empty, composites fall back to \p ixex (i.e.
///       \c k=1), reproducing a single grid-averaged extent per composite,
///       which under-counts multi-composite tensors: for heavy-tailed domains
///       the mean of d^k far exceeds the k-th power of the mean of d.
template <typename Tot, typename Ixex, typename InnerPow>
double inner_aware_volume(Tot const& tot_idxs, Ixex const& ixex,
                          InnerPow const& inner_pow) {
  double mem = ranges::accumulate(tot_idxs.outer, 1., std::multiplies{}, ixex);
  if (option_engaged(inner_pow)) {
    for (auto const& c : tot_idxs.inner) {
      std::size_t k = 0;
      for (auto const& o : tot_idxs.inner)
        if (o.proto_indices() == c.proto_indices()) ++k;
      mem *= inner_pow(c, k);
    }
  } else {
    // No inner_pow, but this tensor HAS composite (CSV/PNO tensor-of-tensor)
    // indices: sizing them by the base extent silently mis-sizes the tensor
    // (each composite counted at its full base-space extent instead of its
    // per-proto domain), which has repeatedly inverted factorization choices
    // (e.g. picking a 4-PAO integral). An empty inner_pow is only valid for a
    // network with NO composites; refuse to guess here.
    if (!ranges::empty(tot_idxs.inner))
      throw std::invalid_argument(
          "inner_aware_volume: composite (CSV/PNO) indices present but no "
          "inner_pow provided -- sizing composites by base extent is a bug. "
          "Pass a real inner_pow (e.g. SizeRegime::inner_pow_fn()).");
    mem = ranges::accumulate(tot_idxs.inner, mem, std::multiplies{}, ixex);
  }
  return mem;
}

/// \brief Cost function returning the flop count of a binary tensor
/// contraction.
///
/// The returned callable computes the product of extents over the union of
/// indices on \c lhs, \c rhs, and \c result. A product of 1 (scalar
/// contraction) is reported as zero flops.
///
/// \param ixex Invocable mapping an Index to its extent.
/// \param inner_pow Optional k-aware CSV/PNO composite extent; see
///        \ref inner_aware_volume. Empty (default) sizes composites by \p ixex.
/// \return A callable <tt>(lhs, rhs, result) -> double</tt> yielding the
/// flop count of the contraction.
template <typename InnerPow = std::function<double(Index const&, std::size_t)>>
auto flops_counter(has_index_extent auto&& ixex, InnerPow inner_pow = {}) {
  return [ixex = std::forward<decltype(ixex)>(ixex),
          inner_pow = std::move(inner_pow)](
             meta::range_of<Index> auto const& lhs,
             meta::range_of<Index> auto const& rhs,
             meta::range_of<Index> auto const& result) -> double {
    using ranges::views::concat;
    // ToT contractibility guard. A contracted index -- shared by lhs & rhs but
    // absent from result -- that carries NO proto-indices is an outer/batch
    // axis in the tensor-of-tensor layout, which a ToT x ToT einsum cannot SUM
    // (it batches a shared outer axis). It is contractible only when at least
    // one operand is flat (proto-free) and thus supplies the index as an
    // ordinary (summable) mode. If BOTH operands are proto-bearing (ToT) and
    // such a contracted non-proto index exists -- the bare C.C overlap over the
    // CSV expansion index mu -- the binary contraction is unevaluable; cost it
    // +inf so the DP routes around it (contract C.t first, de-nesting the
    // amplitude to a flat carrier of mu, then contract mu flat x ToT, exactly
    // how the non-relativistic g.C path handles it). Self-gating: never fires
    // for flat-only networks (no proto) or when either operand is flat (g.C,
    // C.t).
    {
      auto const has_proto = [](Index const& i) {
        return i.has_proto_indices();
      };
      if (ranges::any_of(lhs, has_proto) && ranges::any_of(rhs, has_proto)) {
        for (auto const& k : lhs)
          if (!k.has_proto_indices() && ranges::contains(rhs, k) &&
              !ranges::contains(result, k))
            return std::numeric_limits<double>::max();
      }
      // Proto-value contraction guard (the occ analogue of the C.C case above,
      // and the flat x ToT case the both-ToT test does not reach). A contracted
      // index k -- shared by lhs & rhs, absent from result -- that appears as a
      // PROTO-INDEX of some index in either operand PARAMETRIZES that operand's
      // per-proto inner (PNS/PNO) basis: summing k would mix distinct inner
      // bases. The off-diagonal occ Fock f^i_k t̄^{kj} sums the occ k that is
      // the amplitude's proto-parent, so the pair-(k,j) PNS basis of t̄'s
      // virtuals a<k,j> depends on k; TA's ToT einsum cannot express that. Cost
      // +inf so the DP first de-nests the amplitude against its coefficients
      // (t̄.C.C -> flat carrier of k), after which k is an ordinary mode the
      // flat Fock can contract. NB: this fires for a flat operand too (f is
      // flat), but NOT for f.C / g.C (their contracted mu is a fresh CSV mode,
      // never a proto-value) nor for the intra-pair vv-Fock t̄.f (contracts the
      // shared virtual, an inner carrier -- itself proto-BEARING, hence
      // excluded below).
      {
        auto const proto_value_in = [](auto const& operand, Index const& k) {
          for (auto const& j : operand)
            if (j.has_proto_indices() && ranges::contains(j.proto_indices(), k))
              return true;
          return false;
        };
        for (auto const& k : lhs)
          if (!k.has_proto_indices() && ranges::contains(rhs, k) &&
              !ranges::contains(result, k) &&
              (proto_value_in(lhs, k) || proto_value_in(rhs, k)))
            return std::numeric_limits<double>::max();
      }
    }
    // <IndexSet> is required here: concatenating the three operands repeats
    // every contracted/shared index, so it must be deduplicated before taking
    // the extent product (cf. memsize_counter, which processes each operand
    // separately and so can use the default vector container).
    auto tot_idxs = tot_indices<IndexSet>(concat(lhs, rhs, result));
    double total_flops = inner_aware_volume(tot_idxs, ixex, inner_pow);
    // A product of exactly 1. means the index set was empty (the accumulation
    // init value), i.e. a scalar contraction => zero flops. Extents are
    // integer-valued, so this equality is exact.
    return total_flops == 1. ? 0. : total_flops;
  };
}

/// \brief Cost function returning the total memory footprint of a binary
/// tensor contraction.
///
/// The returned callable sums, over \c lhs, \c rhs, and \c result, the
/// product of extents of each operand's indices. Operands whose extent
/// product is 1 contribute zero.
///
/// \param ixex Invocable mapping an Index to its extent.
/// \param inner_pow Optional k-aware CSV/PNO composite extent; see
///        \ref inner_aware_volume. Empty (default) sizes composites by \p ixex.
/// \return A callable <tt>(lhs, rhs, result) -> double</tt> yielding the
/// summed memory size of the three operands.
template <typename InnerPow = std::function<double(Index const&, std::size_t)>>
auto memsize_counter(has_index_extent auto&& ixex, InnerPow inner_pow = {}) {
  return [ixex = std::forward<decltype(ixex)>(ixex),
          inner_pow = std::move(inner_pow)](
             meta::range_of<Index> auto const& lhs,
             meta::range_of<Index> auto const& rhs,
             meta::range_of<Index> auto const& result) -> double {
    double total_mem{0};
    // Each operand is sized independently, so the default (vector) container of
    // tot_indices suffices -- a single operand's index list has no duplicates,
    // unlike the concatenated set flops_counter must dedup.
    for (auto&& tot_idxs :
         {tot_indices(lhs), tot_indices(rhs), tot_indices(result)}) {
      double mem = inner_aware_volume(tot_idxs, ixex, inner_pow);
      // mem == 1. means this operand had no indices (the accumulation init
      // value), i.e. a scalar; it contributes no memory. Same exact-equality
      // convention as flops_counter above.
      if (mem != 1.) total_mem += mem;
    }
    return total_mem;
  };
}

/// \brief Cost function returning the storage footprint (element count) of a
/// single result tensor.
///
/// Sizes the result's outer indices by \p ixex and its inner (CSV/PNO
/// tensor-of-tensor) composites by \p inner_pow (see \ref inner_aware_volume).
/// Used both as the peak objectives' per-subset footprint (\ref
/// subset_footprints) and as the additive objectives' per-intermediate
/// footprint penalty (see OptimizeOptions::footprint_weight). A scalar (no
/// indices) contributes zero.
///
/// \param ixex Invocable mapping an Index to its extent.
/// \param inner_pow Optional k-aware CSV/PNO composite extent; see
///        \ref inner_aware_volume. Empty (default) sizes composites by \p ixex,
///        which under-counts multi-composite tensors.
/// \return A callable <tt>(result) -> double</tt> yielding the element count of
/// the result tensor.
template <typename InnerPow = std::function<double(Index const&, std::size_t)>>
auto footprint_counter(has_index_extent auto&& ixex, InnerPow inner_pow = {}) {
  return [ixex = std::forward<decltype(ixex)>(ixex),
          inner_pow = std::move(inner_pow)](
             meta::range_of<Index> auto const& result) -> double {
    double mem = inner_aware_volume(tot_indices(result), ixex, inner_pow);
    return mem == 1. ? 0. : mem;
  };
}

///
/// Represents a result of optimization on a range of expressions
/// for a binary evaluation
///
struct OptRes {
  /// Free indices remaining upon evaluation
  IndexSet indices;

  /// The operations count of evaluation
  double ops;

  /// The evaluation sequence
  EvalSequence sequence;

  /// Bitmask splits that resulted into this OptRes
  size_t lp = 0;
  size_t rp = 0;

  /// unique canonical subnets in the optimal tree for this bitmask
  container::vector<size_t> subnets;
};

/// \brief Footprint (dense element count) of every subset's result tensor.
///
/// \c S[n] is the product of extents of the open indices of subset \c n
/// (those remaining after contracting the tensors in \c n, given \c tidxs as
/// the final target indices). \c S[0] (empty subset) and any scalar result are
/// 0. Shared by the peak DP and its tests so both agree on per-subset sizes.
///
/// \param network The TensorNetwork.
/// \param tidxs   Target (open) indices of the network.
/// \param idxsz   Callable mapping an Index to its extent.
/// \param inner_pow Optional k-aware CSV/PNO composite extent forwarded to
///        \ref footprint_counter; see \ref inner_aware_volume.
/// \return \c S[n] = footprint of subset \c n's result tensor.
template <typename TIdxs, typename IdxToSz>
container::vector<double> subset_footprints(
    TensorNetwork const& network, TIdxs const& tidxs, IdxToSz&& idxsz,
    std::function<double(Index const&, std::size_t)> const& inner_pow = {},
    container::vector<char> const* connected = nullptr) {
  container::vector<OptRes> results((size_t{1} << network.tensors().size()));
  init_results(network, tidxs, results, connected);
  auto fp = footprint_counter(std::forward<IdxToSz>(idxsz), inner_pow);
  container::vector<double> S(results.size(), 0.0);
  for (size_t n = 0; n < results.size(); ++n)
    S[n] = (n == 0 || (connected && !(*connected)[n])) ? 0.0
                                                       : fp(results[n].indices);
  return S;
}

/// \brief Collects the distinct batchable indices (in appearance order) across
/// all tensors in \p network.
///
/// Iterates over every tensor slot (bra, ket, and aux) and appends an index to
/// the result the first time it is seen and \p is_batchable returns true for
/// it.  The returned list assigns each index a stable bit position: index at
/// position \c k is bit \c k of a sliced-set bitmask \c B.
///
/// In addition to the top-level slots, two more passes admit pure-occupied
/// indices that \p is_batchable never sees:
///
///  - the pure-occupied protoindices of composite (CSV/PNO/OSV
///    tensor-of-tensor) legs. A composite leg carries its external occupied
///    indices ONLY as protoindices -- they never appear as a top-level
///    bra/ket/aux slot -- so the slot scan alone drops them.
///  - an explicit pure-occupied index that is open (external) on the network
///    root, i.e. a member of \c network.ext_indices(). Such an index is a
///    genuine top-level slot, but \p is_batchable is typically scoped to a
///    non-occupied space (e.g. DF/RI aux), so it would otherwise never be
///    admitted as a batching candidate. Contracted (internal) occupied
///    indices -- those that connect two or more tensors -- are NOT open on
///    the root and so are never admitted by this pass.
///
/// Admitting either lets the batched DP slice that external-occ external
/// mode (mirrored in \ref subset_open_aux). Both passes are guarded by index
/// space (only pure-occupied indices are admitted) so PAO/aux/internal-occ
/// recognition is unchanged.
///
/// \param network  The TensorNetwork to scan.
/// \param is_batchable  Predicate returning true for indices in a batchable
///        space (e.g. a DF/RI auxiliary space).
/// \return Ordered, deduplicated list of batchable indices.
/// \brief Candidate batchable modes: every index (and protoindex) whose space
/// is batchable in EITHER role.
///
/// Batchability is role-based and caller-defined, keeping this layer
/// domain-generic (no index-space kind is named here):
/// - \p is_batchable admits a space batchable when the mode is CONTRACTED
///   (summed at some node);
/// - \p is_batchable_external admits a space batchable when the mode is
///   EXTERNAL (open on the term root -- a spectator carried to the result).
///
/// This returns the UNION of both roles. Each mode's actual role is resolved by
/// \ref PeakBatchedModel::build_context from the root open set, which then
/// drops any mode its role's predicate rejects -- e.g. a mode admitted only as
/// external but appearing contracted is not batchable, which keeps the 2^m
/// search space free of modes that can never be batched in the role they occur
/// in. Protoindices are candidates too: they become plain outer modes in the
/// array view.
inline container::vector<Index> batchable_mode_list(
    TensorNetwork const& network,
    std::function<bool(Index const&)> const& is_batchable,
    std::function<bool(Index const&)> const& is_batchable_external = {}) {
  container::vector<Index> aux;
  // This free function keeps empty-defaulted predicate params for direct
  // single-arg callers, so each is null-guarded below; an all-declining (or
  // all-empty) input yields an empty result naturally -- no early return.
  auto match = [&](Index const& ix) {
    return (is_batchable && is_batchable(ix)) ||
           (is_batchable_external && is_batchable_external(ix));
  };
  for (auto&& t : network.tensors()) {
    auto tp = std::dynamic_pointer_cast<Tensor>(t);
    for (auto&& ix : ranges::views::concat(tp->bra(), tp->ket(), tp->aux())) {
      if (match(ix) && ranges::find(aux, ix) == ranges::end(aux))
        aux.push_back(ix);
      for (auto&& p : ix.proto_indices())
        if (match(p) && ranges::find(aux, p) == ranges::end(aux))
          aux.push_back(p);
    }
  }
  return aux;
}

/// \brief Footprint tables for every sliced-set B of batchable indices.
///
/// Returns a vector of \c 2^m tables, where \c m = \c aux_list.size(). Table
/// \c B is the result of \ref subset_footprints evaluated with an extent
/// function that replaces the full extent of \c aux_list[k] with
/// \c min(full_extent, batch) whenever bit \c k is set in \c B.
///
/// \param network    The TensorNetwork.
/// \param tidxs      Target (open) indices of the network.
/// \param idxsz      Callable mapping an Index to its full extent.
/// \param is_batchable  Predicate identifying batchable indices.
/// \param batch_target_size  Per-index slice-size function (an upper bound): a
///        sliced batchable index \c ix contributes
///        min(full_extent, batch_target_size(aux_list[k])). This is a
///        conservative (over-)estimate of the realized whole-tile batch, which
///        the backend rounds *down* to a tile multiple (never above the
///        target).
/// \param aux_list   Ordered list of distinct batchable indices (as returned
///        by \ref batchable_mode_list).
/// \param inner_pow Optional k-aware CSV/PNO composite extent forwarded to each
///        per-\c B \ref subset_footprints call; see \ref inner_aware_volume.
///        Orthogonal to slicing (composites are not the batchable aux indices).
/// \return \c tables[B][n] = footprint of subset \c n under sliced-set \c B.
template <typename TIdxs, typename IdxToSz>
container::vector<container::vector<double>> sliced_footprints(
    TensorNetwork const& network, TIdxs const& tidxs, IdxToSz&& idxsz,
    std::function<bool(Index const&)> const& is_batchable,
    std::function<std::size_t(Index const&)> const& batch_target_size,
    container::vector<Index> const& aux_list,
    std::function<double(Index const&, std::size_t)> const& inner_pow = {},
    container::vector<char> const* connected = nullptr) {
  // Retained for API compatibility; shrinkability is decided by aux_list
  // membership (which spans all batchability roles), not by this predicate.
  (void)is_batchable;
  std::size_t const m = aux_list.size();
  container::vector<container::vector<double>> tables(std::size_t{1} << m);
  for (std::size_t B = 0; B < tables.size(); ++B) {
    auto extent = [&, B](Index const& ix) -> std::size_t {
      std::size_t e = idxsz(ix);
      // Membership in aux_list IS the authoritative "this is a batchable mode"
      // test: that list spans ALL batchability roles (contracted and external).
      // Gating additionally on the contracted-role predicate silently makes a
      // slice of any mode admitted outside it a NO-OP -- the mode sits in the
      // sliced set B yet keeps its full extent, so the DP sees no benefit and
      // never batches it.
      auto it = ranges::find(aux_list, ix);
      if (it != ranges::end(aux_list)) {
        std::size_t k = static_cast<std::size_t>(it - ranges::begin(aux_list));
        if (B & (std::size_t{1} << k))
          return std::min(e, std::max<std::size_t>(batch_target_size(ix), 1));
      }
      return e;
    };
    tables[B] = subset_footprints(network, tidxs, extent, inner_pow, connected);
  }
  return tables;
}

/// \brief Bitmask of volatile leaf tensors in \p network.
///
/// Bit \c i is set if tensor \c i (in \c network.tensors() order) satisfies
/// \p is_volatile_leaf.  Returns 0 when the predicate is empty (no tensor is
/// volatile, weighting disabled).
///
/// \param network           The TensorNetwork.
/// \param is_volatile_leaf  Predicate identifying volatile leaf tensors; may
///        be empty.
/// \return Bitmask with bit i set iff tensor i is a volatile leaf.
inline std::size_t leaf_volatile_mask(
    TensorNetwork const& network,
    std::function<bool(Tensor const&)> const& is_volatile_leaf) {
  std::size_t mask = 0;
  if (!is_volatile_leaf) return mask;
  std::size_t i = 0;
  for (auto&& t : network.tensors()) {
    auto tp = std::dynamic_pointer_cast<Tensor>(t);
    if (tp && is_volatile_leaf(*tp)) mask |= (std::size_t{1} << i);
    ++i;
  }
  return mask;
}

struct SubNetHash {
  size_t operator()(
      TensorNetwork::SlotCanonicalizationMetadata const& data) const noexcept {
    return data.hash_value();
  }
};

struct SubNetEqual {
  bool operator()(
      TensorNetwork::SlotCanonicalizationMetadata const& left,
      TensorNetwork::SlotCanonicalizationMetadata const& right) const {
    return bliss::ConstGraphCmp::cmp(*left.graph, *right.graph) == 0;
  }
};

/// \brief Seeds the DP table with per-subset open indices and singleton/empty
/// ops counts.
///
/// For each subset \c i of the input tensors, computes the indices that remain
/// open after evaluating that subset and stores them in \c results[i].indices.
/// Initializes \c results[i].ops to 0 for empty and singleton subsets, and
/// to \c max as a sentinel for subsets that will later be filled in by the DP.
/// Singleton subsets also get their one-element \c sequence pre-populated.
template <typename TIdxs>
void init_results(TensorNetwork const& network, TIdxs const& tidxs,
                  container::vector<OptRes>& results,
                  container::vector<char> const* connected = nullptr) {
  using IndexContainer = decltype(OptRes::indices);
  auto tensor_indices = network.tensors()          //
                        | ranges::views::indirect  //
                        | ranges::views::transform(slots);
  auto imed_indices = subset_target_indices(tensor_indices, tidxs, connected);
  SEQUANT_ASSERT(ranges::distance(imed_indices) == ranges::distance(results));
  for (size_t i = 0; i < results.size(); ++i) {
    // Outer-product pruning: a disconnected subset is never an intermediate the
    // DP forms; leave the sentinel. connected[i]==1 for singletons/empty, so
    // leaves are always computed. See outer_product_connectivity.
    if (connected && !(*connected)[i]) {
      results[i].ops = std::numeric_limits<decltype(OptRes::ops)>::max();
      continue;
    }
    results[i].indices =
        imed_indices[i] | ranges::views::move | ranges::to<IndexContainer>;
    results[i].ops = std::popcount(i) > 1
                         ? std::numeric_limits<decltype(OptRes::ops)>::max()
                         : 0;
    if (std::popcount(i) == 1)
      results[i].sequence.emplace_back(std::countr_zero(i));
    // else results[i].sequence is left uninitialized
  }
}

/// \brief Per-tensor adjacency bitmask over "contractible" shared indices.
///
/// `adj[i]` has bit `j` set iff tensors `i` and `j` share at least one
/// top-level (bra/ket/aux) index that is NOT a target index (`tidxs`) -- i.e.
/// an index that is summed somewhere in the term. Protoindices are never
/// compared, so two composites `a<i,j>`, `b<i,j>` that share only occupied
/// protos create no edge. A hyperedge (a contractible index on three-plus
/// tensors) makes all its carriers mutually adjacent. A tensor is never
/// adjacent to itself.
template <typename TIdxs>
inline container::vector<std::size_t> contractible_adjacency(
    TensorNetwork const& network, TIdxs const& tidxs) {
  std::size_t const nt = network.tensors().size();
  container::vector<std::size_t> adj(nt, 0);
  // carrier bitmask per contractible (non-target) top-level index
  container::map<Index, std::size_t, Index::FullLabelCompare> carriers;
  std::size_t i = 0;
  for (auto&& t : network.tensors()) {
    // Iterate the abstract-tensor slots (bra/ket/aux) directly rather than
    // casting to Tensor: TensorNetwork stores AbstractTensorPtr (not
    // necessarily Tensor), and slot bundles may carry null (empty) placeholder
    // indices -- a null carrier would be shared across every tensor with an
    // empty slot and create a spurious edge that needlessly disables pruning.
    // Mirrors init_results, which also enumerates indices via slots().
    for (auto&& ix : slots(*t)) {
      if (!ix.nonnull()) continue;                // skip empty slots
      if (ranges::contains(tidxs, ix)) continue;  // target/output: not summed
      carriers[ix] |= (std::size_t{1} << i);
    }
    ++i;
  }
  for (auto&& [ix, cm] : carriers) {
    if (std::popcount(cm) < 2) continue;  // appears on < 2 tensors
    std::size_t rest = cm;
    while (rest) {
      std::size_t const b = static_cast<std::size_t>(std::countr_zero(rest));
      adj[b] |= (cm & ~(std::size_t{1} << b));  // neighbors, excluding self
      rest &= rest - 1;
    }
  }
  return adj;
}

/// \brief `connected[n] == 1` iff the subgraph induced on the set bits of
/// subset mask `n` is connected under `adjacency` (from \ref
/// contractible_adjacency). Empty and singleton subsets are connected by
/// definition. Size `1 << nt`.
inline container::vector<char> connected_subsets(
    container::vector<std::size_t> const& adjacency, std::size_t nt) {
  container::vector<char> connected(std::size_t{1} << nt, 0);
  for (std::size_t n = 0; n < connected.size(); ++n) {
    if (std::popcount(n) <= 1) {
      connected[n] = 1;
      continue;
    }
    std::size_t const start = n & (~n + 1);  // lowest set bit
    std::size_t reached = start, frontier = start;
    while (frontier) {
      std::size_t next = 0, f = frontier;
      while (f) {
        next |= adjacency[static_cast<std::size_t>(std::countr_zero(f))];
        f &= f - 1;
      }
      next &= n & ~reached;
      reached |= next;
      frontier = next;
    }
    connected[n] = (reached == n) ? 1 : 0;
  }
  return connected;
}

/// \brief Outer-product pruning mask for a term. Returns \ref
/// connected_subsets over \ref contractible_adjacency, EXCEPT it returns an
/// all-connected mask when pruning is disabled (\p prune == false, or the env
/// `SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING` force-override) or when the full
/// network is itself disconnected (a genuine product term -- never a CC
/// residual summand, by the linked-cluster theorem). Callers then apply one
/// uniform `if (!mask[n]) skip` test and get unpruned behavior in those cases.
/// \param prune Whether to prune (OptimizeOptions::prune_outer_products). When
///        false, the returned mask is all-connected (unpruned DP).
template <typename TIdxs>
inline container::vector<char> outer_product_connectivity(
    TensorNetwork const& network, TIdxs const& tidxs, bool prune = true) {
  std::size_t const nt = network.tensors().size();
  std::size_t const sz = std::size_t{1} << nt;
  if (!prune || std::getenv("SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING"))
    return container::vector<char>(sz, 1);
  auto mask = connected_subsets(contractible_adjacency(network, tidxs), nt);
  if (!mask.empty() && !mask.back())  // full network disconnected: disable
    return container::vector<char>(sz, 1);
  return mask;
}

struct SubnetMetadata {
  /// meta_ids[n] is the canonical-subnet id of subset n, or
  /// numeric_limits<size_t>::max() for subsets with popcount < 2.
  container::vector<size_t> meta_ids;
  /// Cost of evaluating one representative of each canonical subnet id,
  /// indexed by id. Populated lazily during the DP.
  container::vector<double> unique_meta_costs;
};

/// \brief Precomputes canonical-subnet identifiers for every subset of size
/// >= 2 so that structurally equivalent subnetworks share a CSE id.
///
/// Builds a `TensorNetwork` for each subset, canonicalizes it, and assigns a
/// dense integer id to each distinct canonical form. The returned
/// `unique_meta_costs` is sized to the number of distinct ids and zero-filled;
/// it is populated during the DP as each canonical subnet's optimal cost
/// becomes known.
///
/// Side effect: `results[n].indices` may be reordered by `canonicalize_slots`.
inline SubnetMetadata build_subnet_metadata(
    TensorNetwork const& network, container::vector<OptRes>& results,
    container::vector<char> const& connected) {
  SubnetMetadata out;
  // Use max as sentinel for entries with popcount < 2 (singletons/empty),
  // which are skipped below and never assigned a real meta ID.
  out.meta_ids.resize(results.size(), std::numeric_limits<size_t>::max());
  container::unordered_map<TensorNetwork::SlotCanonicalizationMetadata, size_t,
                           SubNetHash, SubNetEqual>
      meta_to_id;

  for (size_t n = 0; n < results.size(); ++n) {
    if (std::popcount(n) < 2) continue;
    if (!connected[n]) continue;  // outer-product subset, never an intermediate
    auto ts = bits::on_bits_index(n) | bits::sieve(network.tensors());
    container::vector<ExprPtr> ts_expr;
    for (auto&& t : ts)
      ts_expr.emplace_back(std::dynamic_pointer_cast<Tensor>(t)->clone());

    auto tn = TensorNetwork{ts_expr};
    auto meta = tn.canonicalize_slots(
        TensorCanonicalizer::cardinal_tensor_labels(), &results[n].indices);

    auto [it, inserted] = meta_to_id.try_emplace(std::move(meta), 0);
    if (inserted) it->second = meta_to_id.size() - 1;

    out.meta_ids[n] = it->second;
  }
  out.unique_meta_costs.resize(meta_to_id.size(), 0.0);
  return out;
}

/// \brief Per-subset bitmask of batchable indices that are OPEN in that subset.
///
/// For each subset \c n of the input tensors, bit \c k of \c open_modes[n] is
/// set iff \c aux_list[k] is among the open (external) indices of subset \c n
/// (those that remain after contracting \c n's tensors, with \c tidxs as the
/// final targets). Used by the multi-mode batched DP and oracle to restrict the
/// sliced-set context to indices actually open in a sized subset, so that table
/// lookups (\ref sliced_footprints) are indexed consistently on both sides.
///
/// \param network   The TensorNetwork.
/// \param tidxs     Target (open) indices of the network.
/// \param aux_list  Ordered list of distinct batchable indices (as returned by
///        \ref batchable_mode_list); index \c k maps to bit \c k.
/// \return \c open_modes[n] for every subset \c n.
template <typename TIdxs>
container::vector<std::size_t> subset_open_aux(
    TensorNetwork const& network, TIdxs const& tidxs,
    container::vector<Index> const& aux_list) {
  container::vector<OptRes> results(
      (std::size_t{1} << network.tensors().size()));
  // NOT pruned: is_external_mode consumes open_modes over the full subset
  // lattice (including disconnected subsets), so every entry must be real.
  init_results(network, tidxs, results);
  // A batchable mode may be open either DIRECTLY (a top-level open index) or as
  // a PROTOINDEX of an open index: protoindices become plain outer modes in the
  // array view, so a mode carried as a proto of an open index is open too. Both
  // are checked structurally -- no index-space kind is consulted, keeping this
  // layer domain-generic.
  container::vector<std::size_t> open_modes(results.size(), 0);
  for (std::size_t n = 0; n < results.size(); ++n) {
    for (std::size_t k = 0; k < aux_list.size(); ++k) {
      Index const& ax = aux_list[k];
      bool open = ranges::find(results[n].indices, ax) !=
                  ranges::end(results[n].indices);
      if (!open) {
        for (auto const& ix : results[n].indices) {
          auto const& pr = ix.proto_indices();
          if (ranges::find(pr, ax) != ranges::end(pr)) {
            open = true;
            break;
          }
        }
      }
      if (open) open_modes[n] |= (std::size_t{1} << k);
    }
  }
  return open_modes;
}

/// Per-(subset, sliced-set) state for the multi-mode batched peak DP.
/// Indexed in the flat table by \c n*(2^m)+B (subset \c n, sliced-set \c B).
struct BatchedRes {
  double peak = std::numeric_limits<double>::max();  // min peak for (n, B)
  std::size_t lp = 0, rp = 0;  // winning bipartition (0 for singletons)
  bool lp_first = true;        // winning evaluation order
  std::size_t aprime = 0;  // winning subset of batchable indices sliced here
  double flops = std::numeric_limits<double>::max();  // tie-break (min work
                                                      // among equal-peak)
};

}  // namespace detail
}  // namespace sequant::opt

#endif  // SEQUANT_CORE_OPTIMIZE_SINGLE_TERM_DETAIL_HPP
