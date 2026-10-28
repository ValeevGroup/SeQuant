#ifndef SEQUANT_EVAL_BACKENDS_TILEDARRAY_RESULT_HPP
#define SEQUANT_EVAL_BACKENDS_TILEDARRAY_RESULT_HPP

#ifdef SEQUANT_HAS_TILEDARRAY

#include <sstream>

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/math.hpp>
#include <SeQuant/core/utility/exception.hpp>

#include <TiledArray/einsum/tiledarray.h>
#include <tiledarray.h>

#include <range/v3/view/iota.hpp>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>

namespace sequant {

// implementation details of the TiledArray result backend; prefer
// sequant::detail over an unnamed namespace in a header (see CppCoreGuidelines
// SF.21 / "Use unnamed namespaces in headers ... no" guidance)
namespace detail {

/// @return true if SEQUANT_EVAL_WARN_TOT_SUM_INNER is set in the environment
/// (diagnostic for nested additions whose inner index order disagrees)
inline bool warn_tot_sum_inner() {
  static const bool on =
      std::getenv("SEQUANT_EVAL_WARN_TOT_SUM_INNER") != nullptr;
  return on;
}

// Instrumentation (SEQUANT_SYNC_STATS, analysis-only): count gop.fence() calls
// so the ordered-executor's synchronization overhead can be localized.
inline std::atomic<std::size_t>& fence_counter() {
  static std::atomic<std::size_t> c{0};
  return c;
}
inline void note_fence() {
  fence_counter().fetch_add(1, std::memory_order_relaxed);
}
inline std::atomic<std::size_t>& wait_counter() {
  static std::atomic<std::size_t> c{0};
  return c;
}
inline std::atomic<std::size_t>& slice_counter() {
  static std::atomic<std::size_t> c{0};
  return c;
}
inline void note_slice() {
  slice_counter().fetch_add(1, std::memory_order_relaxed);
}
// Drain policy: every backend op used to end with
// DistArray::wait_for_lazy_cleanup(world), which blocks until EVERY array
// scheduled for lazy deletion so far has finished its pending tasks. On one
// rank that is a cheap local wait; on many ranks it is a pipeline drain after
// each of the (tens of) thousands of ops an iteration issues, so no op's tail
// (remote sends, straggler tiles) overlaps the next op's start. The consumer
// chooses: `true` (default) keeps the drain, `false` lets ops pipeline and
// leaves the cleanup to run alongside the next op (temporaries are released a
// little later, so the peak can grow by the in-flight set).
inline std::atomic<bool>& drain_after_op_flag() {
  static std::atomic<bool> f{true};
  return f;
}
template <typename ArrayT>
inline void drain_lazy_cleanup(ArrayT const& arr) {
  if (drain_after_op_flag().load(std::memory_order_relaxed)) {
    ArrayT::wait_for_lazy_cleanup(arr.world());
    wait_counter().fetch_add(1, std::memory_order_relaxed);
  }
}
inline std::atomic<long long>& slice_ns() {
  static std::atomic<long long> c{0};
  return c;
}
struct FenceReporter {
  ~FenceReporter() {
    if (std::getenv("SEQUANT_SYNC_STATS"))
      std::cerr << "TOTAL gop.fence() calls = " << fence_counter().load()
                << " ; wait_for_lazy_cleanup calls = " << wait_counter().load()
                << " ; slice_mode calls = " << slice_counter().load()
                << " ; slice_mode total = " << (slice_ns().load() / 1e9) << " s"
                << "\n";
  }
};
inline FenceReporter fence_reporter_{};

}  // namespace detail

/// Whether every TA-backend op ends with DistArray::wait_for_lazy_cleanup
/// (see detail::drain_lazy_cleanup). Default true (the historical behavior).
inline void set_ta_drain_after_op(bool on) {
  detail::drain_after_op_flag().store(on, std::memory_order_relaxed);
}
inline bool ta_drain_after_op() {
  return detail::drain_after_op_flag().load(std::memory_order_relaxed);
}

namespace detail {
// Wire PhaseTimer's boundary barrier to a TA world fence so per-region phase
// timers cannot misattribute async (deferred) work across regions. Runs only
// under SEQUANT_UT_PHASE (PhaseTimer::Scope calls barrier() only when enabled).
inline const bool phase_fence_installed_ = [] {
  ::sequant::eval::PhaseTimer::fence_hook() = [] {
    TA::get_default_world().gop.fence();
    ::sequant::detail::note_fence();
  };
  return true;
}();

/// Inner-tensor mode count of a tensor-of-tensor DistArray (0 for a regular,
/// non-nested array). The outer trange carries no inner information, so the
/// inner rank is read from the first local non-empty inner tile and reduced
/// (max) across the world -- this makes it well-defined on ranks holding no
/// local data (or only empty inner tiles), as long as some rank holds a
/// non-empty inner tile. Needed to build a ToT-valid annotation ("outer;inner")
/// in the annotation-free array operations (add_inplace,
/// slice_array_over_mode), since a flat annotation trips DistArray's
/// is_tot_index() check. Being a collective, it is not for per-operation use:
/// the backend's results record their inner rank once and serve it locally
/// (ResultTensorOfTensorTA::inner_rank()).
template <typename ArrayT>
[[nodiscard]] std::size_t tot_inner_rank(ArrayT const& arr) {
  std::size_t r = 0;
  if constexpr (TA::detail::is_tensor_of_tensor_v<
                    typename ArrayT::value_type>) {
    // Inner tensors carry the rank; an outer tile may hold null/empty inner
    // views (e.g. a screened pair), so skip those -- calling range() on an
    // empty view asserts (ArenaTensor) or is UB (TA::Tensor).
    for (auto it = arr.begin(); it != arr.end() && r == 0; ++it) {
      // NOTE: own the tile by value. it->get() (TileReference::get) returns the
      // tile by value; binding it to a reference does not lifetime-extend the
      // temporary through this proxy chain, so a reference would dangle and the
      // inner loop below would read freed stack (ASan: stack-use-after-scope).
      // The copy is shallow over ArenaTensor's Cell pointers, which alias the
      // arena owned by the DistArray-resident tile (valid for arr's lifetime).
      auto const outer_tile = it->get();
      for (auto const& inner : outer_tile) {
        if (inner.empty()) continue;
        if (inner.range().rank() > 0) {
          r = inner.range().rank();
          break;
        }
      }
    }
    arr.world().gop.max(r);
  }
  return r;
}

///
/// \brief Particle-symmetrize a TA::DistArray (tensor-of-scalar or
///        tensor-of-tensor).
///
/// \param arr The array to be symmetrized.
///
/// \pre ToS: rank is even. ToT: outer rank == inner rank.
///
/// \return The symmetrized TA::DistArray.
///
template <typename... Args>
auto column_symmetrize_ta(TA::DistArray<Args...> const& arr) {
  using ranges::views::iota;
  constexpr bool is_tot = TA::detail::is_tensor_of_tensor_v<
      typename TA::DistArray<Args...>::value_type>;

  size_t const outer_rank = arr.trange().rank();

  if constexpr (is_tot) {
    size_t const inner_rank = tot_inner_rank(arr);
    // All-empty-inner ToT (e.g. every CSV pair screened out): no populated tile
    // to read a rank from, so tot_inner_rank() is 0. It represents zero, and
    // symmetrizing zero is zero -- return unchanged instead of failing the
    // equal-rank check below.
    if (inner_rank == 0) return arr;
    if (outer_rank != inner_rank)
      throw Exception(
          "ToT symmetrization requires equal outer and inner rank (outer=" +
          std::to_string(outer_rank) + ", inner=" + std::to_string(inner_rank) +
          ")");
  }

  // ToT (equal-rank, validated above): total rank = outer + inner = 2*outer.
  size_t const total_rank = is_tot ? 2 * outer_rank : outer_rank;

  if (total_rank % 2 != 0)
    throw Exception("This function only supports even-ranked tensors");

  size_t const nparticles = total_rank / 2;

  perm_t perm = iota(size_t{0}, total_rank) | ranges::to<perm_t>;

  auto make_annot = [nparticles](perm_t const& p) {
    if constexpr (is_tot) {
      perm_t const outer(p.begin(), p.begin() + nparticles);
      perm_t const inner(p.begin() + nparticles, p.end());
      return ords_to_annot(outer) + ";" + ords_to_annot(inner);
    } else {
      return ords_to_annot(p);
    }
  };

  auto const lannot = make_annot(perm);

  TA::DistArray<Args...> result;
  auto call_back = [&result, &lannot, &arr, &perm = std::as_const(perm),
                    &make_annot]() {
    auto const rannot = make_annot(perm);
    if (result.is_initialized()) {
      result(lannot) += arr(rannot);
    } else {
      result(lannot) = arr(rannot);
    }
  };

  symmetric_permutation(SymmetricParticleRange{perm.begin(),               //
                                               perm.begin() + nparticles,  //
                                               nparticles},
                        call_back);

  auto const nf = static_cast<double>(rational{1, factorial(nparticles)});
  result(lannot) = nf * result(lannot);

  ::sequant::detail::drain_lazy_cleanup(result);

  return result;
}

///
/// \brief This function implements the antisymmetrization of TA::DistArray.
///
/// \param arr The array to be antisymmetrized.
///
/// \param bra_rank The rank of the bra indices
///
/// \return The antisymmetrized TA::DistArray.
///
template <typename... Args>
auto particle_antisymmetrize_ta(TA::DistArray<Args...> const& arr,
                                size_t bra_rank) {
  using ranges::views::iota;
  size_t const rank = arr.trange().rank();
  SEQUANT_ASSERT(bra_rank <= rank);
  size_t const ket_rank = rank - bra_rank;

  if (bra_rank <= 1 && ket_rank <= 1) {
    // nothing to do
    return arr;
  }

  perm_t perm = iota(size_t{0}, rank) | ranges::to<perm_t>;
  perm_t bra_perm = iota(size_t{0}, bra_rank) | ranges::to<perm_t>;
  perm_t ket_perm = iota(bra_rank, rank) | ranges::to<perm_t>;

  const auto lannot = ords_to_annot(perm);

  auto process_permutations = [&lannot](const TA::DistArray<Args...>& input_arr,
                                        size_t range_rank, perm_t range_perm,
                                        const std::string& other_annot,
                                        bool is_bra) -> TA::DistArray<Args...> {
    if (range_rank <= 1) return input_arr;
    TA::DistArray<Args...> result;

    auto callback = [&](int parity) {
      const auto range_annot = ords_to_annot(range_perm);
      const auto annot = other_annot.empty()
                             ? range_annot
                             : (is_bra ? range_annot + "," + other_annot
                                       : other_annot + "," + range_annot);

      typename decltype(result)::numeric_type p_ = parity == 0 ? 1 : -1;
      if (result.is_initialized()) {
        result(lannot) += p_ * input_arr(annot);
      } else {
        result(lannot) = p_ * input_arr(annot);
      }
    };
    antisymmetric_permutation(ParticleRange{range_perm.begin(), range_rank},
                              callback);
    return result;
  };

  // Process bra permutations first
  const auto ket_annot = ket_rank == 0 ? "" : ords_to_annot(ket_perm);
  auto result = process_permutations(arr, bra_rank, bra_perm, ket_annot, true);

  // Process ket permutations
  const auto bra_annot = bra_rank == 0 ? "" : ords_to_annot(bra_perm);
  result = process_permutations(result, ket_rank, ket_perm, bra_annot, false);

  auto const nf = static_cast<double>(
      rational{1, factorial(bra_rank) * factorial(ket_rank)});
  result(lannot) = nf * result(lannot);

  ::sequant::detail::drain_lazy_cleanup(result);
  return result;
}

template <typename... Args>
inline void log_ta(Args const&... args) noexcept {
  log_result("[TA] ", args...);
}

/// Batch-step data-movement record, emitted at eval log level > 0 (the same
/// gate as the executor's `Eval | ...` records):
///
///   Batch | <Slice|Scatter|Accumulate> | <time>ns | bytes=<B> | <annot>
///
/// Slice: an operand block gathered out of a mode for one batch step;
/// Scatter: a per-batch block written into its pre-sized destination;
/// Accumulate: a per-batch block added into its accumulator. `bytes` is the
/// host size of the moved block, summed over ranks (collective, so every rank
/// must reach this call whenever the level is nonzero -- the level is
/// identical across ranks). These movements are not eval-tree ops, so the
/// `Eval` records never see them; this is the only accounting of their cost.
template <typename... Args>
inline void log_batch_op(char const* kind, std::chrono::nanoseconds elapsed,
                         TA::DistArray<Args...> const& moved,
                         std::string const& annot) noexcept {
  // Gated at run time by the trace level (like the per-op `Eval |` lines a
  // Trace::On evaluation emits), not by SEQUANT_EVAL_TRACE: the batched
  // executor's data movement (Slice / Scatter / Accumulate / Clone) is
  // exactly what a wet trace on a release build must account for.
  auto& l = Logger::instance();
  if (l.eval.level == 0) return;
  // this rank's share of the moved array: no collective here, so a batch op
  // a rank runs on its own never waits on the others (and the line costs a
  // tile walk, not a reduction)
  auto const bytes = TA::size_of<TA::MemorySpace::Host>(moved);
  write_log(l, "Batch | ", kind, " | ", elapsed.count(),
            "ns | local_bytes=", bytes, "B | ", annot, '\n');
}

/// Convert sequant::DeNest to TA::DeNest
inline constexpr TA::DeNest to_ta_denest(DeNest d) noexcept {
  return d == DeNest::True ? TA::DeNest::True : TA::DeNest::False;
}

/// Split a TA tensor-of-tensor (or flat) annotation string into its OUTER
/// labels. The annotation format is `outer0,outer1,...[;inner0,inner1,...]`;
/// only the comma-separated outer labels (the part before the optional `;`)
/// are returned, in order. Used to map a result's outer modes back onto the
/// operands' TiledRange1's when building a result outer TiledRange for an
/// imposed SparseShape.
[[nodiscard]] inline container::svector<std::string> outer_annot_labels(
    std::string const& annot) {
  container::svector<std::string> labels;
  auto const semi = annot.find(';');
  std::string const outer =
      semi == std::string::npos ? annot : annot.substr(0, semi);
  std::size_t pos = 0;
  while (pos <= outer.size()) {
    auto const comma = outer.find(',', pos);
    auto const end = comma == std::string::npos ? outer.size() : comma;
    if (end > pos) labels.emplace_back(outer.substr(pos, end - pos));
    if (comma == std::string::npos) break;
    pos = comma + 1;
  }
  return labels;
}

/// Split a TA ToT annotation string into its INNER labels: the comma-separated
/// labels AFTER the optional `;` (empty if the annotation has no inner part).
[[nodiscard]] inline container::svector<std::string> inner_annot_labels(
    std::string const& annot) {
  container::svector<std::string> labels;
  auto const semi = annot.find(';');
  if (semi == std::string::npos) return labels;
  std::string const inner = annot.substr(semi + 1);
  std::size_t pos = 0;
  while (pos <= inner.size()) {
    auto const comma = inner.find(',', pos);
    auto const end = comma == std::string::npos ? inner.size() : comma;
    if (end > pos) labels.emplace_back(inner.substr(pos, end - pos));
    if (comma == std::string::npos) break;
    pos = comma + 1;
  }
  return labels;
}

/// Does emitting a ToT*ToT general product `result(this) = (l(lan) * r(ran))`
/// require a NON-IDENTITY inner result permutation, which TiledArray's
/// cont_engine does not yet support (it throws "a non-identity inner result
/// permutation is not yet supported")?  The contraction produces the result's
/// inner modes in the natural order (left's non-contracted inner labels in left
/// order, then right's); if \p this_annot's inner labels are in a different
/// order, TA cannot emit it.  Detecting this lets the shaped-product path
/// decline gracefully (the caller falls through to the unshaped einsum prod(),
/// which DOES handle the reorder) instead of throwing.
[[nodiscard]] inline container::svector<std::string> tot_product_natural_inner(
    std::string const& lannot, std::string const& rannot,
    std::string const& this_annot) {
  auto const li = inner_annot_labels(lannot);
  auto const ri = inner_annot_labels(rannot);
  auto const ti = inner_annot_labels(this_annot);
  auto contains = [](container::svector<std::string> const& v,
                     std::string const& s) {
    return std::find(v.begin(), v.end(), s) != v.end();
  };
  // Natural order = left's inner labels that survive (appear in the result),
  // then right's inner labels that survive and are not already taken from left.
  container::svector<std::string> natural;
  for (auto const& s : li)
    if (contains(ti, s)) natural.push_back(s);
  for (auto const& s : ri)
    if (contains(ti, s) && !contains(natural, s)) natural.push_back(s);
  return natural;
}

[[nodiscard]] inline bool tot_product_needs_inner_reorder(
    std::string const& lannot, std::string const& rannot,
    std::string const& this_annot) {
  auto const ti = inner_annot_labels(this_annot);
  if (ti.empty()) return false;  // denest / no inner result modes
  return tot_product_natural_inner(lannot, rannot, this_annot) != ti;
}

/// \p this_annot with its inner labels in the natural order of the product
/// (see tot_product_needs_inner_reorder): the layout the contraction engine
/// can emit directly
[[nodiscard]] inline std::string tot_product_natural_annot(
    std::string const& lannot, std::string const& rannot,
    std::string const& this_annot) {
  auto const semi = this_annot.find(';');
  std::string out =
      semi == std::string::npos ? this_annot : this_annot.substr(0, semi);
  auto const natural = tot_product_natural_inner(lannot, rannot, this_annot);
  if (natural.empty()) return out;
  out += ';';
  for (std::size_t k = 0; k < natural.size(); ++k) {
    if (k) out += ',';
    out += natural[k];
  }
  return out;
}

/// Build the result's OUTER TiledRange for a binary product, given the two
/// operand arrays and the three (left, right, result) annotations. Each result
/// outer label is matched, in order, against the left operand's outer labels
/// (then the right's), and the corresponding TiledRange1 is taken from that
/// operand. A result outer label must appear on at least one operand (true for
/// any well-formed contraction). This is the trange handed to the result-shape
/// provider so it can assemble a SparseShape over the result's outer tiles.
[[nodiscard]] inline TA::TiledRange result_outer_trange(
    TA::TiledRange const& ltr, std::string const& lannot,
    TA::TiledRange const& rtr, std::string const& rannot,
    std::string const& cannot) {
  auto const lo = outer_annot_labels(lannot);
  auto const ro = outer_annot_labels(rannot);
  auto const co = outer_annot_labels(cannot);

  std::vector<TA::TiledRange1> dims;
  dims.reserve(co.size());
  for (auto const& lbl : co) {
    if (auto const i = find_position(lo, lbl)) {
      dims.emplace_back(ltr.dim(*i));
    } else if (auto const j = find_position(ro, lbl)) {
      dims.emplace_back(rtr.dim(*j));
    } else {
      throw Exception("result_outer_trange: result outer label '" + lbl +
                      "' is carried by neither operand of the binary product");
    }
  }
  return TA::TiledRange(dims.begin(), dims.end());
}

/// Map a contiguous element range `[elem_lo, elem_hi)` on a mode's TiledRange1
/// to the tile range `[tile_lo, tile_hi)` it must coincide with. A tiled
/// backend can only cut or scatter whole tiles, so the element bounds must be
/// in-range and fall on tile boundaries; this asserts both (mode_batches()
/// yields exactly such tile-aligned ranges). Shared by slice_mode() (gather a
/// block out) and write_into_slice() (scatter a block in) so both agree on the
/// element-to-tile contract and its alignment preconditions.
[[nodiscard]] inline std::pair<std::size_t, std::size_t> slice_bounds_to_tiles(
    TA::TiledRange1 const& tr1, std::size_t elem_lo, std::size_t elem_hi) {
  SEQUANT_ASSERT(elem_lo >= tr1.elements_range().first && elem_lo < elem_hi &&
                 elem_hi <= tr1.elements_range().second);
  std::size_t const tile_lo = tr1.element_to_tile(elem_lo);
  SEQUANT_ASSERT(tr1.tile(tile_lo).first == elem_lo);  // lo on a tile boundary
  std::size_t const tile_hi = (elem_hi >= tr1.elements_range().second)
                                  ? tr1.tile_extent()
                                  : tr1.element_to_tile(elem_hi);
  SEQUANT_ASSERT(elem_hi >= tr1.elements_range().second ||
                 tr1.tile(tile_hi).first == elem_hi);  // hi on a tile boundary
  return {tile_lo, tile_hi};
}

/// \return the TiledRange of \p arr, as printed by TiledArray
template <typename ArrayT>
[[nodiscard]] std::string ta_trange_string(ArrayT const& arr) {
  std::ostringstream os;
  os << arr.trange();
  return os.str();
}

/// \return the Frobenius norm of \p arr, or -1 if its tile type does not
///         support TA::squared_norm
template <typename ArrayT>
[[nodiscard]] double ta_norm2(ArrayT const& arr) {
  if constexpr (requires(ArrayT const& a) {
                  { TA::squared_norm(a) } -> std::convertible_to<double>;
                })
    return std::sqrt(static_cast<double>(TA::squared_norm(arr)));
  else
    return -1.0;
}

/// diagnostic tile-by-tile comparison of \p a and \p b: counts of tiles
/// nonzero in both, only in one, or in neither, and the norms of the
/// one-sided tiles and of the difference on common tiles (local tiles only)
template <typename ArrayT>
[[nodiscard]] std::string ta_tile_diff(ArrayT const& a, ArrayT const& b) {
  if (a.trange() != b.trange()) return "(trange mismatch)";
  std::size_t only_a = 0, only_b = 0, both = 0, none = 0;
  double n_only_a = 0.0, n_only_b = 0.0, n_both_diff = 0.0;
  for (auto const& ord : a.tiles_range()) {
    bool const za = a.is_zero(ord), zb = b.is_zero(ord);
    if (za && zb) {
      ++none;
      continue;
    }
    if (!za && zb) {
      ++only_a;
      if (a.is_local(ord)) {
        auto const t = a.find_local(ord).get();
        if constexpr (requires { TA::norm(t); })
          n_only_a += std::pow(static_cast<double>(TA::norm(t)), 2);
      }
      continue;
    }
    if (za && !zb) {
      ++only_b;
      if (b.is_local(ord)) {
        auto const t = b.find_local(ord).get();
        if constexpr (requires { TA::norm(t); })
          n_only_b += std::pow(static_cast<double>(TA::norm(t)), 2);
      }
      continue;
    }
    ++both;
    if (a.is_local(ord) && b.is_local(ord)) {
      auto const ta = a.find_local(ord).get();
      auto const tb = b.find_local(ord).get();
      if constexpr (requires { TA::norm(ta.subt(tb)); })
        n_both_diff += std::pow(static_cast<double>(TA::norm(ta.subt(tb))), 2);
    }
  }
  std::ostringstream os;
  os << "tiles: both=" << both << " only_this=" << only_a
     << " only_other=" << only_b << " none=" << none
     << " |only_this|=" << std::sqrt(n_only_a)
     << " |only_other|=" << std::sqrt(n_only_b)
     << " |diff on both|=" << std::sqrt(n_both_diff);
  return os.str();
}

/// \return the TiledRange \p tr, as printed by TiledArray
[[nodiscard]] inline std::string ta_trange_string(TA::TiledRange const& tr) {
  std::ostringstream os;
  os << tr;
  return os.str();
}

}  // namespace detail

/// TA::Tensor memory use logger
/// If TiledArray was configured with TA_TENSOR_MEM_PROFILE set this
/// prints the current use of memory by TA::Tensor objects in host memory space
/// to \p os .
/// \param world the world object to use for logging
/// \param label string to prepend to the profile
void log_ta_tensor_host_memory_use(madness::World& world,
                                   std::string_view label = "");

inline void log_ta_tensor_host_memory_use() {
#if defined(SEQUANT_EVAL_TRACE)
  log_ta_tensor_host_memory_use(TA::get_default_world(), "[TA]");
#endif
}

// defined below; declared here so the result classes' slice_mode() overrides
// can call it.
template <typename... Args>
[[nodiscard]] TA::DistArray<Args...> slice_array_over_mode(
    TA::DistArray<Args...> const& arr, std::size_t mode, std::size_t tile_lo,
    std::size_t tile_hi, std::optional<std::size_t> inner_rank_hint = {});

// defined below; declared here so the result classes' write_into_slice()
// overrides can call it. The scatter inverse of slice_array_over_mode().
template <typename... Args>
void write_array_into_mode(TA::DistArray<Args...>& dest,
                           TA::DistArray<Args...> const& block,
                           std::size_t mode, std::size_t tile_lo,
                           std::size_t tile_hi,
                           std::optional<std::size_t> inner_rank_hint = {});

/// Partition a TiledRange1 into contiguous, tile-aligned element-range batches,
/// each covering at most \p target_batch_size elements: whole tiles are
/// appended to a batch until the next tile would push it over the target, so \p
/// target_batch_size is an UPPER BOUND, not a goal. Tiles are not uniformly
/// sized, so batches are uneven; the last batch may be smaller, and any single
/// tile larger than the target forms its own batch (the one-tile floor -- the
/// only way a batch may exceed the target).
/// Element ranges `[lo, hi)` are in the TiledRange1's element coordinate system
/// (honoring a nonzero element lobound, e.g. a frozen-core offset). This is the
/// TA realization of Result::mode_batches(): the caller requests a target batch
/// size in *elements*, which we convert to whole-tile groups. Returns at least
/// one batch for a non-empty mode.
[[nodiscard]] inline container::svector<std::pair<std::size_t, std::size_t>>
mode_batches_of_trange1(TA::TiledRange1 const& tr1,
                        std::size_t target_batch_size) {
  container::svector<std::pair<std::size_t, std::size_t>> batches;
  auto const& er = tr1.elements_range();
  if (er.second <= er.first) return batches;  // empty mode
  std::size_t const target = std::max<std::size_t>(target_batch_size, 1);
  std::size_t grp_lo = er.first;
  std::size_t acc = 0;
  std::size_t const ntiles = tr1.tile_extent();
  for (std::size_t t = 0; t < ntiles; ++t) {
    auto const& tr = tr1.tile(t);
    std::size_t const tsz = tr.second - tr.first;
    // target_batch_size is an UPPER BOUND on the realized (whole-tile) batch:
    // if appending this tile to a non-empty batch would push it over the
    // target, close the batch at the previous tile boundary first. A batch thus
    // never exceeds the target, except via the one-tile floor (a lone tile
    // larger than the target still forms its own batch).
    if (acc > 0 && acc + tsz > target) {
      batches.emplace_back(grp_lo, tr.first);
      grp_lo = tr.first;
      acc = 0;
    }
    acc += tsz;
    if (t + 1 == ntiles) batches.emplace_back(grp_lo, tr.second);
  }
  return batches;
}

///
/// \brief Result for a tensor value of TA::DistArray type.
/// \tparam ArrayT TA::DistArray type. Tile type of ArrayT is regular tensor of
///                scalars (not a tensor of tensors)
///
namespace detail {

/// the element coordinate at which the ⇓ half of a Kramers-union mode begins:
/// the first element of the tile in the middle of the mode's tile range (a
/// union axis has an even tile count, the ⇑ tiles then the ⇓ tiles)
[[nodiscard]] inline std::int64_t kramers_half_boundary(
    TA::TiledRange1 const& dim) {
  auto const n = dim.tile_extent();
  SEQUANT_ASSERT(n % 2 == 0 && "a Kramers-union mode has an even tile count");
  return static_cast<std::int64_t>(dim.tile(n / 2).first);
}

/// Fills \p bounds, one element boundary per entry of \p modes, from the
/// tiling \p tr when the view's builder did not give them (the modes are a
/// full union axis of \p tr then). Given bounds are kept: they are coordinates
/// of the index space and stay valid on a slice, a relabel, or after the sign
/// moved onto another operand of that space.
inline void fill_signed_bounds(TA::TiledRange const& tr,
                               container::svector<std::size_t> const& modes,
                               container::svector<std::int64_t>& bounds) {
  if (bounds.size() == modes.size()) return;
  bounds.clear();
  for (auto const m : modes) {
    SEQUANT_ASSERT(m < tr.rank());
    bounds.push_back(kramers_half_boundary(tr.dim(m)));
  }
}

/// A copy of @p arr with a pending transform applied: @p phase times (the
/// elementwise conjugate when @p conj of) every tile, and the factor -1 on
/// every tile whose element lobound along mode @p modes[k] is at or past
/// @p bounds[k] -- the ⇓ half of a Kramers-union axis, the sign of the
/// time-reversal flip (see ResultTensorTA::View::signed_halves). Tile by
/// tile, same tiling, shape and process map (the transform leaves every tile
/// norm unchanged); always a deep copy, never a tile shared with @p arr. The
/// copy's tiles are futures of local tasks: no fence is needed.
template <typename ArrayT>
ArrayT signed_copy(ArrayT const& arr, std::int8_t phase, bool conj,
                   container::svector<std::size_t> const& modes,
                   container::svector<std::int64_t> const& bounds) {
  using value_type = typename ArrayT::value_type;
  using numeric_type = typename ArrayT::numeric_type;
  SEQUANT_ASSERT(modes.size() == bounds.size());
  auto& world = arr.world();
  auto const& tr = arr.trange();
  auto const& tiles = tr.tiles_range();
  for (auto const m : modes) SEQUANT_ASSERT(m < tr.rank());
  ArrayT dst(world, tr, arr.shape(), arr.pmap());
  for (auto it = dst.pmap()->begin(); it != dst.pmap()->end(); ++it) {
    auto const ord = *it;
    if (dst.is_zero(ord)) continue;
    auto const idx = tiles.idx(ord);
    double sign = static_cast<double>(phase);
    for (std::size_t k = 0; k < modes.size(); ++k)
      if (static_cast<std::int64_t>(
              tr.dim(modes[k]).tile(idx[modes[k]]).first) >= bounds[k])
        sign = -sign;
    numeric_type const f = numeric_type(1) * sign;
    dst.set(ord, world.taskq.add(
                     [f, conj](value_type const& src) -> value_type {
                       if constexpr (TA::detail::is_complex_v<numeric_type>)
                         if (conj) return src.conj(f);
                       return src.scale(f);
                     },
                     arr.find(ord)));
  }
  return dst;
}

/// the stored mode that holds logical mode \p logical under \p perm (stored
/// mode s holds logical mode perm[s]; empty = identity)
[[nodiscard]] inline std::size_t stored_mode_of(
    container::svector<std::size_t> const& perm, std::size_t logical) {
  if (perm.empty()) return logical;
  for (std::size_t s = 0; s < perm.size(); ++s)
    if (perm[s] == logical) return s;
  SEQUANT_ASSERT(false && "stored_mode_of: no such logical mode");
  return logical;
}

/// splits a comma-separated label list (spaces ignored)
[[nodiscard]] inline container::svector<std::string> split_labels(
    std::string const& s) {
  container::svector<std::string> out;
  std::string cur;
  for (char c : s) {
    if (c == ',') {
      out.push_back(cur);
      cur.clear();
    } else if (c != ' ')
      cur += c;
  }
  if (!cur.empty() || !s.empty()) out.push_back(cur);
  return out;
}

[[nodiscard]] inline std::string join_labels(
    container::svector<std::string> const& toks) {
  std::string out;
  for (std::size_t i = 0; i < toks.size(); ++i) {
    if (i) out += ',';
    out += toks[i];
  }
  return out;
}

/// \p toks (logical order) reordered to the stored order under \p perm; a
/// list whose length is not \p perm's (the synthetic single inner label of a
/// zero nested array) is left as it is
inline void to_stored_order(container::svector<std::string>& toks,
                            container::svector<std::size_t> const& perm) {
  if (perm.empty() || toks.size() != perm.size()) return;
  container::svector<std::string> stored(toks.size());
  for (std::size_t s = 0; s < toks.size(); ++s) stored[s] = toks[perm[s]];
  toks = std::move(stored);
}

/// the annotation of the STORED array for the flat logical annotation
/// \p annot of a view relabeled by \p perm
[[nodiscard]] inline std::string stored_annot(
    std::string const& annot, container::svector<std::size_t> const& perm) {
  if (perm.empty()) return annot;
  auto toks = split_labels(annot);
  to_stored_order(toks, perm);
  return join_labels(toks);
}

/// the annotation of the STORED array for the nested logical annotation
/// \p annot (`outer;inner`) of a view relabeled by \p operm (outer) and
/// \p iperm (inner)
[[nodiscard]] inline std::string stored_annot_nested(
    std::string const& annot, container::svector<std::size_t> const& operm,
    container::svector<std::size_t> const& iperm) {
  auto const semi = annot.find(';');
  if (semi == std::string::npos) return stored_annot(annot, operm);
  auto outer = split_labels(annot.substr(0, semi));
  auto inner = split_labels(annot.substr(semi + 1));
  to_stored_order(outer, operm);
  to_stored_order(inner, iperm);
  return join_labels(outer) + ";" + join_labels(inner);
}

/// the TiledRange a view stands for: \p stored's dims cut to the view's tile
/// block (every mode keeps its element coordinates, as TA's
/// `preserve_lobound` block does) and reordered to the view's logical modes
[[nodiscard]] inline TA::TiledRange logical_trange(
    TA::TiledRange const& stored, container::svector<std::size_t> const& perm,
    container::svector<std::size_t> const& blk_lo,
    container::svector<std::size_t> const& blk_hi) {
  auto const rank = stored.rank();
  std::vector<TA::TiledRange1> dims;
  dims.reserve(rank);
  for (std::size_t m = 0; m < rank; ++m) {
    auto const s = stored_mode_of(perm, m);
    auto const& dim = stored.dim(s);
    if (blk_lo.empty()) {
      dims.push_back(dim);
      continue;
    }
    std::vector<std::size_t> bounds{
        static_cast<std::size_t>(dim.tile(blk_lo[s]).first)};
    for (std::size_t t = blk_lo[s]; t < blk_hi[s]; ++t)
      bounds.push_back(static_cast<std::size_t>(dim.tile(t).second));
    dims.emplace_back(bounds.begin(), bounds.end());
  }
  return TA::TiledRange(dims.begin(), dims.end());
}

/// the outer part of a possibly nested annotation
[[nodiscard]] inline std::string outer_part(std::string const& annot) {
  auto const semi = annot.find(';');
  return semi == std::string::npos ? annot : annot.substr(0, semi);
}

/// Where the half-signs of a product operand go. A sign on a mode that the
/// other operand shares (contracted or Hadamard) moves onto that operand:
/// sum_k (s_k A_k) B_k = sum_k A_k (s_k B_k); a sign on a mode only the
/// result carries moves onto the result; a mode only the signed operand
/// carries cannot move (nothing is planned then).
struct SignMovePlan {
  bool moved = false;
  /// STORED modes of the other operand that take a sign, and its boundary
  container::svector<std::size_t> y_modes;
  container::svector<std::int64_t> y_bounds;
  /// logical modes of the result that take a sign, and its boundary
  container::svector<std::size_t> result_modes;
  container::svector<std::int64_t> result_bounds;
};

/// \param x_modes,x_bounds the signed STORED modes of the signed operand
/// \param x_stored its (outer) labels in stored order
/// \param y_logical the other operand's (outer) labels in logical order
/// \param y_perm the other operand's relabeling (stored mode s holds logical
///        mode y_perm[s]); empty = identity
/// \param c_logical the result's (outer) labels
[[nodiscard]] inline SignMovePlan plan_sign_move(
    container::svector<std::size_t> const& x_modes,
    container::svector<std::int64_t> const& x_bounds,
    container::svector<std::string> const& x_stored,
    container::svector<std::string> const& y_logical,
    container::svector<std::size_t> const& y_perm,
    container::svector<std::string> const& c_logical) {
  SEQUANT_ASSERT(x_modes.size() == x_bounds.size());
  SignMovePlan p;
  for (std::size_t k = 0; k < x_modes.size(); ++k) {
    auto const s = x_modes[k];
    if (s >= x_stored.size()) return SignMovePlan{};
    auto const& label = x_stored[s];
    if (auto it = std::find(y_logical.begin(), y_logical.end(), label);
        it != y_logical.end()) {
      auto const ylog = static_cast<std::size_t>(it - y_logical.begin());
      p.y_modes.push_back(stored_mode_of(y_perm, ylog));
      p.y_bounds.push_back(x_bounds[k]);
    } else if (auto jt = std::find(c_logical.begin(), c_logical.end(), label);
               jt != c_logical.end()) {
      p.result_modes.push_back(
          static_cast<std::size_t>(jt - c_logical.begin()));
      p.result_bounds.push_back(x_bounds[k]);
    } else {
      return SignMovePlan{};
    }
  }
  p.moved = true;
  return p;
}

/// The two operands of a product with the half-signs of the larger signed
/// one moved (see plan_sign_move): onto the smaller operand, which then
/// materializes with them (detail::signed_copy) in place of the larger one,
/// and onto the result. Null pointers = nothing moved (each operand applies
/// its own signs). \p X and \p Y are result classes exposing outer_signed(),
/// outer_stored_labels(annot), outer_perm(), logical_volume() and
/// with_signs()/without_signs().
template <typename X, typename Y>
struct SignMove {
  ResultPtr x, y;
  container::svector<std::size_t> result_modes;
  container::svector<std::int64_t> result_bounds;
};

template <typename X, typename Y>
[[nodiscard]] SignMove<X, Y> move_signs(X const& x, std::string const& xannot,
                                        Y const& y, std::string const& yannot,
                                        std::string const& cannot) {
  SignMove<X, Y> out;
  auto try_move = [&](auto const& big, std::string const& bannot,
                      auto const& small, std::string const& sannot,
                      bool big_is_x) -> bool {
    if (!big.outer_signed() || small.outer_signed()) return false;
    if (big.logical_volume() <= small.logical_volume()) return false;
    auto const plan = plan_sign_move(
        big.signed_modes(), big.signed_bounds(),
        big.outer_stored_labels(bannot), split_labels(outer_part(sannot)),
        small.outer_perm(), split_labels(outer_part(cannot)));
    if (!plan.moved) return false;
    auto bp = big.without_signs();
    auto sp = small.with_signs(plan.y_modes, plan.y_bounds);
    out.x = big_is_x ? bp : sp;
    out.y = big_is_x ? sp : bp;
    out.result_modes = plan.result_modes;
    out.result_bounds = plan.result_bounds;
    return true;
  };
  if (!try_move(x, xannot, y, yannot, true))
    try_move(y, yannot, x, xannot, false);
  return out;
}

}  // namespace detail

template <typename ArrayT, typename = std::enable_if_t<TA::detail::is_tensor_v<
                               typename ArrayT::value_type>>>
class ResultTensorTA final : public Result {
 public:
  using Result::id_t;
  using numeric_type = typename ArrayT::numeric_type;

  /// A pending, zero-copy transform of the parent array the result shares:
  /// the view a result stands for instead of holding a value. Composable
  /// parts, all optional: a tile block (a batch slice), a mode relabeling, an
  /// elementwise conjugation, a phase, and the sign -1 on the ⇓ half of
  /// Kramers-union modes (the time-reversal flip, F = phase · conj ∘ swap ∘
  /// sign, whose swap MPQC serves as an alias array and whose sign,
  /// conjugation and phase ride here). The block and the relabeling are in
  /// the STORED array's mode order; the result's own (logical) modes are the
  /// stored ones permuted by \c perm. A consumer that is itself a TiledArray
  /// expression folds the view in through with_expr(), as `phase *
  /// arr(stored annot).block(blk_lo, blk_hi, preserve_lobound).conj()`; no
  /// expression spells a sign per half, so a signed view materializes
  /// (detail::signed_copy) when consumed, unless a product moves its signs
  /// onto a smaller operand (detail::move_signs). Materialization happens
  /// from a const read path and is unsynchronized: a Result is read and
  /// mutated by one thread at a time, which the evaluation engine guarantees.
  struct View {
    std::int8_t phase = 1;
    bool conj = false;
    /// stored mode \c s holds logical mode \c perm[s]; empty = identity
    container::svector<std::size_t> perm;
    /// tile bounds `[blk_lo[s], blk_hi[s])` per STORED mode; empty = the
    /// whole array
    container::svector<std::size_t> blk_lo, blk_hi;
    /// STORED modes along which the ⇓ half carries the factor -1
    container::svector<std::size_t> signed_halves;
    /// the element coordinate at which the ⇓ half begins, per entry of
    /// signed_halves; filled from the parent's tiling when left empty
    container::svector<std::int64_t> signed_bounds;

    [[nodiscard]] bool trivial() const noexcept {
      return phase == 1 && !conj && perm.empty() && blk_lo.empty() &&
             signed_halves.empty();
    }
    [[nodiscard]] std::size_t stored_mode(std::size_t logical) const {
      return detail::stored_mode_of(perm, logical);
    }
  };

  explicit ResultTensorTA(ArrayT arr) : Result{std::move(arr)} {}

  /// a view of \p parent with \p view pending; an alias of \p parent's tiles
  /// (also when the view is trivial: it reads the tiles as they are and still
  /// copies before it is written to, so no two results share tiles one of
  /// them accumulates into)
  ResultTensorTA(ArrayT parent, View view)
      : Result{}, parent_(std::move(parent)), view_(std::move(view)) {
    SEQUANT_ASSERT(parent_.is_initialized());
    detail::fill_signed_bounds(parent_.trange(), view_.signed_halves,
                               view_.signed_bounds);
  }

  /// @return whether this result is a pending view of another array (see
  ///         is_buffer_alias()); get<>() materializes it on first read
  [[nodiscard]] bool is_view() const noexcept { return !has_value(); }

  [[nodiscard]] bool is_buffer_alias() const noexcept override {
    return !has_value();
  }

  /// the pending view (meaningful while is_view())
  [[nodiscard]] View const& view() const noexcept { return view_; }

  /// the tiled range of the value this result stands for, with the pending
  /// block cut and the modes relabeled; never materializes
  [[nodiscard]] TA::TiledRange logical_trange() const {
    if (has_value()) return Result::get<ArrayT>().trange();
    return detail::logical_trange(parent_.trange(), view_.perm, view_.blk_lo,
                                  view_.blk_hi);
  }

  [[nodiscard]] std::string trange_annot() const override {
    return detail::ta_trange_string(logical_trange());
  }

  /// Hands \p f the TiledArray expression of this result's value for the
  /// logical annotation \p annot: `value(annot)` when the value is held, else
  /// the parent array with the pending view folded in -- its modes relabeled,
  /// the block cut by block(), the conjugation by conj(), the phase as a
  /// scalar factor. A view with half-signs is materialized first. \p f must
  /// return the same type whatever expression it receives.
  template <typename F>
  decltype(auto) with_expr(std::string const& annot, F&& f) const {
    if (!has_value() && !view_.signed_halves.empty()) ensure_materialized();
    if (has_value()) return std::forward<F>(f)(Result::get<ArrayT>()(annot));
    auto const stored = detail::stored_annot(annot, view_.perm);
    auto const ph = view_.phase;
    bool const conj = TA::detail::is_complex_v<numeric_type> && view_.conj;
    auto apply = [&](auto&& e) -> decltype(auto) {
      if constexpr (TA::detail::is_complex_v<numeric_type>) {
        if (conj && ph != 1)
          return std::forward<F>(f)(numeric_type(ph) * e.conj());
        if (conj) return std::forward<F>(f)(e.conj());
      }
      if (ph != 1) return std::forward<F>(f)(numeric_type(ph) * e);
      return std::forward<F>(f)(e);
    };
    if (view_.blk_lo.empty()) return apply(parent_(stored));
    return apply(parent_(stored).block(view_.blk_lo, view_.blk_hi,
                                       TA::preserve_lobound));
  }

  // --- the half-sign move protocol (detail::move_signs) ---
  [[nodiscard]] bool outer_signed() const noexcept {
    return !has_value() && !view_.signed_halves.empty();
  }
  [[nodiscard]] container::svector<std::size_t> const& signed_modes()
      const noexcept {
    return view_.signed_halves;
  }
  [[nodiscard]] container::svector<std::int64_t> const& signed_bounds()
      const noexcept {
    return view_.signed_bounds;
  }
  [[nodiscard]] container::svector<std::string> outer_stored_labels(
      std::string const& annot) const {
    return detail::split_labels(detail::stored_annot(annot, view_.perm));
  }
  [[nodiscard]] container::svector<std::size_t> outer_perm() const {
    return base_view().perm;
  }
  [[nodiscard]] std::size_t logical_volume() const {
    return logical_trange().elements_range().volume();
  }
  /// this result's view with the signs \p modes (STORED) added, as a result
  [[nodiscard]] ResultPtr with_signs(
      container::svector<std::size_t> const& modes,
      container::svector<std::int64_t> const& bounds) const {
    View v = base_view();
    for (std::size_t k = 0; k < modes.size(); ++k) {
      v.signed_halves.push_back(modes[k]);
      v.signed_bounds.push_back(bounds[k]);
    }
    return with_view(std::move(v));
  }
  /// this result's view with its signs dropped, as a result
  [[nodiscard]] ResultPtr without_signs() const {
    View v = base_view();
    v.signed_halves.clear();
    v.signed_bounds.clear();
    return with_view(std::move(v));
  }

 private:
  using this_type = ResultTensorTA<ArrayT>;
  using annot_wrap = Annot<std::string>;

  /// the parent array of a pending view: an empty handle once the value is
  /// held (mutable: released by the const materialization)
  mutable ArrayT parent_;
  View view_;

  [[nodiscard]] id_t type_id() const noexcept override {
    return id_for_type<this_type>();
  }

  // Diagnostic (see Result::fence): force all pending async work in this
  // array's world to complete so a lazily-executed op's timer captures
  // execution, not just dispatch.
  void fence() const noexcept override {
    view_base().world().gop.fence();
    ::sequant::detail::note_fence();
  }

  /// the array a further view of this result is a view of (the parent while
  /// this one is pending, else the held value), and that view
  [[nodiscard]] ArrayT const& view_base() const {
    return has_value() ? Result::get<ArrayT>() : parent_;
  }
  [[nodiscard]] View base_view() const { return has_value() ? View{} : view_; }
  [[nodiscard]] ResultPtr with_view(View v) const {
    return eval_result<this_type>(view_base(), std::move(v));
  }

  [[nodiscard]] std::size_t logical_rank() const {
    return view_base().trange().rank();
  }
  [[nodiscard]] std::string dummy_annot() const {
    return TA::detail::dummy_annotation(
        static_cast<unsigned int>(logical_rank()));
  }

  /// the value this result stands for, as a fresh array: the pending view
  /// applied in one pass (or a deep copy of the held value). A signed view
  /// takes the block, the signed copy, and the relabel, in turn.
  [[nodiscard]] ArrayT materialize() const {
    if (has_value()) return TA::clone(Result::get<ArrayT>());
    auto const ann = dummy_annot();
    auto const t0 = std::chrono::steady_clock::now();
    ArrayT out;
    if (!view_.signed_halves.empty()) {
      auto const stored = detail::stored_annot(ann, view_.perm);
      ArrayT blk;
      if (view_.blk_lo.empty())
        blk = parent_;
      else
        blk(stored) = parent_(stored).block(view_.blk_lo, view_.blk_hi,
                                            TA::preserve_lobound);
      auto const s =
          detail::signed_copy(blk, view_.phase, view_.conj, view_.signed_halves,
                              view_.signed_bounds);
      if (view_.perm.empty())
        out = s;
      else
        out(ann) = s(stored);
    } else {
      with_expr(ann, [&](auto&& e) { out(ann) = e; });
    }
    ::sequant::detail::drain_lazy_cleanup(out);
    detail::log_batch_op("Materialize", std::chrono::steady_clock::now() - t0,
                         out, ann);
    log_ta_tensor_host_memory_use();
    return out;
  }

  /// Result::get<ArrayT>() on a pending view: the view applied, once; the
  /// parent's tiles are released
  void ensure_materialized() const override {
    if (has_value()) return;
    reset_value(materialize());
    parent_ = ArrayT();
  }

  /// the logical labels of \p pre in \p post's order, composed onto this
  /// result's relabeling
  [[nodiscard]] View relabeled(std::string const& pre,
                               std::string const& post) const {
    View v = base_view();
    auto const pre_t = detail::split_labels(pre);
    auto const post_t = detail::split_labels(post);
    if (pre_t == post_t) return v;
    SEQUANT_ASSERT(pre_t.size() == post_t.size());
    container::svector<std::size_t> np(pre_t.size());
    bool identity = true;
    for (std::size_t s = 0; s < pre_t.size(); ++s) {
      auto const& label = pre_t[v.perm.empty() ? s : v.perm[s]];
      auto it = std::find(post_t.begin(), post_t.end(), label);
      SEQUANT_ASSERT(it != post_t.end());
      np[s] = static_cast<std::size_t>(it - post_t.begin());
      if (np[s] != s) identity = false;
    }
    if (identity)
      v.perm.clear();
    else
      v.perm = std::move(np);
    return v;
  }

  /// the extents of the value, for a failure message (no materialization)
  [[nodiscard]] std::string extents() const {
    std::string s = "[";
    for (auto const& e : logical_trange().elements_range().extent())
      s += std::to_string(e) + ",";
    return s + "]";
  }

  [[nodiscard]] ResultPtr sum(
      Result const& other,
      std::array<std::any, 3> const& annot) const override {
    SEQUANT_ASSERT(other.is<this_type>());
    auto const& o = static_cast<this_type const&>(other);
    auto const a = annot_wrap{annot};

    detail::log_ta(a.lannot, " + ", a.rannot, " = ", a.this_annot, "\n");

    ArrayT result;
    // a TA failure names the operands (as the product's wrapper does): the
    // raw TA message carries no annotation
    try {
      with_expr(a.lannot, [&](auto&& le) {
        o.with_expr(a.rannot,
                    [&](auto&& re) { result(a.this_annot) = le + re; });
      });
    } catch (std::exception const& ex) {
      throw std::runtime_error(std::string("TA sum failed: ") + ex.what() +
                               " | " + a.lannot + " + " + a.rannot + " -> " +
                               a.this_annot);
    }
    ::sequant::detail::drain_lazy_cleanup(result);
    log_ta_tensor_host_memory_use();
    return eval_result<this_type>(std::move(result));
  }

  /// A batch slice is a pending view: the tile block of \p mode that covers
  /// `[elem_lo, elem_hi)` (tile-aligned, see slice_bounds_to_tiles),
  /// intersected with the block this result already is. No tile is copied;
  /// the signs, the phase, the conjugation and the relabeling stay pending.
  [[nodiscard]] ResultPtr slice_mode(std::size_t mode, std::size_t elem_lo,
                                     std::size_t elem_hi) const override {
    ::sequant::detail::note_slice();
    auto const& base = view_base();
    auto view = base_view();
    auto const s = view.stored_mode(mode);
    auto const [tile_lo, tile_hi] =
        detail::slice_bounds_to_tiles(base.trange().dim(s), elem_lo, elem_hi);
    if (view.blk_lo.empty()) {
      auto const rank = base.trange().rank();
      view.blk_lo.assign(rank, 0);
      view.blk_hi.resize(rank);
      for (std::size_t d = 0; d < rank; ++d)
        view.blk_hi[d] = base.trange().dim(d).tile_extent();
    }
    SEQUANT_ASSERT(tile_lo >= view.blk_lo[s] && tile_hi <= view.blk_hi[s] &&
                   "slice_mode: the slice must lie within the block this "
                   "result already is");
    view.blk_lo[s] = tile_lo;
    view.blk_hi[s] = tile_hi;
    return eval_result<this_type>(base, std::move(view));
  }

  [[nodiscard]] double norm2() const override {
    if (has_value()) return detail::ta_norm2(Result::get<ArrayT>());
    if constexpr (requires(ArrayT const& a) {
                    { TA::squared_norm(a) } -> std::convertible_to<double>;
                  })
      return with_expr(dummy_annot(), [](auto&& e) {
        return static_cast<double>(e.norm().get());
      });
    else
      return -1.0;
  }
  [[nodiscard]] std::string layout_desc() const override {
    return detail::ta_trange_string(logical_trange());
  }
  [[nodiscard]] std::string tile_diff(Result const& other) const override {
    if (!other.is<this_type>()) return "(kind mismatch)";
    return detail::ta_tile_diff(get<ArrayT>(), other.get<ArrayT>());
  }

  /// Scatter \p block into this result's value (materialized first if this
  /// result is a pending view: a scatter writes, and a view must not write
  /// through to its parent). The block may itself be a pending view: it is
  /// folded into the assignment.
  void write_into_slice(Result const& block, std::size_t mode,
                        std::size_t block_lo, std::size_t block_hi) override {
    SEQUANT_ASSERT(block.is<this_type>());
    auto const& b = static_cast<this_type const&>(block);
    auto& dest = get<ArrayT>();
    auto const rank = dest.trange().rank();
    SEQUANT_ASSERT(mode < rank);
    auto const [tile_lo, tile_hi] = detail::slice_bounds_to_tiles(
        dest.trange().dim(mode), block_lo, block_hi);
    container::svector<std::size_t> lo(rank, 0), hi(rank);
    for (std::size_t d = 0; d < rank; ++d)
      hi[d] = dest.trange().dim(d).tile_extent();
    lo[mode] = tile_lo;
    hi[mode] = tile_hi;
    SEQUANT_ASSERT(b.logical_trange() ==
                       detail::logical_trange(dest.trange(), {}, lo, hi) &&
                   "write_into_slice: the block's layout is not the "
                   "destination sub-block's");
    auto const annot = dummy_annot();
    // preserve_lobound on the destination sub-block: the source block shares
    // its element coordinates (see write_array_into_mode)
    auto const t0 = std::chrono::steady_clock::now();
    b.with_expr(annot, [&](auto&& be) {
      dest(annot).block(lo, hi, TA::preserve_lobound) = be;
    });
    ::sequant::detail::drain_lazy_cleanup(dest);
    detail::log_batch_op("Scatter", std::chrono::steady_clock::now() - t0, dest,
                         annot + " mode=" + std::to_string(mode) + " tiles=[" +
                             std::to_string(tile_lo) + "," +
                             std::to_string(tile_hi) + ")");
    log_ta_tensor_host_memory_use();
  }

  [[nodiscard]] ResultPtr prod(Result const& other,
                               std::array<std::any, 3> const& annot,
                               DeNest DeNestFlag) const override {
    auto const a = annot_wrap{annot};

    if (other.is<ResultScalar<numeric_type>>()) {
      auto scalar = other.get<numeric_type>();

      detail::log_ta(a.lannot, " * ", scalar, " = ", a.this_annot, "\n");

      ArrayT result;
      with_expr(a.lannot,
                [&](auto&& le) { result(a.this_annot) = scalar * le; });
      ::sequant::detail::drain_lazy_cleanup(result);
      log_ta_tensor_host_memory_use();
      return eval_result<this_type>(std::move(result));
    }

    if (a.this_annot.empty()) {
      // the dot product reduces the two TA expressions, which carry the
      // pending transforms: neither operand's view is materialized
      SEQUANT_ASSERT(other.is<this_type>());
      auto const& od = static_cast<this_type const&>(other);
      numeric_type const d =
          with_expr(a.lannot, [&](auto&& le) -> numeric_type {
            return od.with_expr(a.rannot, [&](auto&& re) -> numeric_type {
              return TA::dot(le, re);
            });
          });
      ::sequant::detail::drain_lazy_cleanup(view_base());
      ::sequant::detail::drain_lazy_cleanup(od.view_base());

      detail::log_ta(a.lannot, " * ", a.rannot, " = ", d, "\n");

      log_ta_tensor_host_memory_use();
      return eval_result<ResultScalar<numeric_type>>(d);
    }

    if (!other.is<this_type>()) {
      // a nested operand takes the product (nested x flat, in either order)
      auto annot_swap = annot;
      std::swap(annot_swap[0], annot_swap[1]);
      return other.prod(*this, annot_swap, DeNestFlag);
    }

    auto const& o = static_cast<this_type const&>(other);
    detail::log_ta(a.lannot, " * ", a.rannot, " = ", a.this_annot, "\n");

    // The product is a TiledArray expression -- the expression layer
    // handles contractions, Hadamard and general products natively -- so an
    // operand that is a pending view (block, relabel, conjugation, phase) is
    // taken as is. The half-signs of the larger signed operand move onto the
    // smaller one and the result (detail::move_signs).
    auto const mv =
        detail::move_signs(*this, a.lannot, o, a.rannot, a.this_annot);
    this_type const& lx = mv.x ? static_cast<this_type const&>(*mv.x) : *this;
    this_type const& ly = mv.y ? static_cast<this_type const&>(*mv.y) : o;
    ArrayT result;
    try {
      lx.with_expr(a.lannot, [&](auto&& le) {
        ly.with_expr(a.rannot,
                     [&](auto&& re) { result(a.this_annot) = le * re; });
      });
    } catch (std::exception const& ex) {
      throw std::runtime_error("TA product failed: " + std::string(ex.what()) +
                               " | " + a.lannot + " * " + a.rannot + " -> " +
                               a.this_annot + " | A extents " + extents() +
                               " B extents " + o.extents());
    }
    if (!mv.result_modes.empty())
      result = detail::signed_copy(result, 1, false, mv.result_modes,
                                   mv.result_bounds);
    ::sequant::detail::drain_lazy_cleanup(result);
    log_ta_tensor_host_memory_use();
    return eval_result<this_type>(std::move(result));
  }

  [[nodiscard]] ResultPtr mult_by_phase(std::int8_t factor) const override {
    View v = base_view();
    v.phase = static_cast<std::int8_t>(v.phase * factor);
    return with_view(std::move(v));
  }

  /// Deep copy: \c TA::DistArray's own copy is a shallow (reference-counted)
  /// handle onto the same tiles, so an in-place accumulation into the copy
  /// would be seen by every other holder -- \c TA::clone allocates and copies
  /// the tiles. A pending view is applied into the copy directly, which is
  /// the one allocation either way.
  [[nodiscard]] ResultPtr clone() const override {
    auto const t0 = std::chrono::steady_clock::now();
    auto res = eval_result<this_type>(materialize());
    detail::log_batch_op("Clone", std::chrono::steady_clock::now() - t0,
                         res->template get<ArrayT>(), "");
    return res;
  }

  [[nodiscard]] ResultPtr permute(
      std::array<std::any, 2> const& ann) const override {
    auto const pre_annot = std::any_cast<std::string>(ann[0]);
    auto const post_annot = std::any_cast<std::string>(ann[1]);

    detail::log_ta(pre_annot, " = ", post_annot, " (view)\n");

    return with_view(relabeled(pre_annot, post_annot));
  }

  [[nodiscard]] ResultPtr apply_transform(
      CanonTransform t, std::array<std::any, 2> const& ann) const override {
    auto const pre_annot = std::any_cast<std::string>(ann[0]);
    auto const post_annot = std::any_cast<std::string>(ann[1]);

    detail::log_ta(post_annot, " = apply_transform(", pre_annot, ") (view)\n");

    View v = relabeled(pre_annot, post_annot);
    v.phase = static_cast<std::int8_t>(v.phase * t.phase);
    if constexpr (TA::detail::is_complex_v<numeric_type>)
      if (t.conj) v.conj = !v.conj;
    return with_view(std::move(v));
  }

  /// Accumulate \p other into this result's value (materialized first if
  /// this result is a pending view: no view is written through). The operand
  /// may be a pending view: it is folded into the accumulation.
  void add_inplace(Result const& other) override {
    SEQUANT_ASSERT(other.is<this_type>());
    auto const& o = static_cast<this_type const&>(other);

    auto& t = get<ArrayT>();
    SEQUANT_ASSERT(t.trange() == o.logical_trange());
    auto ann = TA::detail::dummy_annotation(t.trange().rank());

    detail::log_ta(ann, " += ", ann, "\n");

    auto const t0 = std::chrono::steady_clock::now();
    o.with_expr(ann, [&](auto&& oe) { t(ann) += oe; });
    ::sequant::detail::drain_lazy_cleanup(t);
    detail::log_batch_op("Accumulate", std::chrono::steady_clock::now() - t0, t,
                         ann);
    log_ta_tensor_host_memory_use();
  }

  [[nodiscard]] ResultPtr symmetrize() const override {
    return eval_result<this_type>(detail::column_symmetrize_ta(get<ArrayT>()));
  }

  [[nodiscard]] ResultPtr antisymmetrize(size_t bra_rank) const override {
    return eval_result<this_type>(
        detail::particle_antisymmetrize_ta(get<ArrayT>(), bra_rank));
  }

 private:
  [[nodiscard]] std::size_t size_in_bytes() const final {
    if (Logger::instance().eval.level == 0)
      return 0;  // size_of disabled untraced: it is collective
    return measured_size_in_bytes();
  }

  /// the logical size of the value: the held array's, or the parent's scaled
  /// to the block a pending view cuts (exact without a block; a dense
  /// estimate of the block otherwise). A tracer charges a view 0 allocated
  /// bytes through is_buffer_alias(); this is what it stands for.
  [[nodiscard]] std::size_t measured_size_in_bytes() const final {
    auto const& v = view_base();
    auto local_size = TA::size_of<TA::MemorySpace::Host>(v);
    v.world().gop.sum(local_size);
    if (has_value() || view_.blk_lo.empty()) return local_size;
    auto const full = v.trange().elements_range().volume();
    if (full == 0) return 0;
    return static_cast<std::size_t>(
        static_cast<double>(local_size) *
        (static_cast<double>(logical_volume()) / static_cast<double>(full)));
  }
};

template <typename ArrayT,
          typename = std::enable_if_t<
              TA::detail::is_tensor_of_tensor_v<typename ArrayT::value_type>>>
class ResultTensorOfTensorTA final : public Result {
 public:
  using Result::id_t;
  using numeric_type = typename ArrayT::numeric_type;

  /// The nested result's pending view (see ResultTensorTA::View): a tile
  /// block over the OUTER modes, an outer and an inner relabeling (stored
  /// mode s holds logical mode operm[s] / iperm[s]), a conjugation, a phase,
  /// and the half-signs on outer modes.
  struct View {
    std::int8_t phase = 1;
    bool conj = false;
    container::svector<std::size_t> operm;
    container::svector<std::size_t> iperm;
    container::svector<std::size_t> blk_lo, blk_hi;
    container::svector<std::size_t> signed_halves;
    container::svector<std::int64_t> signed_bounds;

    [[nodiscard]] bool trivial() const noexcept {
      return phase == 1 && !conj && operm.empty() && iperm.empty() &&
             blk_lo.empty() && signed_halves.empty();
    }
    [[nodiscard]] bool relabeled() const noexcept {
      return !operm.empty() || !iperm.empty();
    }
    [[nodiscard]] std::size_t stored_mode(std::size_t logical) const {
      return detail::stored_mode_of(operm, logical);
    }
  };

  explicit ResultTensorOfTensorTA(ArrayT arr) : Result{std::move(arr)} {}

  /// a view of \p parent with \p view pending (see the flat result)
  ResultTensorOfTensorTA(ArrayT parent, View view)
      : Result{}, parent_(std::move(parent)), view_(std::move(view)) {
    SEQUANT_ASSERT(parent_.is_initialized());
    detail::fill_signed_bounds(parent_.trange(), view_.signed_halves,
                               view_.signed_bounds);
  }

  [[nodiscard]] bool is_view() const noexcept { return !has_value(); }

  [[nodiscard]] bool is_buffer_alias() const noexcept override {
    return !has_value();
  }

  [[nodiscard]] View const& view() const noexcept { return view_; }

  /// the rank of the inner tensors, read from a populated inner tile of the
  /// array this result stands for (detail::tot_inner_rank(), a
  /// world-collective reduction) and memoized once known; 0 while every inner
  /// tile is empty
  [[nodiscard]] std::size_t inner_rank() const {
    if (inner_rank_) return *inner_rank_;
    auto const r = detail::tot_inner_rank(view_base());
    if (r > 0) inner_rank_ = r;
    return r;
  }
  void note_inner_rank(std::size_t r) const noexcept { inner_rank_ = r; }

  /// the OUTER tiled range of the value this result stands for, with the
  /// pending block cut and the outer modes relabeled; never materializes
  [[nodiscard]] TA::TiledRange logical_trange() const {
    if (has_value()) return Result::get<ArrayT>().trange();
    return detail::logical_trange(parent_.trange(), view_.operm, view_.blk_lo,
                                  view_.blk_hi);
  }

  [[nodiscard]] std::string trange_annot() const override {
    return detail::ta_trange_string(logical_trange());
  }

  /// the STORED annotation for the logical nested annotation \p logical
  [[nodiscard]] std::string translate(std::string const& logical) const {
    if (has_value()) return logical;
    return detail::stored_annot_nested(logical, view_.operm, view_.iperm);
  }

  /// Hands \p f the TiledArray expression of this result's value for the
  /// nested logical annotation \p annot (`outer;inner`); see
  /// ResultTensorTA::with_expr. A view with half-signs is materialized first.
  template <typename F>
  decltype(auto) with_expr(std::string const& annot, F&& f) const {
    if (!has_value() && !view_.signed_halves.empty()) ensure_materialized();
    if (has_value()) return std::forward<F>(f)(Result::get<ArrayT>()(annot));
    auto const stored = translate(annot);
    auto const ph = view_.phase;
    bool const conj = TA::detail::is_complex_v<numeric_type> && view_.conj;
    auto apply = [&](auto&& e) -> decltype(auto) {
      if constexpr (TA::detail::is_complex_v<numeric_type>) {
        if (conj && ph != 1)
          return std::forward<F>(f)(numeric_type(ph) * e.conj());
        if (conj) return std::forward<F>(f)(e.conj());
      }
      if (ph != 1) return std::forward<F>(f)(numeric_type(ph) * e);
      return std::forward<F>(f)(e);
    };
    if (view_.blk_lo.empty()) return apply(parent_(stored));
    return apply(parent_(stored).block(view_.blk_lo, view_.blk_hi,
                                       TA::preserve_lobound));
  }

  // --- the half-sign move protocol (detail::move_signs) ---
  [[nodiscard]] bool outer_signed() const noexcept {
    return !has_value() && !view_.signed_halves.empty();
  }
  [[nodiscard]] container::svector<std::size_t> const& signed_modes()
      const noexcept {
    return view_.signed_halves;
  }
  [[nodiscard]] container::svector<std::int64_t> const& signed_bounds()
      const noexcept {
    return view_.signed_bounds;
  }
  [[nodiscard]] container::svector<std::string> outer_stored_labels(
      std::string const& annot) const {
    return detail::split_labels(detail::outer_part(translate(annot)));
  }
  [[nodiscard]] container::svector<std::size_t> outer_perm() const {
    return base_view().operm;
  }
  [[nodiscard]] std::size_t logical_volume() const {
    return logical_trange().elements_range().volume();
  }
  [[nodiscard]] ResultPtr with_signs(
      container::svector<std::size_t> const& modes,
      container::svector<std::int64_t> const& bounds) const {
    View v = base_view();
    for (std::size_t k = 0; k < modes.size(); ++k) {
      v.signed_halves.push_back(modes[k]);
      v.signed_bounds.push_back(bounds[k]);
    }
    return with_view(std::move(v));
  }
  [[nodiscard]] ResultPtr without_signs() const {
    View v = base_view();
    v.signed_halves.clear();
    v.signed_bounds.clear();
    return with_view(std::move(v));
  }

 private:
  using this_type = ResultTensorOfTensorTA<ArrayT>;
  using annot_wrap = Annot<std::string>;

  using _inner_tensor_type = typename ArrayT::value_type::value_type;

  // "Regular" (non-nested) companion array for ToT * T products. The OUTER
  // tile type must be a TA::Tensor -- inner tile types like btas::Tensor are
  // only valid as the *innermost* tile.
  using compatible_regular_distarray_type =
      TA::DistArray<TA::Tensor<numeric_type>, typename ArrayT::policy_type>;

  // Only @c that_type type is allowed for ToT * T computation
  using that_type = ResultTensorTA<compatible_regular_distarray_type>;

  mutable ArrayT parent_;
  View view_;
  mutable std::optional<std::size_t> inner_rank_;

  static std::size_t inner_count(std::string const& annot) {
    auto const semi = annot.find(';');
    if (semi == std::string::npos) return 0;
    return detail::split_labels(annot.substr(semi + 1)).size();
  }

  static ResultPtr noted(ResultPtr r, std::optional<std::size_t> rank) {
    static_cast<this_type&>(*r).inner_rank_ = rank;
    return r;
  }

  [[nodiscard]] id_t type_id() const noexcept override {
    return id_for_type<this_type>();
  }

  void fence() const noexcept override {
    view_base().world().gop.fence();
    ::sequant::detail::note_fence();
  }

  [[nodiscard]] ArrayT const& view_base() const {
    return has_value() ? Result::get<ArrayT>() : parent_;
  }
  [[nodiscard]] View base_view() const { return has_value() ? View{} : view_; }
  [[nodiscard]] ResultPtr with_view(View v) const {
    return noted(eval_result<this_type>(view_base(), std::move(v)),
                 inner_rank_);
  }

  [[nodiscard]] std::string dummy_annot(std::size_t inner) const {
    return TA::detail::dummy_annotation(
        static_cast<unsigned int>(view_base().trange().rank()),
        static_cast<unsigned int>(inner));
  }

  /// A zero nested array over \p tr: the value of a view of an array whose
  /// inner tiles are all empty (there is no inner rank to annotate a block
  /// with, and the slice of zero is zero). Dense even under a sparse policy
  /// -- harmless, the array is identically zero. Every rank sets its own
  /// tiles with ready values: no fence.
  [[nodiscard]] ArrayT zero_over(TA::TiledRange const& tr) const {
    using value_type = typename ArrayT::value_type;
    ArrayT out{view_base().world(), tr};
    for (auto it = out.begin(); it != out.end(); ++it)
      if (out.is_local(it.index())) *it = value_type{it.make_range()};
    return out;
  }

  [[nodiscard]] ArrayT materialize() const {
    if (has_value()) return TA::clone(Result::get<ArrayT>());
    auto const inner = inner_rank();
    if (inner == 0) return zero_over(logical_trange());
    auto const ann = dummy_annot(inner);
    auto const t0 = std::chrono::steady_clock::now();
    ArrayT out;
    if (!view_.signed_halves.empty()) {
      auto const stored = translate(ann);
      ArrayT blk;
      if (view_.blk_lo.empty())
        blk = parent_;
      else
        blk(stored) = parent_(stored).block(view_.blk_lo, view_.blk_hi,
                                            TA::preserve_lobound);
      auto const s =
          detail::signed_copy(blk, view_.phase, view_.conj, view_.signed_halves,
                              view_.signed_bounds);
      if (view_.relabeled())
        out(ann) = s(stored);
      else
        out = s;
    } else {
      with_expr(ann, [&](auto&& e) { out(ann) = e; });
    }
    ::sequant::detail::drain_lazy_cleanup(out);
    detail::log_batch_op("Materialize", std::chrono::steady_clock::now() - t0,
                         out, ann);
    log_ta_tensor_host_memory_use();
    return out;
  }

  void ensure_materialized() const override {
    if (has_value()) return;
    reset_value(materialize());
    parent_ = ArrayT();
  }

  /// the value as an array (materialized if pending)
  [[nodiscard]] ArrayT const& logical_array() const { return get<ArrayT>(); }

  static void compose_perm(container::svector<std::string> const& pre,
                           container::svector<std::string> const& post,
                           container::svector<std::size_t>& cur) {
    if (pre == post) return;
    SEQUANT_ASSERT(pre.size() == post.size());
    container::svector<std::size_t> np(pre.size());
    bool identity = true;
    for (std::size_t s = 0; s < pre.size(); ++s) {
      auto const& label = pre[cur.empty() ? s : cur[s]];
      auto it = std::find(post.begin(), post.end(), label);
      SEQUANT_ASSERT(it != post.end());
      np[s] = static_cast<std::size_t>(it - post.begin());
      if (np[s] != s) identity = false;
    }
    if (identity)
      cur.clear();
    else
      cur = std::move(np);
  }

  [[nodiscard]] View relabeled(std::string const& pre,
                               std::string const& post) const {
    View v = base_view();
    auto split = [](std::string const& s) {
      auto const semi = s.find(';');
      return std::pair{
          detail::split_labels(semi == std::string::npos ? s
                                                         : s.substr(0, semi)),
          semi == std::string::npos ? container::svector<std::string>{}
                                    : detail::split_labels(s.substr(semi + 1))};
    };
    auto const [pre_o, pre_i] = split(pre);
    auto const [post_o, post_i] = split(post);
    compose_perm(pre_o, post_o, v.operm);
    compose_perm(pre_i, post_i, v.iperm);
    return v;
  }

  [[nodiscard]] std::string extents() const {
    std::string s = "[";
    for (auto const& e : logical_trange().elements_range().extent())
      s += std::to_string(e) + ",";
    return s + "]";
  }

  /// Allocating addition. For nested arrays an addend whose inner tiles are
  /// all empty acts as the additive identity (the arena nested-tensor
  /// addition the pinned TiledArray provides). A relabeled inner view is
  /// materialized first: the sum takes no inner permutation.
  [[nodiscard]] ResultPtr sum(
      Result const& other,
      std::array<std::any, 3> const& annot) const override {
    SEQUANT_ASSERT(other.is<this_type>());
    auto const& o = static_cast<this_type const&>(other);
    auto const a = annot_wrap{annot};

    detail::log_ta(a.lannot, " + ", a.rannot, " = ", a.this_annot, "\n");

    if (detail::warn_tot_sum_inner()) {
      auto inner = [](std::string const& s) {
        auto p = s.find(';');
        return p == std::string::npos ? std::string{} : s.substr(p + 1);
      };
      const auto li = inner(translate(a.lannot)),
                 ri = inner(o.translate(a.rannot)), ti = inner(a.this_annot);
      if (li != ri || li != ti)
        std::cerr << "[sequant-eval] WARNING: nested sum with mismatched inner "
                     "annotations: "
                  << a.lannot << " + " << a.rannot << " = " << a.this_annot
                  << "\n";
    }
    if (!has_value() && !view_.iperm.empty()) ensure_materialized();
    if (!o.has_value() && !o.view_.iperm.empty()) o.ensure_materialized();

    ArrayT result;
    try {
      with_expr(a.lannot, [&](auto&& le) {
        o.with_expr(a.rannot,
                    [&](auto&& re) { result(a.this_annot) = le + re; });
      });
    } catch (std::exception const& ex) {
      throw std::runtime_error(std::string("TA sum failed: ") + ex.what() +
                               " | " + a.lannot + " + " + a.rannot + " -> " +
                               a.this_annot);
    }
    ::sequant::detail::drain_lazy_cleanup(result);
    log_ta_tensor_host_memory_use();
    return noted(eval_result<this_type>(std::move(result)),
                 inner_count(a.this_annot));
  }

  /// A batch slice is a pending view (see ResultTensorTA::slice_mode); the
  /// inner rank is the parent's
  [[nodiscard]] ResultPtr slice_mode(std::size_t mode, std::size_t elem_lo,
                                     std::size_t elem_hi) const override {
    ::sequant::detail::note_slice();
    auto const& base = view_base();
    auto view = base_view();
    auto const s = view.stored_mode(mode);
    auto const [tile_lo, tile_hi] =
        detail::slice_bounds_to_tiles(base.trange().dim(s), elem_lo, elem_hi);
    if (view.blk_lo.empty()) {
      auto const rank = base.trange().rank();
      view.blk_lo.assign(rank, 0);
      view.blk_hi.resize(rank);
      for (std::size_t d = 0; d < rank; ++d)
        view.blk_hi[d] = base.trange().dim(d).tile_extent();
    }
    SEQUANT_ASSERT(tile_lo >= view.blk_lo[s] && tile_hi <= view.blk_hi[s] &&
                   "slice_mode: the slice must lie within the block this "
                   "result already is");
    view.blk_lo[s] = tile_lo;
    view.blk_hi[s] = tile_hi;
    return noted(eval_result<this_type>(base, std::move(view)), inner_rank_);
  }

  [[nodiscard]] double norm2() const override {
    return detail::ta_norm2(get<ArrayT>());
  }
  [[nodiscard]] std::string layout_desc() const override {
    return detail::ta_trange_string(logical_trange());
  }
  [[nodiscard]] std::string tile_diff(Result const& other) const override {
    if (!other.is<this_type>()) return "(kind mismatch)";
    return detail::ta_tile_diff(get<ArrayT>(), other.get<ArrayT>());
  }

  /// Scatter \p block into this result's value (materialized first if this
  /// result is a pending view). The block may itself be a pending view. A
  /// block whose inner tiles are all empty is zero and leaves the
  /// destination slice as it was.
  void write_into_slice(Result const& block, std::size_t mode,
                        std::size_t block_lo, std::size_t block_hi) override {
    SEQUANT_ASSERT(block.is<this_type>());
    auto const& b = static_cast<this_type const&>(block);
    auto& dest = get<ArrayT>();
    auto const rank = dest.trange().rank();
    SEQUANT_ASSERT(mode < rank);
    auto const [tile_lo, tile_hi] = detail::slice_bounds_to_tiles(
        dest.trange().dim(mode), block_lo, block_hi);
    auto const inner = b.inner_rank();
    if (inner == 0) return;
    container::svector<std::size_t> lo(rank, 0), hi(rank);
    for (std::size_t d = 0; d < rank; ++d)
      hi[d] = dest.trange().dim(d).tile_extent();
    lo[mode] = tile_lo;
    hi[mode] = tile_hi;
    SEQUANT_ASSERT(b.logical_trange() ==
                       detail::logical_trange(dest.trange(), {}, lo, hi) &&
                   "write_into_slice: the block's layout is not the "
                   "destination sub-block's");
    // a relabeled inner view is materialized first (see sum)
    if (!b.has_value() && !b.view_.iperm.empty()) b.ensure_materialized();
    auto const annot = TA::detail::dummy_annotation(
        static_cast<unsigned int>(rank), static_cast<unsigned int>(inner));
    auto const t0 = std::chrono::steady_clock::now();
    b.with_expr(annot, [&](auto&& be) {
      dest(annot).block(lo, hi, TA::preserve_lobound) = be;
    });
    ::sequant::detail::drain_lazy_cleanup(dest);
    inner_rank_ = inner;
    detail::log_batch_op("Scatter", std::chrono::steady_clock::now() - t0, dest,
                         annot + " mode=" + std::to_string(mode) + " tiles=[" +
                             std::to_string(tile_lo) + "," +
                             std::to_string(tile_hi) + ")");
    log_ta_tensor_host_memory_use();
  }

  /// A zero nested array shaped like this result's value with \p mode
  /// widened to the full tiling that \p axis_src carries along
  /// \p axis_src_mode: the pre-sized destination of a batched scatter. The
  /// tilings are read from the results' layouts; nothing is materialized.
  [[nodiscard]] ResultPtr pre_sized_zeros_over_mode(
      std::size_t mode, Result const& axis_src,
      std::size_t axis_src_mode) const override {
    auto const self_tr = logical_trange();
    auto const rank = self_tr.rank();
    SEQUANT_ASSERT(mode < rank);
    TA::TiledRange1 const axis_dim = [&]() -> TA::TiledRange1 {
      if (axis_src.is<this_type>()) {
        auto const tr =
            static_cast<this_type const&>(axis_src).logical_trange();
        SEQUANT_ASSERT(axis_src_mode < tr.rank());
        return tr.dim(axis_src_mode);
      }
      SEQUANT_ASSERT(axis_src.is<that_type>());
      auto const tr = static_cast<that_type const&>(axis_src).logical_trange();
      SEQUANT_ASSERT(axis_src_mode < tr.rank());
      return tr.dim(axis_src_mode);
    }();
    std::vector<TA::TiledRange1> dims;
    dims.reserve(rank);
    for (std::size_t d = 0; d < rank; ++d) dims.push_back(self_tr.dim(d));
    dims[mode] = axis_dim;
    auto dest = zero_over(TA::TiledRange(dims.begin(), dims.end()));
    log_ta_tensor_host_memory_use();
    return noted(eval_result<this_type>(std::move(dest)), std::size_t{0});
  }

  /// Emit the general product `le * re` into \p result for the annotations
  /// \p a. TiledArray's contraction engine produces the result's inner modes
  /// in their natural order (the left operand's surviving inner labels in
  /// left order, then the right's); a result wanting another inner order is
  /// produced in the natural one and permuted into place.
  template <typename ResultArrayT, typename L, typename R>
  static void emit_product(ResultArrayT& result, L&& le, R&& re,
                           annot_wrap const& a) {
    if (!detail::tot_product_needs_inner_reorder(a.lannot, a.rannot,
                                                 a.this_annot)) {
      result(a.this_annot) = le * re;
      return;
    }
    auto const natural =
        detail::tot_product_natural_annot(a.lannot, a.rannot, a.this_annot);
    ResultArrayT tmp;
    tmp(natural) = le * re;
    ::sequant::detail::drain_lazy_cleanup(tmp);
    result(a.this_annot) = tmp(natural);
  }

  [[nodiscard]] ResultPtr prod(Result const& other,
                               std::array<std::any, 3> const& annot,
                               DeNest DeNestFlag) const override {
    auto const a = annot_wrap{annot};

    if (other.is<ResultScalar<numeric_type>>()) {
      auto scalar = other.get<numeric_type>();

      detail::log_ta(a.lannot, " * ", scalar, " = ", a.this_annot, "\n");

      ArrayT result;
      with_expr(a.lannot,
                [&](auto&& le) { result(a.this_annot) = scalar * le; });
      ::sequant::detail::drain_lazy_cleanup(result);
      log_ta_tensor_host_memory_use();
      return noted(eval_result<this_type>(std::move(result)),
                   inner_count(a.this_annot));
    } else if (a.this_annot.empty()) {
      // DOT product: the two expressions are reduced, views taken as is
      SEQUANT_ASSERT(other.is<this_type>());
      auto const& o = static_cast<this_type const&>(other);
      numeric_type const d =
          with_expr(a.lannot, [&](auto&& le) -> numeric_type {
            return o.with_expr(a.rannot, [&](auto&& re) -> numeric_type {
              return TA::dot(le, re);
            });
          });
      ::sequant::detail::drain_lazy_cleanup(view_base());
      ::sequant::detail::drain_lazy_cleanup(o.view_base());

      detail::log_ta(a.lannot, " * ", a.rannot, " = ", d, "\n");

      log_ta_tensor_host_memory_use();
      return eval_result<ResultScalar<numeric_type>>(d);
    }

    detail::log_ta(a.lannot, " * ", a.rannot, " = ", a.this_annot, "\n");

    auto const describe = [&](std::string const& what) {
      return "TA ToT product failed: " + what + " | " + a.lannot + " * " +
             a.rannot + " -> " + a.this_annot + " | left outer extents " +
             extents();
    };

    // Every product is a TiledArray expression, so an operand that is a
    // pending view (a block, a relabel, a conjugation, a phase) is taken as
    // is; the half-signs of the larger signed operand move onto the smaller
    // one and the result (detail::move_signs).
    if (other.is<that_type>()) {
      // ToT * T -> ToT
      auto const& o = static_cast<that_type const&>(other);
      auto const mv =
          detail::move_signs(*this, a.lannot, o, a.rannot, a.this_annot);
      this_type const& lx = mv.x ? static_cast<this_type const&>(*mv.x) : *this;
      that_type const& ly = mv.y ? static_cast<that_type const&>(*mv.y) : o;
      ArrayT result;
      try {
        lx.with_expr(a.lannot, [&](auto&& le) {
          ly.with_expr(a.rannot,
                       [&](auto&& re) { emit_product(result, le, re, a); });
        });
      } catch (std::exception const& ex) {
        throw std::runtime_error(describe(ex.what()));
      }
      if (!mv.result_modes.empty())
        result = detail::signed_copy(result, 1, false, mv.result_modes,
                                     mv.result_bounds);
      ::sequant::detail::drain_lazy_cleanup(result);
      log_ta_tensor_host_memory_use();
      return noted(eval_result<this_type>(std::move(result)),
                   inner_count(a.this_annot));

    } else if (other.is<this_type>() && DeNestFlag == DeNest::True) {
      // ToT * ToT -> T: the inner tensors contracted away
      auto const& o = static_cast<this_type const&>(other);
      auto const mv =
          detail::move_signs(*this, a.lannot, o, a.rannot, a.this_annot);
      this_type const& lx = mv.x ? static_cast<this_type const&>(*mv.x) : *this;
      this_type const& ly = mv.y ? static_cast<this_type const&>(*mv.y) : o;
      compatible_regular_distarray_type result;
      try {
        lx.with_expr(a.lannot, [&](auto&& le) {
          ly.with_expr(a.rannot, [&](auto&& re) {
            result(a.this_annot) = le.dot_inner(re);
          });
        });
      } catch (std::exception const& ex) {
        throw std::runtime_error(describe(ex.what()));
      }
      if (!mv.result_modes.empty())
        result = detail::signed_copy(result, 1, false, mv.result_modes,
                                     mv.result_bounds);
      ::sequant::detail::drain_lazy_cleanup(result);
      log_ta_tensor_host_memory_use();
      return eval_result<that_type>(std::move(result));

    } else if (other.is<this_type>() && DeNestFlag == DeNest::False) {
      // ToT * ToT -> ToT
      auto const& o = static_cast<this_type const&>(other);
      auto const mv =
          detail::move_signs(*this, a.lannot, o, a.rannot, a.this_annot);
      this_type const& lx = mv.x ? static_cast<this_type const&>(*mv.x) : *this;
      this_type const& ly = mv.y ? static_cast<this_type const&>(*mv.y) : o;
      ArrayT result;
      try {
        lx.with_expr(a.lannot, [&](auto&& le) {
          ly.with_expr(a.rannot,
                       [&](auto&& re) { emit_product(result, le, re, a); });
        });
      } catch (std::exception const& ex) {
        throw std::runtime_error(describe(ex.what()));
      }
      if (!mv.result_modes.empty())
        result = detail::signed_copy(result, 1, false, mv.result_modes,
                                     mv.result_bounds);
      ::sequant::detail::drain_lazy_cleanup(result);
      log_ta_tensor_host_memory_use();
      return noted(eval_result<this_type>(std::move(result)),
                   inner_count(a.this_annot));
    } else {
      throw detail::invalid_operand();
    }
  }

  [[nodiscard]] ResultPtr mult_by_phase(std::int8_t factor) const override {
    View v = base_view();
    v.phase = static_cast<std::int8_t>(v.phase * factor);
    return with_view(std::move(v));
  }

  /// Deep copy; see \c ResultTensorTA::clone. \c TA::clone deep-copies the
  /// (nested) tiles too, so the copy shares no inner tensor with this one. A
  /// pending view is applied into the copy directly.
  [[nodiscard]] ResultPtr clone() const override {
    auto const t0 = std::chrono::steady_clock::now();
    auto res = noted(eval_result<this_type>(materialize()), inner_rank_);
    detail::log_batch_op("Clone", std::chrono::steady_clock::now() - t0,
                         res->template get<ArrayT>(), "");
    return res;
  }

  [[nodiscard]] ResultPtr permute(
      std::array<std::any, 2> const& ann) const override {
    auto const pre_annot = std::any_cast<std::string>(ann[0]);
    auto const post_annot = std::any_cast<std::string>(ann[1]);

    detail::log_ta(pre_annot, " = ", post_annot, " (view)\n");

    return with_view(relabeled(pre_annot, post_annot));
  }

  [[nodiscard]] ResultPtr apply_transform(
      CanonTransform t, std::array<std::any, 2> const& ann) const override {
    auto const pre_annot = std::any_cast<std::string>(ann[0]);
    auto const post_annot = std::any_cast<std::string>(ann[1]);

    detail::log_ta(post_annot, " = apply_transform(", pre_annot, ") (view)\n");

    View v = relabeled(pre_annot, post_annot);
    v.phase = static_cast<std::int8_t>(v.phase * t.phase);
    if constexpr (TA::detail::is_complex_v<numeric_type>)
      if (t.conj) v.conj = !v.conj;
    return with_view(std::move(v));
  }

  /// Accumulate \p other into this result's value (materialized first if
  /// this result is a pending view). A relabeled inner view operand is
  /// materialized first (see sum); any other view is folded in.
  void add_inplace(Result const& other) override {
    SEQUANT_ASSERT(other.is<this_type>());
    auto const& o = static_cast<this_type const&>(other);

    auto& t = get<ArrayT>();
    if (!o.has_value() && !o.view_.iperm.empty()) o.ensure_materialized();
    SEQUANT_ASSERT(t.trange() == o.logical_trange());
    // ToT annotation needs an inner block ("outer;inner"); an operand whose
    // inner tiles are all empty has no rank to give. Take the rank from
    // whichever operand has data; if both are empty the add is an identity
    // no-op (t += 0), nothing to annotate.
    auto const inner = std::max(this->inner_rank(), o.inner_rank());
    if (inner == 0) return;
    auto ann = TA::detail::dummy_annotation(t.trange().rank(), inner);

    detail::log_ta(ann, " += ", ann, "\n");

    auto const t0 = std::chrono::steady_clock::now();
    o.with_expr(ann, [&](auto&& oe) { t(ann) += oe; });
    ::sequant::detail::drain_lazy_cleanup(t);
    inner_rank_ = inner;  // populated now, if it was not already
    detail::log_batch_op("Accumulate", std::chrono::steady_clock::now() - t0, t,
                         ann);
    log_ta_tensor_host_memory_use();
  }

  [[nodiscard]] ResultPtr symmetrize() const override {
    return noted(
        eval_result<this_type>(detail::column_symmetrize_ta(get<ArrayT>())),
        inner_rank_);
  }

  [[nodiscard]] ResultPtr antisymmetrize(size_t /*bra_rank*/) const override {
    // not implemented yet
    return nullptr;
  }

 private:
  [[nodiscard]] std::size_t size_in_bytes() const final {
    if (Logger::instance().eval.level == 0)
      return 0;  // size_of disabled untraced: it is collective
    return measured_size_in_bytes();
  }

  /// the logical size of the value (see ResultTensorTA)
  [[nodiscard]] std::size_t measured_size_in_bytes() const final {
    auto const& v = view_base();
    auto local_size = TA::size_of<TA::MemorySpace::Host>(v);
    v.world().gop.sum(local_size);
    if (has_value() || view_.blk_lo.empty()) return local_size;
    auto const full = v.trange().elements_range().volume();
    if (full == 0) return 0;
    return static_cast<std::size_t>(
        static_cast<double>(local_size) *
        (static_cast<double>(logical_volume()) / static_cast<double>(full)));
  }
};

/// \brief Restrict a TA::DistArray to a contiguous tile range of one mode.
///
/// Keeps tiles `[tile_lo, tile_hi)` of mode \p mode and all tiles of every
/// other mode; the result's `mode`-th TiledRange1 covers only those tiles.
/// Every mode's original element lobound is preserved (via TA's
/// `preserve_lobound` block), so a spectator index that carries a nonzero
/// lobound (e.g. an active-occupied index with a frozen-core offset) keeps it
/// and still matches the unsliced operand it is contracted against. Implemented
/// with TA's `block()`, so block-sparse shape is preserved. Used to evaluate a
/// tensor
/// network in batches over a contracted index (see
/// make_batched_custom_evaluator): slicing every leaf that carries the index
/// to a tile range, evaluating, and summing reproduces the full contraction
/// because `sum_K = sum_{blocks} sum_{K in block}`.
///
/// \param inner_rank_hint for a nested \p arr, its inner rank (0 for an
///        all-empty-inner array) when the caller knows it (see
///        ResultTensorOfTensorTA::inner_rank()); unset, it is read from the
///        tiles with detail::tot_inner_rank(), a world-collective reduction.
template <typename... Args>
[[nodiscard]] TA::DistArray<Args...> slice_array_over_mode(
    TA::DistArray<Args...> const& arr, std::size_t mode, std::size_t tile_lo,
    std::size_t tile_hi, std::optional<std::size_t> inner_rank_hint) {
  using ranges::views::iota;
  auto const rank = arr.trange().rank();
  SEQUANT_ASSERT(mode < rank);
  SEQUANT_ASSERT(tile_lo < tile_hi &&
                 tile_hi <= arr.trange().dim(mode).tile_extent());
  container::svector<std::size_t> lo(rank, 0), hi(rank);
  for (std::size_t d = 0; d < rank; ++d)
    hi[d] = arr.trange().dim(d).tile_extent();
  lo[mode] = tile_lo;
  hi[mode] = tile_hi;
  // For a tensor-of-tensor array the annotation must label an inner block
  // ("outer;inner"); a flat annotation trips DistArray's is_tot_index() check.
  // The block() is over outer modes only, so both sides share one annotation.
  using value_type = typename TA::DistArray<Args...>::value_type;
  std::string annot;
  if constexpr (TA::detail::is_tensor_of_tensor_v<value_type>) {
    auto const inner_rank =
        inner_rank_hint ? *inner_rank_hint : detail::tot_inner_rank(arr);
    if (inner_rank == 0) {
      // All-empty-inner ToT: tot_inner_rank() is 0, so there's no inner rank to
      // build the ToT annotation block() needs. The slice of zero is zero, so
      // build the sliced outer trange by hand (mode cut to [tile_lo,tile_hi),
      // every dim's original element lobound preserved, matching the
      // preserve_lobound block() below) and return a zero ToT over it. The
      // result is dense even if ArrayT's policy is sparse -- harmless, the
      // array is identically zero.
      std::vector<TA::TiledRange1> tr1s;
      tr1s.reserve(rank);
      for (std::size_t d = 0; d < rank; ++d) {
        auto const& dim = arr.trange().dim(d);
        std::vector<std::size_t> bounds{
            static_cast<std::size_t>(dim.tile(lo[d]).first)};
        for (std::size_t t = lo[d]; t < hi[d]; ++t)
          bounds.push_back(static_cast<std::size_t>(dim.tile(t).second));
        tr1s.emplace_back(bounds.begin(), bounds.end());
      }
      TA::DistArray<Args...> out{arr.world(),
                                 TA::TiledRange(tr1s.begin(), tr1s.end())};
      for (auto it = out.begin(); it != out.end(); ++it)
        if (out.is_local(it.index())) *it = value_type{it.make_range()};
      out.world().gop.fence();
      ::sequant::detail::note_fence();
      return out;
    }
    annot = TA::detail::dummy_annotation(static_cast<unsigned int>(rank),
                                         static_cast<unsigned int>(inner_rank));
  } else {
    annot = detail::ords_to_annot(iota(std::size_t{0}, rank));
  }
  TA::DistArray<Args...> out;
  // preserve_lobound: keep every mode's original element lobound instead of
  // rebasing the block to 0. We only sub-range the batch mode; the other modes
  // are taken in full, but plain block() would still rebase them to 0 -- giving
  // a spectator index like the active-occupied i (whose DistArray carries an
  // ncore lobound) a 0 lobound in this sliced operand, which then mismatches
  // the unsliced operand's TiledRange1 for the same index when the two are
  // combined by einsum's general product (Einsum::index::operator| asserts that
  // a shared index has the same TiledRange1 in both operands). The batch mode
  // keeps its real element offset too, consistently across all sliced operands.
  auto const t0 = std::chrono::steady_clock::now();
  out(annot) = arr(annot).block(lo, hi, TA::preserve_lobound);
  ::sequant::detail::drain_lazy_cleanup(arr);
  detail::log_batch_op("Slice", std::chrono::steady_clock::now() - t0, out,
                       annot + " mode=" + std::to_string(mode) + " tiles=[" +
                           std::to_string(tile_lo) + "," +
                           std::to_string(tile_hi) + ")");
  return out;
}

/// \brief Scatter a per-block DistArray into a contiguous tile range of one
///        mode of a pre-sized destination -- the inverse of
///        slice_array_over_mode().
///
/// Writes \p block into tiles `[tile_lo, tile_hi)` of \p dest's mode \p mode,
/// leaving every other tile of \p dest untouched. \p dest must already be
/// allocated over its full TiledRange (the caller sizes the whole shape), and
/// \p block's TiledRange must equal \p dest's sub-block over `[tile_lo,
/// tile_hi)` (as produced by slice_array_over_mode() for the same mode/range).
/// Implemented with TA's block() on the assignment LHS, so only the addressed
/// sub-block is written and block-sparse shape is preserved. Every mode's
/// element lobound is preserved (via TA's `preserve_lobound`), exactly as the
/// lobound-preserving gather in slice_array_over_mode(): the destination
/// sub-block and the source share element coordinates, so a spectator index
/// carrying a nonzero lobound (e.g. a frozen-core offset) lands at its true
/// offset rather than being rebased to 0. Reconstructs a whole result from a
/// disjoint, gap-free tiling of one mode: scattering each block of a partition
/// reproduces the array `slice_array_over_mode()` would gather back out.
///
/// \param inner_rank_hint for a nested \p block, its inner rank (0 for an
///        all-empty-inner block) when the caller knows it; unset, it is read
///        from the tiles with detail::tot_inner_rank() (a world-collective
///        reduction), as in slice_array_over_mode().
template <typename... Args>
void write_array_into_mode(TA::DistArray<Args...>& dest,
                           TA::DistArray<Args...> const& block,
                           std::size_t mode, std::size_t tile_lo,
                           std::size_t tile_hi,
                           std::optional<std::size_t> inner_rank_hint) {
  using ranges::views::iota;
  auto const rank = dest.trange().rank();
  SEQUANT_ASSERT(mode < rank);
  SEQUANT_ASSERT(tile_lo < tile_hi &&
                 tile_hi <= dest.trange().dim(mode).tile_extent());
  container::svector<std::size_t> lo(rank, 0), hi(rank);
  for (std::size_t d = 0; d < rank; ++d)
    hi[d] = dest.trange().dim(d).tile_extent();
  lo[mode] = tile_lo;
  hi[mode] = tile_hi;
  // For a tensor-of-tensor array the annotation must label an inner block
  // ("outer;inner"); a flat annotation trips DistArray's is_tot_index() check.
  // The block() is over outer modes only, so both sides share one annotation.
  using value_type = typename TA::DistArray<Args...>::value_type;
  std::string annot;
  if constexpr (TA::detail::is_tensor_of_tensor_v<value_type>) {
    auto const inner_rank =
        inner_rank_hint ? *inner_rank_hint : detail::tot_inner_rank(block);
    if (inner_rank == 0) {
      // block has all-empty inner tiles (tot_inner_rank() == 0): it represents
      // zero and there is no inner rank to form the ToT annotation block()
      // needs. A zero contribution leaves the pre-sized destination slice as
      // it was, so skip the scatter entirely -- mirroring the zero-ToT early
      // return in slice_array_over_mode().
      return;
    }
    annot = TA::detail::dummy_annotation(static_cast<unsigned int>(rank),
                                         static_cast<unsigned int>(inner_rank));
  } else {
    annot = detail::ords_to_annot(iota(std::size_t{0}, rank));
  }
  // preserve_lobound: address the destination sub-block in its original element
  // coordinates (keeping every mode's lobound) so it matches the source block,
  // which slice_array_over_mode() also gathered with preserve_lobound. Plain
  // block() would rebase the sub-block to 0 and mismatch the source trange.
  auto const t0 = std::chrono::steady_clock::now();
  dest(annot).block(lo, hi, TA::preserve_lobound) = block(annot);
  ::sequant::detail::drain_lazy_cleanup(dest);
  detail::log_batch_op("Scatter", std::chrono::steady_clock::now() - t0, block,
                       annot + " mode=" + std::to_string(mode) + " tiles=[" +
                           std::to_string(tile_lo) + "," +
                           std::to_string(tile_hi) + ")");
}

/// \brief Compute the result's OUTER TiledRange for a binary product from the
///        type-erased operands and the [left, right, result] annotations.
///
/// Dispatches on each operand's concrete result kind (flat \c ResultTensorTA vs
/// nested \c ResultTensorOfTensorTA) to read its array's TiledRange, then maps
/// the result's outer labels back onto the operands (see
/// detail::result_outer_trange). This is the trange handed to the result-shape
/// provider. \tparam NumericT, \tparam PolicyT name the concrete array types.
/// \tparam InnerTileT the inner tile type of a nested (ToT) operand; defaults
/// to the plain \c TA::Tensor<NumericT>, but the CSV/PNO path uses an
/// arena-pinned inner tile (\c TA::ArenaTensor<NumericT>), which is a distinct
/// (type-id) Result kind and must be named here for \c Result::is<> to
/// recognize it.
template <typename NumericT, typename PolicyT,
          typename InnerTileT = TA::Tensor<NumericT>>
[[nodiscard]] TA::TiledRange result_outer_trange_from_results(
    Result const& left, Result const& right,
    std::array<std::any, 3> const& annot) {
  using FlatArray = TA::DistArray<TA::Tensor<NumericT>, PolicyT>;
  using ToTArray = TA::DistArray<TA::Tensor<InnerTileT>, PolicyT>;
  using ToTResult = ResultTensorOfTensorTA<ToTArray>;

  using FlatResult = ResultTensorTA<FlatArray>;

  auto const a = Annot<std::string>{annot};

  // Each operand's TiledRange by kind: a flat result reports the range of the
  // value it stands for (a pending view is not materialized for it), a nested
  // one its array's. Then defer to the label-matching detail helper.
  auto trange_of = [&](Result const& r) -> TA::TiledRange {
    if (r.is<ToTResult>()) return r.as<ToTResult>().logical_trange();
    return r.as<FlatResult>().logical_trange();
  };
  return detail::result_outer_trange(trange_of(left), a.lannot,
                                     trange_of(right), a.rannot, a.this_annot);
}

/// \brief Emit a binary product through the STANDARD expression layer with an
///        imposed result SparseShape, instead of the einsum path used by
///        Result::prod().
///
/// This is the application side of the result-shape-constraint feature: when a
/// method-supplied provider returns a SparseShape for a Product node, the eval
/// emits the product as `out(ca) = (lhs(la) * rhs(ra)).set_shape(s)` (general
/// product) or `out(ca) = lhs(la).dot_inner(rhs(ra)).set_shape(s)` (the
/// DeNest::True ToT*ToT->flat path), so the result is computed/stored only over
/// the kept region.
///
/// Operand/result nesting determines which `ResultTensor*TA` wraps the output:
///   - both operands flat (T * T)           -> flat result
///   - exactly one ToT (T * ToT / ToT * T)  -> ToT result
///   - both ToT, \p de_nest false           -> ToT result (general product)
///   - both ToT, \p de_nest true            -> flat result (dot_inner)
///
/// \param left,right  The operands (TA tensor or ToT results).
/// \param annot       [left, right, result] annotations.
/// \param shape       The SparseShape to impose on the result's outer modes.
///                    Held by the local expression by pointer; this function
///                    keeps it alive until the assignment completes.
/// \param de_nest     The DeNest flag computed at the binary-product site
///                    (left.tot && right.tot && !result.tot).
/// \return The shaped product wrapped in the appropriate ResultTensor*TA, or a
///         null ResultPtr if the operand kinds are not a supported shaped form
///         (the caller then falls through to the unshaped prod()).
template <typename NumericT, typename PolicyT,
          typename InnerTileT = TA::Tensor<NumericT>>
[[nodiscard]] ResultPtr apply_shaped_product(
    Result const& left, Result const& right,
    std::array<std::any, 3> const& annot, TA::SparseShape<float> const& shape,
    bool de_nest) {
  using FlatArray = TA::DistArray<TA::Tensor<NumericT>, PolicyT>;
  using ToTArray = TA::DistArray<TA::Tensor<InnerTileT>, PolicyT>;
  using FlatResult = ResultTensorTA<FlatArray>;
  using ToTResult = ResultTensorOfTensorTA<ToTArray>;

  auto const a = Annot<std::string>{annot};

  bool const l_tot = left.is<ToTResult>();
  bool const r_tot = right.is<ToTResult>();
  bool const l_flat = left.is<FlatResult>();
  bool const r_flat = right.is<FlatResult>();

  // No world fence after the assignments: the expression engine folds the
  // imposed shape into its own shape (a mask) when the expression is
  // initialized, inside the assignment, so the caller's shape object need not
  // outlive it; a per-product fence would only serialize the asynchronous
  // tile work against the executor's bookkeeping.

  // A flat operand may be a pending view (a batch slice, a permutation, a
  // conjugation): with_expr() folds it into the expression, so it is not
  // materialized here.

  // T * T -> T
  if (l_flat && r_flat) {
    FlatArray out;
    left.as<FlatResult>().with_expr(a.lannot, [&](auto&& le) {
      right.as<FlatResult>().with_expr(a.rannot, [&](auto&& re) {
        out(a.this_annot) = (le * re).set_shape(shape);
      });
    });
    ::sequant::detail::drain_lazy_cleanup(out);
    return eval_result<FlatResult>(std::move(out));
  }

  // T * ToT -> ToT  and  ToT * T -> ToT (general product)
  if ((l_flat && r_tot) || (l_tot && r_flat)) {
    // Same inner-reorder limitation as the ToT*ToT branch: decline if the
    // result's inner annotation needs a non-identity permutation TA can't emit.
    if (detail::tot_product_needs_inner_reorder(a.lannot, a.rannot,
                                                a.this_annot))
      return nullptr;
    ToTArray out;
    if (l_flat) {
      left.as<FlatResult>().with_expr(a.lannot, [&](auto&& le) {
        right.as<ToTResult>().with_expr(a.rannot, [&](auto&& re) {
          out(a.this_annot) = (le * re).set_shape(shape);
        });
      });
    } else {
      left.as<ToTResult>().with_expr(a.lannot, [&](auto&& le) {
        right.as<FlatResult>().with_expr(a.rannot, [&](auto&& re) {
          out(a.this_annot) = (le * re).set_shape(shape);
        });
      });
    }
    ::sequant::detail::drain_lazy_cleanup(out);
    return eval_result<ToTResult>(std::move(out));
  }

  // ToT * ToT
  if (l_tot && r_tot) {
    if (de_nest) {
      // ToT * ToT -> flat T (inner contraction / denest)
      FlatArray out;
      left.as<ToTResult>().with_expr(a.lannot, [&](auto&& le) {
        right.as<ToTResult>().with_expr(a.rannot, [&](auto&& re) {
          out(a.this_annot) = le.dot_inner(re).set_shape(shape);
        });
      });
      ::sequant::detail::drain_lazy_cleanup(out);
      return eval_result<FlatResult>(std::move(out));
    } else {
      // ToT * ToT -> ToT (general product). TiledArray's expression-layer
      // general ToT product cannot emit a result whose INNER annotation needs a
      // non-identity permutation (cont_engine throws). When the result requires
      // such an inner reorder, decline so the caller falls through to the
      // unshaped einsum prod() (which handles it); the shape is simply not
      // imposed on this node -- lossless, just less aggressive.
      if (detail::tot_product_needs_inner_reorder(a.lannot, a.rannot,
                                                  a.this_annot))
        return nullptr;
      ToTArray out;
      left.as<ToTResult>().with_expr(a.lannot, [&](auto&& le) {
        right.as<ToTResult>().with_expr(a.rannot, [&](auto&& re) {
          out(a.this_annot) = (le * re).set_shape(shape);
        });
      });
      ::sequant::detail::drain_lazy_cleanup(out);
      return eval_result<ToTResult>(std::move(out));
    }
  }

  // Unsupported operand kinds for a shaped product (e.g. a scalar operand):
  // decline so the caller falls through to the unshaped prod().
  return nullptr;
}

}  // namespace sequant

#endif  // SEQUANT_HAS_TILEDARRAY

#endif  // SEQUANT_EVAL_BACKENDS_TILEDARRAY_RESULT_HPP
