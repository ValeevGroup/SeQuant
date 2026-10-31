#ifndef SEQUANT_EVAL_BACKENDS_TILEDARRAY_RESULT_HPP
#define SEQUANT_EVAL_BACKENDS_TILEDARRAY_RESULT_HPP

#ifdef SEQUANT_HAS_TILEDARRAY

#include <sstream>

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/eval/elementwise.hpp>
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
#include <utility>

namespace sequant {

// implementation details of the TiledArray result backend; prefer
// sequant::detail over an unnamed namespace in a header (see CppCoreGuidelines
// SF.21 / "Use unnamed namespaces in headers ... no" guidance)
namespace detail {

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
inline void note_wait() {
  wait_counter().fetch_add(1, std::memory_order_relaxed);
}
/// waits for the deferred cleanup of \p arr's world (see
/// TA::DistArray::wait_for_lazy_cleanup) and records the wait
template <typename Array>
void wait_for_lazy_cleanup(Array const& arr) {
  Array::wait_for_lazy_cleanup(arr.world());
  note_wait();
}
inline std::atomic<std::size_t>& slice_counter() {
  static std::atomic<std::size_t> c{0};
  return c;
}
inline void note_slice() {
  slice_counter().fetch_add(1, std::memory_order_relaxed);
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
/// is_tot_index() check.
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

/// \return a zero \p ArrayT over \p trange. Under a sparse policy with a
/// positive SparseShape threshold its shape is all zero, so it holds no tiles;
/// otherwise every tile is set, zero-filled for a regular array and holding
/// empty inner tensors (so tot_inner_rank() == 0) for a tensor of tensors.
template <typename ArrayT>
[[nodiscard]] ArrayT make_zero_array(TA::World& world,
                                     TA::TiledRange const& trange) {
  using shape_type = typename ArrayT::shape_type;
  if constexpr (!shape_type::is_dense()) {
    // A tile is zero only if its norm is below the threshold, so under a zero
    // threshold even a zero norm marks a tile nonzero: the shape would claim
    // tiles that are never set.
    if (shape_type::threshold() > 0)
      return ArrayT(world, trange,
                    shape_type(typename shape_type::value_type(0), trange));
  }
  using value_type = typename ArrayT::value_type;
  ArrayT result(world, trange);
  if constexpr (TA::detail::is_tensor_of_tensor_v<value_type>) {
    for (auto it = result.begin(); it != result.end(); ++it)
      if (result.is_local(it.index())) *it = value_type{it.make_range()};
  } else {
    result.fill_local(typename ArrayT::numeric_type(0));
  }
  world.gop.fence();
  note_fence();
  return result;
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

  ::sequant::detail::wait_for_lazy_cleanup(result);

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

  ::sequant::detail::wait_for_lazy_cleanup(result);
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
inline void log_batch_op([[maybe_unused]] char const* kind,
                         [[maybe_unused]] std::chrono::nanoseconds elapsed,
                         [[maybe_unused]] TA::DistArray<Args...> const& moved,
                         [[maybe_unused]] std::string const& annot) noexcept {
#ifdef SEQUANT_EVAL_TRACE
  auto& l = Logger::instance();
  if (l.eval.level == 0) return;
  auto bytes = TA::size_of<TA::MemorySpace::Host>(moved);
  moved.world().gop.sum(bytes);
  write_log(l, "Batch | ", kind, " | ", elapsed.count(), "ns | bytes=", bytes,
            "B | ", annot, '\n');
#endif
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

/// \brief A pending, zero-copy transformation of a stored array: the view a
///        ResultTensorTA may stand for instead of holding a value.
///
/// Composable parts, all optional: a tile block (a batch slice), a mode
/// permutation, a conjugation. The block and the permutation are expressed in
/// the STORED array's mode order; the result's own (logical) modes are the
/// stored ones permuted by \c perm. A view is folded into a TiledArray
/// expression by ResultTensorTA::with_expr(): `arr(stored annot).block(blk_lo,
/// blk_hi, preserve_lobound).conj()`, so a consumer that is itself an
/// expression (a product, a sum, a scatter) never materializes it.
struct ArrayView {
  /// stored mode \c s holds logical mode \c perm[s]; empty = identity
  container::svector<std::size_t> perm;
  /// tile bounds `[blk_lo[s], blk_hi[s])` per STORED mode, in the stored
  /// array's tile coordinates; empty = the whole array
  container::svector<std::size_t> blk_lo, blk_hi;
  /// elementwise conjugation pending
  bool conj = false;

  [[nodiscard]] bool trivial() const noexcept {
    return perm.empty() && blk_lo.empty() && !conj;
  }

  /// the stored mode that holds logical mode \p logical
  [[nodiscard]] std::size_t stored_mode(std::size_t logical) const {
    if (perm.empty()) return logical;
    for (std::size_t s = 0; s < perm.size(); ++s)
      if (perm[s] == logical) return s;
    SEQUANT_ASSERT(false && "ArrayView::stored_mode: no such logical mode");
    return logical;
  }
};

/// the annotation of the STORED array for the logical annotation \p annot of a
/// view permuted by \p perm (see ArrayView::perm): stored mode \c s is
/// labelled with the logical label of mode \c perm[s]. Flat annotations only.
[[nodiscard]] inline std::string stored_annot(
    std::string const& annot, container::svector<std::size_t> const& perm) {
  if (perm.empty()) return annot;
  auto const labels = outer_annot_labels(annot);
  SEQUANT_ASSERT(labels.size() == perm.size());
  std::string out;
  for (std::size_t s = 0; s < perm.size(); ++s) {
    if (s) out += ',';
    out += labels[perm[s]];
  }
  return out;
}

/// the TiledRange a view stands for: \p stored's dims cut to the view's tile
/// block (every mode keeps its element coordinates, as TA's
/// `preserve_lobound` block does) and reordered to the view's logical modes
[[nodiscard]] inline TA::TiledRange logical_trange(TA::TiledRange const& stored,
                                                   ArrayView const& view) {
  auto const rank = stored.rank();
  std::vector<TA::TiledRange1> dims;
  dims.reserve(rank);
  for (std::size_t m = 0; m < rank; ++m) {
    auto const s = view.stored_mode(m);
    auto const& dim = stored.dim(s);
    if (view.blk_lo.empty()) {
      dims.push_back(dim);
      continue;
    }
    std::vector<std::size_t> bounds{
        static_cast<std::size_t>(dim.tile(view.blk_lo[s]).first)};
    for (std::size_t t = view.blk_lo[s]; t < view.blk_hi[s]; ++t)
      bounds.push_back(static_cast<std::size_t>(dim.tile(t).second));
    dims.emplace_back(bounds.begin(), bounds.end());
  }
  return TA::TiledRange(dims.begin(), dims.end());
}

/// \return the TiledRange \p tr, as printed by TiledArray
[[nodiscard]] inline std::string ta_trange_string(TA::TiledRange const& tr) {
  std::ostringstream os;
  os << tr;
  return os.str();
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

/// \return a new array holding every element of \p arr raised to \p exponent
///         (see Result::pow); for a nested array, every inner element
/// \note The tiles are mapped into a fresh array by TA::foreach, so \p arr,
///       which may be shared (e.g. a cached leaf), is not modified. A sparse
///       array's zero tiles stay zero, which is only right for a positive
///       exponent; other exponents require every tile to be present. The same
///       holds for the empty inner tiles of a nested array.
template <typename ArrayT>
[[nodiscard]] ArrayT ta_elementwise_pow(ArrayT const& arr,
                                        rational const& exponent) {
  using tile_type = typename ArrayT::value_type;
  using numeric_type = typename ArrayT::numeric_type;
  if constexpr (TA::detail::is_tensor_of_tensor_v<tile_type>) {
    // checked here rather than in the tile tasks: an exception thrown inside
    // a MADNESS task aborts the process
    if (!(exponent > 0)) {
      int absent = 0;
      for (auto it = arr.begin(); it != arr.end() && absent == 0; ++it) {
        tile_type const tile = it->get();
        for (std::size_t k = 0; k < tile.size(); ++k)
          if (tile.data()[k].empty()) {
            absent = 1;
            break;
          }
      }
      arr.world().gop.max(absent);
      SEQUANT_ENFORCE(absent == 0,
                      "element-wise power with a non-positive exponent of a "
                      "nested array with an absent inner tile: its value "
                      "(inf or 1) is not representable");
    }
  }
  ElementwisePow const op{exponent};
  auto const elem = [op](numeric_type const& x) { return op(x); };
  auto const pow_tile = [elem](tile_type const& in) -> tile_type {
    if constexpr (TA::detail::is_tensor_of_tensor_v<tile_type>) {
      tile_type out(in.range());
      for (std::size_t k = 0; k < in.size(); ++k)
        if (!in.data()[k].empty()) out.data()[k] = in.data()[k].unary(elem);
      return out;
    } else {
      return in.unary(elem);
    }
  };
  if constexpr (TA::is_dense_v<typename ArrayT::policy_type>) {
    return TA::foreach (arr, [pow_tile](tile_type& out, tile_type const& in) {
      out = pow_tile(in);
    });
  } else {
    SEQUANT_ENFORCE(exponent > 0 || arr.shape().sparsity() == 0,
                    "element-wise power with a non-positive exponent of a "
                    "sparse array with zero tiles");
    return TA::foreach (
        arr, [pow_tile](tile_type& out, tile_type const& in) -> float {
          out = pow_tile(in);
          return out.norm();
        });
  }
}

/// \return \p arr expanded over the modes of \p result_annot that \p
///         operand_annot lacks (see Result::broadcast), or, for a null \p arr,
///         an array every element of which is \p scalar
/// \param world the World of the result
/// \param new_dims the TiledRange1 of each mode of \p result_annot absent from
///        \p operand_annot, in the result's mode order; a mode \p arr carries
///        keeps \p arr's TiledRange1
/// \note Each result tile is computed by a task that waits on the operand tile
///       it repeats (the result's tile index without the added modes), so no
///       fence is needed. For a sparse policy, a result tile is zero iff that
///       operand tile is: its norm is the operand's per-element norm spread
///       over the added modes, raised to the shape threshold if it would fall
///       below it.
template <typename ArrayT>
[[nodiscard]] ArrayT ta_broadcast(
    TA::World& world, ArrayT const* arr,
    typename ArrayT::numeric_type const& scalar,
    std::string const& operand_annot, std::string const& result_annot,
    container::svector<TA::TiledRange1> const& new_dims) {
  using tile_type = typename ArrayT::value_type;
  auto const pos = broadcast_operand_modes(outer_annot_labels(operand_annot),
                                           outer_annot_labels(result_annot));
  auto const rank = pos.size();
  std::size_t const operand_rank = static_cast<std::size_t>(
      std::count_if(pos.begin(), pos.end(), [](auto const& q) { return q; }));
  SEQUANT_ASSERT(arr || operand_rank == 0);
  SEQUANT_ASSERT(new_dims.size() == rank - operand_rank);

  std::vector<TA::TiledRange1> dims;
  dims.reserve(rank);
  for (std::size_t p = 0, k = 0; p < rank; ++p)
    dims.push_back(pos[p] ? arr->trange().dim(*pos[p]) : new_dims[k++]);
  TA::TiledRange const trange(dims.begin(), dims.end());

  // the operand tile a result tile repeats
  auto const operand_tile = [&pos, operand_rank](auto const& result_tile) {
    container::svector<std::size_t> out(operand_rank);
    for (std::size_t p = 0; p < pos.size(); ++p)
      if (pos[p]) out[*pos[p]] = result_tile[p];
    return out;
  };

  auto const make_result = [&]() -> ArrayT {
    using policy_type = typename ArrayT::policy_type;
    if constexpr (TA::is_dense_v<policy_type>) {
      return ArrayT(world, trange);
    } else {
      using shape_type = typename ArrayT::shape_type;
      using norm_type = typename shape_type::value_type;
      norm_type const threshold = shape_type::threshold();
      auto const norm_of = [&](TA::Range::index_type const& tile) -> norm_type {
        auto const volume = trange.make_tile_range(tile).volume();
        if (volume == 0) return 0;
        norm_type per_element;
        if (!arr) {
          if (scalar == typename ArrayT::numeric_type{0}) return 0;
          per_element = static_cast<norm_type>(std::abs(scalar) /
                                               std::sqrt(double(volume)));
        } else {
          auto const op_tile = operand_tile(tile);
          if (arr->is_zero(op_tile)) return 0;
          auto const op_volume =
              arr->trange().make_tile_range(op_tile).volume();
          // the per-element norm of op_tile repeated volume/op_volume times
          per_element = static_cast<norm_type>(
              arr->shape()[op_tile] /
              std::sqrt(double(volume) / double(op_volume)));
        }
        return std::max(per_element, threshold);
      };
      return ArrayT(world, trange,
                    shape_type(norm_of, trange, /*do_not_scale=*/true));
    }
  };
  ArrayT result = make_result();

  if (!arr) {
    result.init_tiles(
        [scalar](TA::Range const& r) { return tile_type(r, scalar); });
    return result;
  }
  auto const kernel = [pos](TA::Range const& r,
                            tile_type const& in) -> tile_type {
    tile_type out(r);
    auto const n = pos.size();
    container::svector<std::size_t> extents(n), out_strides(n), in_strides(n);
    for (std::size_t p = 0; p < n; ++p) {
      extents[p] = r.extent()[p];
      out_strides[p] = r.stride()[p];
      in_strides[p] = pos[p] ? in.range().stride()[*pos[p]] : 0;
    }
    broadcast_copy(out.data(), extents, out_strides, in.data(), in_strides);
    return out;
  };
  for (auto const ord : *result.pmap()) {
    if (result.is_zero(ord)) continue;
    auto const tile = trange.tiles_range().idx(ord);
    result.set(ord, world.taskq.add(kernel, trange.make_tile_range(ord),
                                    arr->find(operand_tile(tile))));
  }
  return result;
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

// Eager gather of a block (defined below): a utility, no longer used by the
// result classes, whose slices are pending views.
template <typename... Args>
[[nodiscard]] TA::DistArray<Args...> slice_array_over_mode(
    TA::DistArray<Args...> const& arr, std::size_t mode, std::size_t tile_lo,
    std::size_t tile_hi);

// Eager scatter of a block (defined below): the inverse of
// slice_array_over_mode(); a utility, no longer used by the result classes.
template <typename... Args>
void write_array_into_mode(TA::DistArray<Args...>& dest,
                           TA::DistArray<Args...> const& block,
                           std::size_t mode, std::size_t tile_lo,
                           std::size_t tile_hi);

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
/// A result may be a VIEW of another result's array -- a tile block of it (a
/// batch slice), a mode permutation, a conjugation, or any combination --
/// pending rather than applied (detail::ArrayView). The view shares the
/// parent's tiles, so it costs no allocation and no pass over the data
/// (TiledArray copies the tiles of a block of a MUTABLE array; with_expr()
/// hands the parent over as a const array, whose preserve_lobound block is a
/// view there):
/// consumers that are TiledArray expressions -- the flat products without a
/// Hadamard index, the dot product, scaling, sums, accumulation of the operand,
/// the batch scatter -- fold it into their expression through with_expr(); a
/// consumer that needs an array (TiledArray's einsum, clone, the symmetrizers,
/// Result::get<ArrayT>()) materializes it once, into this result's own value,
/// after which the parent is released.
///
template <typename ArrayT, typename = std::enable_if_t<TA::detail::is_tensor_v<
                               typename ArrayT::value_type>>>
class ResultTensorTA final : public Result {
 public:
  using Result::id_t;
  using numeric_type = typename ArrayT::numeric_type;
  using view_type = detail::ArrayView;

  explicit ResultTensorTA(ArrayT arr) : Result{std::move(arr)} {}

  /// a view of \p parent with \p view pending. A trivial view (an identity
  /// permutation, a slice over whole modes) is an alias too: it reads the
  /// parent's tiles as they are and still copies before it is written to, so
  /// no two results ever share tiles one of them accumulates into.
  ResultTensorTA(ArrayT parent, view_type view)
      : Result{}, parent_(std::move(parent)), view_(std::move(view)) {}

  [[nodiscard]] std::string trange_annot() const override {
    return detail::ta_trange_string(logical_trange());
  }

  /// true while this result is a pending view sharing another array's tiles
  [[nodiscard]] bool is_buffer_alias() const noexcept override {
    return !has_value();
  }

  /// the TiledRange of the value this result stands for, with the pending
  /// view applied (block cut, modes permuted); never materializes
  [[nodiscard]] TA::TiledRange logical_trange() const {
    if (has_value()) return Result::get<ArrayT>().trange();
    return detail::logical_trange(parent_.trange(), view_);
  }

  /// Hands \p f the TiledArray expression of this result's value for the
  /// logical annotation \p annot: `value(annot)` when the value is held, else
  /// the parent array with the pending view folded in -- its modes relabelled
  /// through the permutation, the block cut by block(), the conjugation by
  /// conj(). \p f must return the same type whatever expression it receives.
  ///
  /// The parent is handed to TiledArray as a CONST array. A block of a const
  /// array that keeps the parent's element coordinates (preserve_lobound) is
  /// a view there: the consuming expression reads the parent's own tiles,
  /// whereas a block of a mutable array is copied tile by tile. The const is
  /// the promise that nothing writes the parent's tiles while the expression
  /// is alive, which holds: a view is never written through, and the parent
  /// is released once this result holds a value of its own.
  template <typename F>
  decltype(auto) with_expr(std::string const& annot, F&& f) const {
    if (has_value()) return std::forward<F>(f)(Result::get<ArrayT>()(annot));
    auto const& parent = std::as_const(parent_);
    auto const stored = detail::stored_annot(annot, view_.perm);
    bool const conj = TA::detail::is_complex_v<numeric_type> && view_.conj;
    if (view_.blk_lo.empty()) {
      if constexpr (TA::detail::is_complex_v<numeric_type>) {
        if (conj) return std::forward<F>(f)(parent(stored).conj());
      }
      return std::forward<F>(f)(parent(stored));
    }
    // The block keeps every mode's element coordinates (TA's
    // preserve_lobound), like the eager slice did, so a spectator mode with a
    // nonzero lobound still matches the operand it is contracted against;
    // the expressions derived from the block (conj(), scaling, negation)
    // carry the lobound along.
    auto blk =
        parent(stored).block(view_.blk_lo, view_.blk_hi, TA::preserve_lobound);
    if constexpr (TA::detail::is_complex_v<numeric_type>) {
      if (conj) return std::forward<F>(f)(blk.conj());
    }
    return std::forward<F>(f)(blk);
  }

 private:
  using this_type = ResultTensorTA<ArrayT>;
  using annot_wrap = Annot<std::string>;

  /// the parent array of a pending view: an empty handle once the value is
  /// held (mutable: released by the const materialization in ensure_value())
  mutable ArrayT parent_;
  view_type view_;

  [[nodiscard]] id_t type_id() const noexcept override {
    return id_for_type<this_type>();
  }

  [[nodiscard]] std::size_t logical_rank() const {
    return has_value() ? Result::get<ArrayT>().trange().rank()
                       : parent_.trange().rank();
  }

  /// a dummy annotation over the logical modes
  [[nodiscard]] std::string dummy_annot() const {
    return TA::detail::dummy_annotation(
        static_cast<unsigned int>(logical_rank()));
  }

  /// the value this result stands for, as a fresh array: the pending view
  /// applied in one pass (or a deep copy of the held value)
  [[nodiscard]] ArrayT materialize() const {
    if (has_value()) return TA::clone(Result::get<ArrayT>());
    auto const ann = dummy_annot();
    auto const t0 = std::chrono::steady_clock::now();
    ArrayT out;
    with_expr(ann, [&](auto&& e) { out(ann) = e; });
    ::sequant::detail::wait_for_lazy_cleanup(out);
    detail::log_batch_op("Materialize", std::chrono::steady_clock::now() - t0,
                         out, ann);
    return out;
  }

  /// Result::get<ArrayT>() on a pending view: the view applied, once; the
  /// parent's tiles are released
  void ensure_value() const override {
    if (has_value()) return;
    set_value(std::make_any<ArrayT>(materialize()));
    parent_ = ArrayT();
  }

  /// the value this result stands for, without storing it: the held value
  /// (shallow), or the pending view applied into a temporary. For read-only
  /// consumers that must not turn a shared view into a held value.
  [[nodiscard]] ArrayT applied_value() const {
    return has_value() ? Result::get<ArrayT>() : materialize();
  }

  /// the array a further view of this result is a view of (the parent while
  /// this one is pending, else the held value), and that view's mode map
  [[nodiscard]] ArrayT const& view_base() const {
    return has_value() ? Result::get<ArrayT>() : parent_;
  }
  [[nodiscard]] view_type base_view() const {
    return has_value() ? view_type{} : view_;
  }

  // Diagnostic (see Result::fence): force all pending async work in this
  // array's world to complete so a lazily-executed op's timer captures
  // execution, not just dispatch.
  void fence() const noexcept override {
    view_base().world().gop.fence();
    ::sequant::detail::note_fence();
  }

  [[nodiscard]] ResultPtr sum(
      Result const& other,
      std::array<std::any, 3> const& annot) const override {
    SEQUANT_ASSERT(other.is<this_type>());
    auto const& o = other.as<this_type>();
    auto const a = annot_wrap{annot};

    detail::log_ta(a.lannot, " + ", a.rannot, " = ", a.this_annot, "\n");

    ArrayT result;
    with_expr(a.lannot, [&](auto&& le) {
      o.with_expr(a.rannot, [&](auto&& re) { result(a.this_annot) = le + re; });
    });
    ::sequant::detail::wait_for_lazy_cleanup(result);
    log_ta_tensor_host_memory_use();
    return eval_result<this_type>(std::move(result));
  }

  /// A batch slice is a pending view: the tile block of logical \p mode that
  /// covers `[elem_lo, elem_hi)` (tile-aligned, see slice_bounds_to_tiles),
  /// intersected with the block this result already is. No tile is copied.
  [[nodiscard]] ResultPtr slice_mode(std::size_t mode, std::size_t elem_lo,
                                     std::size_t elem_hi) const override {
    ::sequant::detail::note_slice();
    auto const t0 = std::chrono::steady_clock::now();
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
    auto r = eval_result<this_type>(base, std::move(view));
    ::sequant::detail::slice_ns().fetch_add(
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - t0)
            .count(),
        std::memory_order_relaxed);
    return r;
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
    return detail::ta_tile_diff(applied_value(),
                                other.as<this_type>().applied_value());
  }

  /// Scatter \p block into this result's value (materialized first if this
  /// result is a pending view: a scatter writes, and a view must not write
  /// through to its parent). The block may itself be a pending view: it is
  /// folded into the assignment.
  void write_into_slice(Result const& block, std::size_t mode,
                        std::size_t block_lo, std::size_t block_hi) override {
    SEQUANT_ASSERT(block.is<this_type>());
    auto const& b = block.as<this_type>();
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
                       detail::logical_trange(dest.trange(), {{}, lo, hi}) &&
                   "write_into_slice: the block's layout is not the "
                   "destination sub-block's");
    auto const annot = dummy_annot();
    // preserve_lobound on the destination sub-block: the source block shares
    // its element coordinates (see write_array_into_mode)
    auto const t0 = std::chrono::steady_clock::now();
    b.with_expr(annot, [&](auto&& be) {
      dest(annot).block(lo, hi, TA::preserve_lobound) = be;
    });
    ::sequant::detail::wait_for_lazy_cleanup(dest);
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

      ::sequant::detail::wait_for_lazy_cleanup(result);
      log_ta_tensor_host_memory_use();
      return eval_result<this_type>(std::move(result));
    }

    if (a.this_annot.empty()) {
      // DOT product
      SEQUANT_ASSERT(other.is<this_type>());
      auto const& o = other.as<this_type>();
      numeric_type const d = with_expr(a.lannot, [&](auto&& le) {
        return o.with_expr(a.rannot, [&](auto&& re) {
          return static_cast<numeric_type>(le.dot(re).get());
        });
      });
      ::sequant::detail::wait_for_lazy_cleanup(view_base());
      ::sequant::detail::wait_for_lazy_cleanup(o.view_base());

      detail::log_ta(a.lannot, " * ", a.rannot, " = ", d, "\n");

      log_ta_tensor_host_memory_use();
      return eval_result<ResultScalar<numeric_type>>(d);
    }

    if (!other.is<this_type>()) {
      // potential T * ToT
      auto annot_swap = annot;
      std::swap(annot_swap[0], annot_swap[1]);
      return other.prod(*this, annot_swap, DeNestFlag);
    }

    // confirmed: other.is<this_type>() is true
    auto const& o = other.as<this_type>();

    detail::log_ta(a.lannot, " * ", a.rannot, " = ", a.this_annot, "\n");

    // The product is a TiledArray expression -- the expression layer handles
    // contractions, Hadamard and general (fused + contracted + free) products
    // natively, which is where TiledArray's einsum routes them too -- so an
    // operand that is a pending view (block, permutation, conjugation) is
    // taken as is: nothing is materialized.
    ArrayT result;
    with_expr(a.lannot, [&](auto&& le) {
      o.with_expr(a.rannot, [&](auto&& re) { result(a.this_annot) = le * re; });
    });
    ::sequant::detail::wait_for_lazy_cleanup(result);
    log_ta_tensor_host_memory_use();
    return eval_result<this_type>(std::move(result));
  }

  /// A pending view is folded into the scaling, not materialized first.
  [[nodiscard]] ResultPtr mult_by_phase(std::int8_t factor) const override {
    auto const ann = dummy_annot();
    ArrayT result;
    with_expr(ann, [&](auto&& e) { result(ann) = numeric_type(factor) * e; });
    ::sequant::detail::wait_for_lazy_cleanup(result);
    return eval_result<this_type>(std::move(result));
  }

  /// Deep copy: \c TA::DistArray's own copy is a shallow (reference-counted)
  /// handle onto the same tiles, so an in-place accumulation into the copy
  /// would be seen by every other holder -- \c TA::clone allocates and copies
  /// the tiles. A pending view is applied into the copy directly.
  [[nodiscard]] ResultPtr clone() const override {
    return eval_result<this_type>(materialize());
  }

  /// A permutation is a pending view: the modes relabelled, no tile moved.
  [[nodiscard]] ResultPtr permute(
      std::array<std::any, 2> const& ann) const override {
    auto const pre_annot = std::any_cast<std::string>(ann[0]);
    auto const post_annot = std::any_cast<std::string>(ann[1]);

    detail::log_ta(pre_annot, " = ", post_annot, "\n");

    return eval_result<this_type>(view_base(),
                                  permuted_view(pre_annot, post_annot, false));
  }

  /// T†{a;i} = conj(T{i;a}) -- bra/ket-swapped layout (annotation rewrite
  /// baked into ann by the IR: operand annot in ann[0], adjoint annot in
  /// ann[1]) plus elementwise conjugation: a pending view, like permute(). For
  /// real numeric_type, conj is a no-op.
  [[nodiscard]] ResultPtr adjoint(
      std::array<std::any, 2> const& ann) const override {
    auto const pre_annot = std::any_cast<std::string>(ann[0]);
    auto const post_annot = std::any_cast<std::string>(ann[1]);

    detail::log_ta(post_annot, " = adjoint(", pre_annot, ")\n");

    return eval_result<this_type>(
        view_base(), permuted_view(pre_annot, post_annot,
                                   TA::detail::is_complex_v<numeric_type>));
  }

  [[nodiscard]] ResultPtr pow(rational const& exponent) const override {
    auto result = detail::ta_elementwise_pow(get<ArrayT>(), exponent);
    ::sequant::detail::wait_for_lazy_cleanup(result);
    log_ta_tensor_host_memory_use();
    return eval_result<this_type>(std::move(result));
  }

  [[nodiscard]] ResultPtr broadcast(
      std::array<std::any, 2> const& ann,
      container::svector<ResultPtr> const& new_mode_ranges) const override {
    auto const pre_annot = std::any_cast<std::string>(ann[0]);
    auto const post_annot = std::any_cast<std::string>(ann[1]);
    detail::log_ta(post_annot, " = broadcast(", pre_annot, ")\n");
    auto result = detail::ta_broadcast(get<ArrayT>().world(), &get<ArrayT>(),
                                       numeric_type{}, pre_annot, post_annot,
                                       new_mode_dims(new_mode_ranges));
    ::sequant::detail::wait_for_lazy_cleanup(result);
    log_ta_tensor_host_memory_use();
    return eval_result<this_type>(std::move(result));
  }

  [[nodiscard]] ResultPtr broadcast_scalar(
      Result const& scalar, std::array<std::any, 2> const& ann,
      container::svector<ResultPtr> const& new_mode_ranges) const override {
    if (!scalar.is<ResultScalar<numeric_type>>())
      throw Exception("broadcast: the scalar's type is not this array's");
    auto const post_annot = std::any_cast<std::string>(ann[1]);
    detail::log_ta(post_annot, " = broadcast(", scalar.get<numeric_type>(),
                   ")\n");
    auto result = detail::ta_broadcast<ArrayT>(
        get<ArrayT>().world(), nullptr, scalar.get<numeric_type>(),
        std::string{}, post_annot, new_mode_dims(new_mode_ranges));
    ::sequant::detail::wait_for_lazy_cleanup(result);
    log_ta_tensor_host_memory_use();
    return eval_result<this_type>(std::move(result));
  }

  /// the view of view_base() that is this result relabelled from \p pre_annot
  /// (its logical modes) to \p post_annot, conjugated once more if \p conj
  [[nodiscard]] view_type permuted_view(std::string const& pre_annot,
                                        std::string const& post_annot,
                                        bool conj) const {
    auto const pre = detail::outer_annot_labels(pre_annot);
    auto const post = detail::outer_annot_labels(post_annot);
    SEQUANT_ASSERT(pre.size() == post.size() && pre.size() == logical_rank());
    auto view = base_view();
    auto const rank = pre.size();
    // stored mode s holds logical mode j of this result, whose label pre[j]
    // sits at position k of the result's annotation: perm'[s] = k
    container::svector<std::size_t> perm(rank);
    for (std::size_t s = 0; s < rank; ++s) {
      auto const j = view.perm.empty() ? s : view.perm[s];
      auto const it = std::find(post.begin(), post.end(), pre[j]);
      SEQUANT_ASSERT(it != post.end() && "permute: label not in the result");
      perm[s] = static_cast<std::size_t>(it - post.begin());
    }
    bool identity = true;
    for (std::size_t s = 0; s < rank; ++s)
      if (perm[s] != s) identity = false;
    view.perm = identity ? container::svector<std::size_t>{} : std::move(perm);
    if (conj) view.conj = !view.conj;
    return view;
  }

  /// Accumulate \p other into this result's value (materialized first if this
  /// result is a pending view: an accumulation writes, and a view must not
  /// write through to its parent). The operand may be a pending view.
  void add_inplace(Result const& other) override {
    SEQUANT_ASSERT(other.is<this_type>());
    auto const& o = other.as<this_type>();

    auto& t = get<ArrayT>();

    SEQUANT_ASSERT(t.trange() == o.logical_trange());
    auto ann = TA::detail::dummy_annotation(t.trange().rank());

    detail::log_ta(ann, " += ", ann, "\n");

    auto const t0 = std::chrono::steady_clock::now();
    o.with_expr(ann, [&](auto&& oe) { t(ann) += oe; });
    ::sequant::detail::wait_for_lazy_cleanup(t);
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
  /// the TiledRange1 of each rank-1 range Result::broadcast receives
  [[nodiscard]] static container::svector<TA::TiledRange1> new_mode_dims(
      container::svector<ResultPtr> const& new_mode_ranges) {
    container::svector<TA::TiledRange1> dims;
    dims.reserve(new_mode_ranges.size());
    for (auto const& r : new_mode_ranges) {
      if (!r->is<this_type>())
        throw Exception(
            "broadcast: a mode range is not of this array's type; the "
            "BackendArrayOps must build arrays of the evaluated type");
      SEQUANT_ASSERT(r->get<ArrayT>().trange().rank() == 1);
      dims.push_back(r->get<ArrayT>().trange().dim(0));
    }
    return dims;
  }

  /// a pending view holds no tiles of its own: 0 bytes
  [[nodiscard]] std::size_t size_in_bytes() const final {
    if (!has_value()) return 0;
    if (Logger::instance().eval.level == 0)
      return 0;  // size_of disabled untraced: it is collective
    return measured_size_in_bytes();
  }

  [[nodiscard]] std::size_t measured_size_in_bytes() const final {
    if (!has_value()) return 0;
    auto& v = Result::get<ArrayT>();
    auto local_size = TA::size_of<TA::MemorySpace::Host>(v);
    v.world().gop.sum(local_size);
    return local_size;
  }
};

template <typename ArrayT,
          typename = std::enable_if_t<
              TA::detail::is_tensor_of_tensor_v<typename ArrayT::value_type>>>
class ResultTensorOfTensorTA final : public Result {
 public:
  using Result::id_t;
  using numeric_type = typename ArrayT::numeric_type;
  using view_type = detail::ArrayView;

  explicit ResultTensorOfTensorTA(ArrayT arr) : Result{std::move(arr)} {}

  /// A view of \p parent with \p view pending: a tile block of its OUTER
  /// modes (a batch slice). A nested result carries no pending permutation or
  /// conjugation (permute() and adjoint() apply theirs); otherwise as
  /// ResultTensorTA's views: shared tiles, no allocation, folded into the
  /// consuming expression by with_expr(), applied once by get<>(), never
  /// written through.
  ResultTensorOfTensorTA(ArrayT parent, view_type view)
      : Result{}, parent_(std::move(parent)), view_(std::move(view)) {
    SEQUANT_ASSERT(view_.perm.empty() && !view_.conj);
  }

  [[nodiscard]] std::string trange_annot() const override {
    return detail::ta_trange_string(logical_trange());
  }

  /// true while this result is a pending view sharing another array's tiles
  [[nodiscard]] bool is_buffer_alias() const noexcept override {
    return !has_value();
  }

  /// the OUTER TiledRange of the value this result stands for, with the
  /// pending block cut; never materializes
  [[nodiscard]] TA::TiledRange logical_trange() const {
    if (has_value()) return Result::get<ArrayT>().trange();
    return detail::logical_trange(parent_.trange(), view_);
  }

  /// the rank of the inner tensors of the array this result stands for, read
  /// from a populated inner tile; 0 when every inner tile is empty
  [[nodiscard]] std::size_t inner_rank() const {
    return detail::tot_inner_rank(view_base());
  }

  /// Hands \p f the TiledArray expression of this result's value for the
  /// nested annotation \p annot (`outer;inner`): `value(annot)` when the value
  /// is held, else the parent array's block over the outer modes, element
  /// coordinates kept (TA's preserve_lobound). The parent is handed over as a
  /// const array, which makes the block a view in TiledArray (see
  /// ResultTensorTA::with_expr). \p f must return the same type whatever
  /// expression it receives.
  template <typename F>
  decltype(auto) with_expr(std::string const& annot, F&& f) const {
    if (has_value()) return std::forward<F>(f)(Result::get<ArrayT>()(annot));
    auto const& parent = std::as_const(parent_);
    if (view_.blk_lo.empty()) return std::forward<F>(f)(parent(annot));
    return std::forward<F>(f)(
        parent(annot).block(view_.blk_lo, view_.blk_hi, TA::preserve_lobound));
  }

 private:
  using this_type = ResultTensorOfTensorTA<ArrayT>;
  using annot_wrap = Annot<std::string>;

  using _inner_tensor_type = typename ArrayT::value_type::value_type;

  // "Regular" (non-nested) companion array for ToT * T products. The OUTER
  // tile type must be a TA::Tensor — inner tile types like btas::Tensor are
  // only valid as the *innermost* tile (they don't support permute/reshape/
  // batch and so can't drive the outer kernel). So we wrap the inner's numeric
  // type in TA::Tensor here, rather than re-using the inner tile type as the
  // outer tile.
  using compatible_regular_distarray_type =
      TA::DistArray<TA::Tensor<numeric_type>, typename ArrayT::policy_type>;

  // Only @c that_type type is allowed for ToT * T computation
  using that_type = ResultTensorTA<compatible_regular_distarray_type>;

  /// the parent array of a pending view: an empty handle once the value is
  /// held (mutable: released by the const materialization in ensure_value())
  mutable ArrayT parent_;
  view_type view_;

  [[nodiscard]] id_t type_id() const noexcept override {
    return id_for_type<this_type>();
  }

  /// the array a further view of this result is a view of (the parent while
  /// this one is pending, else the held value), and that view's block
  [[nodiscard]] ArrayT const& view_base() const {
    return has_value() ? Result::get<ArrayT>() : parent_;
  }
  [[nodiscard]] view_type base_view() const {
    return has_value() ? view_type{} : view_;
  }

  /// a dummy nested annotation over the outer modes and \p inner inner modes
  [[nodiscard]] std::string dummy_annot(std::size_t inner) const {
    return TA::detail::dummy_annotation(
        static_cast<unsigned int>(view_base().trange().rank()),
        static_cast<unsigned int>(inner));
  }

  /// the value this result stands for, as a fresh array: the pending block
  /// applied in one pass (or a deep copy of the held value). A view of an
  /// array whose inner tiles are all empty has no inner rank to annotate a
  /// block with; the slice of zero is zero.
  [[nodiscard]] ArrayT materialize() const {
    if (has_value()) return TA::clone(Result::get<ArrayT>());
    auto const inner = inner_rank();
    if (inner == 0)
      return detail::make_zero_array<ArrayT>(view_base().world(),
                                             logical_trange());
    auto const ann = dummy_annot(inner);
    auto const t0 = std::chrono::steady_clock::now();
    ArrayT out;
    with_expr(ann, [&](auto&& e) { out(ann) = e; });
    ::sequant::detail::wait_for_lazy_cleanup(out);
    detail::log_batch_op("Materialize", std::chrono::steady_clock::now() - t0,
                         out, ann);
    return out;
  }

  /// Result::get<ArrayT>() on a pending view: the view applied, once; the
  /// parent's tiles are released
  void ensure_value() const override {
    if (has_value()) return;
    set_value(std::make_any<ArrayT>(materialize()));
    parent_ = ArrayT();
  }

  /// the value this result stands for, without storing it: the held value
  /// (shallow), or the pending view applied into a temporary. For read-only
  /// consumers that must not turn a shared view into a held value.
  [[nodiscard]] ArrayT applied_value() const {
    return has_value() ? Result::get<ArrayT>() : materialize();
  }

  // Diagnostic (see Result::fence): force all pending async work in this
  // array's world to complete so a lazily-executed op's timer captures
  // execution, not just dispatch.
  void fence() const noexcept override {
    view_base().world().gop.fence();
    ::sequant::detail::note_fence();
  }

  /// Allocating addition. For nested arrays an addend whose inner tiles are
  /// all empty acts as the additive identity, which relies on the arena
  /// nested-tensor addition the pinned TiledArray provides (an addition that
  /// dropped populated cells against a null-cell operand would annihilate the
  /// other addend and surface as a not-a-number norm in an iterative
  /// nested-array evaluation). No guard is needed here; \c add_inplace() keeps
  /// its explicit empty-accumulator branch for the in-place path.
  [[nodiscard]] ResultPtr sum(
      Result const& other,
      std::array<std::any, 3> const& annot) const override {
    SEQUANT_ASSERT(other.is<this_type>());
    auto const& o = other.as<this_type>();
    auto const a = annot_wrap{annot};

    detail::log_ta(a.lannot, " + ", a.rannot, " = ", a.this_annot, "\n");

    ArrayT result;
    with_expr(a.lannot, [&](auto&& le) {
      o.with_expr(a.rannot, [&](auto&& re) { result(a.this_annot) = le + re; });
    });
    ::sequant::detail::wait_for_lazy_cleanup(result);
    log_ta_tensor_host_memory_use();
    return eval_result<this_type>(std::move(result));
  }

  /// A batch slice is a pending view: the tile block of outer \p mode that
  /// covers `[elem_lo, elem_hi)` (tile-aligned, see slice_bounds_to_tiles),
  /// intersected with the block this result already is. No tile is copied.
  [[nodiscard]] ResultPtr slice_mode(std::size_t mode, std::size_t elem_lo,
                                     std::size_t elem_hi) const override {
    ::sequant::detail::note_slice();
    auto const& base = view_base();
    auto view = base_view();
    auto const [tile_lo, tile_hi] = detail::slice_bounds_to_tiles(
        base.trange().dim(mode), elem_lo, elem_hi);
    if (view.blk_lo.empty()) {
      auto const rank = base.trange().rank();
      view.blk_lo.assign(rank, 0);
      view.blk_hi.resize(rank);
      for (std::size_t d = 0; d < rank; ++d)
        view.blk_hi[d] = base.trange().dim(d).tile_extent();
    }
    SEQUANT_ASSERT(tile_lo >= view.blk_lo[mode] &&
                   tile_hi <= view.blk_hi[mode] &&
                   "slice_mode: the slice must lie within the block this "
                   "result already is");
    view.blk_lo[mode] = tile_lo;
    view.blk_hi[mode] = tile_hi;
    return eval_result<this_type>(base, std::move(view));
  }

  /// A pending block is reduced in place, not materialized. The norm does
  /// not depend on inner-mode labels, so the parent is annotated over its
  /// outer modes only (TA's make_tsrexpr) and no inner rank is needed.
  [[nodiscard]] double norm2() const override {
    if (has_value() || view_.blk_lo.empty())
      return detail::ta_norm2(view_base());
    if constexpr (requires(ArrayT const& a) {
                    { TA::squared_norm(a) } -> std::convertible_to<double>;
                  }) {
      auto const& parent = std::as_const(parent_);
      return static_cast<double>(
          parent
              .make_tsrexpr(TA::detail::dummy_annotation(
                  static_cast<unsigned int>(parent.trange().rank())))
              .block(view_.blk_lo, view_.blk_hi, TA::preserve_lobound)
              .norm()
              .get());
    } else {
      return -1.0;
    }
  }
  [[nodiscard]] std::string layout_desc() const override {
    return detail::ta_trange_string(logical_trange());
  }
  [[nodiscard]] std::string tile_diff(Result const& other) const override {
    if (!other.is<this_type>()) return "(kind mismatch)";
    return detail::ta_tile_diff(applied_value(),
                                other.as<this_type>().applied_value());
  }

  /// Scatter \p block into this result's value (materialized first if this
  /// result is a pending view: a scatter writes, and a view must not write
  /// through to its parent). The block may itself be a pending view: it is
  /// folded into the assignment. A block whose inner tiles are all empty is
  /// zero, and leaves the destination slice as it was.
  void write_into_slice(Result const& block, std::size_t mode,
                        std::size_t block_lo, std::size_t block_hi) override {
    SEQUANT_ASSERT(block.is<this_type>());
    auto const& b = block.as<this_type>();
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
                       detail::logical_trange(dest.trange(), {{}, lo, hi}) &&
                   "write_into_slice: the block's layout is not the "
                   "destination sub-block's");
    auto const annot = TA::detail::dummy_annotation(
        static_cast<unsigned int>(rank), static_cast<unsigned int>(inner));
    // preserve_lobound on the destination sub-block: the source block shares
    // its element coordinates (see write_array_into_mode)
    auto const t0 = std::chrono::steady_clock::now();
    b.with_expr(annot, [&](auto&& be) {
      dest(annot).block(lo, hi, TA::preserve_lobound) = be;
    });
    ::sequant::detail::wait_for_lazy_cleanup(dest);
    detail::log_batch_op("Scatter", std::chrono::steady_clock::now() - t0, dest,
                         annot + " mode=" + std::to_string(mode) + " tiles=[" +
                             std::to_string(tile_lo) + "," +
                             std::to_string(tile_hi) + ")");
    log_ta_tensor_host_memory_use();
  }

  /// Emit the general product `le * re` into \p result for the annotations
  /// \p a. TiledArray's contraction engine produces the result's inner modes
  /// in their natural order (the left operand's surviving inner labels in
  /// left order, then the right's); a result wanting another inner order is
  /// produced in the natural one and permuted into place, as einsum did.
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
    ::sequant::detail::wait_for_lazy_cleanup(tmp);
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

      ::sequant::detail::wait_for_lazy_cleanup(result);
      log_ta_tensor_host_memory_use();
      return eval_result<this_type>(std::move(result));
    } else if (a.this_annot.empty()) {
      // DOT product
      SEQUANT_ASSERT(other.is<this_type>());
      auto const& o = other.as<this_type>();
      numeric_type const d = with_expr(a.lannot, [&](auto&& le) {
        return o.with_expr(a.rannot, [&](auto&& re) {
          return static_cast<numeric_type>(le.dot(re).get());
        });
      });
      ::sequant::detail::wait_for_lazy_cleanup(view_base());
      ::sequant::detail::wait_for_lazy_cleanup(o.view_base());

      detail::log_ta(a.lannot, " * ", a.rannot, " = ", d, "\n");

      log_ta_tensor_host_memory_use();
      return eval_result<ResultScalar<numeric_type>>(d);
    }

    detail::log_ta(a.lannot, " * ", a.rannot, " = ", a.this_annot, "\n");

    // Every product is a TiledArray expression, so an operand that is a
    // pending view (a batch slice; for the flat operand also a permutation or
    // a conjugation) is taken as is.
    if (other.is<that_type>()) {
      // ToT * T -> ToT
      auto const& o = other.as<that_type>();
      ArrayT result;
      with_expr(a.lannot, [&](auto&& le) {
        o.with_expr(a.rannot,
                    [&](auto&& re) { emit_product(result, le, re, a); });
      });
      ::sequant::detail::wait_for_lazy_cleanup(result);
      log_ta_tensor_host_memory_use();
      return eval_result<this_type>(std::move(result));

    } else if (other.is<this_type>() && DeNestFlag == DeNest::True) {
      // ToT * ToT -> T: the inner tensors contracted away
      auto const& o = other.as<this_type>();
      compatible_regular_distarray_type result;
      with_expr(a.lannot, [&](auto&& le) {
        o.with_expr(a.rannot, [&](auto&& re) {
          result(a.this_annot) = le.dot_inner(re);
        });
      });
      ::sequant::detail::wait_for_lazy_cleanup(result);
      log_ta_tensor_host_memory_use();
      return eval_result<that_type>(std::move(result));

    } else if (other.is<this_type>() && DeNestFlag == DeNest::False) {
      // ToT * ToT -> ToT
      auto const& o = other.as<this_type>();
      ArrayT result;
      with_expr(a.lannot, [&](auto&& le) {
        o.with_expr(a.rannot,
                    [&](auto&& re) { emit_product(result, le, re, a); });
      });
      ::sequant::detail::wait_for_lazy_cleanup(result);
      log_ta_tensor_host_memory_use();
      return eval_result<this_type>(std::move(result));
    } else {
      throw detail::invalid_operand();
    }
  }

  /// A pending view is folded into the scaling, not materialized first. An
  /// array whose inner tiles are all empty is zero, and so is its phase.
  [[nodiscard]] ResultPtr mult_by_phase(std::int8_t factor) const override {
    auto const inner = inner_rank();
    if (inner == 0) return eval_result<this_type>(materialize());
    auto const ann = dummy_annot(inner);
    ArrayT result;
    with_expr(ann, [&](auto&& e) { result(ann) = numeric_type(factor) * e; });
    ::sequant::detail::wait_for_lazy_cleanup(result);
    log_ta_tensor_host_memory_use();
    return eval_result<this_type>(std::move(result));
  }

  /// Deep copy; see \c ResultTensorTA::clone. \c TA::clone deep-copies the
  /// (nested) tiles too, so the copy shares no inner tensor with this one. A
  /// pending view is applied into the copy directly.
  [[nodiscard]] ResultPtr clone() const override {
    return eval_result<this_type>(materialize());
  }

  [[nodiscard]] ResultPtr permute(
      std::array<std::any, 2> const& ann) const override {
    auto const pre_annot = std::any_cast<std::string>(ann[0]);
    auto const post_annot = std::any_cast<std::string>(ann[1]);

    detail::log_ta(pre_annot, " = ", post_annot, "\n");

    ArrayT result;
    with_expr(pre_annot, [&](auto&& e) { result(post_annot) = e; });
    ::sequant::detail::wait_for_lazy_cleanup(result);
    log_ta_tensor_host_memory_use();
    return eval_result<this_type>(std::move(result));
  }

  [[nodiscard]] ResultPtr adjoint(
      std::array<std::any, 2> const& ann) const override {
    // ToT adjoint: bra/ket-swapped layout (operand annot in ann[0], adjoint
    // annot in ann[1]) plus elementwise conj (a no-op for real numeric_type, so
    // it is elided there to avoid a ConjTsrExpr wrapper). TA's `.conj()`
    // recurses into nested tiles.
    auto const pre_annot = std::any_cast<std::string>(ann[0]);
    auto const post_annot = std::any_cast<std::string>(ann[1]);

    detail::log_ta(post_annot, " = adjoint(", pre_annot, ")\n");

    ArrayT result;
    with_expr(pre_annot, [&](auto&& e) {
      if constexpr (TA::detail::is_complex_v<numeric_type>) {
        result(post_annot) = e.conj();
      } else {
        result(post_annot) = e;
      }
    });
    ::sequant::detail::wait_for_lazy_cleanup(result);
    log_ta_tensor_host_memory_use();
    return eval_result<this_type>(std::move(result));
  }

  [[nodiscard]] ResultPtr pow(rational const& exponent) const override {
    auto result = detail::ta_elementwise_pow(get<ArrayT>(), exponent);
    ::sequant::detail::wait_for_lazy_cleanup(result);
    log_ta_tensor_host_memory_use();
    return eval_result<this_type>(std::move(result));
  }

  [[nodiscard]] ResultPtr broadcast(
      std::array<std::any, 2> const&,
      container::svector<ResultPtr> const&) const override {
    throw Exception("broadcast is not supported for tensor-of-tensor results");
  }

  [[nodiscard]] ResultPtr broadcast_scalar(
      Result const&, std::array<std::any, 2> const&,
      container::svector<ResultPtr> const&) const override {
    throw Exception("broadcast is not supported for tensor-of-tensor results");
  }

  /// Accumulate \p other into this result's value (materialized first if this
  /// result is a pending view). The operand may be a pending view.
  void add_inplace(Result const& other) override {
    SEQUANT_ASSERT(other.is<this_type>());
    auto const& o = other.as<this_type>();

    auto& t = get<ArrayT>();

    SEQUANT_ASSERT(t.trange() == o.logical_trange());
    // ToT annotation needs an inner block ("outer;inner"); tot_inner_rank()
    // reads it from a populated inner tile and is 0 when an operand's inner
    // tiles are all empty. Take the rank from whichever operand has data; if
    // both are empty the add is an identity no-op (t += 0), nothing to
    // annotate.
    auto const inner_rank = std::max(detail::tot_inner_rank(t), o.inner_rank());
    if (inner_rank == 0) return;
    auto ann = TA::detail::dummy_annotation(t.trange().rank(), inner_rank);

    detail::log_ta(ann, " += ", ann, "\n");

    auto const t0 = std::chrono::steady_clock::now();
    o.with_expr(ann, [&](auto&& oe) { t(ann) += oe; });
    ::sequant::detail::wait_for_lazy_cleanup(t);
    detail::log_batch_op("Accumulate", std::chrono::steady_clock::now() - t0, t,
                         ann);
    log_ta_tensor_host_memory_use();
  }

  [[nodiscard]] ResultPtr symmetrize() const override {
    return eval_result<this_type>(detail::column_symmetrize_ta(get<ArrayT>()));
  }

  [[nodiscard]] ResultPtr antisymmetrize(size_t /*bra_rank*/) const override {
    // not implemented yet
    return nullptr;
  }

 private:
  /// a pending view holds no tiles of its own: 0 bytes
  [[nodiscard]] std::size_t size_in_bytes() const final {
    if (!has_value()) return 0;
    if (Logger::instance().eval.level == 0)
      return 0;  // size_of disabled untraced: it is collective
    return measured_size_in_bytes();
  }

  [[nodiscard]] std::size_t measured_size_in_bytes() const final {
    if (!has_value()) return 0;
    auto& v = Result::get<ArrayT>();
    auto local_size = TA::size_of<TA::MemorySpace::Host>(v);
    v.world().gop.sum(local_size);
    return local_size;
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
template <typename... Args>
[[nodiscard]] TA::DistArray<Args...> slice_array_over_mode(
    TA::DistArray<Args...> const& arr, std::size_t mode, std::size_t tile_lo,
    std::size_t tile_hi) {
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
    auto const inner_rank = detail::tot_inner_rank(arr);
    if (inner_rank == 0) {
      // All-empty-inner ToT: tot_inner_rank() is 0, so there's no inner rank to
      // build the ToT annotation block() needs. The slice of zero is zero, so
      // build the sliced outer trange by hand (mode cut to [tile_lo,tile_hi),
      // every dim's original element lobound preserved, matching the
      // preserve_lobound block() below) and return a zero ToT over it.
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
      return detail::make_zero_array<TA::DistArray<Args...>>(
          arr.world(), TA::TiledRange(tr1s.begin(), tr1s.end()));
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
  ::sequant::detail::wait_for_lazy_cleanup(arr);
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
template <typename... Args>
void write_array_into_mode(TA::DistArray<Args...>& dest,
                           TA::DistArray<Args...> const& block,
                           std::size_t mode, std::size_t tile_lo,
                           std::size_t tile_hi) {
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
    auto const inner_rank = detail::tot_inner_rank(block);
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
  ::sequant::detail::wait_for_lazy_cleanup(dest);
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
    ::sequant::detail::wait_for_lazy_cleanup(out);
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
    ::sequant::detail::wait_for_lazy_cleanup(out);
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
      ::sequant::detail::wait_for_lazy_cleanup(out);
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
      ::sequant::detail::wait_for_lazy_cleanup(out);
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
