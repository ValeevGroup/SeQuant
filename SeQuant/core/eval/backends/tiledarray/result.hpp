#ifndef SEQUANT_EVAL_BACKENDS_TILEDARRAY_RESULT_HPP
#define SEQUANT_EVAL_BACKENDS_TILEDARRAY_RESULT_HPP

#ifdef SEQUANT_HAS_TILEDARRAY

#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/math.hpp>
#include <SeQuant/core/utility/exception.hpp>

#include <TiledArray/einsum/tiledarray.h>
#include <tiledarray.h>

#include <range/v3/view/iota.hpp>

#include <algorithm>

namespace sequant {

// implementation details of the TiledArray result backend; prefer
// sequant::detail over an unnamed namespace in a header (see CppCoreGuidelines
// SF.21 / "Use unnamed namespaces in headers ... no" guidance)
namespace detail {

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
      auto const& outer_tile = it->get();
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

  TA::DistArray<Args...>::wait_for_lazy_cleanup(result.world());

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

  TA::DistArray<Args...>::wait_for_lazy_cleanup(result.world());
  return result;
}

template <typename... Args>
inline void log_ta(Args const&... args) noexcept {
  log_result("[TA] ", args...);
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
[[nodiscard]] inline bool tot_product_needs_inner_reorder(
    std::string const& lannot, std::string const& rannot,
    std::string const& this_annot) {
  auto const li = inner_annot_labels(lannot);
  auto const ri = inner_annot_labels(rannot);
  auto const ti = inner_annot_labels(this_annot);
  if (ti.empty()) return false;  // denest / no inner result modes

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
  return natural != ti;  // identity inner perm iff the orders match
}

/// Build the result's OUTER TiledRange for a binary product, given the two
/// operand arrays and the three (left, right, result) annotations. Each result
/// outer label is matched, in order, against the left operand's outer labels
/// (then the right's), and the corresponding TiledRange1 is taken from that
/// operand. A result outer label must appear on at least one operand (true for
/// any well-formed contraction). This is the trange handed to the result-shape
/// provider so it can assemble a SparseShape over the result's outer tiles.
template <typename LArrayT, typename RArrayT>
[[nodiscard]] TA::TiledRange result_outer_trange(LArrayT const& larr,
                                                 std::string const& lannot,
                                                 RArrayT const& rarr,
                                                 std::string const& rannot,
                                                 std::string const& cannot) {
  auto const lo = outer_annot_labels(lannot);
  auto const ro = outer_annot_labels(rannot);
  auto const co = outer_annot_labels(cannot);

  auto find_dim = [](container::svector<std::string> const& labels,
                     std::string const& lbl) -> std::optional<std::size_t> {
    for (std::size_t i = 0; i < labels.size(); ++i)
      if (labels[i] == lbl) return i;
    return std::nullopt;
  };

  std::vector<TA::TiledRange1> dims;
  dims.reserve(co.size());
  for (auto const& lbl : co) {
    if (auto const i = find_dim(lo, lbl)) {
      dims.emplace_back(larr.trange().dim(*i));
    } else if (auto const j = find_dim(ro, lbl)) {
      dims.emplace_back(rarr.trange().dim(*j));
    } else {
      throw Exception("result_outer_trange: result outer label '" + lbl +
                      "' is carried by neither operand of the binary product");
    }
  }
  return TA::TiledRange(dims.begin(), dims.end());
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
    std::size_t tile_hi);

/// Partition a TiledRange1 into contiguous, tile-aligned element-range batches,
/// each covering at least \p target_batch_size elements where possible. Tiles
/// are not uniformly sized, so batches are uneven; the last batch may be
/// smaller, and any single tile larger than the target forms its own batch.
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
    acc += tr.second - tr.first;
    if (acc >= target || t + 1 == ntiles) {
      batches.emplace_back(grp_lo, tr.second);
      grp_lo = tr.second;
      acc = 0;
    }
  }
  return batches;
}

///
/// \brief Result for a tensor value of TA::DistArray type.
/// \tparam ArrayT TA::DistArray type. Tile type of ArrayT is regular tensor of
///                scalars (not a tensor of tensors)
///
template <typename ArrayT, typename = std::enable_if_t<TA::detail::is_tensor_v<
                               typename ArrayT::value_type>>>
class ResultTensorTA final : public Result {
 public:
  using Result::id_t;
  using numeric_type = typename ArrayT::numeric_type;

  explicit ResultTensorTA(ArrayT arr) : Result{std::move(arr)} {}

 private:
  using this_type = ResultTensorTA<ArrayT>;
  using annot_wrap = Annot<std::string>;

  [[nodiscard]] id_t type_id() const noexcept override {
    return id_for_type<this_type>();
  }

  [[nodiscard]] ResultPtr sum(
      Result const& other,
      std::array<std::any, 3> const& annot) const override {
    SEQUANT_ASSERT(other.is<this_type>());
    auto const a = annot_wrap{annot};

    detail::log_ta(a.lannot, " + ", a.rannot, " = ", a.this_annot, "\n");

    ArrayT result;
    result(a.this_annot) =
        get<ArrayT>()(a.lannot) + other.get<ArrayT>()(a.rannot);
    decltype(result)::wait_for_lazy_cleanup(result.world());
    log_ta_tensor_host_memory_use();
    return eval_result<this_type>(std::move(result));
  }

  [[nodiscard]] ResultPtr slice_mode(std::size_t mode, std::size_t elem_lo,
                                     std::size_t elem_hi) const override {
    auto const& tr1 = get<ArrayT>().trange().dim(mode);
    // slice_mode takes element bounds, but a tiled backend can only cut on tile
    // boundaries; mode_batches() returns exactly such (tile-aligned, in-range)
    // bounds. Assert the precondition so misuse is caught rather than silently
    // producing an over- or under-sized slice (which would break batched sums).
    SEQUANT_ASSERT(elem_lo >= tr1.elements_range().first && elem_lo < elem_hi &&
                   elem_hi <= tr1.elements_range().second);
    std::size_t const tile_lo = tr1.element_to_tile(elem_lo);
    SEQUANT_ASSERT(tr1.tile(tile_lo).first ==
                   elem_lo);  // lo on a tile boundary
    std::size_t const tile_hi = (elem_hi >= tr1.elements_range().second)
                                    ? tr1.tile_extent()
                                    : tr1.element_to_tile(elem_hi);
    SEQUANT_ASSERT(elem_hi >= tr1.elements_range().second ||
                   tr1.tile(tile_hi).first ==
                       elem_hi);  // hi on a tile boundary
    return eval_result<this_type>(
        slice_array_over_mode(get<ArrayT>(), mode, tile_lo, tile_hi));
  }

  [[nodiscard]] container::svector<std::pair<std::size_t, std::size_t>>
  mode_batches(std::size_t mode, std::size_t target_batch_size) const override {
    return mode_batches_of_trange1(get<ArrayT>().trange().dim(mode),
                                   target_batch_size);
  }

  [[nodiscard]] ResultPtr prod(Result const& other,
                               std::array<std::any, 3> const& annot,
                               DeNest DeNestFlag) const override {
    auto const a = annot_wrap{annot};

    if (other.is<ResultScalar<numeric_type>>()) {
      auto result = get<ArrayT>();
      auto scalar = other.get<numeric_type>();

      detail::log_ta(a.lannot, " * ", scalar, " = ", a.this_annot, "\n");

      result(a.this_annot) = scalar * result(a.lannot);

      decltype(result)::wait_for_lazy_cleanup(result.world());
      log_ta_tensor_host_memory_use();
      return eval_result<this_type>(std::move(result));
    }

    if (a.this_annot.empty()) {
      // DOT product
      SEQUANT_ASSERT(other.is<this_type>());
      numeric_type d =
          TA::dot(get<ArrayT>()(a.lannot), other.get<ArrayT>()(a.rannot));
      ArrayT::wait_for_lazy_cleanup(get<ArrayT>().world());
      ArrayT::wait_for_lazy_cleanup(other.get<ArrayT>().world());

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

    detail::log_ta(a.lannot, " * ", a.rannot, " = ", a.this_annot, "\n");

    ArrayT result;

    result = TA::einsum(get<ArrayT>()(a.lannot), other.get<ArrayT>()(a.rannot),
                        a.this_annot);
    decltype(result)::wait_for_lazy_cleanup(result.world());
    log_ta_tensor_host_memory_use();
    return eval_result<this_type>(std::move(result));
  }

  [[nodiscard]] ResultPtr mult_by_phase(std::int8_t factor) const override {
    auto pre = get<ArrayT>();
    TA::scale(pre, numeric_type(factor));
    return eval_result<this_type>(std::move(pre));
  }

  [[nodiscard]] ResultPtr permute(
      std::array<std::any, 2> const& ann) const override {
    auto const pre_annot = std::any_cast<std::string>(ann[0]);
    auto const post_annot = std::any_cast<std::string>(ann[1]);

    detail::log_ta(pre_annot, " = ", post_annot, "\n");

    ArrayT result;
    result(post_annot) = get<ArrayT>()(pre_annot);
    ArrayT::wait_for_lazy_cleanup(result.world());
    log_ta_tensor_host_memory_use();
    return eval_result<this_type>(std::move(result));
  }

  [[nodiscard]] ResultPtr adjoint(
      std::array<std::any, 2> const& ann) const override {
    // T†{a;i} = conj(T{i;a}) — bra/ket-swapped layout (annotation rewrite
    // baked into ann by the IR: operand annot in ann[0], adjoint annot in
    // ann[1]) plus elementwise conjugation. For real numeric_type, conj is
    // a no-op; elide it so the TA expression doesn't carry a ConjTsrExpr
    // wrapper unnecessarily.
    auto const pre_annot = std::any_cast<std::string>(ann[0]);
    auto const post_annot = std::any_cast<std::string>(ann[1]);

    detail::log_ta(post_annot, " = adjoint(", pre_annot, ")\n");

    ArrayT result;
    if constexpr (TA::detail::is_complex_v<numeric_type>) {
      result(post_annot) = get<ArrayT>()(pre_annot).conj();
    } else {
      result(post_annot) = get<ArrayT>()(pre_annot);
    }
    ArrayT::wait_for_lazy_cleanup(result.world());
    log_ta_tensor_host_memory_use();
    return eval_result<this_type>(std::move(result));
  }

  void add_inplace(Result const& other) override {
    SEQUANT_ASSERT(other.is<this_type>());

    auto& t = get<ArrayT>();
    auto const& o = other.get<ArrayT>();

    SEQUANT_ASSERT(t.trange() == o.trange());
    auto ann = TA::detail::dummy_annotation(t.trange().rank());

    detail::log_ta(ann, " += ", ann, "\n");

    t(ann) += o(ann);
    ArrayT::wait_for_lazy_cleanup(t.world());
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
    auto& v = get<ArrayT>();
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

  explicit ResultTensorOfTensorTA(ArrayT arr) : Result{std::move(arr)} {}

 private:
  using this_type = ResultTensorOfTensorTA<ArrayT>;
  using annot_wrap = Annot<std::string>;

  using _inner_tensor_type = typename ArrayT::value_type::value_type;

  // "Regular" (non-nested) companion array for ToT * T einsum. The OUTER tile
  // type must be a TA::Tensor — inner tile types like btas::Tensor are only
  // valid as the *innermost* tile (they don't support permute/reshape/batch
  // and so can't drive einsum's outer kernel). So we wrap the inner's numeric
  // type in TA::Tensor here, rather than re-using the inner tile type as the
  // outer tile.
  using compatible_regular_distarray_type =
      TA::DistArray<TA::Tensor<numeric_type>, typename ArrayT::policy_type>;

  // Only @c that_type type is allowed for ToT * T computation
  using that_type = ResultTensorTA<compatible_regular_distarray_type>;

  [[nodiscard]] id_t type_id() const noexcept override {
    return id_for_type<this_type>();
  }

  [[nodiscard]] ResultPtr sum(
      Result const& other,
      std::array<std::any, 3> const& annot) const override {
    SEQUANT_ASSERT(other.is<this_type>());
    auto const a = annot_wrap{annot};

    detail::log_ta(a.lannot, " + ", a.rannot, " = ", a.this_annot, "\n");

    ArrayT result;
    result(a.this_annot) =
        get<ArrayT>()(a.lannot) + other.get<ArrayT>()(a.rannot);
    decltype(result)::wait_for_lazy_cleanup(result.world());
    log_ta_tensor_host_memory_use();
    return eval_result<this_type>(std::move(result));
  }

  [[nodiscard]] ResultPtr slice_mode(std::size_t mode, std::size_t elem_lo,
                                     std::size_t elem_hi) const override {
    auto const& tr1 = get<ArrayT>().trange().dim(mode);
    // slice_mode takes element bounds, but a tiled backend can only cut on tile
    // boundaries; mode_batches() returns exactly such (tile-aligned, in-range)
    // bounds. Assert the precondition so misuse is caught rather than silently
    // producing an over- or under-sized slice (which would break batched sums).
    SEQUANT_ASSERT(elem_lo >= tr1.elements_range().first && elem_lo < elem_hi &&
                   elem_hi <= tr1.elements_range().second);
    std::size_t const tile_lo = tr1.element_to_tile(elem_lo);
    SEQUANT_ASSERT(tr1.tile(tile_lo).first ==
                   elem_lo);  // lo on a tile boundary
    std::size_t const tile_hi = (elem_hi >= tr1.elements_range().second)
                                    ? tr1.tile_extent()
                                    : tr1.element_to_tile(elem_hi);
    SEQUANT_ASSERT(elem_hi >= tr1.elements_range().second ||
                   tr1.tile(tile_hi).first ==
                       elem_hi);  // hi on a tile boundary
    return eval_result<this_type>(
        slice_array_over_mode(get<ArrayT>(), mode, tile_lo, tile_hi));
  }

  [[nodiscard]] container::svector<std::pair<std::size_t, std::size_t>>
  mode_batches(std::size_t mode, std::size_t target_batch_size) const override {
    return mode_batches_of_trange1(get<ArrayT>().trange().dim(mode),
                                   target_batch_size);
  }

  [[nodiscard]] ResultPtr prod(Result const& other,
                               std::array<std::any, 3> const& annot,
                               DeNest DeNestFlag) const override {
    auto const a = annot_wrap{annot};

    if (other.is<ResultScalar<numeric_type>>()) {
      auto result = get<ArrayT>();
      auto scalar = other.get<numeric_type>();

      detail::log_ta(a.lannot, " * ", scalar, " = ", a.this_annot, "\n");

      result(a.this_annot) = scalar * result(a.lannot);

      decltype(result)::wait_for_lazy_cleanup(result.world());
      log_ta_tensor_host_memory_use();
      return eval_result<this_type>(std::move(result));
    } else if (a.this_annot.empty()) {
      // DOT product
      SEQUANT_ASSERT(other.is<this_type>());
      numeric_type d =
          TA::dot(get<ArrayT>()(a.lannot), other.get<ArrayT>()(a.rannot));
      ArrayT::wait_for_lazy_cleanup(get<ArrayT>().world());
      ArrayT::wait_for_lazy_cleanup(other.get<ArrayT>().world());

      detail::log_ta(a.lannot, " * ", a.rannot, " = ", d, "\n");

      log_ta_tensor_host_memory_use();
      return eval_result<ResultScalar<numeric_type>>(d);
    }

    detail::log_ta(a.lannot, " * ", a.rannot, " = ", a.this_annot, "\n");

    if (other.is<that_type>()) {
      // ToT * T -> ToT
      auto result =
          TA::einsum(get<ArrayT>()(a.lannot),
                     other.get<compatible_regular_distarray_type>()(a.rannot),
                     a.this_annot);
      log_ta_tensor_host_memory_use();
      return eval_result<this_type>(std::move(result));

    } else if (other.is<this_type>() && DeNestFlag == DeNest::True) {
      // ToT * ToT -> T
      auto result = TA::einsum<TA::DeNest::True>(
          get<ArrayT>()(a.lannot), other.get<ArrayT>()(a.rannot), a.this_annot);
      log_ta_tensor_host_memory_use();
      return eval_result<that_type>(std::move(result));

    } else if (other.is<this_type>() && DeNestFlag == DeNest::False) {
      // ToT * ToT -> ToT
      auto result = TA::einsum(get<ArrayT>()(a.lannot),
                               other.get<ArrayT>()(a.rannot), a.this_annot);
      log_ta_tensor_host_memory_use();
      return eval_result<this_type>(std::move(result));
    } else {
      throw detail::invalid_operand();
    }
  }

  [[nodiscard]] ResultPtr mult_by_phase(std::int8_t factor) const override {
    auto pre = get<ArrayT>();
    TA::scale(pre, numeric_type(factor));
    log_ta_tensor_host_memory_use();
    return eval_result<this_type>(std::move(pre));
  }

  [[nodiscard]] ResultPtr permute(
      std::array<std::any, 2> const& ann) const override {
    auto const pre_annot = std::any_cast<std::string>(ann[0]);
    auto const post_annot = std::any_cast<std::string>(ann[1]);

    detail::log_ta(pre_annot, " = ", post_annot, "\n");

    ArrayT result;
    result(post_annot) = get<ArrayT>()(pre_annot);
    ArrayT::wait_for_lazy_cleanup(result.world());
    log_ta_tensor_host_memory_use();
    return eval_result<this_type>(std::move(result));
  }

  [[nodiscard]] ResultPtr adjoint(
      std::array<std::any, 2> const& ann) const override {
    // ToT adjoint: bra/ket-swapped layout (operand annot in ann[0], adjoint
    // annot in ann[1]) plus elementwise conj (a no-op for real numeric_type, so
    // it is elided there to avoid a ConjTsrExpr wrapper). Identical to the
    // regular-tensor branch above: TA's `.conj()` recurses into nested tiles.
    auto const pre_annot = std::any_cast<std::string>(ann[0]);
    auto const post_annot = std::any_cast<std::string>(ann[1]);

    detail::log_ta(post_annot, " = adjoint(", pre_annot, ")\n");

    ArrayT result;
    if constexpr (TA::detail::is_complex_v<numeric_type>) {
      result(post_annot) = get<ArrayT>()(pre_annot).conj();
    } else {
      result(post_annot) = get<ArrayT>()(pre_annot);
    }
    ArrayT::wait_for_lazy_cleanup(result.world());
    log_ta_tensor_host_memory_use();
    return eval_result<this_type>(std::move(result));
  }

  void add_inplace(Result const& other) override {
    SEQUANT_ASSERT(other.is<this_type>());

    auto& t = get<ArrayT>();
    auto const& o = other.get<ArrayT>();

    SEQUANT_ASSERT(t.trange() == o.trange());
    // ToT annotation needs an inner block ("outer;inner"); tot_inner_rank()
    // reads it from a populated inner tile and is 0 when an operand's inner
    // tiles are all empty. Take the rank from whichever operand has data; if
    // both are empty the add is an identity no-op (t += 0), nothing to
    // annotate.
    auto const inner_rank =
        std::max(detail::tot_inner_rank(t), detail::tot_inner_rank(o));
    if (inner_rank == 0) return;
    auto ann = TA::detail::dummy_annotation(t.trange().rank(), inner_rank);

    detail::log_ta(ann, " += ", ann, "\n");

    t(ann) += o(ann);
    ArrayT::wait_for_lazy_cleanup(t.world());
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
  [[nodiscard]] std::size_t size_in_bytes() const final {
    auto& v = get<ArrayT>();
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
  out(annot) = arr(annot).block(lo, hi, TA::preserve_lobound);
  TA::DistArray<Args...>::wait_for_lazy_cleanup(arr.world());
  return out;
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
  using FlatResult = ResultTensorTA<FlatArray>;
  using ToTResult = ResultTensorOfTensorTA<ToTArray>;

  auto const a = Annot<std::string>{annot};

  // Dispatch the left/right array extraction on operand kind, then defer to the
  // label-matching detail helper.  Four (flat/ToT)^2 combinations.
  auto with_left = [&](auto const& larr) -> TA::TiledRange {
    if (right.is<ToTResult>())
      return detail::result_outer_trange(larr, a.lannot, right.get<ToTArray>(),
                                         a.rannot, a.this_annot);
    return detail::result_outer_trange(larr, a.lannot, right.get<FlatArray>(),
                                       a.rannot, a.this_annot);
  };

  if (left.is<ToTResult>()) return with_left(left.get<ToTArray>());
  return with_left(left.get<FlatArray>());
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
/// the kept region. The two forms are the ones de-risked in Task 0.
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

  // Every branch fences the world after the assignment so the imposed shape
  // (held by pointer by the expression) is no longer referenced before this
  // function returns and the caller's shape object goes out of scope.

  // T * T -> T
  if (l_flat && r_flat) {
    FlatArray out;
    out(a.this_annot) =
        (left.get<FlatArray>()(a.lannot) * right.get<FlatArray>()(a.rannot))
            .set_shape(shape);
    out.world().gop.fence();
    FlatArray::wait_for_lazy_cleanup(out.world());
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
      out(a.this_annot) =
          (left.get<FlatArray>()(a.lannot) * right.get<ToTArray>()(a.rannot))
              .set_shape(shape);
    } else {
      out(a.this_annot) =
          (left.get<ToTArray>()(a.lannot) * right.get<FlatArray>()(a.rannot))
              .set_shape(shape);
    }
    out.world().gop.fence();
    ToTArray::wait_for_lazy_cleanup(out.world());
    return eval_result<ToTResult>(std::move(out));
  }

  // ToT * ToT
  if (l_tot && r_tot) {
    if (de_nest) {
      // ToT * ToT -> flat T (inner contraction / denest)
      FlatArray out;
      out(a.this_annot) = left.get<ToTArray>()(a.lannot)
                              .dot_inner(right.get<ToTArray>()(a.rannot))
                              .set_shape(shape);
      out.world().gop.fence();
      FlatArray::wait_for_lazy_cleanup(out.world());
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
      out(a.this_annot) =
          (left.get<ToTArray>()(a.lannot) * right.get<ToTArray>()(a.rannot))
              .set_shape(shape);
      out.world().gop.fence();
      ToTArray::wait_for_lazy_cleanup(out.world());
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
