#ifndef SEQUANT_EVAL_BACKENDS_TILEDARRAY_RESULT_HPP
#define SEQUANT_EVAL_BACKENDS_TILEDARRAY_RESULT_HPP

#ifdef SEQUANT_HAS_TILEDARRAY

#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/math.hpp>
#include <SeQuant/core/utility/exception.hpp>

#include <TiledArray/einsum/tiledarray.h>
#include <tiledarray.h>

#include <range/v3/view/iota.hpp>

#include <complex>

namespace sequant {

namespace {

///
/// \brief This function implements the symmetrization of TA::DistArray.
///
/// \param arr The array to be symmetrized
///
/// \pre The rank of the array must be even
///
/// \return The symmetrized TA::DistArray.
///
template <typename... Args>
auto column_symmetrize_ta(TA::DistArray<Args...> const& arr) {
  using ranges::views::iota;

  size_t const rank = arr.trange().rank();
  if (rank % 2 != 0)
    throw Exception("This function only supports even-ranked tensors");

  TA::DistArray<Args...> result;

  perm_t perm = iota(size_t{0}, rank) | ranges::to<perm_t>;

  auto const lannot = ords_to_annot(perm);

  auto call_back = [&result, &lannot, &arr, &perm = std::as_const(perm)]() {
    auto const rannot = ords_to_annot(perm);
    if (result.is_initialized()) {
      result(lannot) += arr(rannot);
    } else {
      result(lannot) = arr(rannot);
    }
  };

  auto const nparticles = rank / 2;
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

}  // namespace

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

    log_ta(a.lannot, " + ", a.rannot, " = ", a.this_annot, "\n");

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

      log_ta(a.lannot, " * ", scalar, " = ", a.this_annot, "\n");

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

      log_ta(a.lannot, " * ", a.rannot, " = ", d, "\n");

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

    log_ta(a.lannot, " * ", a.rannot, " = ", a.this_annot, "\n");

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

    log_ta(pre_annot, " = ", post_annot, "\n");

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

    log_ta(post_annot, " = adjoint(", pre_annot, ")\n");

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

    log_ta(ann, " += ", ann, "\n");

    t(ann) += o(ann);
    ArrayT::wait_for_lazy_cleanup(t.world());
    log_ta_tensor_host_memory_use();
  }

  [[nodiscard]] ResultPtr symmetrize() const override {
    return eval_result<this_type>(column_symmetrize_ta(get<ArrayT>()));
  }

  [[nodiscard]] ResultPtr antisymmetrize(size_t bra_rank) const override {
    return eval_result<this_type>(
        particle_antisymmetrize_ta(get<ArrayT>(), bra_rank));
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

    log_ta(a.lannot, " + ", a.rannot, " = ", a.this_annot, "\n");

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

      log_ta(a.lannot, " * ", scalar, " = ", a.this_annot, "\n");

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

      log_ta(a.lannot, " * ", a.rannot, " = ", d, "\n");

      log_ta_tensor_host_memory_use();
      return eval_result<ResultScalar<numeric_type>>(d);
    }

    log_ta(a.lannot, " * ", a.rannot, " = ", a.this_annot, "\n");

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
      throw invalid_operand();
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

    log_ta(pre_annot, " = ", post_annot, "\n");

    ArrayT result;
    result(post_annot) = get<ArrayT>()(pre_annot);
    ArrayT::wait_for_lazy_cleanup(result.world());
    log_ta_tensor_host_memory_use();
    return eval_result<this_type>(std::move(result));
  }

  [[nodiscard]] ResultPtr adjoint(
      std::array<std::any, 2> const& ann) const override {
    // ToT adjoint: bra/ket-swapped layout (operand annot in ann[0], adjoint
    // annot in ann[1]) plus elementwise conj (no-op for real numeric_type).
    // TA's `.conj()` expression has no overload for a tensor-of-tensors, so for
    // a complex numeric_type we permute first and then conjugate every inner
    // scalar in place via a local-tile walk (the same tile-iteration idiom as
    // tot_inner_rank() above). `it->get()` blocks until each permuted tile is
    // ready, so no extra fence is needed before mutating.
    auto const pre_annot = std::any_cast<std::string>(ann[0]);
    auto const post_annot = std::any_cast<std::string>(ann[1]);

    log_ta(post_annot, " = adjoint(", pre_annot, ")\n");

    ArrayT result;
    result(post_annot) = get<ArrayT>()(pre_annot);
    if constexpr (TA::detail::is_complex_v<numeric_type>) {
      for (auto it = result.begin(); it != result.end(); ++it) {
        auto& outer_tile = it->get();
        for (auto& inner : outer_tile) {
          if (inner.empty())
            continue;  // skip null inner views (screened pairs)
          for (auto& x : inner) x = std::conj(x);
        }
      }
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
    // ToT array: the annotation must carry an inner block ("outer;inner") or
    // DistArray::operator() rejects it (is_tot_index). dummy_annotation's
    // second argument is the inner-mode count, read from the array's inner
    // tiles.
    auto ann =
        TA::detail::dummy_annotation(t.trange().rank(), tot_inner_rank(t));

    log_ta(ann, " += ", ann, "\n");

    t(ann) += o(ann);
    ArrayT::wait_for_lazy_cleanup(t.world());
    log_ta_tensor_host_memory_use();
  }

  [[nodiscard]] ResultPtr symmetrize() const override {
    // not implemented yet
    return nullptr;
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
/// other mode; the result's `mode`-th TiledRange1 covers only those tiles
/// (its element range is shifted to start at 0). Implemented with TA's
/// `block()`, so block-sparse shape is preserved. Used to evaluate a tensor
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
    annot = TA::detail::dummy_annotation(
        static_cast<unsigned int>(rank),
        static_cast<unsigned int>(tot_inner_rank(arr)));
  } else {
    annot = ords_to_annot(iota(std::size_t{0}, rank));
  }
  TA::DistArray<Args...> out;
  out(annot) = arr(annot).block(lo, hi);
  TA::DistArray<Args...>::wait_for_lazy_cleanup(arr.world());
  return out;
}

}  // namespace sequant

#endif  // SEQUANT_HAS_TILEDARRAY

#endif  // SEQUANT_EVAL_BACKENDS_TILEDARRAY_RESULT_HPP
