#ifndef SEQUANT_CORE_EVAL_BACKENDS_TILEDARRAY_KRAMERS_FLIP_HPP
#define SEQUANT_CORE_EVAL_BACKENDS_TILEDARRAY_KRAMERS_FLIP_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <tiledarray.h>

#include <cstddef>
#include <cstdint>
#include <string>

namespace sequant::eval {

/// \brief Whether \p tr1 is a Kramers-union axis: `TA::concat(tr_up, tr_dn)`
///        of two congruently tiled halves (the ⇑ tiles first, then the ⇓
///        tiles, tile k of the ⇓ half as wide as tile k of the ⇑ half).
[[nodiscard]] inline bool union_halves_congruent(TA::TiledRange1 const& tr1) {
  auto const nt = static_cast<std::size_t>(tr1.tile_extent());
  if (nt % 2 != 0) return false;
  auto const h = nt / 2;
  for (std::size_t k = 0; k < h; ++k)
    if (tr1.tile(k).second - tr1.tile(k).first !=
        tr1.tile(k + h).second - tr1.tile(k + h).first)
      return false;
  return true;
}

/// \brief Number of tiles in one half of the Kramers-union axis \p tr1
///        (asserts union_halves_congruent).
[[nodiscard]] inline std::size_t union_half_tiles(TA::TiledRange1 const& tr1) {
  SEQUANT_ASSERT(union_halves_congruent(tr1) &&
                 "a Kramers-union axis has congruently tiled ⇑ and ⇓ halves");
  return static_cast<std::size_t>(tr1.tile_extent()) / 2;
}

namespace detail {

/// The annotation that names every mode of \p arr ("0,1,..." for a flat
/// array, "0,1,...;i0,i1,..." for a tensor of tensors, whose inner rank is
/// read from the first non-empty inner tile and agreed across the world).
/// Empty when the array is a tensor of tensors with no non-empty inner tile
/// anywhere (then it is identically zero).
template <typename ArrayT>
[[nodiscard]] std::string kramers_flip_annotation(ArrayT const& arr) {
  auto const rank = static_cast<unsigned int>(arr.trange().rank());
  if constexpr (TA::detail::is_tensor_of_tensor_v<
                    typename ArrayT::value_type>) {
    std::size_t inner = 0;
    for (auto it = arr.begin(); it != arr.end() && inner == 0; ++it) {
      // own the tile by value: it->get() returns by value and a reference
      // would dangle (see tot_inner_rank in result.hpp)
      auto const outer = it->get();
      if (outer.empty()) continue;
      for (auto const& in : outer)
        if (!in.empty()) {
          inner = in.range().rank();
          break;
        }
    }
    arr.world().gop.max(inner);
    if (inner == 0) return {};
    return TA::detail::dummy_annotation(rank, static_cast<unsigned int>(inner));
  } else {
    return TA::detail::dummy_annotation(rank);
  }
}

}  // namespace detail

/// \brief The time-reversal flip F of \p arr over the OUTER modes \p modes
///        (each a Kramers-union axis), scaled by \p phase.
///
/// Per mode, with the axis split as [⇑ | ⇓]:
///     out[⇑ half] = +conj in[⇓ half],   out[⇓ half] = −conj in[⇑ half]
/// i.e. the sign is that of the TARGET half; over several modes the signs
/// multiply and the conjugation is applied once. F∘F = −1 per mode (the
/// flip is not an involution, which is why it is an eval node and not a
/// retrieval transform). A folded intermediate is (−1)^{n_down} · F(partner):
/// see doc/dev/specs/2026-09-18-union-axis-time-reversal-fold.md in MPQC.
///
/// Every tile of the result is written (two block assignments per mode);
/// no contraction, O(size) data movement. Works for flat and nested
/// (tensor-of-tensor) arrays, dense and sparse policies (the sparse shape's
/// halves are swapped by TA's block assignment).
template <typename ArrayT>
[[nodiscard]] ArrayT kramers_flip_array(
    ArrayT const& arr, container::svector<std::size_t> const& modes,
    std::int8_t phase) {
  using numeric_type = typename ArrayT::numeric_type;
  using value_type = typename ArrayT::value_type;
  auto const& tr = arr.trange();
  auto const rank = tr.rank();
  auto const annot = detail::kramers_flip_annotation(arr);
  if (annot.empty()) return TA::clone(arr);  // identically zero
  if (modes.empty()) {
    // no union axis to swap: F is the conjugation alone (a flipped spelling
    // without union legs is phase * conj of its partner)
    ArrayT out;
    out(annot) = numeric_type(phase) * arr(annot).conj();
    return out;
  }
  ArrayT src = arr;
  bool first = true;
  for (auto const m : modes) {
    SEQUANT_ASSERT(m < rank);
    auto const h = union_half_tiles(tr.dim(m));
    container::svector<std::size_t> lo(rank, 0), hi(rank);
    for (std::size_t d = 0; d < rank; ++d)
      hi[d] = static_cast<std::size_t>(tr.dim(d).tile_extent());
    auto lo_up = lo, hi_up = hi, lo_dn = lo, hi_dn = hi;
    hi_up[m] = h;  // ⇑ tiles [0, h)
    lo_dn[m] = h;  // ⇓ tiles [h, 2h)
    // both halves of dst are overwritten below; the tiles only need to exist
    ArrayT dst(arr.world(), tr, src.shape(), src.pmap());
    dst.init_tiles([](TA::Range const& r) { return value_type(r); });
    // the phase rides on the first pass; the conjugation is applied once
    numeric_type const f_up = first ? numeric_type(phase) : numeric_type(1);
    numeric_type const f_dn = first ? numeric_type(-phase) : numeric_type(-1);
    if (first) {
      dst(annot).block(lo_up, hi_up) =
          f_up * src(annot).block(lo_dn, hi_dn).conj();
      dst(annot).block(lo_dn, hi_dn) =
          f_dn * src(annot).block(lo_up, hi_up).conj();
    } else {
      dst(annot).block(lo_up, hi_up) = f_up * src(annot).block(lo_dn, hi_dn);
      dst(annot).block(lo_dn, hi_dn) = f_dn * src(annot).block(lo_up, hi_up);
    }
    first = false;
    src = dst;
  }
  return src;
}

}  // namespace sequant::eval

#endif  // SEQUANT_CORE_EVAL_BACKENDS_TILEDARRAY_KRAMERS_FLIP_HPP
