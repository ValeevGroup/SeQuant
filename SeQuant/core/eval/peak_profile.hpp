#ifndef SEQUANT_CORE_EVAL_PEAK_PROFILE_HPP
#define SEQUANT_CORE_EVAL_PEAK_PROFILE_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/index.hpp>

#include <algorithm>
#include <cstddef>
#include <utility>

namespace sequant::eval {

///
/// \brief NEW, purely-additive static peak-profile analysis (spec section 9,
/// `doc/dev/specs/2026-08-03-meet-based-home-scope-phase3-design.md`) that
/// consumes the Phase-3a `home_scope` seed. This header lands its two sizing
/// primitives (Phase 3b T1) only -- ZERO production callers, runtime
/// untouched. T2/T3 (a separate task) wire these into an actual sweep.
///
namespace detail {

///
/// \brief The enclosing-batch-context type consulted below.
///
/// \details An ordered stack (outermost-first) of the enclosing realized
/// batch loops, one entry per loop, `{axis mode, {block_lo, block_hi}}`
/// (element range). This is verbatim the same underlying type as \c
/// CacheManager<TreeNode, force_hash_collisions>::BatchContext (see \c
/// cache_manager.hpp:95) -- that alias does not itself depend on \c
/// TreeNode, so this is not a look-alike but literally the same type, and
/// the real runtime \c BatchContext can be passed to \c home_depth_of /
/// \c cell_footprint as-is once a later task wires them in.
///
using BatchContext =
    container::svector<std::pair<Index, std::pair<std::size_t, std::size_t>>>;

///
/// \brief Resolve a residency mode-set to an enclosing-batch-context loop
/// depth.
///
/// \details Returns the deepest (innermost) level \p i of \p ectx whose loop
/// mode `ectx[i].first` is a member of \p home_modes; -1 if no level matches
/// (the cell is invariant to the whole nest -- the chain root). Mirrors the
/// runtime rl-walk at `eval.hpp:1776-1782` verbatim (same innermost-to-
/// outermost scan, same membership test), just against a caller-supplied
/// mode-set instead of the sliced/contracted-modes union computed there.
///
[[nodiscard]] inline int home_depth_of(
    container::svector<Index> const& home_modes,
    BatchContext const& ectx) noexcept {
  for (int i = static_cast<int>(ectx.size()) - 1; i >= 0; --i)
    if (std::find(home_modes.begin(), home_modes.end(), ectx[i].first) !=
        home_modes.end())
      return i;
  return -1;
}

///
/// \brief Home-relative footprint of a cell via the existing \c
/// dryrun::CostModel::memsize.
///
/// \details A cell homed at \p home_modes (depth `d = home_depth_of(
/// home_modes, ectx)`) sizes each CARRIED mode `m` at BLOCK extent (via \p
/// block_of) if `m`'s loop ENCLOSES the home -- i.e. `m` appears in \p ectx
/// at a level `<= d` -- else at FULL (nominal regime) extent. Builds the \c
/// dryrun::ExtentOverrides mapping exactly those block-sized modes to their
/// block extent, then delegates the actual extent-product / composite-moment
/// math to \p cm.memsize() verbatim: no sizing logic is reimplemented here.
///
/// \p block_of is any `Index -> std::size_t` callable giving the block
/// (sliced) element count for a mode.
///
template <typename BlockOfFn>
[[nodiscard]] inline std::size_t cell_footprint(
    container::svector<Index> const& carried,
    container::svector<Index> const& home_modes, BatchContext const& ectx,
    dryrun::CostModel const& cm, BlockOfFn const& block_of) {
  int const d = home_depth_of(home_modes, ectx);
  dryrun::ExtentOverrides ov;
  for (int i = 0; i <= d; ++i) {  // loops at level <= d ENCLOSE the home
    Index const& m = ectx[i].first;
    if (std::find(carried.begin(), carried.end(), m) != carried.end())
      ov[m] = block_of(m);  // sliced inside an enclosing loop
  }
  return cm.memsize(carried, ov);  // non-overridden carried modes size FULL
}

}  // namespace detail

}  // namespace sequant::eval

#endif  // SEQUANT_CORE_EVAL_PEAK_PROFILE_HPP
