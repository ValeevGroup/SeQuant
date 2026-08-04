// Phase 3b T1: tests for the static peak-profile sizing primitives
// (SeQuant/core/eval/peak_profile.hpp). Two free functions are pinned here:
//   - home_depth_of: resolve a residency mode-set to an enclosing-batch-
//     context loop depth (mirrors the runtime rl-walk at eval.hpp:1776-1782).
//   - cell_footprint: home-relative footprint of a carried-index set via the
//     existing dryrun::CostModel::memsize (no sizing logic reimplemented).
// Neither function has a production caller yet (T2/T3 wire them in); this
// task only pins the two primitives' contracts.

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/index.hpp>

#include <catch2/catch_test_macros.hpp>

#include <cstddef>

namespace {

using sequant::Index;
using sequant::container::svector;
using sequant::eval::detail::BatchContext;
using sequant::eval::detail::cell_footprint;
using sequant::eval::detail::home_depth_of;
using sequant::eval::dryrun::CostModel;
using sequant::eval::dryrun::ExtentOverrides;
using sequant::eval::dryrun::SizeRegime;

}  // namespace

TEST_CASE("home_depth_of resolves the deepest enclosing loop in home_modes",
          "[peak_profile]") {
  // Three enclosing loops, outermost-first (matches BatchContext's documented
  // order): o (level 0), i (level 1), a (level 2, innermost).
  Index const o1{L"o_1"}, i1{L"i_1"}, a1{L"a_1"}, x1{L"x_1"};
  BatchContext const ectx{
      {o1, {0, 1}},
      {i1, {0, 1}},
      {a1, {0, 1}},
  };

  // Deepest (innermost) match.
  CHECK(home_depth_of(svector<Index>{a1}, ectx) == 2);
  // Outermost match.
  CHECK(home_depth_of(svector<Index>{o1}, ectx) == 0);
  // No match at all -> chain root.
  CHECK(home_depth_of(svector<Index>{}, ectx) == -1);
  // Mixed set: one member present (i1, level 1), one absent (x1) -> the
  // present member's level, not disturbed by the absent one.
  CHECK(home_depth_of(svector<Index>{i1, x1}, ectx) == 1);
}

TEST_CASE("cell_footprint sizes enclosing-home modes at BLOCK, others FULL",
          "[peak_profile]") {
  // full(i) = full(a) = 10; block(i) = 2 (a fake per-mode slice extent).
  SizeRegime r;
  r.space_extent = {{L"i", 10}, {L"a", 10}};
  CostModel const cm{r};

  Index const p{L"i_1"};  // the mode with an enclosing loop
  Index const q{L"a_1"};  // unbatched: no enclosing loop, no override
  svector<Index> const carried{p, q};

  // p's own loop is the sole enclosing context; block(p) narrows it to 2.
  BatchContext const ectx{{p, {0, 2}}};
  auto const block_of = [](Index const& ix) -> std::size_t {
    return ix == Index{L"i_1"} ? 2 : 10;
  };

  // Cell homed AT p (home_modes = {p}) -> d = 0 -> p's loop (level 0 <= d)
  // encloses the home -> p sized at BLOCK (2), q sized at FULL (10).
  {
    svector<Index> const home_modes{p};
    auto const got = cell_footprint(carried, home_modes, ectx, cm, block_of);
    ExtentOverrides ov;
    ov[p] = 2;
    CHECK(got == cm.memsize(carried, ov));
    // memsize() reports BYTES (elems * CostModel's numeric_size_ = 8.0), not
    // a raw element count: 2 (block p) * 10 (full q) * 8 bytes/elem.
    CHECK(got == 160);
  }

  // Same cell homed ABOVE p's loop (home_modes = {}) -> d = -1 -> no level
  // encloses the home -> both modes sized FULL.
  {
    svector<Index> const home_modes{};
    auto const got = cell_footprint(carried, home_modes, ectx, cm, block_of);
    CHECK(got == cm.memsize(carried));
    // 10 (full p) * 10 (full q) * 8 bytes/elem.
    CHECK(got == 800);
  }
}
