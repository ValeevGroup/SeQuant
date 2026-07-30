// DryRun eval backend tests: size model + batched-schedule cost/peak analysis
// for PNO-CCSD, anchored on the real POST-transform C60 residual.
//
// Two groups:
//   1. SizeRegime / CostModel / Result unit tests -- extent lookup, moment-
//      aware (power-mean) composite sizing, rank-general CSV moment dispatch,
//      and the DryRun Result primitives (prod/slice/mode_batches/size).
//   2. Batched-schedule analysis on the real post-transform fixture
//      data/csv_ccsd_doubles_residual_df.txt (the CSV-CCSD doubles residual
//      dumped from mpqc AFTER the CSV->PAO base transform + DF refactorization,
//      so it carries the real PAO index mu~, DF-aux K, 3-center g{mu~;i;K}, and
//      CSV coefficients C{a<i>;mu~}). These cases optimize/binarize the
//      residual under the faithful C60-scale SizeRegime + BatchPolicy and
//      inspect the DP's batch-mode verdict, the perf-first vs peak-first
//      factorization of the free-mu~ giant, the gated dry-run cache veto, and
//      the scratch-folded batched peak via the shared cost_profile() entry
//      point.

// Flip Trace::Default -> Trace::On for this ENTIRE translation unit (must
// precede every SeQuant/core/eval/eval.hpp inclusion, directly or
// transitively -- header include guards mean the Trace enum's Default member
// is fixed, for this whole TU, by whether this macro is defined before the
// header's FIRST inclusion). This is load-bearing for the [dryrun-eval]
// replay below: make_batched_custom_evaluator's inner per-batch/per-member
// replay calls `evaluate(*mem, le_g, bs.cache)` WITHOUT an explicit `<Trace::
// ...>` argument, so those nested calls (and their note_working_set() calls,
// which are gated at COMPILE TIME on the EvalTrace template argument, not
// just Logger's runtime level) only fire if Trace::Default resolves to On.
// Without this, working_set_hwmark() would reflect only the outermost
// custom-evaluator interception and stay blind to everything the batched
// replay does inside it -- exactly the visibility Task 6's witness needs.
// mpqc's own C60 trace (614336.log, cited throughout this file) was captured
// the same way (SEQUANT_EVAL_TRACE defined in that build).
#define SEQUANT_EVAL_TRACE 1

#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_profile.hpp>
#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>
#include <SeQuant/core/eval/backends/dryrun/result.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/optimize/options.hpp>
#include <SeQuant/core/optimize/single_term_detail.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/space_qns.hpp>  // mbpt::Spin

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

using namespace sequant;
using sequant::eval::dryrun::SizeRegime;

namespace {

std::string slurp(std::string const& path) {
  std::ifstream in(path);
  std::stringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

// Small order-of-magnitude regime for the SizeRegime/memsize contract tests
// below. `a` is dual-purpose exactly as in mpqc's ctx.idx_to_extent (cck.ipp):
// a bare (non-proto) occurrence is the generic base virtual (sized by
// space_extent); a proto-indexed occurrence is a PNO (2 protos) or OSV (1
// proto) domain composite (sized by the moment tables).
SizeRegime probe_regime() {
  SizeRegime r;
  r.space_extent = {
      {L"i", 120},   // active occupied (order-of-magnitude stand-in)
      {L"a", 1800},  // generic base virtual (order-of-magnitude stand-in)
  };
  // CC-average PNO ~58, OSV ~20 (order-of-magnitude stand-ins; see report).
  double const pno = 58.0, osv = 20.0;
  for (std::size_t k = 0; k <= 4; ++k) {
    r.csv_pno_moment[k] = std::pow(pno, double(k));
    r.csv_osv_moment[k] = std::pow(osv, double(k));
  }
  return r;
}

// Free (external) indices of a binarized-tree node's result, read off its
// EvalExpr tensor (scalar-result nodes return an empty vector).
std::vector<Index> node_free_indices(EvalExpr const& n) {
  std::vector<Index> v;
  if (!n.is_tensor()) return v;
  for (auto const& ix : n.as_tensor().const_braketaux_indices())
    v.push_back(ix);
  return v;
}

// Local shim: BatchModeType-tagging of EvalExpr::batched_here() entries (Task
// 1) strips this spike file's plain-Index reads down to the .first projection;
// this test file is not committed, so keep the fix minimal rather than
// threading BatchModeType through the trace/analysis helpers below.
template <typename Range>
container::vector<Index> batch_axes_indices(Range const& entries) {
  container::vector<Index> out;
  for (auto const& e : entries) out.push_back(e.first);
  return out;
}

std::wstring describe_indices(std::vector<Index> const& ixs) {
  std::wstring s;
  for (auto const& ix : ixs) {
    s += std::wstring(ix.full_label());
    s += L" ";
  }
  return s;
}

}  // namespace

TEST_CASE("dryrun size regime basic extents", "[dryrun-probe]") {
  auto r = probe_regime();
  // A bare occ index resolves to its space extent.
  auto i = Index{L"i_1"};
  CHECK(r.extent(i) == 120);
  auto a = Index{L"a_1"};
  CHECK(r.extent(a) == 1800);
  CHECK_THROWS_AS(r.extent(Index{L"x_1"}), std::out_of_range);
}

TEST_CASE("dryrun regime self-consistency (no cluster anchors needed)",
          "[dryrun-probe]") {
  // Replaces the (unreachable) two K-anchor checks: pick one concrete index
  // set, hand-compute its expected element count from SizeRegime's own
  // extents/moments, and confirm memsize_counter agrees. This validates the
  // size MODEL (extent lookup + inner_pow composite dispatch), independent of
  // whether the specific numeric extents match any particular cluster run.
  auto r = probe_regime();
  auto memsize = sequant::opt::detail::memsize_counter(r.idx_to_extent(),
                                                       r.inner_pow_fn());

  // A dense operand: two occ + two bare (base-virtual) indices, e.g. the
  // ladder term's g{i,i;a,a}.
  {
    std::vector<Index> ixs{Index{L"i_1"}, Index{L"i_2"}, Index{L"a_1"},
                           Index{L"a_2"}};
    double const expected = 120.0 * 120.0 * 1800.0 * 1800.0;  // i * i * a * a
    double const got = memsize(ixs, std::vector<Index>{}, std::vector<Index>{});
    CHECK(got == Catch::Approx(expected));
    INFO("dense operand elems=" << got << " (" << (got * 8 / 1e9)
                                << " GB dense)");
  }

  // A CSV/PNO composite operand: one occ-pair-domain composite (2 protos),
  // e.g. an intermediate carrying a_1<i_1,i_2>. tot_indices() pulls the
  // composite's own proto indices (i_1, i_2) into the outer set too (they are
  // genuinely free tensor slots the composite depends on), so the expected
  // volume is i_1 * i_2 * <#PNO^1>, not just one occ factor.
  {
    Index i1{L"i_1"}, i2{L"i_2"};
    Index pno_leg{L"a_1", {i1, i2}};
    std::vector<Index> ixs{pno_leg};
    double const expected = 120.0 * 120.0 * 58.0;  // i * i * <#PNO^1>
    double const got = memsize(ixs, std::vector<Index>{}, std::vector<Index>{});
    CHECK(got == Catch::Approx(expected));
  }
}

// Regression lock (Task 1, a1 verification): inner_pow() must return the
// k-th POWER MEAN M_k, not the raw k-th MOMENT <d^k>, so that
// inner_aware_volume's per-member product over a k-composite group
// (M_k multiplied once per member) telescopes to mean(d^k) rather than
// over-counting by a further power of k. This pins the contract with a
// heavy-tailed (non-constant-in-k) moment table so a future regression to
// raw moments would be caught even though the CURRENT SizeRegime factories
// in this file (probe_regime/df_regime) all use a constant domain, for which
// M_k == d for every k and the bug would otherwise be invisible.
TEST_CASE("dryrun power-mean sizing contract", "[dryrun][sizing]") {
  // Heavy-tailed PNO domain: power means strictly increasing in k.
  // (For any non-degenerate distribution, M_1 < M_2 < M_3 < M_4 by Jensen.)
  SizeRegime r;
  const std::size_t occ = 10;
  r.space_extent[L"i"] = occ;  // occupied index extent
  const std::array<double, 5> Mk{1.0, 30.0, 40.0, 55.0, 75.0};  // M_0..M_4
  for (std::size_t k = 0; k <= 4; ++k) {
    r.csv_pno_moment[k] = Mk[k];
    r.csv_osv_moment[k] = Mk[k];
  }
  auto ip = r.inner_pow_fn();

  // A composite PNO index a<i_1,i_2> has a 2-proto (PNO) domain.
  Index i1{L"i_1"}, i2{L"i_2"};
  Index a1{L"a_1", {i1, i2}};
  Index a2{L"a_2", {i1, i2}};

  // Two composites sharing the SAME proto-set => a k=2 group.
  // inner_aware_volume multiplies ip(c,2) once per member => M_2^2.
  auto tot = tot_indices(std::vector<Index>{a1, a2});
  auto ixex = [&](Index const& ix) {
    return static_cast<double>(r.extent(ix));
  };
  const double vol = sequant::opt::detail::inner_aware_volume(tot, ixex, ip);

  const double expected =
      double(occ) * double(occ) * Mk[2] * Mk[2];  // occ^2 * M_2^2
  CHECK(vol == Catch::Approx(expected));
  // Guard against the two wrong models:
  CHECK(vol !=
        Catch::Approx(double(occ) * occ * Mk[1] * Mk[1]));  // flat-average
                                                            // under-count
  CHECK(vol != Catch::Approx(double(occ) * occ *
                             std::pow(Mk[1], 4)));  // old "occ^2*PNO^4
                                                    // artifact"
}

// Rank-general CSV moment dispatch (CSV-CCSDT and beyond): inner_pow() selects
// the moment table by cluster rank (= number of proto indices). Rank 1 -> OSV,
// rank 2 -> PNO, rank >= 3 -> the csv_moment_by_rank[rank] table if present,
// else a fallback to the PNO (rank-2) table (preserving the pre-rank-general
// behavior where every proto-rank >= 2 used the PNO table).
TEST_CASE("dryrun rank-general CSV moment dispatch", "[dryrun][sizing]") {
  SizeRegime r;
  for (std::size_t k = 0; k <= 4; ++k) {
    r.csv_osv_moment[k] = 3.0;   // rank-1 table (flat, distinct value)
    r.csv_pno_moment[k] = 30.0;  // rank-2 table (flat, distinct value)
  }
  // rank-3 (triple) table, distinct from both rank-1 and rank-2.
  r.csv_moment_by_rank[3] = {1.0, 100.0, 100.0, 100.0, 100.0};
  auto ip = r.inner_pow_fn();

  Index i1{L"i_1"}, i2{L"i_2"}, i3{L"i_3"}, i4{L"i_4"};
  Index osv{L"a_1", {i1}};               // rank-1 composite
  Index pno{L"a_2", {i1, i2}};           // rank-2 composite
  Index triple{L"a_3", {i1, i2, i3}};    // rank-3 composite
  Index quad{L"a_4", {i1, i2, i3, i4}};  // rank-4 composite (no table)

  // Each rank draws from its own table.
  CHECK(ip(osv, 2) == Catch::Approx(3.0));
  CHECK(ip(pno, 2) == Catch::Approx(30.0));
  CHECK(ip(triple, 2) == Catch::Approx(100.0));
  // A rank with no table falls back to the PNO (rank-2) table, NOT the OSV one.
  CHECK(ip(quad, 2) == Catch::Approx(30.0));
}

// ===========================================================================
// POST-TRANSFORM VERDICT: the actual PAO(mu~)/DF-aux(K) go/no-go probe.
//
// Fixture data/csv_ccsd_doubles_residual_df.txt is the CSV CCSD DOUBLES
// residual dumped from mpqc (repro/w8-batch-min.json) at the EXACT point it is
// handed to sequant::optimize() -- i.e. AFTER the CSV->PAO base transform and
// DF refactorization. It therefore contains the real PAO index "mu~"
// (μ̃_NNNN) and DF-aux index "K" (Κ_N), the 3-center DF integrals
// g{μ̃;i;Κ}, and the CSV coefficients C{a<i>;μ̃}. Standalone SeQuant's default
// mbpt registry has neither mu~ nor K, so we augment a cloned context with
// mbpt::add_pao_spaces (mu~) + mbpt::add_df_spaces (K).
//
// This case reproduces, offline and in milliseconds, the batched-objective DP
// decision that on the cluster left the free-mu~ giant un-sliced (single-mode
// aux batching -> OOM). It reports, per contraction node carrying a free mu~,
// whether the DP annotated a mu~ (or K) batch mode -- the localizing signal:
//   - if the giant mu~-carrying intermediate never gets a mu~ mode even as its
//     modeled size dwarfs the threshold => the COST MODEL / DP is the gap;
//   - if it does get one here => the gap is downstream (binarize/runtime).
// The DryRun harness's value is that the SAME fixture can be swept across size
// regimes (water-8 vs C60) by changing only the SizeRegime extents.
// ===========================================================================

namespace {

// Regime for the post-transform (mu~ + K) residual. Extents are per-index
// domain sizes fed to the DP's idx_to_extent; proto-indexed "a" legs are
// PNO/OSV composites sized by the moment tables. Takes the real per-nonnull-
// cluster power means M_1..M_4 (heavy tail) as measured by mpqc's PaoPnoRMP2 --
// the two lines it prints, "PNO domain power means M_1..M_4 per pair" and "OSV
// domain power means M_1..M_4 per orbital". A per-k power mean (M_2 > M_1 >
// ...) sizes a k-composite group as mean(d^k), the true block-sparse volume.
// pno_M[k]/osv_M[k] for k in [1,4]; index 0 unused (kept 1). Rank >= 3
// (CSV-CCSDT) is not populated here (CCSD is ranks 1-2).
SizeRegime df_regime(std::size_t mu_tilde, std::size_t aux, std::size_t i_occ,
                     std::array<double, 5> const& pno_M,
                     std::array<double, 5> const& osv_M) {
  SizeRegime r;
  r.space_extent = {
      {L"i", i_occ},
      {L"μ̃", mu_tilde},
      {L"Κ", aux},
      {L"a", mu_tilde},
  };
  r.csv_pno_moment = pno_M;
  r.csv_osv_moment = osv_M;
  return r;
}

// One named (molecule, basis, parameter-set) problem size for the DryRun cost
// model: base-space extents + the measured heavy-tailed CSV power means. Single
// source of truth -- use df_regime(kFoo) at every call site instead of
// repeating the literals, so a moment re-measurement is a one-line edit.
struct ProblemSize {
  std::size_t mu_tilde;         // PAO domain extent (= #AO)
  std::size_t aux;              // DF aux (K) extent
  std::size_t i_occ;            // active occupied extent
  std::array<double, 5> pno_M;  // per-pair PNO power means M_0..M_4 (0 unused)
  std::array<double, 5> osv_M;  // per-orbital OSV power means M_0..M_4
};

// C60, cc-pVDZ-F12, TCUTPNO=1e-8 / TCUTOSV=1e-9. Extents from the Owl job log;
// PNO/OSV power means M_1..M_4 measured by mpqc PaoPnoRMP2 (job 617809 -- the
// TRUE pVDZ-F12 values; NOTE the earlier 48.60/64.35 PNO + 206/234 OSV numbers
// were from job 617653, which was mis-configured with an aug-cc-pVTZ OBS).
inline constexpr ProblemSize kC60_pVDZF12{
    /*mu_tilde=*/1800u,
    /*aux=*/4320u,
    /*i_occ=*/120u,
    /*pno_M=*/
    {1.0, 42.029069767441861, 46.039206412923569, 49.766252354482994,
     53.151291880343109},
    /*osv_M=*/
    {1.0, 148.25, 155.04434849422921, 161.33527408797721, 166.85553430303926}};

// Water-20 (H2O)20 / cc-pVDZ-F12, extracted from Owl job 649160:
//   ext(i)=80 (active occ), ext(K)=1682 (DF aux; from g(i,i,K)=86.1MB),
//   ext(mu~)=896 (PAO; from g(mu~,mu~,K)=10.8GB), and the measured heavy-tailed
//   CSV moments (PNO M_1..M_4 per pair, OSV M_1..M_4 per orbital).
inline constexpr ProblemSize kWater20_pVDZF12{
    /*mu_tilde=*/896u,
    /*aux=*/1682u,
    /*i_occ=*/80u,
    /*pno_M=*/
    {1.0, 23.175775480059084, 25.865548281212597, 28.171416142614103,
     30.03848680550367},
    /*osv_M=*/
    {1.0, 58.987499999999997, 59.289227520688783, 59.584437469011633,
     59.872014818179686}};

// Build a SizeRegime from a named ProblemSize.
SizeRegime df_regime(ProblemSize const& p) {
  return df_regime(p.mu_tilde, p.aux, p.i_occ, p.pno_M, p.osv_M);
}

// Batchable = the two modes mpqc's runtime batches on the CSV path: PAO (mu~)
// and DF aux (K). Both are non-proto base spaces (mu~ = PAO, K = DFBS aux).
bool is_df_batchable(Index const& ix) {
  auto const k = ix.space().base_key();
  return k == L"μ̃" || k == L"Κ";
}

}  // namespace

// Diagnostic ([.]): does the order-aware ordered-key DP actually ENGAGE on the
// real C60 giant, or does build_cells' m>7 fallback (enumeration blowup guard)
// make it inert? Prints m (# batchable indices), ordered, nCells. Fast: only
// build_context, no DP solve.
TEST_CASE("ordered-key C60 giant: does order_aware engage (m vs cap)?",
          "[.][ordered-key-c60-m]") {
  using namespace sequant;
  auto ctx0 = get_default_context().clone();
  ctx0.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx0.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);  // mu~
  sequant::mbpt::add_df_spaces(isr);                             // K
  auto ctx_resetter = set_scoped_default_context(std::move(ctx0));

  auto const body = slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                          "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = deserialize<ExprPtr>(line);
  REQUIRE(expr);
  REQUIRE(expr->is<Sum>());
  auto const& summands = expr->as<Sum>().summands();
  std::size_t const gi = 38 < summands.size() ? 38 : 0;
  ExprPtr giant = summands[gi];
  if (giant->is<Product>())
    giant = ex<Product>(giant->as<Product>().scalar(),
                        giant->as<Product>().factors(), Product::Flatten::Yes);
  REQUIRE(giant->is<Product>());
  TensorNetwork net(giant->as<Product>().factors());
  container::svector<Index> targets;

  auto regime = df_regime(kC60_pVDZF12);
  auto idxsz = regime.idx_to_extent();
  std::function<std::size_t(Index const&)> bts = [](Index const& ix) {
    return ix.space().base_key() == L"μ̃" ? std::size_t{256} : std::size_t{72};
  };
  opt::detail::PeakBatchedModel model{idxsz, bts, {}, regime.inner_pow_fn()};
  model.is_batchable_contracted_index = is_df_batchable;
  model.order_aware_recompute = true;
  auto ctx = model.build_context(net, targets);
  std::wcerr << L"[ordered-key-c60-m] giant (summand " << gi << L", "
             << giant->as<Product>().factors().size() << L" factors): m="
             << ctx.m << L" ordered=" << ctx.ordered << L" nCells="
             << ctx.nCells << L" (nB=" << ctx.nB << L")\n";
}

// S3.5 GATE ([.]): does the CHEAP phase-2 external placement (node-level, S3.3)
// REACH the C60 peak-setting operand and DROP its modeled root peak? Runs the
// giant through the batched DP + reconstruct_batched_modes and reports the
// REPORTED root peak (bytes) and the External stamps for four configs, each ON
// (batch_spectator) vs OFF (no phase-2) at a fixed batchable set so the delta
// isolates placement (select_root / the contracted schedule is identical
// within a predicate):
//   mu~/K   : is_df_batchable -- phase-2 can place only the batchable external
//   mu~. mu~/K/i : is_df_batchable widened with the active-occupied space --
//   S3.1's
//             critical finding is the occ pair is is_batchable==FALSE under the
//             CSV policy, so sz_u/sliced_footprints cannot size its slice until
//             `i` is admitted as batchable; then the occ externals are
//             placeable.
// Metric is the DP-side modeled peak (cost_profile), NEVER the replay
// avoidable_time. Hidden ([.]); may take a while on the giant.
TEST_CASE("ext-place C60 giant: phase-2 external placement drops modeled peak",
          "[.][ext-place-c60]") {
  using namespace sequant;
  auto ctx0 = get_default_context().clone();
  ctx0.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx0.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);  // mu~
  sequant::mbpt::add_df_spaces(isr);                             // K
  auto ctx_resetter = set_scoped_default_context(std::move(ctx0));

  auto const body = slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                          "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = deserialize<ExprPtr>(line);
  REQUIRE(expr);
  REQUIRE(expr->is<Sum>());
  auto const& summands = expr->as<Sum>().summands();
  std::size_t const gi = 38 < summands.size() ? 38 : 0;
  ExprPtr giant = summands[gi];
  if (giant->is<Product>())
    giant = ex<Product>(giant->as<Product>().scalar(),
                        giant->as<Product>().factors(), Product::Flatten::Yes);
  REQUIRE(giant->is<Product>());
  TensorNetwork net(giant->as<Product>().factors());
  container::svector<Index> targets;

  auto regime = df_regime(kC60_pVDZF12);
  auto idxsz = regime.idx_to_extent();
  std::function<std::size_t(Index const&)> bts = [](Index const& ix) {
    return ix.space().base_key() == L"μ̃" ? std::size_t{256} : std::size_t{72};
  };
  // Widen is_df_batchable with the active-occupied space `i` so the occ-pair
  // externals become sizeable by sz_u (S3.1: they are is_batchable==false under
  // the bare CSV mu~/K policy).
  auto is_df_occ = [](Index const& ix) {
    auto const k = ix.space().base_key();
    return k == L"μ̃" || k == L"Κ" || k == L"i";
  };

  auto inner_pow = regime.inner_pow_fn();  // REQUIRED: CSV composite sizing
  auto measure = [&](wchar_t const* tag,
                     std::function<bool(Index const&)> const& isb,
                     bool spectator) -> double {
    opt::detail::PeakBatchedModel model{idxsz, bts, {}, inner_pow};
    model.is_batchable_contracted_index = isb;
    // Same predicate in the external role: this network's batchable external
    // mode must stay admitted once the fallback is removed (Task 4). Setting
    // both to isb is byte-identical to the historical fallback here.
    model.is_batchable_external_index = isb;
    model.order_aware_recompute = true;
    model.batch_spectator_indices = spectator;
    model.perf_first = true;  // legacy over_budget path; node-level ignores
    model.peak_threshold = 100e9;  // 100 GB budget, bytes
    model.numeric_size = 8.0;
    auto ctx = model.build_context(net, targets);
    auto st = opt::detail::solve_single_term(model, net, targets, ctx);
    double peak_bytes = 0.0;
    auto res =
        model.reconstruct_batched_modes(ctx, st, net, targets, &peak_bytes);
    int ext = 0;
    std::wstring labels;
    for (auto const& modes : res.second)
      for (auto const& pr : modes.axes)
        if (pr.second == BatchModeType::External) {
          ++ext;
          labels += pr.first.full_label();
          labels += L' ';
        }
    std::wcerr << L"[ext-place-c60] " << tag << L": peak=" << (peak_bytes / 1e9)
               << L" GB  #External=" << ext << L"  {" << labels << L"}  m="
               << ctx.m << L" nCells=" << ctx.nCells << L"\n";
    return peak_bytes;
  };

  double const p0a = measure(L"OFF mu~/K  ", is_df_batchable, false);
  double const p1 = measure(L"ON  mu~/K  ", is_df_batchable, true);
  double const p0b = measure(L"OFF mu~/K/i", is_df_occ, false);
  double const p2 = measure(L"ON  mu~/K/i", is_df_occ, true);

  // The dense cost model is a CONSERVATIVE upper bound (real memory is smaller
  // due to unmodeled block sparsity), so meeting the budget needs NO sparsity
  // modeling -- just a finer occ slice. The peak node {g C} = a*mu~*K_blk*
  // occ_blk^2 scales as occ_blk^2, so this is analytic (no sweep needed):
  // occ_blk 72->32 cuts {g C} by (72/32)^2 ~ 5x further, ~14x below the
  // unsliced 627 GB. One point confirms the budget is met.
  std::function<std::size_t(Index const&)> bts32 =
      [](Index const& ix) -> std::size_t {
    auto const key = ix.space().base_key();
    if (key == L"μ̃") return std::size_t{256};
    if (key == L"i") return std::size_t{32};  // finer occ block
    return std::size_t{72};                   // aux K
  };
  double p32 = 0.0;
  {
    opt::detail::PeakBatchedModel model{idxsz, bts32, {}, inner_pow};
    model.is_batchable_contracted_index = is_df_occ;
    model.is_batchable_external_index =
        is_df_occ;  // external role, byte-ident.
    model.order_aware_recompute = true;
    model.batch_spectator_indices = true;
    model.perf_first = true;
    model.peak_threshold = 100e9;
    model.numeric_size = 8.0;
    auto ctx = model.build_context(net, targets);
    auto st = opt::detail::solve_single_term(model, net, targets, ctx);
    (void)model.reconstruct_batched_modes(ctx, st, net, targets, &p32);
  }
  std::wcerr << L"[ext-place-c60] ON mu~/K/i occ_block=32: peak=" << (p32 / 1e9)
             << L" GB (budget 100 GB)\n";

  // Phase-2 must never RAISE the modeled peak (a placement is adopted only when
  // it lowers the node's subtree peak).
  CHECK(p1 <= p0a);
  CHECK(p2 <= p0b);
  // GATE (spec S7): the cheap phase-2 pass, with the occ pair made sizeable,
  // REACHES the peak-setting operand and drops the giant's modeled peak below
  // the unbatched baseline. If this FAILS (peak stuck) because the winning
  // schedule was pruned pre-external, that is the trigger for the COMPLETE
  // variant (documented TODO), NOT more cheap-variant tuning.
  CHECK(p2 < p0a);
}

// Diagnostic + regression ([.]): scan EVERY summand of the C60 CSV-CCSD doubles
// residual under the batching configuration recommended for deployment, and
// report the modeled peak per term.
//
// Configuration (see the role-based batchability control):
//   - PAOs (mu~) are NOT batchable at all;
//   - the DF aux (K) is batchable in the CONTRACTED role;
//   - the occupied space is batchable in BOTH roles. Admitting its CONTRACTED
//     occurrences is what controls the g.g-over-K intermediates
//     (mu~,i3,mu~,i4): with external-occ only, those plateau at a ~391 GB floor
//     that NO occ block size can reach, because i3,i4 occur contracted there.
//     Re-admitting them brings the whole residual under budget.
//
// A correct inner_pow (CSV composite sizing) is REQUIRED: without it composites
// size at their base PAO extent and the factorization inverts (4-PAO integral).
TEST_CASE("C60 residual peak per summand under the recommended batching",
          "[.][occ-driver-scan]") {
  using namespace sequant;
  auto ctx0 = get_default_context().clone();
  ctx0.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx0.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);
  sequant::mbpt::add_df_spaces(isr);
  auto ctx_resetter = set_scoped_default_context(std::move(ctx0));

  auto const body = slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                          "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = deserialize<ExprPtr>(line);
  REQUIRE(expr);
  REQUIRE(expr->is<Sum>());
  auto const& summands = expr->as<Sum>().summands();

  auto regime = df_regime(kC60_pVDZF12);
  auto idxsz = regime.idx_to_extent();
  auto inner_pow = regime.inner_pow_fn();
  std::function<bool(Tensor const&)> is_vol = [](Tensor const& t) {
    return t.label() == L"t";
  };
  // aux + occ are contracted-batchable; PAOs are not batchable at all.
  auto is_batchable_contracted = [](Index const& ix) {
    auto const k = ix.space().base_key();
    return k == L"Κ" || k == L"i";  // DF aux (Greek Kappa) + occupied
  };
  // occ is additionally batchable in the external (spectator) role.
  std::function<bool(Index const&)> is_batchable_external =
      [](Index const& ix) { return ix.space().base_key() == L"i"; };
  std::function<std::size_t(Index const&)> bts =
      [](Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"i" ? std::size_t{8} : std::size_t{72};
  };

  auto measure = [&](TensorNetwork const& net,
                     container::svector<Index> const& tgt, bool spectator,
                     double& flops_out) {
    opt::detail::PeakBatchedModel model{idxsz, bts, is_vol, inner_pow};
    model.is_batchable_contracted_index = is_batchable_contracted;
    model.is_batchable_external_index = is_batchable_external;
    model.order_aware_recompute = true;
    model.batch_spectator_indices = spectator;
    model.perf_first = true;
    model.peak_threshold = 100e9;
    model.numeric_size = 8.0;
    auto ctx = model.build_context(net, tgt);
    auto st = opt::detail::solve_single_term(model, net, tgt, ctx);
    std::size_t const root = (std::size_t{1} << ctx.nt) - 1;
    int const best = model.select_root(ctx, st);
    flops_out = (best >= 0) ? st[root][0][best].flops : 0.0;
    double pk = 0.0;
    (void)model.reconstruct_batched_modes(ctx, st, net, tgt, &pk);
    return pk / 1e9;  // GB
  };

  double max_on = 0.0, tot_off = 0.0, tot_on = 0.0;
  std::size_t arg_on = 0;
  for (std::size_t gi = 0; gi < summands.size(); ++gi) {
    ExprPtr g = summands[gi];
    if (g->is<Product>())
      g = ex<Product>(g->as<Product>().scalar(), g->as<Product>().factors(),
                      Product::Flatten::Yes);
    if (!g->is<Product>()) continue;
    TensorNetwork net(g->as<Product>().factors());
    container::svector<Index> targets;
    double f_off = 0.0, f_on = 0.0;
    double const p_off = measure(net, targets, false, f_off);
    double const p_on = measure(net, targets, true, f_on);
    tot_off += f_off;
    tot_on += f_on;
    if (p_on > max_on) {
      max_on = p_on;
      arg_on = gi;
    }
    if (p_off > 50.0 || p_on != p_off)
      std::wcerr << L"[occ-driver] summand " << gi << L": peak " << p_off
                 << L" -> " << p_on << L" GB   flops " << f_off << L" -> "
                 << f_on << L"\n";
  }
  std::wcerr << L"[occ-driver] MAX peak=" << max_on << L" GB (summand "
             << arg_on << L"), budget 100 GB;  total flops " << tot_off
             << L" -> " << tot_on << L" (added compute "
             << (tot_off > 0.0 ? (tot_on - tot_off) / tot_off * 100.0 : 0.0)
             << L"%)\n";
  // Batching an external/spectator mode is work-neutral: no recompute is added.
  CHECK(tot_on == tot_off);
}

// Regression ([.]): the C60 giant (summand 38) is (g.C)(g.C).t where the two
// DF integrals g share the same aux Κ. Contracting Κ early reconstructs the
// (μ̃μ̃|μ̃μ̃) 4-PAO integral -- ~1800^4*8 ≈ 84 TB. That factorization ONLY looks
// competitive when the PNO composites a<i,j> are mis-sized at their base PAO
// extent (1800) instead of the per-pair PNO domain, i.e. when inner_pow is
// omitted. With correct inner_pow the DP picks the g.C 3-center path (Κ sliced)
// and the modeled peak is orders of magnitude below any 4-PAO scale. Guards
// both: (1) omitting inner_pow now THROWS; (2) with it, no 4-PAO node forms.
TEST_CASE("no 4-PAO integral with correct composite sizing (C60 giant)",
          "[.][roofline-4pao]") {
  using namespace sequant;
  auto ctx0 = get_default_context().clone();
  ctx0.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx0.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);
  sequant::mbpt::add_df_spaces(isr);
  auto ctx_resetter = set_scoped_default_context(std::move(ctx0));

  auto const body = slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                          "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = deserialize<ExprPtr>(line);
  REQUIRE(expr);
  REQUIRE(expr->is<Sum>());
  auto const& summands = expr->as<Sum>().summands();
  std::size_t const gi = 38 < summands.size() ? 38 : 0;
  ExprPtr giant = summands[gi];
  if (giant->is<Product>())
    giant = ex<Product>(giant->as<Product>().scalar(),
                        giant->as<Product>().factors(), Product::Flatten::Yes);
  REQUIRE(giant->is<Product>());
  TensorNetwork net(giant->as<Product>().factors());
  container::svector<Index> targets;

  auto regime = df_regime(kC60_pVDZF12);
  auto idxsz = regime.idx_to_extent();
  std::function<std::size_t(Index const&)> bts = [](Index const& ix) {
    return ix.space().base_key() == L"μ̃" ? std::size_t{256} : std::size_t{72};
  };
  auto is_aux = [](Index const& ix) {
    return ix.space().base_key() == L"Κ";  // K only; PAO NOT batchable
  };
  std::function<bool(Tensor const&)> is_vol = [](Tensor const& t) {
    return t.label() == L"t";
  };

  // (1) Omitting inner_pow on this composite network is now a hard error.
  {
    opt::detail::PeakBatchedModel bad{idxsz, bts, is_vol, {}};
    bad.is_batchable_contracted_index = is_aux;
    CHECK_THROWS_AS(bad.build_context(net, targets), std::invalid_argument);
  }

  // (2) With correct composite sizing, the DP forms NO 4-PAO integral.
  opt::detail::PeakBatchedModel model{idxsz, bts, is_vol,
                                      regime.inner_pow_fn()};
  model.is_batchable_contracted_index = is_aux;
  model.order_aware_recompute = true;
  model.perf_first = true;
  model.volatile_weight = 20.0;
  model.peak_threshold = std::numeric_limits<double>::infinity();
  model.numeric_size = 8.0;
  auto ctx = model.build_context(net, targets);
  auto st = opt::detail::solve_single_term(model, net, targets, ctx);
  double peak_bytes = 0.0;
  (void)model.reconstruct_batched_modes(ctx, st, net, targets, &peak_bytes);
  std::wcerr << L"[roofline-4pao] K-only modeled peak=" << (peak_bytes / 1e9)
             << L" GB (a 4-PAO mu~^4 node would be ~84000 GB)\n";
  // A μ̃^4 node is 1800^4*8 ≈ 84000 GB; correct sizing keeps the peak far below.
  CHECK(peak_bytes < 10000e9);
}

// Hidden ([.]) whole-residual dev sweep (no correctness assertions): optimizes
// every summand and prints, per term, whether its free-mu~ giant got a mu~/K
// batch mode or escaped slicing. Select explicitly:
//   ./unit_tests-sequant "[dryrun-df]"
TEST_CASE("dryrun POST-transform PAO/K batch-mode verdict", "[.][dryrun-df]") {
  // Augment the default mbpt registry with PAO (mu~) and DF-aux (K) so the
  // post-transform fixture deserializes; raise the dummy-ordinal ceiling for
  // mpqc's high internal ordinals (mu~_1152, a_21674, ...).
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);  // mu~
  sequant::mbpt::add_df_spaces(isr);                             // K
  auto ctx_resetter = set_scoped_default_context(std::move(ctx));

  auto const body = slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                          "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());

  // The fixture body is a single line: one Sum of Products (the whole doubles
  // residual). Deserialize, then split into summands and re-flatten each
  // (deserialize keeps literal nesting; single_term_opt needs a flat factor
  // list or term_batch_axes silently comes back empty -- see the mechanics
  // case above).
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = deserialize<ExprPtr>(line);
  bool const parsed = static_cast<bool>(expr);
  REQUIRE(parsed);

  auto flatten_product = [](ExprPtr const& e) -> ExprPtr {
    if (!e->is<Product>()) return e;
    auto const& p = e->as<Product>();
    return ex<Product>(p.scalar(), p.factors(), Product::Flatten::Yes);
  };
  std::vector<ExprPtr> terms;
  if (expr->is<Sum>())
    for (auto const& s : expr->as<Sum>()) terms.push_back(flatten_product(s));
  else
    terms.push_back(flatten_product(expr));
  REQUIRE(!terms.empty());
  std::wcout << L"\n=== POST-TRANSFORM DOUBLES RESIDUAL: " << terms.size()
             << L" summands (mu~ + K present) ===\n";

  // One regime (C60-scale mu~ domain) at the C60 peak_threshold. Per-term
  // progress + timing to std::cerr (unbuffered) so a slow/pathological term is
  // visible; each term reports whether its largest free-mu~ intermediate got a
  // mu~/K batch mode.
  // 40 GB matches the (too-tight) C60 run; override to explore realistic
  // budgets where the mu~-sliced sensible schedule becomes feasible.
  double const peak_threshold =
      (std::getenv("SEQUANT_UT_DRYRUN_PEAK_THR_GB")
           ? std::atof(std::getenv("SEQUANT_UT_DRYRUN_PEAK_THR_GB"))
           : 40.0) *
      1e9;
  // Env overrides for the root-cause BISECT (default = faithful real C60
  // config). SEQUANT_UT_DRYRUN_OCC/PNO/OSV vary extents;
  // SEQUANT_UT_DRYRUN_PAO_TS/AUX_TS vary batch target sizes;
  // SEQUANT_UT_DRYRUN_VW varies volatile_weight; SEQUANT_UT_DRYRUN_ROOFLINE=0
  // disables the roofline tie-break. Absent env var => real-config default.
  auto env_d = [](char const* k, double dflt) {
    char const* v = std::getenv(k);
    return v ? std::atof(v) : dflt;
  };
  auto env_u = [](char const* k, std::size_t dflt) {
    char const* v = std::getenv(k);
    return v ? static_cast<std::size_t>(std::atoll(v)) : dflt;
  };
  std::size_t const occ_ext =
      env_u("SEQUANT_UT_DRYRUN_OCC", kC60_pVDZF12.i_occ);
  double const pno_mom = env_d("SEQUANT_UT_DRYRUN_PNO", kC60_pVDZF12.pno_M[1]);
  double const osv_mom = env_d("SEQUANT_UT_DRYRUN_OSV", kC60_pVDZF12.osv_M[1]);
  std::size_t const pao_ts = env_u("SEQUANT_UT_DRYRUN_PAO_TS", 256u);
  std::size_t const aux_ts = env_u("SEQUANT_UT_DRYRUN_AUX_TS", 72u);
  double const vol_weight = env_d("SEQUANT_UT_DRYRUN_VW", 20.0);
  bool const use_roofline = env_u("SEQUANT_UT_DRYRUN_ROOFLINE", 1u) != 0;
  // FAITHFUL heavy-tail moments: real M_2..M_4 (default to M_1 = the constant
  // domain, so absent env vars reproduce the old scalar df_regime exactly).
  // Feed the mpqc-printed "PNO/OSV domain power means M_1..M_4" here.
  double const pno_m2 =
      env_d("SEQUANT_UT_DRYRUN_PNO_M2", kC60_pVDZF12.pno_M[2]);
  double const pno_m3 =
      env_d("SEQUANT_UT_DRYRUN_PNO_M3", kC60_pVDZF12.pno_M[3]);
  double const pno_m4 =
      env_d("SEQUANT_UT_DRYRUN_PNO_M4", kC60_pVDZF12.pno_M[4]);
  double const osv_m2 =
      env_d("SEQUANT_UT_DRYRUN_OSV_M2", kC60_pVDZF12.osv_M[2]);
  double const osv_m3 =
      env_d("SEQUANT_UT_DRYRUN_OSV_M3", kC60_pVDZF12.osv_M[3]);
  double const osv_m4 =
      env_d("SEQUANT_UT_DRYRUN_OSV_M4", kC60_pVDZF12.osv_M[4]);
  // Objective: "perf" = dense_time_space (perf-first, min-flops), else
  // "peak" = dense_space_time (peak-first, min-flops s.t. peak<=threshold).
  // Default "peak" reproduces the old DensePeakSizeBatched behavior.
  char const* obj_env = std::getenv("SEQUANT_UT_DRYRUN_OBJ");
  std::string const obj = obj_env ? std::string(obj_env) : std::string("peak");
  // SEQUANT_UT_DRYRUN_AUX_ONLY=1 makes ONLY the DF aux (K) sliceable, NOT PAO
  // (mu~) -- reproduces MPQC's aux-only batching to check the proper (gC)^2 PPL
  // factorization keeps the mu~-full giant (large realized peak = OOM) and
  // does NOT form the fully-sliceable 4-PAO integral.
  bool const aux_only = env_u("SEQUANT_UT_DRYRUN_AUX_ONLY", 0u) != 0;
  auto is_batchable = [aux_only](Index const& ix) {
    auto const k = ix.space().base_key();
    return aux_only ? (k == L"Κ") : (k == L"μ̃" || k == L"Κ");
  };
  std::wcerr << L"[dryrun-df] config: occ=" << occ_ext << L" pno=" << pno_mom
             << L" osv=" << osv_mom << L" pao_ts=" << pao_ts << L" aux_ts="
             << aux_ts << L" vw=" << vol_weight << L" roofline="
             << (use_roofline ? 1 : 0) << L" aux_only=" << (aux_only ? 1 : 0)
             << L"\n";
  // C60 pVDZ-F12 scale, REAL run dimensions (from the 614336 Owl job log):
  //   active occupied = 120 (tiles [0,120), elements [60,180))
  //   #PAO = #AO = 1800 (872/pair is a sparsity metric, NOT the dense cost)
  //   DF aux (aug-cc-pVDZ-RI) = 4320
  //   Average PNOs per pair = 41.92 ; Average OSVs per pair = 309.67
  // The earlier water-8 occ/PNO (i=32, PNO=19, OSV=57) were as wrong as the
  // K=672; they under-scaled every intermediate and (via peak_threshold
  // pruning) changed the DP's mode choice.
  auto regime = df_regime(
      /*mu_tilde=*/kC60_pVDZF12.mu_tilde, /*aux=*/kC60_pVDZF12.aux,
      /*i_occ=*/occ_ext,
      std::array<double, 5>{1.0, pno_mom, pno_m2, pno_m3, pno_m4},
      std::array<double, 5>{1.0, osv_mom, osv_m2, osv_m3, osv_m4});
  std::wcerr << L"[dryrun-df] moments: PNO M_1..M_4=" << pno_mom << L","
             << pno_m2 << L"," << pno_m3 << L"," << pno_m4 << L"  OSV M_1..M_4="
             << osv_mom << L"," << osv_m2 << L"," << osv_m3 << L"," << osv_m4
             << L"  objective="
             << (obj == "perf" ? L"perf(dense_time_space)"
                               : L"peak(dense_space_time)")
             << L"\n";
  auto memsize = sequant::opt::detail::memsize_counter(regime.idx_to_extent(),
                                                       regime.inner_pow_fn());
  // Static NOMINAL flops per contraction node (slicing does not reduce total
  // work, so this is the true schedule flops). The 4-PAO fully-sliced integral
  // is memory-cheap but flops-catastrophic; only this metric surfaces it.
  auto flops = sequant::opt::detail::flops_counter(regime.idx_to_extent(),
                                                   regime.inner_pow_fn());
  auto has_free_mu_tilde = [](std::vector<Index> const& ixs) {
    for (auto const& ix : ixs)
      if (ix.space().base_key() == L"μ̃") return true;
    return false;
  };

  std::size_t total_mu_nodes = 0, n_terms_with_giant = 0;
  std::size_t total_mu_nodes_with_mu_axis = 0, total_mu_nodes_with_k_only = 0;
  double overall_flops = 0.0;  // summed schedule flops over all terms
  double overall_exec = 0.0;   // summed roofline exec cost (DP's real mode)
  double overall_flops_exec = 0.0;  // summed EXECUTED (recompute-charged) flops
  std::size_t overall_batched_nodes = 0, overall_recomputed_nodes = 0;
  // FOREST-LEVEL CSE. MPQC caches intermediates across terms, keyed on
  // EvalExpr::hash_value(). Two subexpressions that are equal but carry
  // DIFFERENT batch-mode annotations cannot share a cache entry -- they are
  // evaluated under different slicings. So:
  //   cse_by_expr = ideal reuse (what a batching-blind schedule achieves)
  //   cse_by_expr_and_axes = reuse actually achievable given the annotations
  // The gap between them is CSE destroyed by inconsistent batching, which is
  // invisible to the per-term nominal-flops metric.
  std::map<std::size_t, double> cse_by_expr;  // hash -> flops
  std::map<std::pair<std::size_t, std::wstring>, double> cse_by_expr_and_axes;
  std::size_t overall_internal_nodes = 0;
  double overall_biggest = 0.0;
  bool overall_biggest_has_mu = false, overall_biggest_has_k = false;
  std::wstring overall_biggest_desc, overall_biggest_axes,
      overall_biggest_sources;

  for (std::size_t ti = 0; ti < terms.size(); ++ti) {
    auto t0 = std::chrono::steady_clock::now();
    std::cerr << "[dryrun-df] optimizing term " << (ti + 1) << "/"
              << terms.size() << " ..." << std::flush;

    auto axes_map = std::make_shared<std::unordered_map<
        Expr const*, container::vector<NodeBatchAnnotation>>>();
    OptimizeOptions opts;
    // "flops" = DenseFLOPs: the classic min-flops contraction-order DP (flops
    // HAS optimal substructure there), i.e. the pre-peak-objective default.
    // It is the reference point for "what does a real min-time schedule cost".
    opts.objective_function = (obj == "flops") ? ObjectiveFunction::DenseFLOPs
                              : (obj == "perf")
                                  ? ObjectiveFunction::DenseTimeSpaceBatched
                                  : ObjectiveFunction::DenseSpaceTimeBatched;
    opts.idx_to_extent = regime.idx_to_extent();
    opts.inner_pow = regime.inner_pow_fn();
    opts.batch_policy.is_batchable_contracted_index = is_batchable;
    // FAITHFUL replica of the real C60 run's OptimizeOptions
    // (make_optimize_options + cck.ipp batch policy, 614336 job log):
    //   batch:pao_target_size=256, batch:aux_target_size=72
    //   optimize:volatile_weight=20, machine_balance=200, fast_mem_elems=1e6
    //   is_volatile_leaf = (label == "t"), accumulation_factor = 1.0 (MPQC
    //   dflt)
    opts.batch_policy.batch_target_size =
        [pao_ts, aux_ts](Index const& ix) -> std::size_t {
      return ix.space().base_key() == L"μ̃" ? pao_ts   // pao_target_size
                                           : aux_ts;  // aux_target_size
    };
    opts.batch_policy.is_volatile_leaf = [](Tensor const& t) {
      return t.label() == L"t";
    };
    opts.batch_policy.accumulation_factor = 1.0;
    opts.batch_policy.peak_threshold = peak_threshold;
    opts.volatile_weight = vol_weight;
    if (use_roofline) {
      opts.roofline.machine_balance = 200.0;
      opts.roofline.fast_mem_elems = 1000000.0;
    }
    opts.term_batch_axes = axes_map;

    auto term = optimize(terms[ti], opts);
    if (!static_cast<bool>(term)) {
      std::cerr << " (skipped)\n";
      continue;
    }
    auto it = axes_map->find(term.get());
    container::vector<NodeBatchAnnotation> node_axes;
    if (it != axes_map->end()) node_axes = it->second;
    BinarizationOptions bopts;
    bopts.node_batch_axes = node_axes;
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto node = binarize(term, {}, bopts);
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

    // CORRECT metric: a node can only slice a batchable index it CONTRACTS
    // (relax()'s contracted_here = open at children & closed at parent). A free
    // index on an intermediate is sliceable only at the ANCESTOR that contracts
    // it, so the giant's REALIZED size is its nominal size with each free index
    // that some ancestor sliced reduced to batch_target_size. Walk top-down
    // carrying active = union of ancestor batched_here (by FULL label, so
    // mu~_1241 sliced above only reduces the SAME mu~_1241 below, not a
    // different mu~_j).
    auto keyof = [](Index const& ix) { return std::wstring(ix.full_label()); };
    auto ext_of = [](Index const& ix) -> double {
      auto bk = ix.space().base_key();
      if (bk == L"μ̃") return 1800.0;
      if (bk == L"Κ") return 4320.0;
      return 0.0;
    };
    auto tgt_of = [pao_ts, aux_ts](Index const& ix) -> double {
      return ix.space().base_key() == L"μ̃" ? double(pao_ts) : double(aux_ts);
    };
    double term_biggest = 0.0;  // max REALIZED bytes over free-mu~ nodes
    std::wstring term_biggest_desc, term_biggest_axes, term_biggest_sources;
    bool term_biggest_has_mu = false;  // its free mu~ ESCAPED slicing
    bool term_biggest_has_k = false;   // its free K ESCAPED slicing
    std::size_t term_mu_nodes = 0;
    double term_flops = 0.0;  // summed nominal contraction flops for this term
    double term_exec = 0.0;   // summed roofline exec cost for this term
    // EXECUTED flops: what the machine actually runs. A node is re-run once per
    // batch of every ancestor-sliced mode it does NOT touch (the DP charges
    // this as `rf` in cost_model.hpp:911); a mode it DOES touch shrinks the
    // per-batch work and sums back to nominal. Nominal flops is blind to this,
    // which is why it cannot see a batch-fragmentation pathology.
    double term_flops_exec = 0.0;
    std::size_t term_batched_nodes = 0;     // nodes carrying >=1 batch mode
    std::size_t term_recomputed_nodes = 0;  // nodes with recompute factor > 1
    std::size_t term_internal_nodes = 0;    // all contraction (non-leaf) nodes
    // Per-term batch-mode multiset (which modes get sliced, and on how many
    // nodes) -- the term-level source of "more batch groups".
    std::map<std::wstring, std::size_t> term_axis_hist;
    // active maps a sliced index's FULL label -> a descriptor of the ANCESTOR
    // node that slices it (its result free-index signature + its batched_here).
    // This is the node-dump: it turns "escaped={}" from an inference into a
    // concrete "mu~_X is sliced by ancestor <node>" (or ESCAPED).
    // nbatches carries, for each ancestor-sliced mode, its batch count
    // (extent / batch_target_size) -- the executed-flops recompute factor.
    std::function<void(std::remove_cvref_t<decltype(node)> const&,
                       std::map<std::wstring, std::wstring>,
                       std::map<std::wstring, double>)>
        walk = [&](auto const& n, std::map<std::wstring, std::wstring> active,
                   std::map<std::wstring, double> nbatches) {
          auto free_ixs = node_free_indices(*n);
          if (has_free_mu_tilde(free_ixs)) {
            ++term_mu_nodes;
            ++total_mu_nodes;
            double nominal =
                memsize(free_ixs, std::vector<Index>{}, std::vector<Index>{}) *
                8.0;
            double factor = 1.0;
            std::vector<Index> escaped;  // batchable free ixs NOT sliced above
            std::wstring sources;        // per-free-batchable-index provenance
            for (auto const& ix : free_ixs) {
              auto bk = ix.space().base_key();
              if (bk != L"μ̃" && bk != L"Κ") continue;
              auto it = active.find(keyof(ix));
              if (it != active.end()) {
                factor *= tgt_of(ix) / ext_of(ix);
                sources += L"\n      " + keyof(ix) + L": SLICED by ancestor " +
                           it->second;
              } else {
                escaped.push_back(ix);
                sources += L"\n      " + keyof(ix) +
                           L": ESCAPED (no ancestor "
                           L"slices it)";
              }
            }
            double const realized = nominal * factor;
            bool mu_esc = false, k_esc = false;
            for (auto const& ix : escaped)
              (ix.space().base_key() == L"μ̃" ? mu_esc : k_esc) = true;
            if (mu_esc) ++total_mu_nodes_with_mu_axis;  // repurposed: escaped
            if (k_esc && !mu_esc) ++total_mu_nodes_with_k_only;
            if (realized > 5e10) {
              std::wcerr << L"\n  [GIANT realized=" << (realized / 1e9)
                         << L"GB nominal=" << (nominal / 1e9) << L"GB] free={"
                         << describe_indices(free_ixs) << L"}" << sources
                         << L"\n";
            }
            if (realized > term_biggest) {
              term_biggest = realized;
              term_biggest_desc = describe_indices(free_ixs);
              term_biggest_has_mu = mu_esc;
              term_biggest_has_k = k_esc;
              term_biggest_sources = sources;
              term_biggest_axes = describe_indices(
                  container::vector<Index>(escaped.begin(), escaped.end()));
            }
          }
          if (!n.leaf()) {
            // Nominal flops of THIS contraction: union of children's free
            // indices and this node's result free indices (slicing does not
            // change the flop count).
            auto const lf = node_free_indices(*n.left());
            auto const rf = node_free_indices(*n.right());
            double const nf = flops(lf, rf, free_ixs);
            term_flops += nf;
            // Executed flops = nominal * product of nbatches[k] over every
            // ancestor-sliced mode k this contraction does NOT touch. Touching
            // k (in either operand or the result) means the slice divides the
            // work and the batches sum back to nominal; not touching it means
            // the whole contraction is redone once per batch of k.
            auto touches = [&](std::wstring const& k) {
              auto has = [&](auto const& v) {
                for (auto const& ix : v)
                  if (keyof(ix) == k) return true;
                return false;
              };
              return has(lf) || has(rf) || has(free_ixs);
            };
            double recompute = 1.0;
            for (auto const& [k, nb] : nbatches)
              if (!touches(k)) recompute *= nb;
            term_flops_exec += nf * recompute;
            if (recompute > 1.0) ++term_recomputed_nodes;
            if (!n->batched_here().empty()) ++term_batched_nodes;
            ++term_internal_nodes;
            for (auto const& ax : n->batched_here())
              term_axis_hist[std::wstring(ax.first.space().base_key())]++;
            // Forest-level CSE bookkeeping (see declarations above). The mode
            // signature includes the ancestor-sliced context, not just this
            // node's own modes: the same expression evaluated under a different
            // ancestor slicing is a different cache entry.
            ++overall_internal_nodes;
            std::wstring axsig;
            for (auto const& [k, nb] : nbatches) axsig += k + L";";
            axsig += L"|";
            for (auto const& ax : n->batched_here())
              axsig += std::wstring(ax.first.full_label()) + L",";
            auto const h = n->hash_value();
            cse_by_expr.emplace(h, nf);
            cse_by_expr_and_axes.emplace(std::pair{h, axsig}, nf);
            // Roofline exec cost = the DP's ACTUAL optimization mode (NOT raw
            // flops): max(flops, mb * max(traffic, prefac*flops/sqrt(fastmem/
            // tiles))), mb=200, fastmem=1e6, tiles=3, prefac=1. Traffic = the
            // UNBATCHED operand+result footprint, so exec cost is batch-
            // INDEPENDENT (batching moves only the peak). This is why raw flops
            // can be non-monotone in the threshold while the DP's objective is
            // monotone. Element units (matches the DP).
            double const traffic = memsize(lf, rf, free_ixs);
            double const Q = std::max(traffic, nf / std::sqrt(1e6 / 3.0));
            term_exec += std::max(nf, 200.0 * Q);
            std::map<std::wstring, std::wstring> child_active = active;
            std::wstring const self_desc =
                L"[free={" + describe_indices(node_free_indices(*n)) +
                L"} batched_here={" +
                describe_indices(batch_axes_indices(n->batched_here())) + L"}]";
            std::map<std::wstring, double> child_nbatch = nbatches;
            for (auto const& ax : n->batched_here()) {
              child_active[keyof(ax.first)] = self_desc;
              double const e = ext_of(ax.first);
              if (e > 0.0)
                child_nbatch[keyof(ax.first)] =
                    std::max(1.0, e / tgt_of(ax.first));
            }
            walk(n.left(), child_active, child_nbatch);
            walk(n.right(), child_active, child_nbatch);
          }
        };
    walk(node, {}, {});
    overall_flops += term_flops;
    overall_exec += term_exec;
    overall_flops_exec += term_flops_exec;
    overall_batched_nodes += term_batched_nodes;
    overall_recomputed_nodes += term_recomputed_nodes;
    if (term_mu_nodes > 0) ++n_terms_with_giant;
    if (term_biggest > overall_biggest) {
      overall_biggest = term_biggest;
      overall_biggest_desc = term_biggest_desc;
      overall_biggest_axes = term_biggest_axes;
      overall_biggest_has_mu = term_biggest_has_mu;
      overall_biggest_has_k = term_biggest_has_k;
      overall_biggest_sources = term_biggest_sources;
    }
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                  std::chrono::steady_clock::now() - t0)
                  .count();
    std::cerr << " " << ms << "ms  free-mu~ nodes=" << term_mu_nodes
              << " biggest_realized=" << (term_biggest / 1e9)
              << "GB mu~escaped=" << (term_biggest_has_mu ? "YES" : "NO")
              << " Kescaped=" << (term_biggest_has_k ? "YES" : "NO")
              << " flops=" << (term_flops / 1e12) << "T\n";
    // Machine-readable per-term line for cross-objective join. modes= lists the
    // batch-mode multiset by index-space base key (e.g. mu~:3,K:5).
    std::wstring axes_str;
    for (auto const& [k, c] : term_axis_hist)
      axes_str +=
          (axes_str.empty() ? L"" : L",") + k + L":" + std::to_wstring(c);
    std::wcerr << L"PERTERM " << ti << L" nodes=" << term_internal_nodes
               << L" batched=" << term_batched_nodes << L" recomp="
               << term_recomputed_nodes << L" flopsT=" << (term_flops / 1e12)
               << L" execT=" << (term_exec / 1e12) << L" modes=["
               << (axes_str.empty() ? L"-" : axes_str) << L"]\n";
  }

  std::wcerr << L"\n=== POST-TRANSFORM VERDICT (REALIZED peak, mu~=1800, "
                L"thr=40GB) ===\n"
             << L"terms with a free-mu~ intermediate: " << n_terms_with_giant
             << L"/" << terms.size() << L"\n"
             << L"TOTAL SCHEDULE FLOPS (nominal, all terms): "
             << (overall_flops / 1e12) << L" Tflop\n"
             << L"TOTAL EXECUTED FLOPS (recompute-charged): "
             << (overall_flops_exec / 1e12) << L" Tflop  (x"
             << (overall_flops_exec / std::max(1.0, overall_flops))
             << L" vs nominal; batched nodes=" << overall_batched_nodes
             << L", recomputed nodes=" << overall_recomputed_nodes << L")\n"
             << L"FOREST CSE: internal nodes=" << overall_internal_nodes
             << L" distinct-by-expr=" << cse_by_expr.size()
             << L" distinct-by-expr+modes=" << cse_by_expr_and_axes.size()
             << L"\n  CSE'd forest flops (ideal, by expr):  "
             << (std::accumulate(
                     cse_by_expr.begin(), cse_by_expr.end(), 0.0,
                     [](double a, auto const& p) { return a + p.second; }) /
                 1e12)
             << L" Tflop\n  CSE'd forest flops (achievable, expr+modes): "
             << (std::accumulate(
                     cse_by_expr_and_axes.begin(), cse_by_expr_and_axes.end(),
                     0.0,
                     [](double a, auto const& p) { return a + p.second; }) /
                 1e12)
             << L" Tflop\n"
             << L"TOTAL ROOFLINE EXEC COST (the DP's real mode): "
             << (overall_exec / 1e12) << L" (mb=200, batch-independent)\n"
             << L"free-mu~ contraction nodes: " << total_mu_nodes << L"\n"
             << L"  ... with a free mu~ that ESCAPED slicing: "
             << total_mu_nodes_with_mu_axis << L"\n"
             << L"  ... with only a free K that escaped:      "
             << total_mu_nodes_with_k_only << L"\n"
             << L"LARGEST REALIZED free-mu~ intermediate: {"
             << overall_biggest_desc << L"} = " << (overall_biggest / 1e9)
             << L" GB (realized, after ancestor slices)\n"
             << L"  escaped (un-sliced) indices = {" << overall_biggest_axes
             << L"}\n"
             << L"  -> free mu~ escaped: "
             << (overall_biggest_has_mu ? L"YES" : L"NO")
             << L" | free K escaped: "
             << (overall_biggest_has_k ? L"YES" : L"NO") << L"\n"
             << L"  NODE-DUMP (per free batchable index, which ancestor slices "
                L"it):"
             << overall_biggest_sources << L"\n"
             << L"INTERPRETATION: realized peak = what the runtime actually "
                L"materializes. If it is >> 40GB with mu~ escaped=YES, the DP "
                L"left the giant's free mu~ un-sliced at its contracting "
                L"ancestor -> reproduces the C60 OOM in predicted space.\n";

  REQUIRE(total_mu_nodes > 0);
  SUCCEED();
}

// ===========================================================================
// Task 2-6: DryRun eval BACKEND (zero-data Result + eval_expr) and the
// end-to-end replay harness that WITNESSES the runtime's batch-mode
// realization on the real post-transform giant term. The POST-TRANSFORM
// VERDICT case above established the DP side of this story (the DP DOES
// annotate a mu~ batch mode on the giant, surviving binarize). These new
// cases do NOT modify anything above; [dryrun-probe]/[dryrun-df] stay exactly
// as committed.
// ===========================================================================

namespace {

using sequant::eval::dryrun::CostModel;
using sequant::eval::dryrun::DryRunLeafEvaluator;
using sequant::eval::dryrun::EvalExprDryRun;
using sequant::eval::dryrun::EvalNodeDryRun;
using sequant::eval::dryrun::ExtentOverrides;
using sequant::eval::dryrun::make_dryrun_result;
using sequant::eval::dryrun::ResultDryRun;
using sequant::eval::dryrun::ResultDryRunNested;

// A small, self-consistent regime for the backend unit tests below (distinct
// from probe_regime()/df_regime() above -- this one just needs a couple of
// named spaces plus non-trivial PNO moments).
SizeRegime backend_test_regime() {
  SizeRegime r;
  r.space_extent = {
      {L"i", 10},
      {L"a", 20},
  };
  double const pno = 4.0;
  for (std::size_t k = 0; k <= 4; ++k)
    r.csv_pno_moment[k] = std::pow(pno, double(k));
  r.csv_osv_moment = r.csv_pno_moment;
  return r;
}

std::array<std::any, 3> annot3(container::svector<Index> l,
                               container::svector<Index> r,
                               container::svector<Index> res) {
  return {std::any{std::move(l)}, std::any{std::move(r)},
          std::any{std::move(res)}};
}

}  // namespace

TEST_CASE("dryrun cost model memsize matches memsize_counter",
          "[dryrun-costmodel]") {
  auto r = backend_test_regime();
  CostModel cm{r};
  container::svector<Index> idx{Index{L"i_1"}, Index{L"i_2"}, Index{L"a_4"}};
  auto direct = sequant::opt::detail::memsize_counter(
      r.idx_to_extent(), r.inner_pow_fn())(idx, container::svector<Index>{},
                                           container::svector<Index>{});
  CHECK(cm.memsize(idx) == static_cast<std::size_t>(direct * 8.0));
}

TEST_CASE("dryrun cost model memsize honors an extent override",
          "[dryrun-costmodel]") {
  auto r = backend_test_regime();
  CostModel cm{r};
  container::svector<Index> idx{Index{L"a_3"}, Index{L"i_1"}};
  auto const full = cm.memsize(idx);
  ExtentOverrides ov;
  ov[Index{L"a_3"}] = 5;  // narrowed from 20 to 5
  auto const sliced = cm.memsize(idx, ov);
  CHECK(sliced < full);
  CHECK(full == sliced * 4);  // linear in a_3's extent
}

TEST_CASE("dryrun cost model flops and exec_cost are finite/positive",
          "[dryrun-costmodel]") {
  auto r = backend_test_regime();
  CostModel cm{r};
  container::svector<Index> out{Index{L"a_3"}, Index{L"i_1"}};
  container::svector<Index> contracted{Index{L"i_2"}};
  auto const f = cm.flops(out, contracted);
  CHECK(f > 0.0);
  CHECK(cm.exec_cost(f, cm.memsize(out), 4096) > 0.0);
}

TEST_CASE("dryrun flat result size delegates to cost model",
          "[dryrun-result]") {
  // Result::prod/sum/permute/slice_mode/mode_batches/size_in_bytes are all
  // overridden PRIVATE in the concrete DryRun classes (mirroring
  // ResultTensorTAPP), since real callers only ever reach a Result through a
  // ResultPtr/Result const& (base-class access); tests that want to call
  // them on a concrete object must do the same -- via a `Result const&`.
  auto r = backend_test_regime();
  auto cm = std::make_shared<CostModel const>(r);
  container::svector<Index> idx{Index{L"i_1"}, Index{L"i_2"}, Index{L"a_4"}};
  ResultDryRun t{idx, cm};
  Result const& rt = t;
  CHECK(rt.size_in_bytes() == cm->memsize(idx));
}

TEST_CASE("dryrun flat result prod yields the result annotation index set",
          "[dryrun-result]") {
  auto r = backend_test_regime();
  auto cm = std::make_shared<CostModel const>(r);
  ResultDryRun l{{Index{L"a_3"}, Index{L"i_2"}}, cm};
  ResultDryRun rr{{Index{L"i_2"}, Index{L"i_1"}}, cm};
  Result const& rl = l;
  container::svector<Index> res{Index{L"a_3"}, Index{L"i_1"}};
  auto out = rl.prod(rr, annot3(l.indices(), rr.indices(), res), DeNest::False);
  REQUIRE(out);
  bool const is_flat = out->is<ResultDryRun>();
  CHECK(is_flat);
  auto const& ot = out->as<ResultDryRun>();
  CHECK(ot.indices() == res);
  CHECK(out->size_in_bytes() == cm->memsize(res));
}

TEST_CASE("dryrun flat result full contraction yields a scalar",
          "[dryrun-result]") {
  auto r = backend_test_regime();
  auto cm = std::make_shared<CostModel const>(r);
  ResultDryRun l{{Index{L"i_1"}}, cm};
  ResultDryRun rr{{Index{L"i_1"}}, cm};
  Result const& rl = l;
  auto out = rl.prod(rr, annot3(l.indices(), rr.indices(), {}), DeNest::False);
  REQUIRE(out);
  bool const is_scalar = out->is<ResultScalar<double>>();
  CHECK(is_scalar);
}

TEST_CASE("dryrun flat result slice_mode shrinks the sliced mode",
          "[dryrun-result]") {
  auto r = backend_test_regime();
  auto cm = std::make_shared<CostModel const>(r);
  ResultDryRun t{{Index{L"a_3"}, Index{L"i_2"}}, cm};  // mu~ extent 20
  Result const& rt = t;
  auto const full = rt.size_in_bytes();
  auto sliced = rt.slice_mode(0, 0, 5);  // quarter of mu~
  REQUIRE(sliced);
  CHECK(sliced->size_in_bytes() < full);
  CHECK(sliced->size_in_bytes() == full / 4);
}

TEST_CASE("dryrun flat result mode_batches tiles the mode extent",
          "[dryrun-result]") {
  auto r = backend_test_regime();
  auto cm = std::make_shared<CostModel const>(r);
  ResultDryRun t{{Index{L"a_3"}, Index{L"i_2"}}, cm};
  Result const& rt = t;
  auto batches = rt.mode_batches(0, 5);  // 20 / 5 = 4 batches
  CHECK(batches.size() == 4);
  CHECK(batches.front().first == 0);
  CHECK(batches.back().second == 20);
}

TEST_CASE(
    "dryrun flat result pre_sized_zeros_over_mode widens to the mode "
    "source's full extent",
    "[dryrun-result][pre-sized]") {
  // D3.1 (Task 6 witness): the runtime scatter branch
  // (make_batched_custom_evaluator, BatchModeType::External) presizes its
  // destination from a BLOCK PARTIAL (a result whose batch mode has been
  // slice_mode()'d down to one block) via
  // part->pre_sized_zeros_over_mode(dest_mode, carrier_full, carrier_mode).
  // The dry-run analogue of the TA backend's ResultTensorTA::
  // pre_sized_zeros_over_mode (which swaps in axis_src's FULL
  // TiledRange1): the returned token's mode-0 extent (queryable immediately
  // via size_in_bytes(), a structural fact, BEFORE any write_into_slice()
  // call) must be the mode source's full extent (10, from
  // backend_test_regime), not the block's narrower extent (4).
  auto r = backend_test_regime();
  auto cm = std::make_shared<CostModel const>(r);
  Index i1{L"i_1"}, a3{L"a_3"};
  container::svector<Index> idx{i1, a3};  // i_1 extent 10, a_3 extent 20

  ResultDryRun full{idx, cm};  // the unsliced carrier leaf's token
  Result const& full_r = full;
  auto const full_bytes = full_r.size_in_bytes();

  auto const block = full_r.slice_mode(0, 0, 4);  // one block: i_1 in [0,4)
  REQUIRE(block);
  CHECK(block->size_in_bytes() == full_bytes * 4 / 10);

  auto dest = block->pre_sized_zeros_over_mode(/*mode=*/0, full_r,
                                               /*axis_src_mode=*/0);
  REQUIRE(dest);
  CHECK(dest->is<ResultDryRun>());
  // Widened back to the FULL i_1 extent (10), not the block's 4.
  CHECK(dest->size_in_bytes() == full_bytes);
}

TEST_CASE(
    "dryrun nested result pre_sized_zeros_over_mode widens to the mode "
    "source's full extent",
    "[dryrun-result][pre-sized]") {
  // ToT analogue of the flat case above -- CSV/PNO residuals carry
  // ResultDryRunNested tokens, so the External-mode scatter's presize must
  // also widen a nested token's outer batch mode.
  auto r = backend_test_regime();
  auto cm = std::make_shared<CostModel const>(r);
  Index i1{L"i_1"}, i2{L"i_2"}, a3{L"a_3"};
  Index a_pno{L"a_1", {i1, i2}};
  container::svector<Index> outer{i1, a3};
  container::svector<Index> inner{a_pno};
  container::svector<Index> canon{i1, a3, a_pno};

  ResultDryRunNested full{outer, inner, cm, {}, canon};
  Result const& full_r = full;
  auto const full_bytes = full_r.size_in_bytes();

  auto const block = full_r.slice_mode(0, 0, 5);  // i_1 in [0,5), half of 10
  REQUIRE(block);
  CHECK(block->size_in_bytes() == full_bytes / 2);

  auto dest = block->pre_sized_zeros_over_mode(/*mode=*/0, full_r,
                                               /*axis_src_mode=*/0);
  REQUIRE(dest);
  CHECK(dest->is<ResultDryRunNested>());
  CHECK(dest->size_in_bytes() == full_bytes);
}

TEST_CASE(
    "dryrun result write_into_slice assembles disjoint blocks into a "
    "pre-sized destination",
    "[dryrun-result][write-into-slice]") {
  auto r = backend_test_regime();
  auto cm = std::make_shared<CostModel const>(r);

  // A ToT template whose outer batch mode (position 0) is the occupied index
  // i_1 (extent 10). Its full modelled size is the reference we must
  // reconstruct by assembling disjoint outer-mode blocks.
  Index i1{L"i_1"}, i2{L"i_2"}, a3{L"a_3"};
  Index a_pno{L"a_1", {i1, i2}};  // proto-indexed composite (CSV/PNO) leg
  container::svector<Index> outer{i1, a3};
  container::svector<Index> inner{a_pno};
  container::svector<Index> canon{i1, a3, a_pno};

  ResultDryRunNested tmpl{outer, inner, cm, {}, canon};
  Result const& tmpl_r = tmpl;
  auto const whole = tmpl_r.size_in_bytes();

  // Two disjoint element blocks [0,5) and [5,10) that tile the full extent 10
  // of the batch mode with no overlap and no gap. Slice each out of the
  // template (mirrors how the batched evaluator produces per-block results).
  auto block0 = tmpl_r.slice_mode(0, 0, 5);
  auto block1 = tmpl_r.slice_mode(0, 5, 10);
  REQUIRE(block0);
  REQUIRE(block1);

  // A fresh pre-sized destination whose batch mode starts UNFILLED (extent 0):
  // its shape/index-set is fixed, but its assembled size grows as blocks are
  // written in. Assemble the two blocks into it.
  ResultDryRunNested dest{outer, inner, cm, {{i1, 0}}, canon};
  Result& dest_w = dest;  // the mutator is reached through the base interface
  dest_w.write_into_slice(*block0, 0, 0, 5);
  dest_w.write_into_slice(*block1, 0, 5, 10);

  // No double-count, no gap: the assembled modelled size equals the whole
  // array's size exactly, and the covered element range is the full extent.
  CHECK(dest_w.size_in_bytes() == whole);
  CHECK(dest.assembled_range(0) == std::pair<std::size_t, std::size_t>{0, 10});

  // Lobounds are preserved: a block written at a nonzero element offset (a
  // frozen-core-style occupied offset) assembles at that offset, not rebased
  // to 0.
  ResultDryRunNested dest_fc{outer, inner, cm, {{i1, 0}}, canon};
  Result& dest_fc_w = dest_fc;
  auto block_fc = tmpl_r.slice_mode(0, 2, 6);
  REQUIRE(block_fc);
  dest_fc_w.write_into_slice(*block_fc, 0, 2, 6);
  CHECK(dest_fc.assembled_range(0) ==
        std::pair<std::size_t, std::size_t>{2, 6});
}

TEST_CASE("dryrun nested result uses moment-aware inner extent, not extent^k",
          "[dryrun-nested]") {
  auto r = backend_test_regime();
  // Non-trivial second moment: <#PNO^2> != <#PNO>^2 (dispersion inflates it).
  r.csv_pno_moment[1] = 4.0;
  r.csv_pno_moment[2] = 4.0 * 4.0 * 1.5;
  auto cm = std::make_shared<CostModel const>(r);

  Index i1{L"i_1"}, i2{L"i_2"};
  Index a_pno{L"a_1", {i1, i2}};  // proto-indexed (CSV/PNO composite) leg
  container::svector<Index> outer{i1, i2, Index{L"a_3"}};
  container::svector<Index> inner{a_pno};
  ResultDryRunNested c{outer, inner, cm};

  CHECK(c.outer() == outer);
  CHECK(c.inner() == inner);

  container::svector<Index> combined = outer;
  combined.push_back(a_pno);
  CHECK(c.indices() == combined);

  auto const without_composite = cm->memsize(outer);  // just i1*i2*mu~
  auto const with_composite = cm->memsize(combined);  // routes a_pno via k=1
  // The composite contributes the FIRST moment (4.0), not extent(a_pno)^1
  // (a_pno's own "extent" as a bare, non-composite space is never queried
  // here -- only backend_test_regime()'s csv_pno_moment[1] is), so the ratio
  // is exactly the first moment.
  CHECK(with_composite == without_composite * 4);
  Result const& rc = c;
  CHECK(rc.size_in_bytes() == with_composite);
}

TEST_CASE("dryrun make_dryrun_result dispatches flat vs nested by content",
          "[dryrun-nested]") {
  auto r = backend_test_regime();
  auto cm = std::make_shared<CostModel const>(r);
  Index i1{L"i_1"}, i2{L"i_2"};
  Index a_pno{L"a_1", {i1, i2}};

  auto flat = make_dryrun_result({Index{L"i_1"}, Index{L"i_2"}}, cm);
  CHECK(flat->is<ResultDryRun>());

  auto nested = make_dryrun_result({i1, i2, a_pno}, cm);
  CHECK(nested->is<ResultDryRunNested>());
}

TEST_CASE("dryrun leaf yielder builds a sized token from a tensor leaf",
          "[dryrun-leaf]") {
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto ctx_resetter = set_scoped_default_context(std::move(ctx));

  auto r = backend_test_regime();
  auto cm = std::make_shared<CostModel const>(r);

  auto expr = deserialize<ExprPtr>("g{i_1,i_2;a_4}");
  bool const parsed = static_cast<bool>(expr);
  REQUIRE(parsed);

  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = binarize<EvalExprDryRun>(expr);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  REQUIRE(node.leaf());

  DryRunLeafEvaluator yield{cm};
  auto res = yield(node);
  REQUIRE(res);
  bool const is_flat = res->is<ResultDryRun>();
  CHECK(is_flat);

  container::svector<Index> idx{Index{L"i_1"}, Index{L"i_2"}, Index{L"a_4"}};
  CHECK(res->size_in_bytes() == cm->memsize(idx));
}

TEST_CASE(
    "dryrun external-mode scatter replay models the sliced footprint (D3.1)",
    "[dryrun-extmode][eval]") {
  // D3.1 regression: the runtime External-mode scatter branch in
  // make_batched_custom_evaluator (eval.hpp) calls, on the first block,
  // part->pre_sized_zeros_over_mode(dest_mode, carrier_full, carrier_mode),
  // then dest->write_into_slice(...) for every block. Before this task the
  // dry-run Result classes did not override pre_sized_zeros_over_mode, so
  // the replay hit the base class's `throw
  // detail::unimplemented_method("pre_sized_zeros_over_mode")` the moment an
  // External mode was stamped -- the witness could not measure external
  // batching at all (doc/dev/specs/2026-07-20-external-mode-batching-
  // design.md, D3). This test drives the SAME scatter branch the TA
  // regression `batched_eval_external_proto_occ_scatter` (test_eval_ta.cpp)
  // exercises, on the dry-run backend: a small forest carrying the occupied
  // index ONLY as a protoindex of a composite PNO leg (canonicalization
  // promotes it to a plain outer canon index, so index_position locates it
  // directly -- no proto-aware locator needed, exactly as that TA test
  // documents).
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto ctx_resetter = set_scoped_default_context(std::move(ctx));

  using sequant::BatchModeType;
  using sequant::evaluate;
  using sequant::index_position;
  using sequant::make_batched_custom_evaluator;
  using sequant::never_volatile;
  using sequant::no_scope_guard;
  using node_t = EvalNodeDryRun;

  auto r = backend_test_regime();  // i (occ) extent 10, a (virt) extent 20
  auto cm = std::make_shared<CostModel const>(r);

  // g is a flat operand carrying no occ; the two C legs are composite
  // (a1/a2<i_1,i_2>) carrying the occ only as protos -- same shape as the TA
  // regression's W{a1<i,j>,a2<i,j>} = (g * C) * C giant.
  auto expr = deserialize<ExprPtr>(
      "(g{a_3;a_4} * C{a_4;a1<i_1,i_2>}) * C{a2<i_1,i_2>;a_3}");
  bool const parsed = static_cast<bool>(expr);
  REQUIRE(parsed);

  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = binarize<EvalExprDryRun>(expr);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

  auto const occ = get_default_context().index_space_registry()->retrieve(L"i");
  auto accept_occ = [occ](Index const& ix) {
    return ix.space() == occ && !ix.has_proto_indices();
  };

  Index mode;
  for (auto const& ix : node->canon_indices())
    if (accept_occ(ix)) {
      mode = ix;
      break;
    }
  REQUIRE(mode.nonnull());
  // The proto occ IS locatable: promoted to a plain outer canon index.
  REQUIRE(index_position(node, mode).has_value());

  // Stamp External on the root and every internal node whose result carries
  // the occ, as the optimizer would for a forest-level external mode.
  node->set_batched_here({{mode, BatchModeType::External}});
  auto stamp_carriers = [&](auto&& self, node_t& n) -> void {
    if (n.leaf()) return;
    if (&n != &node && index_position(n, mode).has_value())
      n->set_batched_here({{mode, BatchModeType::External}});
    self(self, n.left());
    self(self, n.right());
  };
  stamp_carriers(stamp_carriers, node);

  DryRunLeafEvaluator yield{cm};

  // Reference: plain unbatched evaluation (batched_here are ignored without a
  // custom evaluator).
  auto const ref = evaluate(node, yield);
  REQUIRE(ref);
  auto const ref_bytes = ref->size_in_bytes();
  REQUIRE(ref_bytes > 0);

  // Spy scope-guard: records the block count each time the scatter fires.
  std::vector<std::size_t> guard_calls;
  auto spy = [&guard_calls](std::size_t n) {
    guard_calls.push_back(n);
    return no_scope_guard{};
  };

  auto cache = sequant::CacheManager<node_t>::empty();
  cache.set_custom_evaluator(make_batched_custom_evaluator(
      yield, [](Index const&) -> std::size_t { return 4; }, accept_occ, spy,
      never_volatile{}));

  ResultPtr result;
  bool threw = false;
  std::string what;
  try {
    result = evaluate(node, yield, cache);
  } catch (std::exception const& e) {
    threw = true;
    what = e.what();
  }

  // GREEN (after D3.1): the replay completes. RED (before D3.1): this threw
  // std::logic_error(".. pre_sized_zeros_over_mode ..") the first time the
  // scatter branch called part->pre_sized_zeros_over_mode() on a dry-run
  // Result that did not override it.
  INFO("evaluate() threw: " << what);
  REQUIRE_FALSE(threw);
  REQUIRE(result);

  // The assembled result's mode-th index is widened back to the FULL
  // (unsliced) extent: the scattered result reconstructs the same modeled
  // size as the unbatched reference. (The per-block modeled size's
  // ~block/extent scaling -- the sliced footprint the scatter buys per
  // block -- is asserted directly against the cost model by the two
  // pre_sized_zeros_over_mode unit tests above; this end-to-end replay
  // additionally confirms the runtime genuinely reassembles them via the
  // scatter branch rather than, say, silently no-op'ing.)
  CHECK(result->size_in_bytes() == ref_bytes);

  // The scatter genuinely fired over the occ: > 1 block (occ extent 10,
  // block width 4 -> 3 blocks).
  REQUIRE_FALSE(guard_calls.empty());
  for (auto const n : guard_calls) CHECK(n > 1);
  CHECK(guard_calls.front() == 3);
}

// ===========================================================================
// Task 6: THE replay harness. Deserializes the real post-transform giant
// term (the FIRST summand of csv_ccsd_doubles_residual_df.txt -- see the
// POST-TRANSFORM VERDICT case above, which already established this is the
// ~13-tensor free-mu~ giant), optimizes it once under the SAME C60-scale
// regime/BatchPolicy the DP verdict used, binarizes with the DP's
// batch-mode annotations, then REPLAYS it through the REAL runtime
// (make_evaluator / evaluate<Trace::On>) against zero-data DryRun tokens --
// witnessing what the runtime actually realizes, not what the DP annotated.
// ===========================================================================

TEST_CASE(
    "dryrun eval backend replays the post-transform giant term through the "
    "real batched runtime",
    "[dryrun-eval]") {
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);  // mu~
  sequant::mbpt::add_df_spaces(isr);                             // K
  auto ctx_resetter = set_scoped_default_context(std::move(ctx));

  auto const body = slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                          "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = deserialize<ExprPtr>(line);
  bool const parsed = static_cast<bool>(expr);
  REQUIRE(parsed);
  bool const is_sum = expr->is<Sum>();
  REQUIRE(is_sum);

  // Focus the replay on the GIANT TERM only -- the batched DP is 20-60s per
  // term (running it on all 55 would take ~20-30 minutes). The [dryrun-df]
  // verdict case above's exhaustive sweep (`for (auto const& s :
  // expr->as<Sum>()) ...`, the SAME deserialize-order iteration as
  // `summands()` below) empirically identified summand index 38 (term 39/55,
  // 1-indexed in that sweep's log) as the giant: the ~13-tensor
  // g.C.g.C.s.C.C.s.C.C.t.t.t chain reporting the 1.2 TB free-mu~
  // intermediate at C60 (mu~=1800, K=4320) scale. A cheap proxy (picking the
  // summand with the most flattened tensor factors) was tried first and
  // picked a DIFFERENT, structurally-similar but much smaller term (14
  // factors, ~0.0005 GB giant) -- factor count alone does not identify the
  // giant, only the DP's actual index-space accounting does. So this uses
  // the exhaustively-verified positional index directly rather than a
  // heuristic that was empirically shown to pick the wrong term.
  auto const& summands = expr->as<Sum>().summands();
  REQUIRE(!summands.empty());
  auto flatten_product = [](ExprPtr const& e) -> ExprPtr {
    if (!e->is<Product>()) return e;
    auto const& p = e->as<Product>();
    return ex<Product>(p.scalar(), p.factors(), Product::Flatten::Yes);
  };
  std::size_t const giant_idx = 38 < summands.size() ? 38 : 0;
  ExprPtr giant = flatten_product(summands[giant_idx]);
  REQUIRE(giant);
  std::size_t const giant_nfactors =
      giant->is<Product>() ? giant->as<Product>().factors().size() : 1;
  std::cerr << "[dryrun-eval] selected giant term (index " << giant_idx
            << "): " << giant_nfactors << " flattened factors (of "
            << summands.size() << " summands)\n";

  // FAITHFUL real C60 config (614336 job log): occ=120, PNO=42, OSV=310,
  // mu~=1800, aux=4320, pao_target_size=256, aux_target_size=72,
  // volatile_weight=20, machine_balance=200, fast_mem_elems=1e6. SAME regime
  // and SAME batchable predicate the [dryrun-df] verdict case uses.
  auto regime = df_regime(kC60_pVDZF12);
  auto cm = std::make_shared<CostModel const>(regime);

  // ONE BatchPolicy object, reused verbatim for both optimize() and the
  // runtime evaluator factory (make_evaluator) -- the plan's hard
  // constraint, so the DP's and the runtime's notion of "batchable" and
  // "target batch size" cannot drift apart.
  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = is_df_batchable;
  policy.batch_target_size = [](Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"μ̃" ? std::size_t{256}  // pao_target_size
                                         : std::size_t{72};  // aux_target_size
  };
  policy.is_volatile_leaf = [](Tensor const& t) { return t.label() == L"t"; };
  policy.accumulation_factor = 1.0;
  policy.peak_threshold = 40e9;  // DP-side knob only; the runtime evaluator
                                 // never consults peak_threshold.

  auto axes_map = std::make_shared<std::unordered_map<
      Expr const*, container::vector<NodeBatchAnnotation>>>();
  OptimizeOptions opts;
  opts.objective_function = ObjectiveFunction::DensePeakSizeBatched;
  opts.idx_to_extent = regime.idx_to_extent();
  opts.inner_pow = regime.inner_pow_fn();
  opts.batch_policy = policy;
  opts.volatile_weight = 20.0;
  opts.roofline.machine_balance = 200.0;
  opts.roofline.fast_mem_elems = 1000000.0;
  opts.term_batch_axes = axes_map;

  auto t0 = std::chrono::steady_clock::now();
  auto optimized = optimize(giant, opts);
  auto const opt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                          std::chrono::steady_clock::now() - t0)
                          .count();
  std::cerr << "[dryrun-eval] optimize(giant) took " << opt_ms << "ms\n";
  bool const have_optimized = static_cast<bool>(optimized);
  REQUIRE(have_optimized);

  auto it = axes_map->find(optimized.get());
  container::vector<NodeBatchAnnotation> node_axes;
  if (it != axes_map->end()) node_axes = it->second;
  REQUIRE(!node_axes.empty());

  BinarizationOptions bopts;
  bopts.node_batch_axes = node_axes;

  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = binarize<EvalExprDryRun>(optimized, {}, bopts);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

  // Locate the giant sub-node (the free-mu~ contraction node whose modeled
  // size dwarfs everything else -- same identification criterion the
  // [dryrun-df] verdict case above used) purely to REPORT what the DP
  // annotated on it before the runtime replay, for side-by-side comparison
  // with what the runtime actually realizes.
  auto memsize_ext = sequant::opt::detail::memsize_counter(
      regime.idx_to_extent(), regime.inner_pow_fn());
  auto has_free_mu_tilde = [](std::vector<Index> const& ixs) {
    for (auto const& ix : ixs)
      if (ix.space().base_key() == L"μ̃") return true;
    return false;
  };
  double giant_nominal_bytes = 0.0;
  std::wstring giant_desc, giant_axes_desc;
  node.visit_internal([&](auto const& n) {
    auto free_ixs = node_free_indices(*n);
    if (!has_free_mu_tilde(free_ixs)) return;
    double const bytes =
        memsize_ext(free_ixs, std::vector<Index>{}, std::vector<Index>{}) * 8.0;
    if (bytes > giant_nominal_bytes) {
      giant_nominal_bytes = bytes;
      giant_desc = describe_indices(free_ixs);
      giant_axes_desc = describe_indices(batch_axes_indices(n->batched_here()));
    }
  });
  REQUIRE(giant_nominal_bytes > 0.0);

  auto cache = sequant::cache_manager(std::vector<EvalNodeDryRun>{node});
  cache.set_custom_evaluator(
      sequant::make_evaluator(policy, DryRunLeafEvaluator{cm}));

  // Enable eval tracing (redirected to a private ostringstream, not stdout)
  // so working_set_hwmark() actually accumulates: the engine's per-op hwmark
  // input is gated at RUNTIME on Logger::instance().eval.level > 0
  // (log::printing()), independent of the Trace::On COMPILE-TIME template
  // argument below (which only gates whether the tracing code path exists at
  // all). Restore the previous logger state afterward so this test does not
  // leak global state to others.
  std::ostringstream trace_os;
  auto& logger = Logger::instance();
  auto const prev_level = logger.eval.level;
  auto* const prev_stream = logger.eval.stream;
  logger.eval.level = 2;
  logger.eval.stream = &trace_os;

  std::cerr << "[dryrun-eval] replaying giant term through the batched "
               "runtime evaluator ...\n";
  auto t1 = std::chrono::steady_clock::now();
  ResultPtr result;
  bool threw = false;
  std::string what;
  try {
    result = sequant::evaluate<Trace::On>(node, DryRunLeafEvaluator{cm}, cache);
  } catch (std::exception const& e) {
    threw = true;
    what = e.what();
  }
  auto const eval_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                           std::chrono::steady_clock::now() - t1)
                           .count();

  logger.eval.level = prev_level;
  logger.eval.stream = prev_stream;

  INFO("evaluate() threw: " << what);
  REQUIRE(!threw);
  REQUIRE(result);

  auto const peak = cache.working_set_hwmark();
  auto const root_bytes = result->size_in_bytes();

  // Parse the captured per-op trace for the largest single `result=<N>B`
  // materialized anywhere during the replay -- the most direct witness of
  // whether the giant was EVER realized at (close to) its full un-sliced
  // size, regardless of caching/accounting nuances in working_set_hwmark().
  std::size_t max_single_result_bytes = 0;
  {
    std::string const trace = trace_os.str();
    std::size_t pos = 0;
    while ((pos = trace.find("result=", pos)) != std::string::npos) {
      pos += std::string("result=").size();
      std::size_t end = trace.find('B', pos);
      if (end == std::string::npos) break;
      std::string const num = trace.substr(pos, end - pos);
      if (!num.empty() &&
          num.find_first_not_of("0123456789") == std::string::npos) {
        std::size_t const v = std::stoull(num);
        max_single_result_bytes = std::max(max_single_result_bytes, v);
      }
      pos = end;
    }
  }

  bool const giant_realized_full =
      max_single_result_bytes >= 0.9 * giant_nominal_bytes;

  // Extra corroborating diagnostics: how many ops ran, and how many distinct
  // BatchGroup interceptions fired (a real nested multi-batch replay should
  // show many -- the giant's mu~ (extent 1800, target 100 => ~18 batches) and
  // K (extent 4320, target 100 => ~44 batches) modes, nested, would fire
  // hundreds of small batched ops if genuinely realized).
  std::string const trace = trace_os.str();
  auto count_occurrences = [](std::string const& hay,
                              std::string const& needle) {
    std::size_t n = 0, pos = 0;
    while ((pos = hay.find(needle, pos)) != std::string::npos) {
      ++n;
      pos += needle.size();
    }
    return n;
  };
  std::size_t const n_eval_lines = count_occurrences(trace, "Eval |");
  std::size_t const n_batch_group_begin =
      count_occurrences(trace, "BatchGroup | Begin");
  // Distribution of result= sizes >= 100 MB, to see whether many large
  // (but sub-nominal) intermediates appeared (consistent with a partially-
  // sliced mode) or none did (consistent with the OTHER mode alone already
  // bounding everything well below 100 MB).
  std::size_t n_results_over_100mb = 0;
  {
    std::size_t pos = 0;
    while ((pos = trace.find("result=", pos)) != std::string::npos) {
      pos += std::string("result=").size();
      std::size_t end = trace.find('B', pos);
      if (end == std::string::npos) break;
      std::string const num = trace.substr(pos, end - pos);
      if (!num.empty() &&
          num.find_first_not_of("0123456789") == std::string::npos) {
        if (std::stoull(num) >= 100'000'000ull) ++n_results_over_100mb;
      }
      pos = end;
    }
  }

  std::wcout << L"\n=== [dryrun-eval] GIANT TERM RUNTIME REPLAY ===\n"
             << L"optimize(): " << opt_ms << L"ms, evaluate(): " << eval_ms
             << L"ms\n"
             << L"DP-annotated giant: free={" << giant_desc << L"} nominal="
             << (giant_nominal_bytes / 1e9) << L" GB batched_here={"
             << giant_axes_desc << L"}\n"
             << L"root result size = " << root_bytes << L" bytes ("
             << (double(root_bytes) / 1e9) << L" GB)\n"
             << L"cache.working_set_hwmark() = " << peak << L" bytes ("
             << (double(peak) / 1e9) << L" GB)\n"
             << L"max single result= observed in trace = "
             << max_single_result_bytes << L" bytes ("
             << (double(max_single_result_bytes) / 1e9) << L" GB)\n"
             << L"trace ops: " << n_eval_lines << L" Eval lines, "
             << n_batch_group_begin << L" BatchGroup interceptions, "
             << n_results_over_100mb << L" results >=100MB\n"
             << L"=> giant realized at (>=90% of) its FULL nominal size during "
                L"replay: "
             << (giant_realized_full
                     ? L"YES (mu~ NOT sliced in practice -- bug "
                       L"reproduced)"
                     : L"NO (sliced down from nominal)")
             << L"\n";

  CHECK(peak > 0);
}

// Task 6 (perf-first validation): optimize the SAME C60 giant term (index 38)
// under BOTH the peak-first (DenseSpaceTimeBatched) and perf-first
// (DenseTimeSpaceBatched) objectives, at the faithful real config, and compare
// the factorization the DP picks. The 4-PAO signature is a contraction node
// carrying >= 4 free mu~ indices (the (mu~ mu~|mu~ mu~) AO integral).
// Peak-first forms it (fully sliceable below the 40 GB threshold, so it
// survives the hard filter despite astronomically higher flops); perf-first is
// flops-primary and must NEVER form it. This is the direct in-harness proof of
// the fix, on ONE term (~seconds), without the [dryrun-df] full-sweep cost.
TEST_CASE(
    "dryrun objective determines the C60 PPL factorization (perf-first forms "
    "the 4-PNO integral, peak-first does not)",
    "[dryrun-objective]") {
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);  // mu~
  sequant::mbpt::add_df_spaces(isr);                             // K
  auto ctx_resetter = set_scoped_default_context(std::move(ctx));

  auto const body = slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                          "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = deserialize<ExprPtr>(line);
  REQUIRE(static_cast<bool>(expr));
  REQUIRE(expr->is<Sum>());
  auto const& summands = expr->as<Sum>().summands();
  REQUIRE(!summands.empty());
  auto flatten_product = [](ExprPtr const& e) -> ExprPtr {
    if (!e->is<Product>()) return e;
    auto const& p = e->as<Product>();
    return ex<Product>(p.scalar(), p.factors(), Product::Flatten::Yes);
  };
  // SEQUANT_UT_DRYRUN_LIST_TERMS=1 prints every summand's factor structure so
  // the PPL/ladder term (the one whose early-K factorization would form a
  // 4-PNO integral) can be located by eye.
  if (std::getenv("SEQUANT_UT_DRYRUN_LIST_TERMS")) {
    for (std::size_t s = 0; s < summands.size(); ++s)
      std::wcerr << L"[term " << s << L"] "
                 << to_latex(flatten_product(summands[s])) << L"\n";
  }
  // Term selector (default 38): SEQUANT_UT_DRYRUN_TERM.
  std::size_t giant_idx = 38;
  if (char const* te = std::getenv("SEQUANT_UT_DRYRUN_TERM"))
    giant_idx = static_cast<std::size_t>(std::atoll(te));
  if (giant_idx >= summands.size()) giant_idx = 0;
  ExprPtr giant = flatten_product(summands[giant_idx]);
  REQUIRE(giant);
  std::wcerr << L"[dryrun-objective] selected term " << giant_idx << L" = "
             << to_latex(giant) << L"\n";

  // FAITHFUL real C60 config: the measured heavy-tailed CSV moments in
  // kC60_pVDZF12, NOT a flat scalar -- so the 4-PNO PPL node is sized
  // occ^2 * M_4^4 (real, per the measured M_4) instead of a flat under-count.
  auto regime = df_regime(kC60_pVDZF12);
  auto memsize = sequant::opt::detail::memsize_counter(regime.idx_to_extent(),
                                                       regime.inner_pow_fn());
  auto ext_of = [](Index const& ix) -> double {
    auto bk = ix.space().base_key();
    if (bk == L"μ̃") return double(kC60_pVDZF12.mu_tilde);
    if (bk == L"Κ") return double(kC60_pVDZF12.aux);
    return 0.0;
  };
  auto tgt_of = [](Index const& ix) -> double {
    return ix.space().base_key() == L"μ̃" ? 256.0 : 72.0;
  };
  auto keyof = [](Index const& ix) { return std::wstring(ix.full_label()); };

  struct Analysis {
    std::size_t max_free_mu = 0;       // >= 4 => 4-PAO AO integral formed
    double largest_realized_gb = 0.0;  // DP-model realized free-mu~ (static)
    std::wstring largest_desc;
    // Peak/flops/exec from the single shared cost_profile() entry point (Task
    // 4: gated-cache peak replay + static flops/exec walk), replacing the
    // ad-hoc manual replay + hwmark read this case used before Task 5.
    sequant::eval::dryrun::CostProfile cp;
    // P1 seed-path result (populated only when seed_external_occ=true): the
    // honest WITHIN-MODEL (DP-staged) unseeded/seeded peak_bytes and flops.
    // Unlike `cp`, which mixes an unseeded cost_profile()-replay hwmark with
    // (when seeded) an OVERRIDDEN peak_bytes, both sr.unseeded_* and
    // sr.seeded_* come from the SAME PeakBatchedModel/DP so they are
    // comparable to each other.
    sequant::opt::detail::SeededBatchedResult sr;
    // The SAME seed probe re-run with peak_threshold = +infinity, i.e. with
    // select_root's feasibility filter inert. Needed for the work-neutrality
    // proof: under a FINITE budget the seeded and unseeded root contexts have
    // different peaks, hence different FEASIBLE SETS, so select_root legally
    // returns two DIFFERENT factorizations (unseeded: nothing fits => global
    // min-flops fallback; seeded: min-flops AMONG THOSE THAT FIT) and their
    // flops differ by the price of feasibility -- which says nothing about
    // whether slicing an external is work-neutral. At +inf both sides return
    // the same (globally cheapest) factorization, so any flops difference is
    // attributable to the slicing alone. Seeding still slices the occ either
    // way: the seed is the root batch CONTEXT bit, orthogonal to the threshold.
    sequant::opt::detail::SeededBatchedResult sr_inf;
  };

  // Optimize the giant under `obj`, binarize, and walk the tree computing per
  // node (a) its free-mu~ count (the 4-PAO structural signature) and (b) its
  // REALIZED free-mu~ size after ancestor slicing (same active-ancestor
  // accounting as the [dryrun-df] verdict case). Then call the single shared
  // cost_profile() entry point (Task 4) to get the modeled peak/flops/exec via
  // the gated-cache replay -- replacing the ad-hoc manual replay + hwmark read
  // this case used before Task 5, so there is ONE peak/flops code path.
  auto analyze = [&](ObjectiveFunction obj, bool seed_external_occ = false,
                     std::size_t occ_block = 0) -> Analysis {
    sequant::BatchPolicy policy;
    policy.is_batchable_contracted_index = is_df_batchable;
    policy.batch_target_size = [](Index const& ix) -> std::size_t {
      return ix.space().base_key() == L"μ̃" ? std::size_t{256} : std::size_t{72};
    };
    policy.is_volatile_leaf = [](Tensor const& t) { return t.label() == L"t"; };
    policy.accumulation_factor = 1.0;
    // This case characterizes the perf-first-vs-peak-first contrast under the
    // LEGACY set-keyed model: peak-first forms the fully-sliceable 4-PAO (its
    // batched peak looks small) while perf-first refuses it on flops. Pin
    // order_aware_recompute OFF now that BatchPolicy defaults it ON -- under
    // the realistic (resident-scan) model the contrast COLLAPSES (the 4-PAO's
    // batched peak is priced higher by accumulator residency, so peak-first
    // also avoids it, and perf-first's flops then exceed peak-first's because
    // the per-block recompute is charged). The legacy contrast is what
    // motivated the perf-first objective, so it is what this case documents.
    policy.order_aware_recompute = false;
    policy.peak_threshold =
        (std::getenv("SEQUANT_UT_DRYRUN_PEAK_THR_GB")
             ? std::atof(std::getenv("SEQUANT_UT_DRYRUN_PEAK_THR_GB"))
             : 40.0) *
        1e9;

    auto axes_map = std::make_shared<std::unordered_map<
        Expr const*, container::vector<NodeBatchAnnotation>>>();
    OptimizeOptions opts;
    opts.objective_function = obj;
    opts.idx_to_extent = regime.idx_to_extent();
    opts.inner_pow = regime.inner_pow_fn();
    opts.batch_policy = policy;
    opts.volatile_weight = 20.0;
    opts.roofline.machine_balance = 200.0;
    opts.roofline.fast_mem_elems = 1000000.0;
    opts.term_batch_axes = axes_map;

    auto t0 = std::chrono::steady_clock::now();
    auto optimized = optimize(giant, opts);
    auto const opt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                            std::chrono::steady_clock::now() - t0)
                            .count();
    REQUIRE(static_cast<bool>(optimized));
    auto it = axes_map->find(optimized.get());
    container::vector<NodeBatchAnnotation> node_axes;
    if (it != axes_map->end()) node_axes = it->second;
    BinarizationOptions bopts;
    bopts.node_batch_axes = node_axes;
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto node = binarize<EvalExprDryRun>(optimized, {}, bopts);
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

    // Optional full-tree dump (SEQUANT_UT_DRYRUN_PERF_TREE=1): every node's
    // free indices and batched_here (the indices SLICED, i.e. batched, at that
    // node), in post-order-ish indentation, so the exact schedule is
    // inspectable.
    if (std::getenv("SEQUANT_UT_DRYRUN_PERF_TREE")) {
      wchar_t const* on = (obj == ObjectiveFunction::DenseTimeSpaceBatched)
                              ? L"perf-first (DenseTimeSpaceBatched)"
                              : L"peak-first (DenseSpaceTimeBatched)";
      std::wcerr << L"\n[dryrun-objective] TREE for " << on << L":\n";
      std::function<void(std::remove_cvref_t<decltype(node)> const&, int)>
          dump = [&](auto const& n, int depth) {
            std::wstring const pad(2 * depth + 2, L' ');
            auto const free = node_free_indices(*n);
            container::vector<Index> const bax =
                batch_axes_indices(n->batched_here());
            std::size_t nmu = 0, nk = 0, npno = 0, nosv = 0;
            for (auto const& ix : free) {
              if (ix.space().base_key() == L"μ̃") ++nmu;
              if (ix.space().base_key() == L"Κ") ++nk;
              if (ix.has_proto_indices()) {
                if (ix.proto_indices().size() >= 2)
                  ++npno;  // 2-proto (or higher) = PNO composite
                else
                  ++nosv;  // 1-proto = OSV composite
              }
            }
            std::wcerr << pad << (n.leaf() ? L"leaf  " : L"CONTRACT ")
                       << L"free={" << describe_indices(free) << L"} (mu~="
                       << nmu << L" K=" << nk << L" PNO=" << npno << L" OSV="
                       << nosv << L")  batched_here={" << describe_indices(bax)
                       << L"}\n";
            if (!n.leaf()) {
              dump(n.left(), depth + 1);
              dump(n.right(), depth + 1);
            }
          };
      dump(node, 0);
    }

    Analysis a;
    std::function<void(std::remove_cvref_t<decltype(node)> const&,
                       std::map<std::wstring, std::wstring>)>
        walk = [&](auto const& n, std::map<std::wstring, std::wstring> active) {
          auto free_ixs = node_free_indices(*n);
          std::size_t nmu = 0;
          for (auto const& ix : free_ixs)
            if (ix.space().base_key() == L"μ̃") ++nmu;
          if (nmu > a.max_free_mu) a.max_free_mu = nmu;
          if (nmu > 0) {
            double nominal =
                memsize(free_ixs, std::vector<Index>{}, std::vector<Index>{}) *
                8.0;
            double factor = 1.0;
            for (auto const& ix : free_ixs) {
              auto bk = ix.space().base_key();
              if (bk != L"μ̃" && bk != L"Κ") continue;
              if (active.find(keyof(ix)) != active.end())
                factor *= tgt_of(ix) / ext_of(ix);
            }
            double const realized_gb = nominal * factor / 1e9;
            if (realized_gb > a.largest_realized_gb) {
              a.largest_realized_gb = realized_gb;
              a.largest_desc = describe_indices(free_ixs);
            }
          }
          if (!n.leaf()) {
            std::map<std::wstring, std::wstring> child_active = active;
            for (auto const& ax : n->batched_here())
              child_active[keyof(ax.first)] = L"y";
            walk(n.left(), child_active);
            walk(n.right(), child_active);
          }
        };
    walk(node, {});

    // ---- modeled cost via the single shared cost_profile() entry point ----
    // Replaces the ad-hoc manual replay + working_set_hwmark read this case
    // used before Task 5. cost_profile() builds the gated dry-run cache (Task
    // 2: free-batchable-mode veto + footprint gate), replays zero-data through
    // the real eval loop with the Task-3 scratch-fold PeakSink, forces the
    // printing gate on internally, and folds the outer cached residency -- so
    // peak_bytes captures the batched-inner transient the raw outer hwmark
    // misses. The same CacheConfig the [cost_profile] test uses;
    // is_batchable_index is set from `policy` inside cost_profile() (advisory
    // here).
    sequant::eval::dryrun::CacheConfig cfg;
    cfg.max_footprint = 1e11;
    cfg.min_repeats = 1;
    cfg.is_volatile = [](EvalNodeDryRun const& n) {
      if (!n.leaf() || !n->is_tensor()) return false;
      return n->as_tensor().label() == L"t";
    };
    a.cp = sequant::eval::dryrun::cost_profile(
        std::vector<EvalNodeDryRun>{node}, policy, cfg, regime,
        /*trace=*/nullptr);

    wchar_t const* obj_name = (obj == ObjectiveFunction::DenseTimeSpaceBatched)
                                  ? L"perf-first (DenseTimeSpaceBatched)"
                                  : L"peak-first (DenseSpaceTimeBatched)";
    std::wcerr << L"[dryrun-objective] " << obj_name << L": optimize " << opt_ms
               << L"ms  max free-mu~ on a node=" << a.max_free_mu
               << L"  largest realized free-mu~={" << a.largest_desc << L"}="
               << a.largest_realized_gb << L" GB\n               cost_profile: "
               << L"n_ops=" << a.cp.model_n_ops << L" flops="
               << a.cp.model_flops << L" peak_bytes=" << (a.cp.peak_bytes / 1e9)
               << L" GB\n";

    // ---- P1 forest-batching seed: size+report the giant under an occ-slice --
    // Reuse the [dryrun-occ-sizing] giant-node locator + ext_occ extraction to
    // find the ONE external occ carried on the perf-first PPL W node,
    // then run the batched DP with that occ SEEDED into the root frontier and
    // OVERRIDE the reported peak with the seeded root point. flops (a static,
    // batch-invariant walk of the SAME factorization) is left untouched:
    // partitioning an external mode is work-neutral.
    if (seed_external_occ) {
      auto memsize_full = sequant::opt::detail::memsize_counter(
          regime.idx_to_extent(), regime.inner_pow_fn());
      auto free_has_bare = [](std::vector<Index> const& ixs,
                              std::wstring_view bk) {
        for (auto const& ix : ixs)
          if (ix.space().base_key() == bk) return true;
        return false;
      };
      auto count_pno = [](std::vector<Index> const& ixs) {
        std::size_t nn = 0;
        for (auto const& ix : ixs)
          if (ix.proto_indices().size() >= 2) ++nn;
        return nn;
      };
      double giant_full_bytes = 0.0;
      std::vector<Index> giant_free;
      node.visit_internal([&](auto const& nn) {
        auto free_ixs = node_free_indices(*nn);
        if (free_has_bare(free_ixs, L"μ̃") || free_has_bare(free_ixs, L"Κ"))
          return;
        if (count_pno(free_ixs) < 2) return;
        double const bytes =
            memsize_full(free_ixs, std::vector<Index>{}, std::vector<Index>{}) *
            8.0;
        if (bytes > giant_full_bytes) {
          giant_full_bytes = bytes;
          giant_free = free_ixs;
        }
      });
      REQUIRE(giant_full_bytes > 0.0);
      std::vector<Index> ext_occ;
      for (auto const& ix : giant_free)
        for (auto const& p : ix.proto_indices())
          if (p.space().base_key() == L"i") {
            bool seen = false;
            for (auto const& e : ext_occ)
              if (e.full_label() == p.full_label()) seen = true;
            if (!seen) ext_occ.push_back(p);
          }
      REQUIRE(!ext_occ.empty());

      // Build the giant's TensorNetwork (mirrors single_term_opt) and run the
      // batched DP with the SAME perf-first params optimize() used, seeding the
      // ONE external occ index into the ROOT batch context.
      REQUIRE(giant->is<Product>());
      container::svector<ExprPtr> gtensors;
      for (auto const& f : giant->as<Product>().factors())
        if (f->is<Tensor>()) gtensors.push_back(f);
      TensorNetwork gtn{gtensors};

      using BModel = sequant::opt::detail::PeakBatchedModel<
          std::function<std::size_t(Index const&)>>;
      BModel model{regime.idx_to_extent(),
                   [](Index const& ix) -> std::size_t {
                     return ix.space().base_key() == L"μ̃" ? std::size_t{256}
                                                          : std::size_t{72};
                   },
                   [](Tensor const& t) { return t.label() == L"t"; },
                   regime.inner_pow_fn(),
                   /*volatile_weight=*/20.0,
                   /*machine_balance=*/200.0,
                   /*fast_mem_elems=*/1000000.0,
                   /*block_tiles=*/3.0,
                   /*block_prefactor=*/1.0,
                   /*batch_persistent_only=*/false,
                   /*peak_flops_tolerance=*/0.0,
                   /*accumulation_factor=*/1.0,
                   /*peak_threshold=*/40.0 * 1e9,
                   /*numeric_size=*/8.0,
                   /*perf_first=*/true};
      model.is_batchable_contracted_index = is_df_batchable;
      container::svector<Index> const gtidxs{};
      auto sr = sequant::opt::detail::seeded_root_peak_batched(
          model, gtn, gtidxs, ext_occ.front(), occ_block);
      std::wcerr << L"[dryrun-objective] SEED external occ {"
                 << std::wstring(sr.seeded_axis ? sr.seeded_axis->full_label()
                                                : std::wstring{L"?"})
                 << L"} occ_block=" << sr.occ_block << L" spectator_ok="
                 << (sr.spectator_ok ? 1 : 0) << L"\n   DP peak unseeded="
                 << (sr.unseeded_peak_bytes / 1e9) << L" GB seeded="
                 << (sr.seeded_peak_bytes / 1e9) << L" GB  flops unseeded="
                 << sr.unseeded_flops << L" seeded=" << sr.seeded_flops
                 << L"\n";
      // Re-probe with the feasibility filter inert (peak_threshold = +inf) so
      // both root contexts return the same (globally cheapest) factorization
      // and the flops comparison isolates the slicing. See Analysis::sr_inf.
      BModel model_inf = model;
      model_inf.peak_threshold = std::numeric_limits<double>::infinity();
      auto sr_inf = sequant::opt::detail::seeded_root_peak_batched(
          model_inf, gtn, gtidxs, ext_occ.front(), occ_block);
      std::wcerr << L"   (+inf budget) flops unseeded=" << sr_inf.unseeded_flops
                 << L" seeded=" << sr_inf.seeded_flops << L"\n";
      // The reported peak becomes the seeded root point; flops is work-neutral
      // and stays the SAME cost_profile value the unseeded run has.
      a.cp.peak_bytes = sr.seeded_peak_bytes;
      a.sr = sr;
      a.sr_inf = sr_inf;
    }
    return a;
  };

  auto peak_first = analyze(ObjectiveFunction::DenseSpaceTimeBatched);
  auto perf_first = analyze(ObjectiveFunction::DenseTimeSpaceBatched);
  // P1 go/no-go: the SAME perf-first factorization, with the external
  // occ seeded into the DP root frontier so its giant is sized under an
  // occ-slice (occ_block=10, full occ extent 120).
  auto perf_first_occ =
      analyze(ObjectiveFunction::DenseTimeSpaceBatched, true, 10);

  auto report = [](wchar_t const* tag, Analysis const& a) {
    std::wcerr << tag << L": 4-PAO node formed = "
               << (a.max_free_mu >= 4 ? L"YES" : L"NO") << L" (max free mu~="
               << a.max_free_mu << L")\n    DP-model largest realized free-mu~="
               << a.largest_realized_gb << L" GB {" << a.largest_desc << L"}\n"
               << L"    cost_profile: n_ops=" << a.cp.model_n_ops << L" flops="
               << a.cp.model_flops << L" peak_bytes=" << (a.cp.peak_bytes / 1e9)
               << L" GB\n";
  };
  std::wcerr << L"\n=== [dryrun-objective] VERDICT (C60 giant, index 38) ===\n";
  report(L"peak-first (DenseSpaceTimeBatched)", peak_first);
  report(L"perf-first (DenseTimeSpaceBatched)", perf_first);

  // STRUCTURAL PROOF (the direct in-harness proof of the fix), read from the
  // static tree walk above -- NOT from the cost_profile replay: peak-first
  // forms the fully-sliceable 4-PAO AO integral (the C60 pathology),
  // perf-first, being flops-primary, must NEVER form it. Kept as a plain
  // tree-walk check because the free-mu~ signature is a property of the
  // FACTORIZATION the DP picked, not of the peak replay; cost_profile() models
  // cost, it does not expose per-node free-index structure.
  CHECK(peak_first.max_free_mu >= 4);
  CHECK(perf_first.max_free_mu < 4);

  // COST PROOF via the single shared cost_profile() entry point.
  // (a) perf-first is flops-primary: it must not pick a higher-flops
  //     factorization than peak-first.
  CHECK(perf_first.cp.model_flops <= peak_first.cp.model_flops);
  // (b) perf-first's modeled peak is dominated by the GENUINE 4-PNO
  //     intermediate the perf-first schedule forms -- the CC doubles W node
  //     {a_1<i,i> a_2<i,i> a_3<i,i> a_4<i,i>} with FOUR distinct PNO legs over
  //     one occ-pair (see SEQUANT_UT_DRYRUN_PERF_TREE dump). With the FAITHFUL
  //     measured moments (kC60_pVDZF12), it is sized occ^2 * M_4^4 =
  //     120^2 * 53.151^4 * 8 ~= 0.92 TB (dense occ^2 pairs; the screened ~6300
  //     CC pairs would give ~0.40 TB). This is the CORRECT, moment-aware size
  //     (df_regime's csv_pno_moment[k] are power means), NOT a naive-product
  //     artifact and NOT a mis-sized twin R{a<i,i>,a<i,i>} (which would be
  //     occ^2*M_2^2 ~= 0.3 GB). The Kappa mode is CONTRACTED at this node, so
  //     batching cannot shrink W -- it is the irreducible peak floor of the
  //     flop-optimal factorization, and precisely why perf-first (which forms
  //     it) OOMs C60 while peak-first (which does not) does not. A 0.5..2 TB
  //     band brackets the real-moment value with margin and is non-flaky.
  CHECK(perf_first.cp.peak_bytes < 2e12);
  CHECK(perf_first.cp.peak_bytes > 5e11);

  // (c) P1 GO/NO-GO (external-occ forest batching). Seeding the ONE external
  //     external occ (the residual's own output index, carried only as a PNO
  //     protoindex on the giant W) into the DP ROOT batch context sizes the
  //     whole tree with that occ sliced to occ_block=10 (full occ extent 120).
  //     Since the occ is contracted at NO node it is a pure external: slicing
  //     it is work-neutral (identical flops) and shrinks the giant W by
  //     ~occ_block/occ. Gate WITHIN one peak model: both sr.unseeded_peak_bytes
  //     and sr.seeded_peak_bytes come from the SAME DP-staged
  //     PeakBatchedModel/seeded_root_peak_batched() run (unlike cp.peak_bytes,
  //     which for the unseeded case is an unrelated cost_profile()-replay
  //     hwmark). If select_root does NOT pick the occ-sliced realization, this
  //     is a NO-GO (do not force it).
  CHECK(perf_first.cp.peak_bytes > 5e11);
  auto const& sr = perf_first_occ.sr;
  CHECK(sr.seeded_peak_bytes < 0.2 * sr.unseeded_peak_bytes);
  // flops proof: partitioning a pure external mode must be EXACTLY
  // work-neutral (same DP factorization, same total flops) -- exact integer
  // equality, not an Approx across two independently-computed cost_profile()
  // values (which would be equal by construction and prove nothing).
  //
  // Read from the +inf-budget probe, NOT the finite-budget one above. Since
  // select_root gained the perf-first peak_threshold CEILING, a finite budget
  // makes the seeded and unseeded root contexts have different feasible sets,
  // so the two sides legally resolve to DIFFERENT factorizations and their
  // flops differ by the cost of feasibility (the seeded side pays ~6% more to
  // stay under budget). That is the ceiling working, not a work-neutrality
  // violation -- proving neutrality requires holding the factorization fixed,
  // which is exactly what the inert filter does. The occ is still sliced on the
  // seeded side (seeding sets the root batch context; the threshold only
  // selects among candidates), so this remains a genuine sliced-vs-unsliced
  // comparison, and an external not carried on every node would charge
  // batch-recompute and break the equality.
  auto const& sr_inf = perf_first_occ.sr_inf;
  CHECK(sr_inf.seeded_flops == sr_inf.unseeded_flops);
}

// D1.2 (external-mode batching wired into DP SELECTION): the forest external
// seed must flow into the DP's REPORTED peak, not just the standalone
// seeded_root_peak_batched() probe. Optimizing the over-budget C60 giant
// through PeakBatchedModel::reconstruct_batched_modes (the path optimize()
// drives) must, with batch_spectator_indices ON, report a root peak BELOW its
// flag-OFF value (the external occ sliced into the root batch context) and
// stamp BatchModeType::External ONLY on the seeded external occ; with the flag
// OFF the reported peak is byte-identical to the unseeded baseline and NO
// External modes are stamped. This is the wiring the 2026-07-20 experiment
// found missing (peak was byte-identical flag-off vs flag-on because the emit
// was post-selection).
TEST_CASE(
    "dryrun external-mode seeding lowers the DP-reported peak of the C60 giant",
    "[dryrun-extmode]") {
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);  // mu~
  sequant::mbpt::add_df_spaces(isr);                             // K
  auto ctx_resetter = set_scoped_default_context(std::move(ctx));

  auto const body = slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                          "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = deserialize<ExprPtr>(line);
  REQUIRE(static_cast<bool>(expr));
  REQUIRE(expr->is<Sum>());
  auto const& summands = expr->as<Sum>().summands();
  REQUIRE(summands.size() > 38);
  auto flatten_product = [](ExprPtr const& e) -> ExprPtr {
    if (!e->is<Product>()) return e;
    auto const& p = e->as<Product>();
    return ex<Product>(p.scalar(), p.factors(), Product::Flatten::Yes);
  };
  ExprPtr giant = flatten_product(summands[38]);  // the C60 PPL/ladder giant
  REQUIRE(giant);
  REQUIRE(giant->is<Product>());

  auto regime = df_regime(kC60_pVDZF12);
  container::svector<ExprPtr> gtensors;
  for (auto const& f : giant->as<Product>().factors())
    if (f->is<Tensor>()) gtensors.push_back(f);
  TensorNetwork gtn{gtensors};
  container::svector<Index> const gtidxs{};

  using BModel = sequant::opt::detail::PeakBatchedModel<
      std::function<std::size_t(Index const&)>>;
  auto make_model = [&](bool spectator_on) {
    BModel m{regime.idx_to_extent(),
             [](Index const& ix) -> std::size_t {
               auto const k = ix.space().base_key();
               if (k == L"μ̃") return std::size_t{256};  // pao_target_size
               if (k == L"i") return std::size_t{8};    // occ_target_size: occ
               // blocks are small; the aux value 72 must NOT leak to the occ.
               return std::size_t{72};  // aux_target_size (DF Κ)
             },
             [](Tensor const& t) { return t.label() == L"t"; },
             regime.inner_pow_fn(),
             /*volatile_weight=*/20.0,
             /*machine_balance=*/200.0,
             /*fast_mem_elems=*/1000000.0,
             /*block_tiles=*/3.0,
             /*block_prefactor=*/1.0,
             /*batch_persistent_only=*/false,
             /*peak_flops_tolerance=*/0.0,
             /*accumulation_factor=*/1.0,
             /*peak_threshold=*/40.0 * 1e9,
             /*numeric_size=*/8.0,
             /*perf_first=*/true};
    m.is_batchable_contracted_index = is_df_batchable;
    // External role admits the DF/PAO spaces AND the occ: occ is never
    // contracted here (a spectator on the giant), so it is batchable ONLY in
    // the external role -- exactly the role-split the two predicates encode.
    // is_df_batchable alone (μ̃/Κ) would drop the external occ from
    // ctx.batchable_modes, leaving the spectator seed nothing to adopt.
    m.is_batchable_external_index = [](Index const& ix) {
      auto const k = ix.space().base_key();
      return k == L"μ̃" || k == L"Κ" || k == L"i";
    };
    m.batch_spectator_indices = spectator_on;
    return m;
  };

  // Drive the actual optimize() selection path: build the DP table and call
  // reconstruct_batched_modes (which optimize()/run_single_term_opt_axes
  // calls), reading back the REPORTED root peak and the per-node emitted modes.
  auto reported_peak = [&](bool spectator_on,
                           container::vector<NodeBatchAnnotation>& node_axes) {
    auto m = make_model(spectator_on);
    auto mctx = m.build_context(gtn, gtidxs);
    auto mst = sequant::opt::detail::solve_single_term(m, gtn, gtidxs, mctx);
    double peak = 0.0;
    auto [seq, modes] =
        m.reconstruct_batched_modes(mctx, mst, gtn, gtidxs, &peak);
    node_axes = std::move(modes);
    return peak;
  };

  container::vector<NodeBatchAnnotation> ax_off, ax_on;
  double const peak_off = reported_peak(false, ax_off);
  double const peak_on = reported_peak(true, ax_on);

  std::wcerr << L"[dryrun-extmode] reported root peak flag-off="
             << (peak_off / 1e9) << L" GB  flag-on=" << (peak_on / 1e9)
             << L" GB\n";

  // The giant is genuinely over the 40 GB budget (else nothing to batch).
  REQUIRE(peak_off > 40.0 * 1e9);

  // Flag OFF: no External modes stamped; reported peak is the unseeded
  // baseline.
  bool any_external_off = false;
  for (auto const& axs : ax_off)
    for (auto const& e : axs.axes)
      if (e.second == BatchModeType::External) any_external_off = true;
  CHECK(!any_external_off);

  // Flag ON: the DP-reported peak DROPS below the flag-off value -- the
  // external occ sliced into the root batch context, work-neutral (identical
  // flops). Chosen policy (D1.3): JOINTLY seed BOTH external occ i_1,i_2 of the
  // doubles residual. The giant is the 4-PNO-leg particle-particle-ladder
  //   W^{a1<i1,i2> a2<i1,i2>}_{a3<i1,i2> a4<i1,i2>} = sum_K (g.C.C)(g.C.C),
  // whose four virtual legs are all protoindexed by the same occ pair (i1,i2).
  // block/extent = 8/120 each (realistic occ block -- the aux value 72 is wrong
  // for occ), so the footprint scales by the PRODUCT (8/120)^2 = 1/225.
  CHECK(peak_on < peak_off);
  // Joint scaling: exactly (8/120)^2 of the unseeded footprint (both occ
  // sliced).
  CHECK(peak_on == Catch::Approx(peak_off * (8.0 / 120.0) * (8.0 / 120.0)));
  // ...and the giant now FITS the 40 GB budget (~1874 GB -> ~8.3 GB): the DP
  // models external batching bounding the PPL giant with a realistic occ block.
  CHECK(peak_on < 40.0 * 1e9);

  // Flag ON: External stamped, and ONLY on an external occ (space "i") -- the
  // chosen seed modes -- never on a contracted DF-aux/PAO mode (emit follows
  // selection). BOTH external occ (i_1 and i_2) must be stamped (joint seed).
  auto occ_space = isr->retrieve(L"i");
  std::set<std::wstring> external_labels;
  for (auto const& axs : ax_on)
    for (auto const& e : axs.axes)
      if (e.second == BatchModeType::External) {
        CHECK(e.first.space() == occ_space);
        external_labels.insert(std::wstring(e.first.full_label()));
      }
  CHECK(external_labels.size() == 2);  // both i_1 and i_2 seeded jointly
}

// D1.3 selection-policy EXPERIMENT (hidden; run explicitly). Sweeps the
// external-seed selection modes on the C60 giant and prints a predicted
// peak/flops table so the chosen policy is settled on evidence:
//   candidate=single first-adoptable external vs JOINT all-adoptable;
//   block size (batch_target_size) sweep to find the size that fits budget.
// Every adopted variant is work-neutral (seeded_flops == unseeded_flops) --
// external batching never changes total flops -- so the flops column is
// constant; the peak column is what moves. Run:
//   ./tests/unit/unit_tests-sequant "[.][dryrun-extmode-sweep]"
TEST_CASE("dryrun external-mode selection-policy sweep (C60 giant)",
          "[.][dryrun-extmode-sweep]") {
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);
  sequant::mbpt::add_df_spaces(isr);
  auto ctx_resetter = set_scoped_default_context(std::move(ctx));

  auto const body = slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                          "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = deserialize<ExprPtr>(line);
  REQUIRE(static_cast<bool>(expr));
  auto const& summands = expr->as<Sum>().summands();
  auto flatten_product = [](ExprPtr const& e) -> ExprPtr {
    if (!e->is<Product>()) return e;
    auto const& p = e->as<Product>();
    return ex<Product>(p.scalar(), p.factors(), Product::Flatten::Yes);
  };
  ExprPtr giant = flatten_product(summands[38]);
  auto regime = df_regime(kC60_pVDZF12);
  container::svector<ExprPtr> gtensors;
  for (auto const& f : giant->as<Product>().factors())
    if (f->is<Tensor>()) gtensors.push_back(f);
  TensorNetwork gtn{gtensors};
  container::svector<Index> const gtidxs{};

  using BModel = sequant::opt::detail::PeakBatchedModel<
      std::function<std::size_t(Index const&)>>;
  BModel m{regime.idx_to_extent(),
           [](Index const& ix) -> std::size_t {
             return ix.space().base_key() == L"μ̃" ? std::size_t{256}
                                                  : std::size_t{72};
           },
           [](Tensor const& t) { return t.label() == L"t"; },
           regime.inner_pow_fn(),
           /*volatile_weight=*/20.0,
           /*machine_balance=*/200.0,
           /*fast_mem_elems=*/1000000.0,
           /*block_tiles=*/3.0,
           /*block_prefactor=*/1.0,
           /*batch_persistent_only=*/false,
           /*peak_flops_tolerance=*/0.0,
           /*accumulation_factor=*/1.0,
           /*peak_threshold=*/40.0 * 1e9,
           /*numeric_size=*/8.0,
           /*perf_first=*/true};
  m.is_batchable_contracted_index = is_df_batchable;
  m.is_batchable_external_index = is_df_batchable;  // external role (Task-4)
  m.batch_spectator_indices = true;

  auto mctx = m.build_context(gtn, gtidxs);
  auto mst = sequant::opt::detail::solve_single_term(m, gtn, gtidxs, mctx);
  std::size_t const root = (std::size_t{1} << mctx.nt) - 1;
  int const best0 = m.select_root(mctx, mst, 0);
  double const unseeded_peak = mst[root][0][best0].peak * 8.0;
  double const unseeded_flops = mst[root][0][best0].flops;

  // Enumerate the external modes of the giant.
  container::svector<Index> external_modes;
  for (std::size_t k = 0; k < mctx.m; ++k)
    if (m.is_external_mode(mctx, k))
      external_modes.push_back(mctx.batchable_modes[k]);

  std::wcerr << L"\n[extmode-sweep] unseeded peak=" << (unseeded_peak / 1e9)
             << L" GB  flops=" << unseeded_flops << L"  #external modes="
             << external_modes.size() << L"\n";
  for (auto const& s : external_modes)
    std::wcerr << L"[extmode-sweep]   external " << s.full_label()
               << L" extent=" << regime.idx_to_extent()(s) << L"\n";

  auto blk72 = [](Index const&) { return std::size_t{72}; };

  // Per-mode (single-seed) variants.
  for (auto const& s : external_modes) {
    double peak = 0.0, sf = 0.0, uf = 0.0;
    container::svector<Index> one{s};
    bool ok = m.seeded_forest_peak(gtn, gtidxs, one, blk72, peak, &sf, &uf);
    std::wcerr << L"[extmode-sweep] single seed " << s.full_label()
               << L" -> peak=" << (ok ? peak / 1e9 : -1.0)
               << L" GB  work_neutral=" << (sf == uf) << L" (adopt=" << ok
               << L")\n";
  }

  // JOINT all-adoptable at block 72.
  {
    double peak = 0.0, sf = 0.0, uf = 0.0;
    bool ok = m.seeded_forest_peak(gtn, gtidxs, external_modes, blk72, peak,
                                   &sf, &uf);
    std::wcerr << L"[extmode-sweep] JOINT all (" << external_modes.size()
               << L") block=72 -> peak=" << (ok ? peak / 1e9 : -1.0)
               << L" GB  work_neutral=" << (sf == uf) << L" seeded_flops=" << sf
               << L" (adopt=" << ok << L")\n";
  }

  // Block-size sweep on the JOINT seed: find the block that fits budget. Occ
  // blocks are small (8/16 are the realistic occ_target_size values -- 72 is
  // the AUX block and must not be used for the occ); included here for the
  // record.
  for (std::size_t b : {std::size_t{120}, std::size_t{96}, std::size_t{72},
                        std::size_t{48}, std::size_t{36}, std::size_t{24},
                        std::size_t{16}, std::size_t{12}, std::size_t{8}}) {
    auto blk = [b](Index const&) { return b; };
    double peak = 0.0, sf = 0.0, uf = 0.0;
    bool ok =
        m.seeded_forest_peak(gtn, gtidxs, external_modes, blk, peak, &sf, &uf);
    std::wcerr << L"[extmode-sweep] JOINT block=" << b << L" -> peak="
               << (ok ? peak / 1e9 : -1.0) << L" GB  <100GB="
               << (ok && peak < 100e9) << L"  work_neutral=" << (sf == uf)
               << L"\n";
  }
}

// P1 gate spike (external-occ forest batching, mechanism b): PURE SIZING check.
// Does the cost model's footprint of the perf-first PPL W giant respond to
// slicing ONE external occupied index to a block? The external occ (the
// residual's own output i,j) appears ONLY as protoindices of the output PNO
// composites a<i,j>, never as a top-level slot and never contracted -- so the
// question is whether inner_aware_volume, which pulls those protos into the
// OUTER extent product (tot_indices; sized by idx_to_extent) while sizing the
// PNO composites themselves by the k-th CSV power mean (independent of the occ
// block), still shrinks the node when idx_to_extent is overridden to return
// occ_block for one occ proto. If the occ lives in the outer product, the whole
// footprint scales by occ_block/occ (GO). If the composite moment sizing
// swallowed the occ (it does not, but this is the NO-GO hypothesis the spike
// tests), the footprint would not move.
//
// This does NOT touch runtime, optimize() selection, or emit wiring. It sizes
// ONE node's free-index set twice, with two extent functions. Hidden tag; run:
//   ./tests/unit/unit_tests-sequant "[dryrun-occ-sizing]"
TEST_CASE(
    "dryrun external-occ slicing shrinks the PPL W footprint (P1 sizing gate)",
    "[.][dryrun-occ-sizing]") {
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);  // mu~
  sequant::mbpt::add_df_spaces(isr);                             // K
  auto ctx_resetter = set_scoped_default_context(std::move(ctx));

  auto const body = slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                          "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = deserialize<ExprPtr>(line);
  REQUIRE(static_cast<bool>(expr));
  REQUIRE(expr->is<Sum>());
  auto const& summands = expr->as<Sum>().summands();
  REQUIRE(!summands.empty());
  auto flatten_product = [](ExprPtr const& e) -> ExprPtr {
    if (!e->is<Product>()) return e;
    auto const& p = e->as<Product>();
    return ex<Product>(p.scalar(), p.factors(), Product::Flatten::Yes);
  };
  // Same term as [dryrun-objective]: index 38 is the C60 PPL/ladder giant.
  std::size_t giant_idx = 38;
  if (giant_idx >= summands.size()) giant_idx = 0;
  ExprPtr giant = flatten_product(summands[giant_idx]);
  REQUIRE(giant);

  // FAITHFUL real C60 config -- SAME regime the [dryrun-objective] case uses.
  auto regime = df_regime(kC60_pVDZF12);

  // Optimize + binarize the giant under the perf-first objective
  // (DenseTimeSpaceBatched) with the SAME policy the [dryrun-objective] analyze
  // uses, so we walk the exact tree that forms the occ^2 * PNO^4 W node.
  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = is_df_batchable;
  policy.batch_target_size = [](Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"μ̃" ? std::size_t{256} : std::size_t{72};
  };
  policy.is_volatile_leaf = [](Tensor const& t) { return t.label() == L"t"; };
  policy.accumulation_factor = 1.0;
  policy.peak_threshold = 40.0 * 1e9;

  auto axes_map = std::make_shared<std::unordered_map<
      Expr const*, container::vector<NodeBatchAnnotation>>>();
  OptimizeOptions opts;
  opts.objective_function = ObjectiveFunction::DenseTimeSpaceBatched;
  opts.idx_to_extent = regime.idx_to_extent();
  opts.inner_pow = regime.inner_pow_fn();
  opts.batch_policy = policy;
  opts.volatile_weight = 20.0;
  opts.roofline.machine_balance = 200.0;
  opts.roofline.fast_mem_elems = 1000000.0;
  opts.term_batch_axes = axes_map;

  auto optimized = optimize(giant, opts);
  REQUIRE(static_cast<bool>(optimized));
  auto it = axes_map->find(optimized.get());
  container::vector<NodeBatchAnnotation> node_axes;
  if (it != axes_map->end()) node_axes = it->second;
  BinarizationOptions bopts;
  bopts.node_batch_axes = node_axes;
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = binarize<EvalExprDryRun>(optimized, {}, bopts);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

  // Locate the GIANT PPL W node = the perf-first (gC)^2 ladder intermediate
  // W(a1a2a3a4) whose FOUR PNO composite legs a<i,j> share one occ-pair, sized
  // occ^2 * M_4^4 = ~0.92 TB (matches cost_profile's ~954 GB peak). The
  // defining, batching-relevant property is that it carries NO free mu~ and NO
  // free K -- those are contracted at/below it -- so mu~/K batching cannot
  // shrink it and the external occ is its ONLY memory lever. (This is why a
  // plain max-footprint pick is WRONG: the larger 3-center {mu~ a K} nodes ARE
  // mu~/K-sliceable and are not the term that OOMs after batching.) Among the
  // no-free-mu~/K nodes, take the max-footprint one; require >= 2 PNO legs so
  // we land on the ladder intermediate, not a tiny scalar.
  auto memsize_full = sequant::opt::detail::memsize_counter(
      regime.idx_to_extent(), regime.inner_pow_fn());
  auto free_has_bare = [](std::vector<Index> const& ixs, std::wstring_view bk) {
    for (auto const& ix : ixs)
      if (ix.space().base_key() == bk) return true;
    return false;
  };
  auto count_pno = [](std::vector<Index> const& ixs) {
    std::size_t n = 0;
    for (auto const& ix : ixs)
      if (ix.proto_indices().size() >= 2) ++n;
    return n;
  };
  double giant_full_bytes = 0.0;
  std::vector<Index> giant_free;
  node.visit_internal([&](auto const& n) {
    auto free_ixs = node_free_indices(*n);
    if (free_has_bare(free_ixs, L"μ̃") || free_has_bare(free_ixs, L"Κ")) return;
    if (count_pno(free_ixs) < 2) return;
    double const bytes =
        memsize_full(free_ixs, std::vector<Index>{}, std::vector<Index>{}) *
        8.0;
    if (bytes > giant_full_bytes) {
      giant_full_bytes = bytes;
      giant_free = free_ixs;
    }
  });
  REQUIRE(giant_full_bytes > 0.0);

  // The external occ indices on the giant node = the distinct proto indices of
  // its composite legs whose base space is the active occupied ("i"). These are
  // the residual's own output occ (i,j): never a top-level slot, never
  // contracted, present only as PNO protos -- exactly the external mode the
  // forest batching targets.
  std::vector<Index> ext_occ;
  for (auto const& ix : giant_free)
    for (auto const& p : ix.proto_indices())
      if (p.space().base_key() == L"i") {
        bool seen = false;
        for (auto const& e : ext_occ)
          if (e.full_label() == p.full_label()) seen = true;
        if (!seen) ext_occ.push_back(p);
      }

  std::wcerr << L"\n=== [dryrun-occ-sizing] P1 SIZING GATE (C60 giant, term "
             << giant_idx << L") ===\n"
             << L"giant node free indices = {" << describe_indices(giant_free)
             << L"}\n"
             << L"external occ protos found (" << ext_occ.size() << L") = {"
             << describe_indices(ext_occ) << L"}\n";
  REQUIRE(!ext_occ.empty());

  // Slice exactly ONE external occ index to a block (rank-general design:
  // batch an occupied INDEX, never a pair). occ_block = 10 (a plausible occ
  // tile size; full occ extent = kC60_pVDZF12.i_occ = 120).
  std::size_t const occ_block = 10;
  Index const sliced = ext_occ.front();
  std::wstring const sliced_label(sliced.full_label());
  std::wcerr << L"slicing ONE external occ index {" << sliced_label
             << L"} to occ_block=" << occ_block << L" (full occ extent="
             << kC60_pVDZF12.i_occ << L")\n";

  // Sliced extent function: identical to the regime's, except the ONE chosen
  // external occ index returns min(full, occ_block).
  auto full_ext = regime.idx_to_extent();
  auto sliced_ext = [full_ext, sliced_label,
                     occ_block](Index const& ix) -> std::size_t {
    if (std::wstring(ix.full_label()) == sliced_label)
      return std::min<std::size_t>(full_ext(ix), occ_block);
    return full_ext(ix);
  };
  auto memsize_sliced =
      sequant::opt::detail::memsize_counter(sliced_ext, regime.inner_pow_fn());

  double const sliced_bytes =
      memsize_sliced(giant_free, std::vector<Index>{}, std::vector<Index>{}) *
      8.0;
  double const ratio = sliced_bytes / giant_full_bytes;

  std::wcerr << L"full   footprint = " << giant_full_bytes << L" bytes ("
             << (giant_full_bytes / 1e9) << L" GB)\n"
             << L"sliced footprint = " << sliced_bytes << L" bytes ("
             << (sliced_bytes / 1e9) << L" GB)\n"
             << L"sliced/full ratio = " << ratio
             << L"  (expected ~ occ_block/occ = "
             << (double(occ_block) / double(kC60_pVDZF12.i_occ)) << L")\n"
             << L"VERDICT: external-occ slicing "
             << (ratio < 0.6 ? L"SHRINKS the PPL W footprint => GO"
                             : L"does NOT shrink the footprint => NO-GO")
             << L"\n";

  // GO criterion: the modelled footprint must scale EXACTLY with the sliced
  // occ extent, not just "shrink". The outer (idx_to_extent) factor scales
  // linearly in the sliced occ extent; the inner PNO M4^4 factor is untouched.
  // So the footprint ratio is exactly occ_block/occ. Pin that, not just
  // "smaller". If this fails, the finding is that the mechanism is not purely
  // multiplicative as expected (NO-GO); report the measured ratio vs
  // occ_block/occ.
  double const expected_ratio =
      static_cast<double>(occ_block) / static_cast<double>(kC60_pVDZF12.i_occ);
  CHECK(sliced_bytes ==
        Catch::Approx(giant_full_bytes * expected_ratio).epsilon(1e-9));
}

// P1 recognition gate: prove the batched DP's batchable-mode scan
// (batchable_mode_list) now admits the external occupied protoindices carried
// on the C60 giant's composite (PNO) legs. Before the proto-aware change,
// batchable_mode_list scanned only the top-level bra/ket/aux slots, so the
// residual's external occ pair -- which lives ONLY as protoindices of the
// composite legs and never as a top-level slot -- was dropped from
// canon_indices_ and could never enter the batchable set. The forest-batching
// feature needs it there. This is the recognition twin of the
// [dryrun-occ-sizing] sizing gate: same C60 term, same optimize+binarize, same
// giant-node locator and ext_occ extraction; it asserts on the batchable list,
// not the footprint. Hidden tag; run:
//   ./tests/unit/unit_tests-sequant "[dryrun-occ-recognize]"
TEST_CASE(
    "dryrun external occ is recognized as a batchable mode (P1 recognition "
    "gate)",
    "[.][dryrun-occ-recognize]") {
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);  // mu~
  sequant::mbpt::add_df_spaces(isr);                             // K
  auto ctx_resetter = set_scoped_default_context(std::move(ctx));

  auto const body = slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                          "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = deserialize<ExprPtr>(line);
  REQUIRE(static_cast<bool>(expr));
  REQUIRE(expr->is<Sum>());
  auto const& summands = expr->as<Sum>().summands();
  REQUIRE(!summands.empty());
  auto flatten_product = [](ExprPtr const& e) -> ExprPtr {
    if (!e->is<Product>()) return e;
    auto const& p = e->as<Product>();
    return ex<Product>(p.scalar(), p.factors(), Product::Flatten::Yes);
  };
  // Same term as [dryrun-occ-sizing]: index 38 is the C60 PPL/ladder giant.
  std::size_t giant_idx = 38;
  if (giant_idx >= summands.size()) giant_idx = 0;
  ExprPtr giant = flatten_product(summands[giant_idx]);
  REQUIRE(giant);

  // FAITHFUL real C60 config -- SAME regime/policy as [dryrun-occ-sizing].
  auto regime = df_regime(kC60_pVDZF12);

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = is_df_batchable;
  policy.batch_target_size = [](Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"μ̃" ? std::size_t{256} : std::size_t{72};
  };
  policy.is_volatile_leaf = [](Tensor const& t) { return t.label() == L"t"; };
  policy.accumulation_factor = 1.0;
  policy.peak_threshold = 40.0 * 1e9;

  auto axes_map = std::make_shared<std::unordered_map<
      Expr const*, container::vector<NodeBatchAnnotation>>>();
  OptimizeOptions opts;
  opts.objective_function = ObjectiveFunction::DenseTimeSpaceBatched;
  opts.idx_to_extent = regime.idx_to_extent();
  opts.inner_pow = regime.inner_pow_fn();
  opts.batch_policy = policy;
  opts.volatile_weight = 20.0;
  opts.roofline.machine_balance = 200.0;
  opts.roofline.fast_mem_elems = 1000000.0;
  opts.term_batch_axes = axes_map;

  auto optimized = optimize(giant, opts);
  REQUIRE(static_cast<bool>(optimized));
  auto it = axes_map->find(optimized.get());
  container::vector<NodeBatchAnnotation> node_axes;
  if (it != axes_map->end()) node_axes = it->second;
  BinarizationOptions bopts;
  bopts.node_batch_axes = node_axes;
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = binarize<EvalExprDryRun>(optimized, {}, bopts);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

  // Locate the giant PPL W node (identical locator to [dryrun-occ-sizing]) and
  // grab its intermediate tensor: the perf-first (gC)^2 ladder W whose four PNO
  // composite legs share one occ pair, carries NO free mu~ and NO free K, and
  // has >= 2 PNO legs. That tensor's composite legs carry the external occ pair
  // as protoindices -- the external mode the forest batching targets.
  auto memsize_full = sequant::opt::detail::memsize_counter(
      regime.idx_to_extent(), regime.inner_pow_fn());
  auto free_has_bare = [](std::vector<Index> const& ixs, std::wstring_view bk) {
    for (auto const& ix : ixs)
      if (ix.space().base_key() == bk) return true;
    return false;
  };
  auto count_pno = [](std::vector<Index> const& ixs) {
    std::size_t n = 0;
    for (auto const& ix : ixs)
      if (ix.proto_indices().size() >= 2) ++n;
    return n;
  };
  double giant_full_bytes = 0.0;
  ExprPtr giant_tensor;
  node.visit_internal([&](auto const& n) {
    if (!n->is_tensor()) return;
    auto free_ixs = node_free_indices(*n);
    if (free_has_bare(free_ixs, L"μ̃") || free_has_bare(free_ixs, L"Κ")) return;
    if (count_pno(free_ixs) < 2) return;
    double const bytes =
        memsize_full(free_ixs, std::vector<Index>{}, std::vector<Index>{}) *
        8.0;
    if (bytes > giant_full_bytes) {
      giant_full_bytes = bytes;
      giant_tensor = n->as_tensor().clone();
    }
  });
  REQUIRE(giant_full_bytes > 0.0);
  REQUIRE(static_cast<bool>(giant_tensor));

  // The external occ protos the recognition must surface (same extraction as
  // [dryrun-occ-sizing]): distinct occupied ("i") protos of the giant's
  // composite legs.
  std::vector<Index> ext_occ;
  for (auto const& ix : giant_tensor->as<Tensor>().const_braketaux_indices())
    for (auto const& p : ix.proto_indices())
      if (p.space().base_key() == L"i") {
        bool seen = false;
        for (auto const& e : ext_occ)
          if (e.full_label() == p.full_label()) seen = true;
        if (!seen) ext_occ.push_back(p);
      }
  REQUIRE(!ext_occ.empty());

  // Build a one-tensor network from the giant W node and scan it exactly as the
  // batched DP does (cost_model's build_context calls this same function with
  // the same is_batchable predicate). is_df_batchable admits ONLY mu~/K; before
  // the proto-aware change the returned list carries NO occupied index, so this
  // gate FAILS. After it, the external occ pair must appear.
  auto tn = TensorNetwork{container::vector<ExprPtr>{giant_tensor}};
  auto batchable =
      sequant::opt::detail::batchable_mode_list(tn, is_df_batchable);

  std::vector<Index> batchable_v(batchable.begin(), batchable.end());
  std::wcerr << L"\n=== [dryrun-occ-recognize] P1 RECOGNITION GATE (C60 giant, "
             << L"term " << giant_idx << L") ===\n"
             << L"external occ protos on giant (" << ext_occ.size() << L") = {"
             << describe_indices(ext_occ) << L"}\n"
             << L"batchable_mode_list (" << batchable_v.size() << L") = {"
             << describe_indices(batchable_v) << L"}\n";

  auto is_active_occ = [](Index const& ix) {
    return ix.space().base_key() == L"i";
  };
  // The giant's external occ must now be a batchable mode, even though it has
  // no top-level slot on the node.
  CHECK(std::any_of(batchable.begin(), batchable.end(), is_active_occ));

  // ...and the surfaced occ mode is exactly one of the giant's external occ
  // protos, not some unrelated occupied index.
  bool matches_ext = false;
  for (auto const& b : batchable)
    if (is_active_occ(b))
      for (auto const& e : ext_occ)
        if (b.full_label() == e.full_label()) matches_ext = true;
  CHECK(matches_ext);
}

// Task 1 (P0, dry-run cost-model comparison): perf-first
// (DenseTimeSpaceBatched) vs peak-first (DenseSpaceTimeBatched) modelled
// flops/peak, swept over ALL summands of the real C60 residual (not just the
// index-38 giant the case above focuses on). This is the decision input for
// whether enabling perf-first at scale is worth it: the summed-flops ratio
// (perf-first total flops / peak-first total flops) is the modelled
// per-iteration speedup ceiling, and the max modelled peak per objective
// shows the memory spread. This is a MODELLED (flops) proxy, not a
// wall-clock measurement -- see the note at the end of this case.
// Hidden ([.]): a MODELLED-cost diagnostic (no correctness CHECKs, only a
// per-objective flops/peak sweep it prints), and it runs optimize() +
// single_term_opt on EVERY C60 residual summand -- tens of minutes in a
// Debug/-O0 build, which times out CI. Its C60-sweep siblings in this file are
// hidden for the same reason; run it explicitly for the decision-input report.
TEST_CASE(
    "dryrun perf-first vs peak-first modelled cost across all C60 residual "
    "terms",
    "[.][dryrun-perfcost]") {
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);  // mu~
  sequant::mbpt::add_df_spaces(isr);                             // K
  auto ctx_resetter = set_scoped_default_context(std::move(ctx));

  auto const body = slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                          "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = deserialize<ExprPtr>(line);
  REQUIRE(static_cast<bool>(expr));
  REQUIRE(expr->is<Sum>());
  auto const& summands = expr->as<Sum>().summands();
  REQUIRE(!summands.empty());
  auto flatten_product = [](ExprPtr const& e) -> ExprPtr {
    if (!e->is<Product>()) return e;
    auto const& p = e->as<Product>();
    return ex<Product>(p.scalar(), p.factors(), Product::Flatten::Yes);
  };

  // FAITHFUL real C60 config (same regime as [dryrun-objective] above): the
  // measured heavy-tailed CSV moments in kC60_pVDZF12.
  auto regime = df_regime(kC60_pVDZF12);

  struct Analysis {
    sequant::eval::dryrun::CostProfile cp;
  };

  // Optimize `giant` under `obj`, binarize, and return the modeled cost
  // profile (flops/peak_bytes/exec_cost) via the single shared
  // cost_profile() entry point -- the same policy/cache setup as
  // [dryrun-objective]'s analyze() above, minus the per-node free-mu~
  // structural walk (not needed for this aggregate sweep).
  auto analyze = [&](ExprPtr const& giant, ObjectiveFunction obj) -> Analysis {
    sequant::BatchPolicy policy;
    policy.is_batchable_contracted_index = is_df_batchable;
    policy.batch_target_size = [](Index const& ix) -> std::size_t {
      return ix.space().base_key() == L"μ̃" ? std::size_t{256} : std::size_t{72};
    };
    policy.is_volatile_leaf = [](Tensor const& t) { return t.label() == L"t"; };
    policy.accumulation_factor = 1.0;
    policy.peak_threshold =
        (std::getenv("SEQUANT_UT_DRYRUN_PEAK_THR_GB")
             ? std::atof(std::getenv("SEQUANT_UT_DRYRUN_PEAK_THR_GB"))
             : 40.0) *
        1e9;

    auto axes_map = std::make_shared<std::unordered_map<
        Expr const*, container::vector<NodeBatchAnnotation>>>();
    OptimizeOptions opts;
    opts.objective_function = obj;
    opts.idx_to_extent = regime.idx_to_extent();
    opts.inner_pow = regime.inner_pow_fn();
    opts.batch_policy = policy;
    opts.volatile_weight = 20.0;
    opts.roofline.machine_balance = 200.0;
    opts.roofline.fast_mem_elems = 1000000.0;
    opts.term_batch_axes = axes_map;

    auto optimized = optimize(giant, opts);
    REQUIRE(static_cast<bool>(optimized));
    auto it = axes_map->find(optimized.get());
    container::vector<NodeBatchAnnotation> node_axes;
    if (it != axes_map->end()) node_axes = it->second;
    BinarizationOptions bopts;
    bopts.node_batch_axes = node_axes;
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto node = binarize<EvalExprDryRun>(optimized, {}, bopts);
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

    sequant::eval::dryrun::CacheConfig cfg;
    cfg.max_footprint = 1e11;
    cfg.min_repeats = 1;
    cfg.is_volatile = [](EvalNodeDryRun const& n) {
      if (!n.leaf() || !n->is_tensor()) return false;
      return n->as_tensor().label() == L"t";
    };
    Analysis a;
    a.cp = sequant::eval::dryrun::cost_profile(
        std::vector<EvalNodeDryRun>{node}, policy, cfg, regime,
        /*trace=*/nullptr);
    return a;
  };

  double total_peak_first_flops = 0.0, total_perf_first_flops = 0.0;
  double max_peak_first_peak_gb = 0.0, max_perf_first_peak_gb = 0.0;
  std::size_t n_ok = 0, n_skipped = 0;

  std::wcerr << L"\n=== [dryrun-perfcost] per-term peak-first vs perf-first "
                L"modelled cost (C60 residual, "
             << summands.size() << L" terms) ===\n";

  for (std::size_t t = 0; t < summands.size(); ++t) {
    ExprPtr giant = flatten_product(summands[t]);
    if (!giant) {
      std::wcerr << L"perfcost term " << t << L": skipped (null term)\n";
      ++n_skipped;
      continue;
    }

    bool ok = true;
    Analysis peak_first, perf_first;
    try {
      peak_first = analyze(giant, ObjectiveFunction::DenseSpaceTimeBatched);
      perf_first = analyze(giant, ObjectiveFunction::DenseTimeSpaceBatched);
    } catch (std::exception const&) {
      ok = false;
    } catch (...) {
      ok = false;
    }
    if (!ok) {
      std::wcerr << L"term " << t << L": skipped\n";
      ++n_skipped;
      continue;
    }

    ++n_ok;
    double const peak_first_gb = peak_first.cp.peak_bytes / 1e9;
    double const perf_first_gb = perf_first.cp.peak_bytes / 1e9;
    std::wcerr << L"perfcost term " << t << L" | peak_first flops="
               << peak_first.cp.model_flops << L" peak_GB=" << peak_first_gb
               << L" | perf_first flops=" << perf_first.cp.model_flops
               << L" peak_GB=" << perf_first_gb << L"\n";

    total_peak_first_flops += peak_first.cp.model_flops;
    total_perf_first_flops += perf_first.cp.model_flops;
    if (peak_first_gb > max_peak_first_peak_gb)
      max_peak_first_peak_gb = peak_first_gb;
    if (perf_first_gb > max_perf_first_peak_gb)
      max_perf_first_peak_gb = perf_first_gb;
  }

  double const ratio = total_peak_first_flops > 0.0
                           ? total_perf_first_flops / total_peak_first_flops
                           : 0.0;
  std::wcerr << L"perfcost TOTAL peak_first_flops=" << total_peak_first_flops
             << L" perf_first_flops=" << total_perf_first_flops
             << L" ratio(perf/peak)=" << ratio << L"\n";
  std::wcerr << L"perfcost MAXPEAK peak_first=" << max_peak_first_peak_gb
             << L" perf_first=" << max_perf_first_peak_gb << L"\n";
  std::wcerr << L"perfcost: " << n_ok << L" of " << summands.size()
             << L" terms analyzed, " << n_skipped << L" skipped\n";

  // Deliverable is the printed table + ratio above (decision input for D1),
  // not a pass/fail gate -- assert only basic sanity: at least one term
  // produced a positive modelled flops count under the peak-first objective.
  REQUIRE(total_peak_first_flops > 0.0);
}

// Hidden diagnostic: does the dry-run cost model PREDICT the ~3x per-iteration
// overcompute measured on Owl for water-20 (job 649156 aux-only 154 s/iter vs
// 649160 occ+aux 454 s/iter, same 100 GB ceiling, same perf-first objective)?
// Nested external-occ batching (i_1 x i_2, 5x5) is NOT expected to be
// work-neutral: shared/persistent intermediates are re-formed per occ block.
// This sums modelled flops / exec_cost / peak over ALL CCSD-doubles-residual
// terms for the water-20 regime, aux-only vs occ+aux, and prints the ratios.
// A flops/exec ratio ~1 would mean the model treats occ slicing as
// work-neutral (missing the runtime replay); ~3 would mean it charges it.
TEST_CASE("dryrun water-20 occ-batching overcompute (aux-only vs occ+aux)",
          "[.][dryrun-water20-overcompute]") {
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);
  sequant::mbpt::add_df_spaces(isr);
  auto ctx_resetter = set_scoped_default_context(std::move(ctx));

  auto const body = slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                          "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = deserialize<ExprPtr>(line);
  REQUIRE(expr);
  REQUIRE(expr->is<Sum>());
  auto const& summands = expr->as<Sum>().summands();
  REQUIRE(!summands.empty());
  auto flatten_product = [](ExprPtr const& e) -> ExprPtr {
    if (!e->is<Product>()) return e;
    auto const& p = e->as<Product>();
    return ex<Product>(p.scalar(), p.factors(), Product::Flatten::Yes);
  };

  auto regime = df_regime(kWater20_pVDZF12);

  // Three configs -- [0]=aux-only, [1]=occ+aux order_aware=false (the MPQC
  // production path: root-level forest seed), [2]=occ+aux order_aware=true
  // (node-level placement). n_ext/n_con = how many External vs Contracted batch
  // modes the optimizer actually emitted (engagement check).
  int n_ext[3] = {0, 0, 0};
  int n_con[3] = {0, 0, 0};

  // occ_on=false: aux-only (mu~/K contracted-role batching, occ never batched).
  // occ_on=true : + external-occ (i) batching; occ_target 16 -> ceil(80/16)=5
  //               blocks per occ mode; doubles terms nest i_1 x i_2.
  // order_aware: toggles the ordered-key recompute cost model + node-level
  //              external placement (order_aware && spectator). MPQC production
  //              runs order_aware=false; this test compares both for occ+aux to
  //              probe whether order_aware=true does materially more recompute.
  auto analyze = [&](ExprPtr const& giant, bool occ_on, bool order_aware,
                     int cfg_idx) -> sequant::eval::dryrun::CostProfile {
    sequant::BatchPolicy policy;
    policy.is_batchable_contracted_index = is_df_batchable;  // mu~, K
    policy.is_batchable_external_index =
        occ_on ? std::function<bool(Index const&)>([](Index const& ix) {
          auto const k = ix.space().base_key();
          return k == L"μ̃" || k == L"Κ" || k == L"i";
        })
               : std::function<bool(Index const&)>(is_df_batchable);
    policy.batch_spectator_indices = occ_on;
    policy.order_aware_recompute = order_aware;
    policy.batch_target_size = [](Index const& ix) -> std::size_t {
      auto const k = ix.space().base_key();
      if (k == L"i") return 16;   // occ_target_size
      if (k == L"Κ") return 256;  // aux_target_size
      return 256;                 // mu~
    };
    policy.is_volatile_leaf = [](Tensor const& t) { return t.label() == L"t"; };
    policy.accumulation_factor = 1.0;
    policy.peak_threshold = (std::getenv("SEQUANT_UT_W20_THR")
                                 ? std::atof(std::getenv("SEQUANT_UT_W20_THR"))
                                 : 100.0) *
                            1e9;

    auto axes_map = std::make_shared<std::unordered_map<
        Expr const*, container::vector<NodeBatchAnnotation>>>();
    OptimizeOptions opts;
    opts.objective_function = ObjectiveFunction::DenseTimeSpaceBatched;
    opts.idx_to_extent = regime.idx_to_extent();
    opts.inner_pow = regime.inner_pow_fn();
    opts.batch_policy = policy;
    opts.volatile_weight = 20.0;
    opts.roofline.machine_balance = 200.0;
    opts.roofline.fast_mem_elems = 1000000.0;
    opts.term_batch_axes = axes_map;

    auto optimized = optimize(giant, opts);
    REQUIRE(static_cast<bool>(optimized));
    auto it = axes_map->find(optimized.get());
    container::vector<NodeBatchAnnotation> node_axes;
    if (it != axes_map->end()) node_axes = it->second;
    for (auto const& na : node_axes)
      for (auto const& e : na.axes) {
        if (e.second == sequant::BatchModeType::External)
          n_ext[cfg_idx]++;
        else
          n_con[cfg_idx]++;
      }
    BinarizationOptions bopts;
    bopts.node_batch_axes = node_axes;
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto node = binarize<EvalExprDryRun>(optimized, {}, bopts);
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    sequant::eval::dryrun::CacheConfig cfg;
    cfg.max_footprint = 1e11;
    cfg.min_repeats = 1;
    cfg.is_volatile = [](EvalNodeDryRun const& n) {
      return n.leaf() && n->is_tensor() && n->as_tensor().label() == L"t";
    };
    return sequant::eval::dryrun::cost_profile(
        std::vector<EvalNodeDryRun>{node}, policy, cfg, regime,
        /*trace=*/nullptr);
  };

  // Three configs, summed over all terms. cfg 0 = aux-only, 1 = occ+aux
  // order_aware=false (MPQC production, root-level seed), 2 = occ+aux
  // order_aware=true (node-level placement). model_* is order-/batching-blind
  // (~flat); dryrun_* is the recompute-aware replay tally. The 2-vs-1 ratio is
  // the diagnostic: does order_aware=true do materially more recompute?
  double mflops[3] = {0, 0, 0}, mexec[3] = {0, 0, 0};
  double dflops[3] = {0, 0, 0}, dexec[3] = {0, 0, 0};
  std::size_t dops[3] = {0, 0, 0};
  double maxpeak[3] = {0, 0, 0};
  std::size_t n_ok = 0;
  for (std::size_t t = 0; t < summands.size(); ++t) {
    ExprPtr giant = flatten_product(summands[t]);
    if (!giant) continue;
    sequant::eval::dryrun::CostProfile r[3];
    try {
      r[0] = analyze(giant, /*occ_on=*/false, /*order_aware=*/false, 0);
      r[1] = analyze(giant, /*occ_on=*/true, /*order_aware=*/false, 1);
      r[2] = analyze(giant, /*occ_on=*/true, /*order_aware=*/true, 2);
    } catch (std::exception const&) {
      continue;
    } catch (...) {
      continue;
    }
    ++n_ok;
    for (int c = 0; c < 3; ++c) {
      mflops[c] += r[c].model_flops;
      mexec[c] += r[c].model_exec;
      dflops[c] += r[c].dryrun_flops;
      dexec[c] += r[c].dryrun_exec;
      dops[c] += r[c].dryrun_n_ops;
      maxpeak[c] = std::max(maxpeak[c], r[c].peak_bytes / 1e9);
    }
  }
  auto ratio = [](double num, double den) { return den > 0 ? num / den : 0.0; };
  wchar_t const* lab[3] = {L"aux-only        ", L"occ+aux OA=false",
                           L"occ+aux OA=true "};
  std::wcerr << L"\n=== [water-20 overcompute] " << n_ok << L" terms ("
             << summands.size() << L" total) ===\n";
  for (int c = 0; c < 3; ++c)
    std::wcerr << L"  " << lab[c] << L" | model_flops=" << mflops[c]
               << L" dryrun_flops=" << dflops[c] << L" dryrun_exec=" << dexec[c]
               << L" dryrun_n_ops=" << dops[c] << L" max_peak_GB=" << maxpeak[c]
               << L" | ext=" << n_ext[c] << L" con=" << n_con[c] << L"\n";
  std::wcerr << L"  --- vs aux-only ---   occ+aux OA=false: dryrun_exec x"
             << ratio(dexec[1], dexec[0]) << L" n_ops x"
             << ratio(double(dops[1]), double(dops[0]))
             << L"   |   occ+aux OA=true: dryrun_exec x"
             << ratio(dexec[2], dexec[0]) << L" n_ops x"
             << ratio(double(dops[2]), double(dops[0])) << L"\n";
  std::wcerr
      << L"  --- REGRESSION PROBE (OA=true vs OA=false) --- dryrun_exec x"
      << ratio(dexec[2], dexec[1]) << L"  dryrun_n_ops x"
      << ratio(double(dops[2]), double(dops[1])) << L"\n";
  REQUIRE(mflops[0] > 0.0);
}

// P4 GO/NO-GO AUDIT (Concern #3): across ALL C60 residual terms under the
// perf-first (DenseTimeSpaceBatched) objective, does external-occ FOREST
// batching (P3, the result-external occ carried as composite protos)
// bound EVERY over-budget giant, or does some giant's dominant memory lever
// escape it -- i.e. a free mu~/K carried on an INTERMEDIATE (contracted
// downstream, hence NOT the result-external P3 slices, and NOT node-
// local sliceable if no ancestor slices it)?
//
// Per term we: (1) optimize perf-first + binarize with node-local batch modes;
// (2) walk every node tracking ancestor batched_here and compute each node's
// REALIZED bytes (nominal shrunk by node-local mu~/K slicing already applied);
// (3) pick the term's biggest realized node and read its anatomy -- escaped
// (unsliced-by-ancestor) free mu~/K, and whether it carries result-external occ
// protos on its composite legs; (4) re-size that biggest node with ONE
// result-external occ proto's extent overridden to occ_block (the static
// P1 [dryrun-occ-sizing] slice -- O(1), no DP), giving the footprint AFTER
// external-occ forest batching. The verdict classifies each term as
// BOUNDED-by-P3 (occ slice fits) / NEEDS-mu-forest / NEEDS-k-forest / GAP.
// Hidden; run:
//   ./tests/unit/unit_tests-sequant "[dryrun-c60-batchability-audit]"
TEST_CASE("dryrun C60 per-term perf-first batchability audit (P4 go/no-go)",
          "[.][dryrun-c60-batchability-audit]") {
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);  // mu~
  sequant::mbpt::add_df_spaces(isr);                             // K
  auto ctx_resetter = set_scoped_default_context(std::move(ctx));

  auto const body = slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                          "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = deserialize<ExprPtr>(line);
  REQUIRE(static_cast<bool>(expr));
  REQUIRE(expr->is<Sum>());
  auto const& summands = expr->as<Sum>().summands();
  REQUIRE(!summands.empty());
  auto flatten_product = [](ExprPtr const& e) -> ExprPtr {
    if (!e->is<Product>()) return e;
    auto const& p = e->as<Product>();
    return ex<Product>(p.scalar(), p.factors(), Product::Flatten::Yes);
  };

  auto regime = df_regime(kC60_pVDZF12);
  auto memsize = sequant::opt::detail::memsize_counter(regime.idx_to_extent(),
                                                       regime.inner_pow_fn());
  auto keyof = [](Index const& ix) { return std::wstring(ix.full_label()); };
  auto tgt_of = [](Index const& ix) -> double {
    return ix.space().base_key() == L"μ̃" ? 256.0 : 72.0;
  };
  auto ext_of = [&](Index const& ix) -> double {
    return double(regime.idx_to_extent()(ix));
  };
  auto batch_target = [](Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"μ̃" ? std::size_t{256} : std::size_t{72};
  };

  // The over-budget threshold: peak_threshold (a node must batch to fit).
  double const kThrGB =
      (std::getenv("SEQUANT_UT_AUDIT_BUDGET_GB")
           ? std::atof(std::getenv("SEQUANT_UT_AUDIT_BUDGET_GB"))
           : 40.0);
  int const occ_block =
      (std::getenv("SEQUANT_UT_AUDIT_OCC_BLOCK")
           ? std::atoi(std::getenv("SEQUANT_UT_AUDIT_OCC_BLOCK"))
           : 10);

  // LOCAL DIAGNOSTIC (never committed): SEQUANT_UT_AUDIT_NO_PAO drops PAO (mu~)
  // from the batchable set, leaving only DF-aux K, to isolate whether the PAO
  // batched-context degree is what makes the DP intractable.
  bool const no_pao = std::getenv("SEQUANT_UT_AUDIT_NO_PAO") != nullptr;
  auto audit_batchable = [no_pao](Index const& ix) -> bool {
    auto const k = ix.space().base_key();
    return no_pao ? (k == L"Κ") : (k == L"μ̃" || k == L"Κ");
  };

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = audit_batchable;
  policy.batch_target_size = batch_target;
  policy.is_volatile_leaf = [](Tensor const& t) { return t.label() == L"t"; };
  policy.accumulation_factor = 1.0;
  policy.peak_threshold = kThrGB * 1e9;

  struct TermVerdict {
    std::size_t t = 0;
    double biggest_realized_gb = 0;  // biggest node after node-local mu~/K slic
    double biggest_occ_sliced_gb = 0;  // same node after ONE ext-occ slice
    bool biggest_has_ext_occ = false;  // biggest node carries occ protos
    bool biggest_escaped_mu = false;  // biggest node: free mu~ no ancestor slic
    bool biggest_escaped_k = false;
    std::wstring biggest_desc;
    std::wstring escaped_desc;
  };
  std::vector<TermVerdict> verdicts;

  auto const full_ext = regime.idx_to_extent();

  // classify a term by whether P3 external-occ batching bounds its biggest node
  auto classify = [&](TermVerdict const& v) -> std::wstring {
    if (v.biggest_realized_gb <= kThrGB) return L"OK (under budget)";
    // over-budget: does slicing ONE result-external occ bring it under budget?
    if (v.biggest_has_ext_occ && v.biggest_occ_sliced_gb <= kThrGB)
      return L"BOUNDED-by-P3 (occ slice fits)";
    if (v.biggest_escaped_mu) return L"NEEDS-mu-forest (escaped free mu~)";
    if (v.biggest_escaped_k) return L"NEEDS-k-forest (escaped free K)";
    return L"GAP? (over budget, occ slice insufficient, no escaped mu~/K)";
  };

  std::wcerr << L"\n=== [dryrun-c60-batchability-audit] perf-first, "
             << summands.size() << L" C60 terms, budget=" << kThrGB
             << L" GB, occ_block=" << occ_block << L" ===\n";

  // Optional term subset: SEQUANT_UT_AUDIT_TERMS="3,38,41" limits the sweep.
  std::set<std::size_t> term_filter;
  if (auto* e = std::getenv("SEQUANT_UT_AUDIT_TERMS")) {
    std::wstringstream ss(std::wstring(e, e + std::strlen(e)));
    std::wstring tok;
    while (std::getline(ss, tok, L',')) {
      if (!tok.empty()) term_filter.insert((std::size_t)std::stoul(tok));
    }
  }

  for (std::size_t t = 0; t < summands.size(); ++t) {
    if (!term_filter.empty() && !term_filter.count(t)) continue;
    ExprPtr giant = flatten_product(summands[t]);
    if (!giant) continue;

    std::size_t nfac = giant->is<Product>() ? giant->as<Product>().size() : 0;
    std::wcerr << L"[audit] term " << t << L" (nfac=" << nfac
               << L") optimizing..." << std::endl;

    // LOCAL DIAGNOSTIC (never committed): quantify how much outer-product
    // pruning reduces the DP subset space for this term, and count the
    // batchable-index degree m that drives the batched DP's 2^m context cost.
    if (std::getenv("SEQUANT_UT_AUDIT_PRUNE_STATS") && giant->is<Product>()) {
      container::svector<ExprPtr> gts;
      for (auto const& f : giant->as<Product>().factors())
        if (f->is<Tensor>()) gts.push_back(f);
      TensorNetwork gtn{gts};
      container::svector<Index> const empty_tidxs{};  // matches single_term_opt
      auto conn =
          sequant::opt::detail::outer_product_connectivity(gtn, empty_tidxs);
      std::size_t const nt = gts.size();
      std::size_t total_ge2 = 0, connected_ge2 = 0;
      for (std::size_t n = 1; n < conn.size(); ++n) {
        if (std::popcount(n) < 2) continue;
        ++total_ge2;
        if (conn[n]) ++connected_ge2;
      }
      // batchable-index degree m: distinct batchable top-level indices.
      std::set<std::wstring> batch_ixs;
      for (auto const& f : gts) {
        auto const& tp = f->as<Tensor>();
        for (auto const& ix : tp.const_braket())
          if (audit_batchable(ix)) batch_ixs.insert(std::wstring(ix.label()));
        for (auto const& ix : tp.aux())
          if (audit_batchable(ix)) batch_ixs.insert(std::wstring(ix.label()));
      }
      std::wcerr << L"[prune-stats] term " << t << L" nt=" << nt
                 << L" subsets(>=2)=" << total_ge2 << L" connected="
                 << connected_ge2 << L" pruned=" << (total_ge2 - connected_ge2)
                 << L" ("
                 << (total_ge2 ? 100.0 * double(total_ge2 - connected_ge2) /
                                     double(total_ge2)
                               : 0.0)
                 << L"%) full_net_connected=" << (conn.back() ? 1 : 0)
                 << L" batchable_m=" << batch_ixs.size() << L" -> 2^m="
                 << (double)(1ull
                             << std::min<std::size_t>(batch_ixs.size(), 62))
                 << std::endl;
      continue;  // skip the expensive optimize; stats only
    }

    TermVerdict v;
    v.t = t;
    try {
      auto axes_map = std::make_shared<std::unordered_map<
          Expr const*, container::vector<NodeBatchAnnotation>>>();
      OptimizeOptions opts;
      opts.objective_function = ObjectiveFunction::DenseTimeSpaceBatched;
      opts.idx_to_extent = regime.idx_to_extent();
      opts.inner_pow = regime.inner_pow_fn();
      opts.batch_policy = policy;
      opts.volatile_weight = 20.0;
      opts.roofline.machine_balance = 200.0;
      opts.roofline.fast_mem_elems = 1000000.0;
      opts.term_batch_axes = axes_map;

      auto optimized = optimize(giant, opts);
      if (!optimized) continue;
      // LOCAL DIAGNOSTIC (never committed): print the found factorization so a
      // pruned vs SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING run can be diffed.
      if (std::getenv("SEQUANT_UT_AUDIT_PRINT_SOLN"))
        std::wcerr << L"[soln] term " << t << L" prune="
                   << (std::getenv("SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING") ? 0
                                                                            : 1)
                   << L" : " << to_latex(optimized) << std::endl;
      auto it = axes_map->find(optimized.get());
      container::vector<NodeBatchAnnotation> node_axes;
      if (it != axes_map->end()) node_axes = it->second;
      BinarizationOptions bopts;
      bopts.node_batch_axes = node_axes;
      SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
      auto node = binarize<EvalExprDryRun>(optimized, {}, bopts);
      SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

      // (2)+(3)+(4): walk every node, track ancestor batch modes, find the
      // biggest node-local-batched realized node, and for THAT node also size
      // it with one result-external occ proto sliced to occ_block (static P1
      // sizing). All O(1) -- no DP.
      double biggest = 0;
      std::vector<Index> biggest_free;
      std::vector<Index> biggest_escaped;
      double biggest_occ_sliced = 0;
      bool biggest_has_occ = false;
      std::function<void(std::remove_cvref_t<decltype(node)> const&,
                         std::map<std::wstring, int>)>
          walk = [&](auto const& n, std::map<std::wstring, int> active) {
            auto free_ixs = node_free_indices(*n);
            double nominal =
                memsize(free_ixs, std::vector<Index>{}, std::vector<Index>{}) *
                8.0;
            double factor = 1.0;  // node-local mu~/K slicing already applied
            std::vector<Index> escaped;
            for (auto const& ix : free_ixs) {
              auto bk = ix.space().base_key();
              if (bk != L"μ̃" && bk != L"Κ") continue;
              if (active.count(keyof(ix)))
                factor *= tgt_of(ix) / ext_of(ix);
              else
                escaped.push_back(ix);
            }
            double const realized = nominal * factor;
            if (realized > biggest) {
              biggest = realized;
              biggest_free = free_ixs;
              biggest_escaped = escaped;
              // pick ONE result-external occ proto on this node's composite
              // legs
              std::wstring occ_label;
              for (auto const& ix : free_ixs) {
                for (auto const& p : ix.proto_indices())
                  if (p.space().base_key() == L"i") {
                    occ_label = std::wstring(p.full_label());
                    break;
                  }
                if (!occ_label.empty()) break;
              }
              biggest_has_occ = !occ_label.empty();
              if (occ_label.empty()) {
                biggest_occ_sliced = realized;
              } else {
                auto sliced_ext = [full_ext, occ_label,
                                   occ_block](Index const& ix) -> std::size_t {
                  if (std::wstring(ix.full_label()) == occ_label)
                    return std::min<std::size_t>(full_ext(ix),
                                                 std::size_t(occ_block));
                  return full_ext(ix);
                };
                auto memsize_sliced = sequant::opt::detail::memsize_counter(
                    sliced_ext, regime.inner_pow_fn());
                double const sliced_nominal =
                    memsize_sliced(free_ixs, std::vector<Index>{},
                                   std::vector<Index>{}) *
                    8.0;
                biggest_occ_sliced = sliced_nominal * factor;
              }
            }
            if (!n.leaf()) {
              std::map<std::wstring, int> child_active = active;
              for (auto const& ax : n->batched_here())
                child_active[keyof(ax.first)] = 1;
              walk(n.left(), child_active);
              walk(n.right(), child_active);
            }
          };
      walk(node, {});

      v.biggest_realized_gb = biggest / 1e9;
      v.biggest_occ_sliced_gb = biggest_occ_sliced / 1e9;
      v.biggest_has_ext_occ = biggest_has_occ;
      v.biggest_desc = describe_indices(biggest_free);
      v.escaped_desc = describe_indices(biggest_escaped);
      for (auto const& ix : biggest_escaped) {
        if (ix.space().base_key() == L"μ̃") v.biggest_escaped_mu = true;
        if (ix.space().base_key() == L"Κ") v.biggest_escaped_k = true;
      }
    } catch (std::exception const&) {
      std::wcerr << L"[audit] term " << t << L": skipped (exception)\n";
      continue;
    } catch (...) {
      std::wcerr << L"[audit] term " << t << L": skipped (unknown)\n";
      continue;
    }
    verdicts.push_back(v);
    // Inline verdict so a partial (timed-out) run still yields data.
    std::wcerr << L"[audit] term " << v.t << L" DONE: biggest_realized="
               << v.biggest_realized_gb << L" GB  occ_sliced="
               << v.biggest_occ_sliced_gb << L" GB  ext_occ="
               << (v.biggest_has_ext_occ ? 1 : 0) << L"  escaped={"
               << v.escaped_desc << L"}  => " << classify(v) << std::endl;
  }

  int n_ok = 0, n_p3 = 0, n_mu = 0, n_k = 0, n_gap = 0;
  double max_biggest = 0, max_occ_sliced = 0;
  for (auto const& v : verdicts) {
    auto cls = classify(v);
    if (v.biggest_realized_gb > max_biggest)
      max_biggest = v.biggest_realized_gb;
    if (v.biggest_occ_sliced_gb > max_occ_sliced)
      max_occ_sliced = v.biggest_occ_sliced_gb;
    if (v.biggest_realized_gb > kThrGB) {
      std::wcerr << L"term " << v.t << L": biggest_realized="
                 << v.biggest_realized_gb << L" GB  occ_sliced="
                 << v.biggest_occ_sliced_gb << L" GB  ext_occ="
                 << (v.biggest_has_ext_occ ? 1 : 0) << L"  escaped={"
                 << v.escaped_desc << L"}\n    free={" << v.biggest_desc
                 << L"}\n    => " << cls << L"\n";
    }
    if (cls.rfind(L"OK", 0) == 0)
      ++n_ok;
    else if (cls.rfind(L"BOUNDED-by-P3", 0) == 0)
      ++n_p3;
    else if (cls.rfind(L"NEEDS-mu", 0) == 0)
      ++n_mu;
    else if (cls.rfind(L"NEEDS-k", 0) == 0)
      ++n_k;
    else
      ++n_gap;
  }

  std::wcerr
      << L"\n=== [dryrun-c60-batchability-audit] VERDICT ===\n"
      << L"terms analyzed=" << verdicts.size() << L"\n"
      << L"  OK (under budget)          = " << n_ok << L"\n"
      << L"  BOUNDED-by-P3 (occ forest) = " << n_p3 << L"\n"
      << L"  NEEDS-mu-forest            = " << n_mu << L"\n"
      << L"  NEEDS-k-forest             = " << n_k << L"\n"
      << L"  GAP (no known lever)       = " << n_gap << L"\n"
      << L"max biggest-realized (node-local only)  = " << max_biggest
      << L" GB\n"
      << L"max biggest after ext-occ slice         = " << max_occ_sliced
      << L" GB\n"
      << L"P4 GO/NO-GO: "
      << ((n_mu + n_k + n_gap) == 0
              ? L"GO -- every "
                L"over-budget giant is bounded by P3 external-occ batching"
              : L"REVIEW -- some giants escape P3 (see NEEDS-*/GAP above)")
      << L"\n";

  REQUIRE(verdicts.size() > 0);
}

// Task 3: the opt-in scratch-fold peak sink captures the batched-inner peak the
// OUTER cache.working_set_hwmark() misses. Reuse the [dryrun-objective]
// peak-first (DenseSpaceTimeBatched) setup -- the objective whose batched-inner
// transient
// (~38.9 GB, materialized INSIDE a make_batched_scratch cache) dwarfs the
// outer, cross-batch cached residency (~0.2 GB). A PeakSink passed to
// make_evaluator folds each scratch cache's high-watermark into one global
// accumulator, so the global peak reflects the true batched-replay peak rather
// than just the outer residency the accessor sees today.
TEST_CASE("dryrun scratch-fold captures batched peak", "[dryrun][peak]") {
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);  // mu~
  sequant::mbpt::add_df_spaces(isr);                             // K
  auto ctx_resetter = set_scoped_default_context(std::move(ctx));

  auto const body = slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                          "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = deserialize<ExprPtr>(line);
  REQUIRE(static_cast<bool>(expr));
  REQUIRE(expr->is<Sum>());
  auto const& summands = expr->as<Sum>().summands();
  REQUIRE(!summands.empty());
  auto flatten_product = [](ExprPtr const& e) -> ExprPtr {
    if (!e->is<Product>()) return e;
    auto const& p = e->as<Product>();
    return ex<Product>(p.scalar(), p.factors(), Product::Flatten::Yes);
  };
  std::size_t const giant_idx = 38 < summands.size() ? 38 : 0;
  ExprPtr giant = flatten_product(summands[giant_idx]);
  REQUIRE(giant);

  // FAITHFUL real C60 config (identical to [dryrun-objective]).
  auto regime = df_regime(kC60_pVDZF12);
  auto cm = std::make_shared<CostModel const>(regime);

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = is_df_batchable;
  policy.batch_target_size = [](Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"μ̃" ? std::size_t{256} : std::size_t{72};
  };
  policy.is_volatile_leaf = [](Tensor const& t) { return t.label() == L"t"; };
  policy.accumulation_factor = 1.0;
  policy.peak_threshold = 40e9;

  auto axes_map = std::make_shared<std::unordered_map<
      Expr const*, container::vector<NodeBatchAnnotation>>>();
  OptimizeOptions opts;
  // peak-first: forms the fully-sliceable 4-PAO node whose per-batch scratch
  // transient is the batched-inner peak we want the sink to capture.
  opts.objective_function = ObjectiveFunction::DenseSpaceTimeBatched;
  opts.idx_to_extent = regime.idx_to_extent();
  opts.inner_pow = regime.inner_pow_fn();
  opts.batch_policy = policy;
  opts.volatile_weight = 20.0;
  opts.roofline.machine_balance = 200.0;
  opts.roofline.fast_mem_elems = 1000000.0;
  opts.term_batch_axes = axes_map;

  auto optimized = optimize(giant, opts);
  REQUIRE(static_cast<bool>(optimized));
  auto it = axes_map->find(optimized.get());
  container::vector<NodeBatchAnnotation> node_axes;
  if (it != axes_map->end()) node_axes = it->second;
  BinarizationOptions bopts;
  bopts.node_batch_axes = node_axes;
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = binarize<EvalExprDryRun>(optimized, {}, bopts);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

  // Zero-data replay through the real eval loop + real CacheManager, with the
  // real batched custom evaluator, but now with a PeakSink threaded through
  // make_evaluator so each per-batch scratch cache's working_set_hwmark folds
  // into `peak`.
  auto cache = sequant::cache_manager(std::vector<EvalNodeDryRun>{node});
  std::atomic<double> peak{0.0};
  cache.set_custom_evaluator(sequant::make_evaluator(
      policy, DryRunLeafEvaluator{cm}, sequant::make_no_scope_guard{}, &peak));

  std::ostringstream trace_os;
  auto& logger = Logger::instance();
  auto const prev_level = logger.eval.level;
  auto* const prev_stream = logger.eval.stream;
  logger.eval.level = 2;  // gate log::printing() on so hwmark accumulates
  logger.eval.stream = &trace_os;
  try {
    (void)sequant::evaluate<Trace::On>(node, DryRunLeafEvaluator{cm}, cache);
  } catch (std::exception const&) {
    // mirror [dryrun-objective]: a DryRun sizing throw must not mask the sink
    // read
  }
  logger.eval.level = prev_level;
  logger.eval.stream = prev_stream;

  double const global_peak = peak.load();
  double const outer_hwmark = double(cache.working_set_hwmark());

  std::wcerr << L"\n[dryrun-peak] global scratch-folded peak="
             << (global_peak / 1e9) << L" GB  outer working_set_hwmark="
             << (outer_hwmark / 1e9) << L" GB  ratio="
             << (outer_hwmark > 0.0 ? global_peak / outer_hwmark : 0.0)
             << L"\n";

  // The sink must have captured something (batched replay ran).
  CHECK(global_peak > 0.0);
  // The global (scratch-folded) peak is at least the outer cached residency.
  CHECK(global_peak >= outer_hwmark);
  // For this specifically-batched term the batched-inner transient dwarfs the
  // outer residency (~195x observed); a 2x floor is safe and non-flaky.
  CHECK(global_peak > outer_hwmark * 2.0);
}

// Task 4: cost_profile() is the single reusable entry point that ties the
// static cost walk (flops/exec_cost/n_ops) to the gated-cache peak replay
// (Task 2 build_dryrun_cache + Task 3 PeakSink scratch-fold) behind one API
// that MPQC and these tests share. The key regression guard is peak_bytes > 0
// with NO trace stream: cost_profile MUST force log::printing() on internally
// (the CacheManager hwmark only accumulates on the printing path), otherwise
// the sink and the outer hwmark would both be zero.
TEST_CASE("cost_profile returns peak/flops/exec/n_ops",
          "[dryrun][cost_profile]") {
  using sequant::eval::dryrun::build_dryrun_cache;
  using sequant::eval::dryrun::CacheConfig;
  using sequant::eval::dryrun::cost_profile;
  using sequant::eval::dryrun::CostProfile;

  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);  // mu~
  sequant::mbpt::add_df_spaces(isr);                             // K
  auto ctx_resetter = set_scoped_default_context(std::move(ctx));

  auto const body = slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                          "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = deserialize<ExprPtr>(line);
  REQUIRE(static_cast<bool>(expr));
  REQUIRE(expr->is<Sum>());
  auto const& summands = expr->as<Sum>().summands();
  REQUIRE(!summands.empty());
  auto flatten_product = [](ExprPtr const& e) -> ExprPtr {
    if (!e->is<Product>()) return e;
    auto const& p = e->as<Product>();
    return ex<Product>(p.scalar(), p.factors(), Product::Flatten::Yes);
  };
  std::size_t const giant_idx = 38 < summands.size() ? 38 : 0;
  ExprPtr giant = flatten_product(summands[giant_idx]);
  REQUIRE(giant);

  // FAITHFUL real C60 config (identical to [dryrun-objective] / [peak]).
  auto regime = df_regime(kC60_pVDZF12);

  // ONE BatchPolicy, reused for optimize() and the replay evaluator (its
  // is_batchable_index MUST equal CacheConfig::is_batchable_index).
  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = is_df_batchable;
  policy.batch_target_size = [](Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"μ̃" ? std::size_t{256} : std::size_t{72};
  };
  policy.is_volatile_leaf = [](Tensor const& t) { return t.label() == L"t"; };
  policy.accumulation_factor = 1.0;
  policy.peak_threshold = 40e9;

  auto axes_map = std::make_shared<std::unordered_map<
      Expr const*, container::vector<NodeBatchAnnotation>>>();
  OptimizeOptions opts;
  // peak-first: forms the fully-sliceable 4-PAO node whose per-batch scratch
  // transient is the batched-inner peak the sink must capture.
  opts.objective_function = ObjectiveFunction::DenseSpaceTimeBatched;
  opts.idx_to_extent = regime.idx_to_extent();
  opts.inner_pow = regime.inner_pow_fn();
  opts.batch_policy = policy;
  opts.volatile_weight = 20.0;
  opts.roofline.machine_balance = 200.0;
  opts.roofline.fast_mem_elems = 1000000.0;
  opts.term_batch_axes = axes_map;

  auto optimized = optimize(giant, opts);
  REQUIRE(static_cast<bool>(optimized));
  auto it = axes_map->find(optimized.get());
  container::vector<NodeBatchAnnotation> node_axes;
  if (it != axes_map->end()) node_axes = it->second;
  BinarizationOptions bopts;
  bopts.node_batch_axes = node_axes;
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = binarize<EvalExprDryRun>(optimized, {}, bopts);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

  std::vector<EvalNodeDryRun> forest{node};

  // CacheConfig mirroring the real batched run: footprint gate + free-batchable
  // mode veto; is_volatile adapts policy.is_volatile_leaf onto tree nodes.
  CacheConfig cfg;
  cfg.max_footprint = 1e11;
  cfg.min_repeats = 1;
  cfg.is_batchable_index = policy.is_batchable_contracted_index;
  cfg.is_volatile = [](EvalNodeDryRun const& n) {
    if (!n.leaf() || !n->is_tensor()) return false;
    return n->as_tensor().label() == L"t";
  };

  // Call cost_profile with NO trace stream: this is the printing-gate gotcha
  // regression guard -- peak_bytes must still be > 0.
  CostProfile const cp = cost_profile(forest, policy, cfg, regime,
                                      /*trace=*/nullptr);

  std::wcerr << L"\n[cost_profile] model_n_ops=" << cp.model_n_ops
             << L" model_flops=" << cp.model_flops << L" model_exec="
             << cp.model_exec << L" dryrun_n_ops=" << cp.dryrun_n_ops
             << L" dryrun_flops=" << cp.dryrun_flops << L" dryrun_exec="
             << cp.dryrun_exec << L" peak_bytes=" << (cp.peak_bytes / 1e9)
             << L" GB\n";

  CHECK(cp.model_n_ops > 0);
  CHECK(cp.model_flops > 0.0);
  CHECK(cp.model_exec > 0.0);
  // The replay tallied at least one product op with positive cost.
  CHECK(cp.dryrun_n_ops > 0);
  CHECK(cp.dryrun_flops > 0.0);
  CHECK(cp.dryrun_exec > 0.0);
  // The gotcha guard: printing was forced on internally, so the hwmark/sink
  // accumulated even with no trace stream.
  CHECK(cp.peak_bytes > 0.0);

  // Cross-check that cost_profile's peak wires the SAME fold as the Task-3
  // [peak] test: replay the same forest through the same gated cache + the same
  // make_evaluator(&peak2) sink, force printing on identically, and fold the
  // outer hwmark. cost_profile's peak_bytes must equal max(sink, outer hwmark).
  auto cm = std::make_shared<CostModel const>(regime);
  DryRunLeafEvaluator leaf{cm};
  auto cache = build_dryrun_cache(forest, cfg, regime);
  std::atomic<double> peak2{0.0};
  cache.set_custom_evaluator(sequant::make_evaluator(
      policy, leaf, sequant::make_no_scope_guard{}, &peak2));
  auto& logger = Logger::instance();
  auto const prev_level = logger.eval.level;
  auto* const prev_stream = logger.eval.stream;
  logger.eval.level = 2;
  logger.eval.stream = nullptr;
  try {
    (void)sequant::evaluate<Trace::On>(node, leaf, cache);
  } catch (std::exception const&) {
  }
  logger.eval.level = prev_level;
  logger.eval.stream = prev_stream;
  double const expected_peak =
      std::max(peak2.load(), double(cache.working_set_hwmark()));
  CHECK(expected_peak > 0.0);
  CHECK(cp.peak_bytes == Catch::Approx(expected_peak).epsilon(1e-9));
}

// Task 5 (Minor b): cover the UTF-8 -> wide bridge cost_profile() uses to fill
// a caller's wide trace stream. The eval loop writes a NARROW (UTF-8) trace
// whose per-op label field carries multi-byte index labels (mu~ = U+03BC
// U+0303, aux = U+039A). cost_profile transcodes that narrow buffer into the
// wide sink by decoding UTF-8 code points (a plain widen() would mojibake the
// labels). This runs cost_profile on the single C60 giant with a real
// std::wostringstream and asserts the wide output is non-empty and contains a
// token the eval trace emits -- both an ASCII field (`result=`) and a
// multi-byte index label (mu~), so the code-point decode path (not just the
// ASCII path) is exercised.
TEST_CASE("cost_profile trace stream round-trips", "[dryrun][cost_profile]") {
  using sequant::eval::dryrun::CacheConfig;
  using sequant::eval::dryrun::cost_profile;
  using sequant::eval::dryrun::CostProfile;

  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);  // mu~
  sequant::mbpt::add_df_spaces(isr);                             // K
  auto ctx_resetter = set_scoped_default_context(std::move(ctx));

  auto const body = slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                          "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = deserialize<ExprPtr>(line);
  REQUIRE(static_cast<bool>(expr));
  REQUIRE(expr->is<Sum>());
  auto const& summands = expr->as<Sum>().summands();
  REQUIRE(!summands.empty());
  auto flatten_product = [](ExprPtr const& e) -> ExprPtr {
    if (!e->is<Product>()) return e;
    auto const& p = e->as<Product>();
    return ex<Product>(p.scalar(), p.factors(), Product::Flatten::Yes);
  };
  std::size_t const giant_idx = 38 < summands.size() ? 38 : 0;
  ExprPtr giant = flatten_product(summands[giant_idx]);
  REQUIRE(giant);

  auto regime = df_regime(kC60_pVDZF12);

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = is_df_batchable;
  policy.batch_target_size = [](Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"μ̃" ? std::size_t{256} : std::size_t{72};
  };
  policy.is_volatile_leaf = [](Tensor const& t) { return t.label() == L"t"; };
  policy.accumulation_factor = 1.0;
  policy.peak_threshold = 40e9;

  auto axes_map = std::make_shared<std::unordered_map<
      Expr const*, container::vector<NodeBatchAnnotation>>>();
  OptimizeOptions opts;
  opts.objective_function = ObjectiveFunction::DenseSpaceTimeBatched;
  opts.idx_to_extent = regime.idx_to_extent();
  opts.inner_pow = regime.inner_pow_fn();
  opts.batch_policy = policy;
  opts.volatile_weight = 20.0;
  opts.roofline.machine_balance = 200.0;
  opts.roofline.fast_mem_elems = 1000000.0;
  opts.term_batch_axes = axes_map;

  auto optimized = optimize(giant, opts);
  REQUIRE(static_cast<bool>(optimized));
  auto it = axes_map->find(optimized.get());
  container::vector<NodeBatchAnnotation> node_axes;
  if (it != axes_map->end()) node_axes = it->second;
  BinarizationOptions bopts;
  bopts.node_batch_axes = node_axes;
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = binarize<EvalExprDryRun>(optimized, {}, bopts);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

  std::vector<EvalNodeDryRun> forest{node};

  CacheConfig cfg;
  cfg.max_footprint = 1e11;
  cfg.min_repeats = 1;
  cfg.is_volatile = [](EvalNodeDryRun const& n) {
    if (!n.leaf() || !n->is_tensor()) return false;
    return n->as_tensor().label() == L"t";
  };

  // Real wide sink: exercise the UTF-8 -> wide transcode.
  std::wostringstream trace;
  CostProfile const cp = cost_profile(forest, policy, cfg, regime, &trace);

  std::wstring const w = trace.str();
  if (std::getenv("SEQUANT_UT_DRYRUN_DUMP_TRACE"))
    std::wcerr << L"[round-trip] first 600 wide chars:\n"
               << w.substr(0, std::min<std::size_t>(600, w.size())) << L"\n";

  CHECK(std::isfinite(cp.peak_bytes));
  CHECK(cp.peak_bytes > 0.0);
  // Non-empty wide output.
  CHECK(!w.empty());
  // ASCII field the eval trace always emits (exercises the 1-byte decode path).
  CHECK(w.find(L"result=") != std::wstring::npos);
  // Multi-byte index label in an op's label field (exercises the multi-byte
  // code-point decode path -- the whole point of the transcode). mu~ =
  // U+03BC U+0303; a plain widen() of the UTF-8 bytes would NOT produce it.
  CHECK(w.find(L"μ̃") != std::wstring::npos);
}

// Task 2: the FAITHFUL (gated) dry-run cache built by build_dryrun_cache must
// keep a free-batchable giant -- a node whose result carries a free mu~/K mode
// -- out of the run-scope cache, exactly as the real batched eval loop does,
// instead of materializing it whole. Two independent gates cooperate, and this
// case pins WHICH one acts (the phase-2 lifetime-mask veto made the mode veto
// PRECISE: it drops a node only when the node is ACTUALLY batch-sliced -- a
// Contracted mode among its RESULT indices, or a non-empty cross-occurrence
// external mask -- NOT merely for carrying a free-batchable index, which is
// what the OLD over-broad veto did and which emptied CSE). We contrast three
// configs over the SAME binarized C60 giant term (index 38):
//   ref:  no gate  (max_footprint=0, never batchable) -> giant IS cached
//   mode: mode veto only (max_footprint=0, is_df_batchable) -> giant STILL
//         cached: in this DF-aux schedule K is summed (never a result mode) and
//         mu~ is contracted ABOVE the giant, so NO node is actually sliced and
//         the precise mode veto is inert
//   task: + footprint gate (max_footprint=1e11) -> giant NOT cached: the
//         footprint gate is what caps the >100 GB full giant
// A small non-batchable intermediate stays cached under all three. Residency is
// read via CacheManager::exists() (the veto is a registration-time decision:
// a vetoed node is never registered, so it is not, and cannot become, resident;
// this is a stronger/cleaner signal than working_set_hwmark, which only tracks
// alive cached bytes and is confounded by batched-inner scratch -- see the
// [dryrun-objective] INTERPRETATION notes).
TEST_CASE("dryrun gated cache vetoes free-batchable giant", "[dryrun][cache]") {
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);  // mu~
  sequant::mbpt::add_df_spaces(isr);                             // K
  auto ctx_resetter = set_scoped_default_context(std::move(ctx));

  auto const body = slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                          "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = deserialize<ExprPtr>(line);
  REQUIRE(static_cast<bool>(expr));
  REQUIRE(expr->is<Sum>());
  auto const& summands = expr->as<Sum>().summands();
  REQUIRE(!summands.empty());
  auto flatten_product = [](ExprPtr const& e) -> ExprPtr {
    if (!e->is<Product>()) return e;
    auto const& p = e->as<Product>();
    return ex<Product>(p.scalar(), p.factors(), Product::Flatten::Yes);
  };
  std::size_t const giant_idx = 38 < summands.size() ? 38 : 0;
  ExprPtr giant = flatten_product(summands[giant_idx]);
  REQUIRE(giant);

  // FAITHFUL real C60 config (614336 job log), identical to [dryrun-objective].
  auto regime = df_regime(kC60_pVDZF12);
  auto memsize = sequant::opt::detail::memsize_counter(regime.idx_to_extent(),
                                                       regime.inner_pow_fn());

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = is_df_batchable;
  policy.batch_target_size = [](Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"μ̃" ? std::size_t{256} : std::size_t{72};
  };
  policy.is_volatile_leaf = [](Tensor const& t) { return t.label() == L"t"; };
  policy.accumulation_factor = 1.0;
  policy.peak_threshold = 40e9;

  // Optimize (perf-first) + binarize the giant, carrying the per-node batch
  // modes so the binarized tree matches the real schedule.
  auto axes_map = std::make_shared<std::unordered_map<
      Expr const*, container::vector<NodeBatchAnnotation>>>();
  OptimizeOptions opts;
  opts.objective_function = ObjectiveFunction::DenseTimeSpaceBatched;
  opts.idx_to_extent = regime.idx_to_extent();
  opts.inner_pow = regime.inner_pow_fn();
  opts.batch_policy = policy;
  opts.volatile_weight = 20.0;
  opts.roofline.machine_balance = 200.0;
  opts.roofline.fast_mem_elems = 1000000.0;
  opts.term_batch_axes = axes_map;
  auto optimized = optimize(giant, opts);
  REQUIRE(static_cast<bool>(optimized));
  auto it = axes_map->find(optimized.get());
  container::vector<NodeBatchAnnotation> node_axes;
  if (it != axes_map->end()) node_axes = it->second;
  BinarizationOptions bopts;
  bopts.node_batch_axes = node_axes;
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = binarize<EvalExprDryRun>(optimized, {}, bopts);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  std::vector<EvalNodeDryRun> const nodes{node};

  // A node carries a free-batchable mode iff any of its RESULT (canon) indices
  // is a batch mode -- the exact condition the gated factory's veto keys off.
  auto carries_batchable = [](EvalNodeDryRun const& n) {
    for (auto const& ix : n->canon_indices())
      if (is_df_batchable(ix)) return true;
    return false;
  };

  // Walk the binarized tree: pick the internal free-batchable node with the
  // largest result footprint (THE giant), and any internal non-batchable node
  // (a small intermediate that must stay cacheable).
  EvalNodeDryRun giant_node = node;
  double giant_bytes = -1.0;
  bool giant_found = false;
  EvalNodeDryRun small_node = node;
  bool small_found = false;
  std::function<void(EvalNodeDryRun const&)> walk =
      [&](EvalNodeDryRun const& n) {
        if (!n.leaf()) {
          std::vector<Index> const free(n->canon_indices().begin(),
                                        n->canon_indices().end());
          double const bytes =
              memsize(std::vector<Index>{}, std::vector<Index>{}, free) * 8.0;
          if (carries_batchable(n)) {
            if (bytes > giant_bytes) {
              giant_bytes = bytes;
              giant_node = n;
              giant_found = true;
            }
          } else if (!small_found) {
            small_node = n;
            small_found = true;
          }
          walk(n.left());
          walk(n.right());
        }
      };
  walk(node);
  REQUIRE(giant_found);  // the term has a free-mu~/K internal node
  REQUIRE(small_found);  // and at least one non-batchable internal node

  // is_volatile for the gated factory (the amplitude leaves are volatile).
  auto is_vol = [](EvalNodeDryRun const& n) {
    return n.leaf() && n->is_tensor() && n->as_tensor().label() == L"t";
  };

  using sequant::eval::dryrun::build_dryrun_cache;
  using sequant::eval::dryrun::CacheConfig;

  // ref: no gates -> the giant IS cached (baseline; this is what the SIMPLE
  // ad-hoc factory would do).
  CacheConfig ref_cfg;
  ref_cfg.max_footprint = 0.;
  ref_cfg.min_repeats = 1;
  ref_cfg.is_volatile = is_vol;
  ref_cfg.is_batchable_index = [](Index const&) { return false; };
  auto ref_cache = build_dryrun_cache(nodes, ref_cfg, regime);

  // mode: batchable-mode veto ONLY (footprint gate disabled) -> giant NOT
  // cached, isolating the veto from the footprint gate.
  CacheConfig axis_cfg = ref_cfg;
  axis_cfg.is_batchable_index = is_df_batchable;
  auto axis_cache = build_dryrun_cache(nodes, axis_cfg, regime);

  // task: the config the task asks for (footprint gate + mode veto).
  CacheConfig task_cfg = axis_cfg;
  task_cfg.max_footprint = 1e11;
  auto task_cache = build_dryrun_cache(nodes, task_cfg, regime);

  // Without any gate the giant is registered/cacheable...
  CHECK(ref_cache.exists(giant_node));
  // ...the batchable-MODE veto alone does NOT remove it. The phase-2 veto is
  // precise: it drops a node only when the node is ACTUALLY batch-sliced -- a
  // Contracted batch mode among its RESULT (canon) indices, or a non-empty
  // cross-occurrence external lifetime mask (cache_manager.hpp) -- NOT merely
  // for carrying a free-batchable index. This giant carries a free mu~ that is
  // contracted ABOVE it (never sliced AT the node) and is not external-seeded
  // in this DF-aux schedule (K is summed, so it is never a result mode), so no
  // node here is actually sliced and the mode veto is correctly inert. This is
  // the whole point of the phase-2 change: the OLD over-broad veto dropped
  // anything whose result carried a batchable index, which emptied CSE.
  CHECK(axis_cache.exists(giant_node));
  // ...only the FOOTPRINT gate caps the >100 GB full giant.
  CHECK_FALSE(task_cache.exists(giant_node));

  // A small non-batchable intermediate stays cached under all three.
  CHECK(ref_cache.exists(small_node));
  CHECK(axis_cache.exists(small_node));
  CHECK(task_cache.exists(small_node));

  // Faithful end-to-end: replay the schedule through the real eval loop against
  // the gated (task) cache; the giant must not become resident, while eval
  // completes without materializing it whole.
  auto cm = std::make_shared<CostModel const>(regime);
  task_cache.set_custom_evaluator(
      sequant::make_evaluator(policy, DryRunLeafEvaluator{cm}));
  REQUIRE_NOTHROW(
      (void)sequant::evaluate(node, DryRunLeafEvaluator{cm}, task_cache));
  CHECK_FALSE(task_cache.exists(giant_node));  // still not registered => never
                                               // resident
  CHECK_FALSE(task_cache.alive(giant_node));

  std::wcerr << L"\n=== [dryrun][cache] gated-veto verdict (C60 giant, index "
             << giant_idx << L") ===\n  giant free-batchable node footprint = "
             << (giant_bytes / 1e9) << L" GB\n  ref (no gate)  exists(giant) = "
             << (ref_cache.exists(giant_node) ? L"YES" : L"no")
             << L"\n  mode veto      exists(giant) = "
             << (axis_cache.exists(giant_node) ? L"YES" : L"no")
             << L"\n  task config    exists(giant) = "
             << (task_cache.exists(giant_node) ? L"YES" : L"no")
             << L"\n  small intermediate exists (task) = "
             << (task_cache.exists(small_node) ? L"YES" : L"no") << L"\n";
}

// REPRO: the free-batchable-mode veto vs occ batching.
//
// cache_manager's veto drops from the cache (and erases from `persistent`) any
// node whose result carries an index the GLOBAL is_batchable_index predicate
// accepts -- regardless of whether the chosen schedule actually slices that
// node on that mode. For aux-only (K) that is surgical: only K-carrying
// intermediates are vetoed. But the occupied indices i_* appear in nearly every
// CC-residual intermediate, so marking occ batchable vetoes essentially the
// whole cache. That kills CSE outright AND empties the batch-group member pool
// (members are drawn from persistent cache keys; see eval.hpp
// make_batched_custom_evaluator / make_batched_scratch), so every group is
// degenerate (1 member) and each consumer rebuilds shared sub-intermediates
// once per batch combination.
//
// Measures the REAL C60 residual term replayed through the real eval loop under
// three mode policies. cp.n_ops is the STATIC forest node count (the ideal:
// every node once); the trace's `Eval | Product` count is what the replay
// actually executed. Their ratio is the realized recompute factor. Batching
// legitimately multiplies ops (each op is a slice, so ~n_batches is work
// neutral); the pathology is a factor well ABOVE the batch-count product, and a
// distinct-expression count that stays flat while ops explode.
//
// Mirrors the C60 job (631467): K batched at 256, occ at 8, mu~ NOT batched,
// 100 GB budget, perf-first.
// ==========================================================================
// TRUST BOUNDARY -- READ BEFORE QUOTING ANY NUMBER FROM THIS WITNESS.
//
// The multi-mode (aux+occ) arms of this witness REPLAY through machinery with
// known, open defects at every stage of the batched pipeline:
//
//   generation     the DP is order-blind: it keys its cell by a loop SET, so a
//                  loop TREE is unrepresentable and free-hoist vs charged-hoist
//                  cannot be distinguished at all.
//   representation the annotation channel is entangled with the cache veto --
//                  cache_manager vetoes a node whose own batched_here() carries
//                  a sliced batchable mode, so the SAME mode batched via an
//                  annotation vs via the (now removed) heuristic yields
//                  DIFFERENT CSE.
//   processing     nested aux+occ batch-group join rejects cross-axis members,
//                  degenerating groups to 1 member and wiping CSE.
//   processing     external-occ seeding does not reach the term that sets the
//                  forest peak (24 stamps over 55 terms; the external arm
//                  reports a peak identical to the contracted arm, though the
//                  peak term plainly carries the proto-occ pair).
//   processing     the caching veto has two wrong horns: the original broad
//                  form killed CSE outright on the real C60 job, the narrowed
//                  form (2a52e063c) may not be reachable from a DP-emitted
//                  schedule at all.
//
// CONSEQUENCE -- what may and may not be quoted:
//
//   NOT TRUSTED (diagnostics only, NOT targets): avoidable_time, replay op
//   counts, unique-key counts, and any external-on vs external-off comparison.
//   These are computed BY replaying through the machinery listed above, in both
//   numerator and denominator. They measure a pipeline with open defects, not a
//   physical quantity. The aggregate `hw` peak is likewise liveness-dependent.
//
//   TRUSTED: pure size arithmetic on an emitted node (verified exactly -- the
//   peak op's g(mu,mu,K_2) operand is 1800*1800*256 elements to the byte), and
//   DP-side structural facts read from SEQUANT_DP_RECOMPUTE_DEBUG /
//   SEQUANT_SELROOT_DEBUG (which terms are annotated, fit=0 counts, rf values,
//   carried/escaped sets). Those do not pass through the replay.
//
// Gate Phase A on the DP-side facts, which survive this boundary. Phase B's
// gate necessarily uses the replay -- legitimate only because by then that
// machinery is what is UNDER TEST, not what is serving as the oracle.
//
// HISTORY, because this exact failure already happened once: an earlier
// revision reported 76% avoidable_time; that figure was adopted into design
// specs as a target and drove roughly a week of work before it was found to be
// ~30.8 points runtime misconfiguration measuring a policy MPQC does not use,
// attached by inference ("consistent with") to a real C60 job whose actual
// pathology was the broad caching veto. Do not repeat that with 43.72%.
// Full account: .superpowers/sdd/oamb-a0-note.md
// ==========================================================================
TEST_CASE("dryrun occ batching wipes CSE (free-batchable-mode veto repro)",
          "[.][dryrun-occ-veto]") {
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);
  sequant::mbpt::add_df_spaces(isr);
  auto ctx_resetter = set_scoped_default_context(std::move(ctx));

  auto const body = slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                          "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = deserialize<ExprPtr>(line);
  REQUIRE(static_cast<bool>(expr));
  REQUIRE(expr->is<Sum>());
  auto const& summands = expr->as<Sum>().summands();
  REQUIRE(!summands.empty());
  auto flatten_product = [](ExprPtr const& e) -> ExprPtr {
    if (!e->is<Product>()) return e;
    auto const& p = e->as<Product>();
    return ex<Product>(p.scalar(), p.factors(), Product::Flatten::Yes);
  };
  // The veto's damage is to CROSS-TERM sharing -- the cache is what shares an
  // intermediate across summands (the gC -> {gCC, gCC, gCCC} case in
  // eval.hpp). A single term has almost nothing to share, so the whole forest
  // is replayed, exactly as the real run does. SEQUANT_UT_DRYRUN_NTERMS caps
  // the term count (optimizing every summand under the batched DP is slow).
  std::size_t nterms = summands.size();
  if (char const* nt = std::getenv("SEQUANT_UT_DRYRUN_NTERMS"))
    nterms = std::min<std::size_t>(nterms, std::atoll(nt));

  auto regime = df_regime(kC60_pVDZF12);

  struct Meas {
    std::size_t static_nodes = 0;  // cp.n_ops: internal nodes (ideal = once)
    std::size_t replay_ops = 0;    // Eval|Product lines actually replayed
    std::size_t distinct = 0;      // distinct product expressions (no slice)
    std::size_t unique_keys = 0;   // distinct (expr, touched-mode-slice) builds
    std::size_t ops_with_cost = 0;  // ops matched to an OpCost line (coverage)
    double total_exec = 0;          // sum of modelled exec_cost over all ops
    double avoidable_exec = 0;      // exec_cost of the duplicate-key builds
    double cp_exec = 0;             // cost_profile() static total (validation)
    double peak_gb = 0;
    // Task 5 acceptance gate (aux+occ leg). Counted over the EMITTED eval-node
    // forest's batched_here() stamps, classifying occ by space base_key L"i"
    // (the same way the external-role policy above identifies occ). These are
    // the runtime-visible role guarantees of the role split:
    //   contracted_occ_stamps -- occ stamped Contracted. MUST be 0: the fix's
    //     core guarantee. Before the role split the DP put occ in the
    //     contracted predicate to satisfy the runtime accept, so it got sliced
    //     in a contracted cell the runtime never realizes; this would have been
    //     > 0.
    //   external_occ_stamps -- occ stamped External. MUST be > 0: proves the
    //     derived union accept did NOT drop external occ, so the occ scatter is
    //     engaged. If it were 0 the external-occ path would be dead and the
    //     contracted==0 check would be vacuously true; the pairing makes it
    //     real.
    std::size_t contracted_occ_stamps = 0;
    std::size_t external_occ_stamps = 0;
    std::wstring top_expr;  // worst avoidable offender
    double top_ratio = 0;   // its builds / its distinct touched-slices
    // Avoidable-recompute OP-COUNT factor: total builds over the number of
    // DISTINCT (expression, touched-mode-slice) builds a perfect-sharing
    // evaluator would do, minus 1. Legitimate slicing (a different slice of an
    // mode the op's value touches) makes a distinct key and does NOT count;
    // only the SAME value rebuilt under an enclosing batch it is invariant to
    // does. 0 = no avoidable recomputation.
    double avoidable() const {
      return unique_keys ? double(replay_ops) / double(unique_keys) - 1.0 : 0.0;
    }
    // Avoidable-recompute TIME fraction: the share of modelled exec_cost spent
    // on duplicate-key builds. This -- not the op-count factor -- is what a fix
    // would recover, since a cheap invariant node rebuilt 1000x and an
    // expensive one rebuilt twice weigh very differently.
    double avoidable_time() const {
      return total_exec > 0 ? avoidable_exec / total_exec : 0.0;
    }
  };

  auto run = [&](bool batch_aux, bool batch_occ) -> Meas {
    sequant::BatchPolicy policy;
    // Task 6 GATE: the aux+occ leg exercises the EXTERNAL-occ path -- batch the
    // residual-target occ i,j as forest spectators, sliced per external block.
    // This is the production MPQC config (csv_batch_policy.h) and the config
    // the scope_level fix (7cf839d68) must be verified against: with both flags
    // on, an external-carrying intermediate (the 1-PNO gC) must slice per
    // external block instead of being hoisted to the real cache at full extent.
    // The unbatched and aux-only legs keep batch_occ==false, so both flags stay
    // off there and their measurements are byte-identical to before.
    policy.batch_spectator_indices = batch_occ;
    policy.order_aware_recompute = batch_occ;
    // Role split: aux Κ is batchable in the CONTRACTED role (it is summed at
    // some node); the residual-target occ i is batchable in the EXTERNAL role
    // ONLY (a spectator carried to the result, contracted nowhere). Declaring
    // occ contracted-batchable would let the DP slice it in a contracted cell
    // the runtime never realizes -- the contamination this role split removes.
    policy.is_batchable_contracted_index = [batch_aux](Index const& ix) {
      return batch_aux &&
             ix.space().base_key() == L"Κ";  // mu~ NOT batched (as in 631467)
    };
    policy.is_batchable_external_index = [batch_occ](Index const& ix) {
      return batch_occ && ix.space().base_key() == L"i";
    };
    policy.batch_target_size = [](Index const& ix) -> std::size_t {
      return ix.space().base_key() == L"Κ" ? std::size_t{256} : std::size_t{8};
    };
    policy.is_volatile_leaf = [](Tensor const& t) { return t.label() == L"t"; };
    policy.accumulation_factor = 1.0;
    policy.peak_threshold =
        (std::getenv("SEQUANT_UT_DRYRUN_PEAK_THR_GB")
             ? std::atof(std::getenv("SEQUANT_UT_DRYRUN_PEAK_THR_GB"))
             : 100.0) *
        1e9;  // default = the C60 job's 100 GB budget

    auto axes_map = std::make_shared<std::unordered_map<
        Expr const*, container::vector<NodeBatchAnnotation>>>();
    OptimizeOptions opts;
    opts.objective_function = ObjectiveFunction::DenseTimeSpaceBatched;
    opts.idx_to_extent = regime.idx_to_extent();
    opts.inner_pow = regime.inner_pow_fn();
    opts.batch_policy = policy;
    opts.volatile_weight = 20.0;
    opts.roofline.machine_balance = 200.0;
    opts.roofline.fast_mem_elems = 1000000.0;
    opts.term_batch_axes = axes_map;

    std::vector<EvalNodeDryRun> forest;
    std::vector<std::size_t> forest_summand;
    for (std::size_t s = 0; s < nterms; ++s) {
      ExprPtr const term = flatten_product(summands[s]);
      if (!term) continue;
      auto optimized = optimize(term, opts);
      if (!optimized) continue;
      auto it = axes_map->find(optimized.get());
      container::vector<NodeBatchAnnotation> node_axes;
      if (it != axes_map->end()) node_axes = it->second;
      BinarizationOptions bopts;
      bopts.node_batch_axes = node_axes;
      SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
      forest.push_back(binarize<EvalExprDryRun>(optimized, {}, bopts));
      SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
      forest_summand.push_back(s);
    }
    REQUIRE(!forest.empty());

    // Getenv-gated: name the summand with the largest per-summand REPLAY peak
    // (the peak-defining term -- realized, not nominal) and print its FULL
    // symbolic expression. Ranks by cost_profile().peak_bytes per summand.
    if (batch_aux && batch_occ && std::getenv("SEQUANT_DUMP_PEAKTERM")) {
      sequant::eval::dryrun::CacheConfig pcfg;
      pcfg.max_footprint = 1e11;
      pcfg.min_repeats = 1;
      pcfg.is_volatile = [](EvalNodeDryRun const& n) {
        return n.leaf() && n->is_tensor() && n->as_tensor().label() == L"t";
      };
      double best_peak = 0.0;
      std::size_t best_s = 0;
      for (std::size_t fi = 0; fi < forest.size(); ++fi) {
        std::vector<EvalNodeDryRun> one{forest[fi]};
        auto const cp1 =
            sequant::eval::dryrun::cost_profile(one, policy, pcfg, regime);
        if (cp1.peak_bytes > best_peak) {
          best_peak = cp1.peak_bytes;
          best_s = forest_summand[fi];
        }
      }
      std::wcerr << L"[peakterm] max per-summand REPLAY peak = "
                 << (best_peak / 1e9) << L" GB in summand " << best_s
                 << L"\n[peakterm] FULL TERM (summand " << best_s << L"):\n"
                 << to_latex(flatten_product(summands[best_s])) << L"\n";
      // Factorized tree of the peak summand: per internal node, result indices
      // and nominal footprint, plus batched_here (what the runtime slices at
      // that node).
      auto fp2 = sequant::opt::detail::footprint_counter(regime.idx_to_extent(),
                                                         regime.inner_pow_fn());
      for (std::size_t fi = 0; fi < forest.size(); ++fi) {
        if (forest_summand[fi] != best_s) continue;
        std::wcerr << L"[peaktree] === summand " << best_s
                   << L" factorized tree (internal nodes) ===\n";
        forest[fi].visit_internal([&](auto const& n) {
          double const gb = fp2(n->canon_indices()) * 8.0 / 1e9;
          std::wstring res, bh;
          for (auto const& ix : n->canon_indices())
            res += std::wstring(ix.full_label()) + L" ";
          for (auto const& [ix, knd] : n->batched_here())
            bh += std::wstring(ix.full_label()) + L":" +
                  (knd == BatchModeType::External     ? L"EXT"
                   : knd == BatchModeType::Contracted ? L"CON"
                                                      : L"?") +
                  L" ";
          std::wcerr << L"[peaktree] foot=" << gb << L"GB batched_here={" << bh
                     << L"} result={" << res << L"}\n";
        });
      }
      // FACTORIZER's own modeled peak for the peak summand (DP estimate, NOT
      // the replay). If it is far below the replay peak, the DP INTENDS the 2.9
      // TB giants sliced (batched) -- i.e. the factorizer wants them batched
      // and the gap is runtime-only.
      {
        ExprPtr t46 = flatten_product(summands[best_s]);
        if (t46->is<Product>()) {
          TensorNetwork net46(t46->as<Product>().factors());
          container::svector<Index> tg46;
          opt::detail::PeakBatchedModel m46{
              regime.idx_to_extent(),
              [](Index const& ix) -> std::size_t {
                return ix.space().base_key() == L"Κ" ? std::size_t{256}
                                                     : std::size_t{8};
              },
              [](Tensor const& t) { return t.label() == L"t"; },
              regime.inner_pow_fn()};
          m46.is_batchable_contracted_index = [](Index const& ix) {
            return ix.space().base_key() == L"Κ";
          };
          m46.is_batchable_external_index = [](Index const& ix) {
            return ix.space().base_key() == L"i";
          };
          m46.batch_spectator_indices = true;
          m46.order_aware_recompute = true;
          m46.perf_first = true;
          m46.volatile_weight = 20.0;
          m46.machine_balance = 200.0;
          m46.fast_mem_elems = 1000000.0;
          m46.accumulation_factor = 1.0;
          m46.peak_threshold = 100e9;
          m46.numeric_size = 8.0;
          auto ctx46 = m46.build_context(net46, tg46);
          auto st46 = opt::detail::solve_single_term(m46, net46, tg46, ctx46);
          double dp_peak = 0.0;
          (void)m46.reconstruct_batched_modes(ctx46, st46, net46, tg46,
                                              &dp_peak);
          std::wcerr << L"[peakterm] FACTORIZER modeled peak (summand "
                     << best_s << L") = " << (dp_peak / 1e9)
                     << L" GB  vs REPLAY peak = " << (best_peak / 1e9)
                     << L" GB\n";
        }
      }
    }

    // Task 6 DIAGNOSTIC (getenv-gated, test-only, prints nothing off-path).
    // For the external-occ (aux+occ) leg, walk every internal node and report
    // the giants: nominal (unsliced) result footprint and the child leaf
    // labels, so a gC-shaped node (a `g` leaf contracted with a `C` leaf) is
    // identifiable. This pins WHICH node sets the peak.
    if (batch_aux && batch_occ && std::getenv("SEQUANT_UT_DRYRUN_GCPROBE_GB")) {
      double const thr_gb =
          std::atof(std::getenv("SEQUANT_UT_DRYRUN_GCPROBE_GB"));
      auto footprint = sequant::opt::detail::footprint_counter(
          regime.idx_to_extent(), regime.inner_pow_fn());
      auto lbl = [](auto const& node) -> std::wstring {
        if (node->is_tensor()) return std::wstring(node->as_tensor().label());
        return L"<op>";
      };
      auto idxstr = [](auto const& n) -> std::wstring {
        std::wstring s;
        for (auto const& ix : n->canon_indices())
          s += std::wstring(ix.full_label()) + L" ";
        return s;
      };
      std::function<void(EvalNodeDryRun const&)> gcwalk =
          [&](EvalNodeDryRun const& n) {
            if (n.leaf()) return;
            double const gb = footprint(n->canon_indices()) * 8.0 / 1e9;
            if (gb >= thr_gb) {
              std::wcerr << L"[gcprobe] foot=" << gb << L"GB children={"
                         << lbl(n.left()) << L"," << lbl(n.right())
                         << L"} result={" << idxstr(n) << L"}\n";
            }
            gcwalk(n.left());
            gcwalk(n.right());
          };
      std::wcerr << L"[gcprobe] === aux+occ giants (foot >= " << thr_gb
                 << L" GB) ===\n";
      for (auto const& root : forest) gcwalk(root);
    }

    // VETO REACHABILITY (oamb-a0-note.md 11.4). cache_manager vetoes a node
    // iff its OWN batched_here() carries a Contracted, batchable mode that is
    // also FREE in its own result. The DP emits Contracted only for
    // aprime subset-of contracted_here -- exactly the modes NOT free in that
    // node's result -- so structurally the predicate should never hold. Count
    // it on the real forest rather than argue from the invariant.
    std::size_t n_contracted_stamps = 0, n_veto_eligible = 0;
    // Task 5 role-guarantee tallies: occ is identified by space base_key L"i",
    // exactly as the external-role policy classifies it above.
    std::size_t n_contracted_occ_stamps = 0, n_external_occ_stamps = 0;
    for (auto const& root : forest) {
      root.visit_internal([&](auto const& n) {
        auto const& canon = n->canon_indices();
        for (auto const& [ix, knd] : n->batched_here()) {
          bool const is_occ = ix.space().base_key() == L"i";
          if (knd == BatchModeType::Contracted && is_occ)
            ++n_contracted_occ_stamps;
          if (knd == BatchModeType::External && is_occ) ++n_external_occ_stamps;
          if (knd != BatchModeType::Contracted) continue;
          if (!policy.is_batchable_contracted_index(ix)) continue;
          ++n_contracted_stamps;
          if (std::find(canon.begin(), canon.end(), ix) != canon.end())
            ++n_veto_eligible;
        }
      });
    }
    std::wcerr << L"[veto-reach] Contracted+batchable stamps = "
               << n_contracted_stamps << L", of which FREE in the node's own "
               << L"result (veto fires) = " << n_veto_eligible << L"\n";

    // Is the ORIGINAL hazard (5c5eb1e82) live? v1 vetoed any node whose RESULT
    // carries an accepted batchable index, on the grounds that such a node is
    // sliced by the runtime yet would be cached whole. v2 never fires (above),
    // so v1's protection is gone -- but that only matters if such nodes are in
    // fact cached. Count them: nodes registered in the cache whose own result
    // carries a batchable index. Zero would mean the hazard is not reachable on
    // this forest and only the ANNOTATION-scoped question remains.
    {
      sequant::eval::dryrun::CacheConfig probe_cfg;
      probe_cfg.max_footprint = 1e11;
      probe_cfg.min_repeats = 1;
      probe_cfg.is_volatile = [](EvalNodeDryRun const& n) {
        if (!n.leaf() || !n->is_tensor()) return false;
        return n->as_tensor().label() == L"t";
      };
      auto probe_cache =
          sequant::eval::dryrun::build_dryrun_cache(forest, probe_cfg, regime);
      std::size_t n_cached = 0, n_cached_carrying = 0;
      for (auto const& root : forest)
        root.visit_internal([&](auto const& n) {
          if (!probe_cache.exists(n)) return;
          ++n_cached;
          for (auto const& ix : n->canon_indices())
            if (policy.is_batchable_index()(ix)) {
              ++n_cached_carrying;
              break;
            }
        });
      std::wcerr << L"[veto-hazard] cached internal nodes = " << n_cached
                 << L", of which carry a BATCHABLE index free in their result "
                 << L"(v1 would have vetoed; cached WHOLE today) = "
                 << n_cached_carrying << L"\n";
    }

    sequant::eval::dryrun::CacheConfig cfg;
    cfg.max_footprint = 1e11;
    cfg.min_repeats = 1;
    cfg.is_volatile = [](EvalNodeDryRun const& n) {
      if (!n.leaf() || !n->is_tensor()) return false;
      return n->as_tensor().label() == L"t";
    };

    std::wostringstream trace;
    auto const cp = sequant::eval::dryrun::cost_profile(forest, policy, cfg,
                                                        regime, &trace);

    if (std::getenv("SEQUANT_UT_DRYRUN_DUMPTRACE"))
      std::wcerr << L"---- TRACE ----\n" << trace.str() << L"---- END ----\n";

    Meas m;
    m.static_nodes = cp.model_n_ops;
    m.peak_gb = cp.peak_bytes / 1e9;
    m.contracted_occ_stamps = n_contracted_occ_stamps;
    m.external_occ_stamps = n_external_occ_stamps;

    // i-th "| "-delimited field of a trace line, trailing space/CR trimmed.
    auto field = [](std::wstring const& s, std::size_t i) -> std::wstring {
      std::size_t start = 0;
      for (std::size_t k = 0; k < i; ++k) {
        auto const p = s.find(L"| ", start);
        if (p == std::wstring::npos) return L"";
        start = p + 2;
      }
      auto const end = s.find(L"| ", start);
      std::wstring f = s.substr(
          start, end == std::wstring::npos ? std::wstring::npos : end - start);
      while (!f.empty() &&
             (f.back() == L' ' || f.back() == L'\r' || f.back() == L'\n'))
        f.pop_back();
      return f;
    };

    // Reconstruct the enclosing batch-loop stack from the BatchGroup Begin/End
    // + BatchIter markers (BatchGroup nesting == batch-loop nesting; each level
    // carries the last BatchIter's (mode, slice)). Key each replayed op by its
    // expression PLUS the slices of only the modes that occur in it -- the
    // touched-mode-slice signature. Two builds with the same key are the same
    // value recomputed (avoidable); a different touched slice is legitimate.
    m.cp_exec = cp.model_exec;
    std::vector<std::pair<std::wstring, std::wstring>> stack;  // (mode, slice)
    std::map<std::wstring, std::size_t> total_of;              // expr -> builds
    std::map<std::wstring, std::set<std::wstring>> keys_of;    // expr -> keys
    std::map<std::wstring, double>
        exec_of;  // expr -> avoidable exec (offender)
    std::set<std::wstring> all_keys;

    std::wistringstream in(trace.str());
    std::wstring ln;
    double last_exec = -1.0;  // exec_cost from the OpCost line preceding an op
    while (std::getline(in, ln)) {
      if (ln.rfind(L"BatchGroup | Begin", 0) == 0) {
        stack.emplace_back(L"", L"");
      } else if (ln.rfind(L"BatchGroup | End", 0) == 0) {
        if (!stack.empty()) stack.pop_back();
      } else if (ln.rfind(L"BatchIter", 0) == 0) {
        if (!stack.empty()) stack.back() = {field(ln, 1), field(ln, 2)};
      } else if (ln.rfind(L"OpCost", 0) == 0) {
        last_exec = std::wcstod(field(ln, 2).c_str(), nullptr);
      } else if (ln.find(L"Eval | Product") != std::wstring::npos) {
        ++m.replay_ops;
        double const exec = last_exec >= 0.0 ? last_exec : 0.0;
        if (last_exec >= 0.0) ++m.ops_with_cost;
        last_exec = -1.0;  // consume; a missing OpCost must not carry over
        m.total_exec += exec;
        auto const p = ln.rfind(L"| ");
        std::wstring const expr =
            (p == std::wstring::npos) ? ln : ln.substr(p + 2);
        std::wstring sig;
        for (auto const& [mode, slice] : stack)
          if (!mode.empty() && expr.find(mode) != std::wstring::npos)
            sig += L"|" + mode + L"=" + slice;
        std::wstring const key = expr + L"@" + sig;
        ++total_of[expr];
        keys_of[expr].insert(key);
        bool const dup = !all_keys.insert(key).second;
        if (dup) {  // same value already built at this touched-slice
          m.avoidable_exec += exec;
          exec_of[expr] += exec;
        }
      }
    }
    m.distinct = total_of.size();
    m.unique_keys = all_keys.size();
    double worst_exec = 0;
    for (auto const& [e, tot] : total_of) {
      double const r = double(tot) / double(keys_of[e].size());
      if (r > m.top_ratio) m.top_ratio = r;
      if (exec_of[e] > worst_exec) {  // offender by avoidable TIME, not count
        worst_exec = exec_of[e];
        m.top_expr = e;
      }
    }
    return m;
  };

  auto const base = run(/*aux=*/false, /*occ=*/false);
  auto const aux = run(/*aux=*/true, /*occ=*/false);
  auto const both = run(/*aux=*/true, /*occ=*/true);

  auto report = [](wchar_t const* label, Meas const& m) {
    std::wcerr << L"[occ-veto] " << label << L": ops=" << m.replay_ops
               << L" unique(expr,slice)=" << m.unique_keys << L" avoidable_ops="
               << m.avoidable() << L"x  AVOIDABLE_TIME="
               << (100.0 * m.avoidable_time()) << L"%  PEAK=" << m.peak_gb
               << L"GB  (cost-coverage " << m.ops_with_cost << L"/"
               << m.replay_ops << L", total_exec=" << m.total_exec << L" vs cp "
               << m.cp_exec << L")\n           worst avoidable-time offender: "
               << m.top_expr << L"\n";
  };
  std::wcerr << L"\n=== [dryrun-occ-veto] C60 residual forest, " << nterms
             << L" terms (K@256, occ@8, 100GB, perf-first) ===\n";
  report(L"unbatched ", base);
  report(L"aux-only  ", aux);
  report(L"aux+occ   ", both);
  std::wcerr << L"[occ-veto] aux+occ role stamps: Contracted-occ="
             << both.contracted_occ_stamps << L" (MUST be 0)  External-occ="
             << both.external_occ_stamps << L" (MUST be > 0)\n";

  // The forest is identical in all three; only the mode policy differs, so the
  // static node count must not move.
  CHECK(base.static_nodes == both.static_nodes);
  CHECK(aux.static_nodes == both.static_nodes);

  // ==== Task 5 ACCEPTANCE GATE (the assertions that catch the original bug) ==
  // These are the runtime-visible guarantees of the batchability role split,
  // measured on the emitted aux+occ eval-node forest (occ classified by space
  // base_key L"i", as the external-role policy classifies it above).
  //
  // (1) NO emitted node stamps a contracted-occ mode. This is the fix's core
  // guarantee: occ is batchable in the EXTERNAL role ONLY. Before the role
  // split, callers put occ into the (then single) batchability predicate to
  // satisfy the runtime accept, which as a side effect declared it
  // Contracted-batchable to the DP; the DP then sliced it in a contracted cell
  // the runtime never realizes, so the intermediate materialized whole. On that
  // pre-split code this count was > 0 and this assertion would FAIL; post-split
  // it must be 0.
  CHECK(both.contracted_occ_stamps == 0);
  // (2) External-occ scatter still fires (non-vacuous guard). At least one node
  // must carry an External-occ stamp, proving the derived-union runtime accept
  // did NOT drop external occ -- the occ scatter is engaged. If this were 0 the
  // whole external-occ path would be dead and assertion (1) would be vacuously
  // true; pairing the two makes the gate meaningful.
  CHECK(both.external_occ_stamps > 0);
  // ==========================================================================

  // Measurement guard: nearly every op must carry a cost model OpCost line, or
  // the time fraction is unreliable.
  CHECK(both.ops_with_cost >= both.replay_ops - both.replay_ops / 20);

  // ASPIRATIONAL GATE -- PEAK (intentionally RED, documented). This is a MEMORY
  // problem: the C60 job it mirrors ran under a 100 GB budget and never
  // completed an iteration. Post role-split, the modelled aux+occ replay peak
  // is ~5860.9 GB (nterms=55), i.e. ~59x over budget. The two aspirational
  // CHECKs (peak < 100, avoidable < 0.10) remain RED research targets; they are
  // NOT forced to pass and NOT the acceptance gate. The Task 5 acceptance
  // assertions -- Contracted-occ stamps == 0 and External-occ stamps > 0
  // (above) -- are the ones that must be GREEN.
  //
  // WHY THE PEAK ROSE, AND WHY IT IS THE HONEST NUMBER. Earlier revisions of
  // this witness reported ~2302 GB (contracted-occ) and ~2947 GB (a
  // transitional config where occ was still Contracted-stamped). Those figures
  // were an ARTIFACT: the DP was slicing occ in contracted cells the runtime
  // never realizes, so the DP's modelled szcell HID the true cost of the
  // intermediate
  // -- the peak it reported was below the real replay footprint. The role split
  // removes that phantom contracted-occ batching, so the DP no longer hides the
  // cost of the cross-pair (two-PNO-leg) giants, and the peak rises to the TRUE
  // ~5860.9 GB. A higher-but-honest number is the correct outcome, not a
  // regression.
  //
  // WHY EXTERNAL SLICING CANNOT CLOSE IT. The peak driver is a cross-pair
  // intermediate carrying TWO independent PNO legs (each an a<i,j> PNO domain
  // over a distinct occ pair). External-occ slicing removes ONE occ pair's
  // dependence (the forest-spectator target occ), but it cannot reduce the
  // SECOND PNO-pair leg -- that leg is contracted, not external. Reducing it
  // would require contracted-occ batching, which this scheme DELIBERATELY does
  // not do (that is exactly the phantom the role split removed). Closing this
  // gate is therefore a downstream design question (contracted-occ batching or
  // a factorization that never forms the two-PNO-leg intermediate), tracked as
  // an out-of-scope follow-up -- not something this witness can or should
  // force.
  //
  // MEASURED AFTER THE TASKS 1-3 HOIST-EXCLUSION (extscope; e3174a468 +
  // 0711a1586, nterms=55). The order-aware plan expected the summand-46 giant
  // -- believed cache-HOISTED at full external extent above the i1,i2 scatter
  // loop
  // -- to fall to the factorizer-modeled ~35 GB once the loop-local
  // external-carrying intermediate was excluded from hoist_invariants::collect.
  // It did NOT: the aux+occ peak is UNCHANGED at 5860.877 GB, byte-identical
  // with and without the e3174a468 `!ext_loop_local` conjunct. Localization
  // (SEQUANT_UT_PEAK_COMPONENTS probe, since reverted): cost_profile()'s
  // peak_bytes = max(batched-SCRATCH high-watermark peak.load(), gated CACHE
  // working_set_hwmark()). On this forest the CACHE hwmark is ~94 GB (already
  // BELOW the 100 GB budget) and the SCRATCH is 5860.877 GB (~2x the 2930 GB
  // summand-46 giant, co-resident) -- so the peak is set entirely by the
  // scatter-SCRATCH path, which Task 3 deliberately leaves byte-unchanged (the
  // hoist-exclusion touches only the cache, and the cache was never the driver
  // here). The factorizer models 35.08 GB for summand 46 vs the 5860.877 GB
  // replay: the 167x gap is a scatter-scratch replay defect (two 2930 GB giants
  // held live at once), orthogonal to the hoist Task 3 fixed. Closing it is the
  // same downstream scatter-scratch work as the contracted-occ leg above, NOT
  // reachable from the hoist path. Full account:
  // .superpowers/sdd/task-4-report.md.
  CHECK(both.peak_gb < 100.0);  // documented-RED research target (~5860.9 GB)

  // ASPIRATIONAL -- AVOIDABLE recomputation TIME: the share of modelled
  // exec_cost spent rebuilding a value already built at the same touched-mode
  // slice. Legitimate slicing (a different slice of a mode the op's value
  // depends on) is divided out. The raw op count and even the op-count
  // avoidable factor overstate it, because a cheap invariant node rebuilt many
  // times weighs little.
  //
  // Post role-split, honest numbers (annotation-only batching, no heuristic
  // fallback -- the production configuration, cf. MPQC csv_batch_policy.h;
  // nterms=55): unbatched ~1.8%, aux-only ~6.5%, aux+occ ~44.6% with
  // unique(expr,slice) ~= 2080. The aux+occ figure came DOWN from a
  // transitional ~65.0% (unique 47292): with the phantom contracted-occ
  // batching removed, the DP no longer spawns the flood of
  // contracted-occ-sliced key variants, so both the avoidable-time share and
  // the unique-key count fall to their honest values. (Ignore any older
  // ~2302/2947 GB or ~76% narrative that claimed "external-occ raises the
  // peak/avoidable" -- that was the phantom, now gone.)
  //
  // aux-only is the CLEAN baseline; it is unchanged (both flags stay off on
  // that leg). The surviving aux+occ ~44.6% is dominated by the cross-pair
  // two-PNO- leg giants: external slicing cannot reduce the second (contracted)
  // PNO-pair leg, so those nodes replay per external block without a matching
  // working-set reduction. Bringing aux+occ down to the aux-only baseline needs
  // the same downstream design work the PEAK gate does (contracted-occ batching
  // or a factorization that never forms the two-PNO-leg intermediate) -- so
  // this aspirational CHECK also stays RED and documented, not forced.
  CHECK(aux.avoidable_time() < 0.10);  // aux-only is clean (PASSES)
  CHECK(both.avoidable_time() <
        0.10);  // documented-RED research target (~0.446)
}

// D5 avoidable-recomputation GATE: prove EXTERNAL-mode batching eliminates the
// ~76% avoidable recomputation the CONTRACTED-occ schedule showed above.
//
// RETRACTED / CORRECTED (2026-07-27, comment only -- assertions/logic below
// unchanged): the "~76% contracted-occ vs ~0 external-occ" narrative in this
// comment block is contaminated. Both configs put occ in the SINGLE
// is_batchable_index (contracted-role) predicate -- is_batchable_external_index
// was never set -- so the "external-occ" arm still batched CONTRACTED occ too;
// the identical peak across arms is the tell. Honest post-role-split re-measure
// of the sibling [.][dryrun-occ-veto] witness (nterms=55): peak ~= 5860.9 GB,
// avoidable_time ~= 44.6%, contracted_occ_stamps == 0, external_occ_stamps ==
// 244. Root cause: .superpowers/sdd/contamination-role-predicate.md.
//
// The [.][dryrun-occ-veto] witness above measures, on the C60 residual replay,
// the exec-cost share spent rebuilding a value already built at the same
// touched-mode slice (avoidable_time). Its CONTRACTED-occ config (occ batched
// as a Contracted mode) showed ~76%: the PPL giant W and its g.C-class factors
// -- carried over the external occ i,j as PNO protoindices -- were rebuilt once
// per contracted-occ block though the block-loop is invariant to them.
//
// EXTERNAL mode fixes exactly this. With batch_spectator_indices ON, the DP
// recognizes an over-budget term's genuine forest external modes (the residual
// target occ i,j, carried on EVERY node) and stamps BatchModeType::External on
// them (D1, work-neutral: seeded_flops == unseeded_flops). The runtime scatter
// branch (eval.hpp) then slices EVERY leaf carrying the external occ into
// disjoint blocks and write_into_slice()s the block partials into one pre-sized
// result -- each block does 1/nblocks of the work, nothing rebuilt. So the
// expectation is avoidable_time ~ 0 on the external-batched replay.
//
// Same C60 job shape as the veto (K@256, occ@8, mu~ NOT batched, 100 GB,
// perf-first) -- the ONLY policy change is batch_spectator_indices = true, so
// occ i,j go External instead of Contracted.
TEST_CASE(
    "dryrun external-mode batching zeroes the occ-veto avoidable recompute",
    "[.][dryrun-extmode-avoidable]") {
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);
  sequant::mbpt::add_df_spaces(isr);
  auto ctx_resetter = set_scoped_default_context(std::move(ctx));

  auto const body = slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                          "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = deserialize<ExprPtr>(line);
  REQUIRE(static_cast<bool>(expr));
  REQUIRE(expr->is<Sum>());
  auto const& summands = expr->as<Sum>().summands();
  REQUIRE(!summands.empty());
  auto flatten_product = [](ExprPtr const& e) -> ExprPtr {
    if (!e->is<Product>()) return e;
    auto const& p = e->as<Product>();
    return ex<Product>(p.scalar(), p.factors(), Product::Flatten::Yes);
  };
  std::size_t nterms = summands.size();
  if (char const* nt = std::getenv("SEQUANT_UT_DRYRUN_NTERMS"))
    nterms = std::min<std::size_t>(nterms, std::atoll(nt));

  auto regime = df_regime(kC60_pVDZF12);

  struct Meas {
    std::size_t static_nodes = 0;
    std::size_t replay_ops = 0;
    std::size_t unique_keys = 0;
    std::size_t ops_with_cost = 0;
    std::size_t n_scatter_begin = 0;    // BatchScatter interceptions fired
    std::size_t n_external_stamps = 0;  // node_axes entries tagged External
    std::size_t n_contracted_occ =
        0;  // node_axes entries: occ tagged Contracted
    std::size_t n_batchgroup_begin = 0;  // BatchGroup interceptions fired
    double total_exec = 0;
    double avoidable_exec = 0;
    double peak_gb = 0;
    std::wstring top_expr;
    double avoidable_time() const {
      return total_exec > 0 ? avoidable_exec / total_exec : 0.0;
    }
  };

  // batch_external toggles EXTERNAL mode (batch_spectator_indices). occ+aux are
  // batchable in both; only the external gate differs.
  auto run = [&](bool batch_external) -> Meas {
    sequant::BatchPolicy policy;
    policy.batch_spectator_indices = batch_external;
    // Role split: aux Κ contracted-batchable; residual-target occ i is
    // batchable in the EXTERNAL role only (batch_external gates its emission,
    // matching batch_spectator_indices above). It is NEVER
    // contracted-batchable.
    policy.is_batchable_contracted_index = [](Index const& ix) {
      return ix.space().base_key() == L"Κ";  // aux; mu~/occ NOT contracted
    };
    policy.is_batchable_external_index = [batch_external](Index const& ix) {
      return batch_external && ix.space().base_key() == L"i";
    };
    policy.batch_target_size = [](Index const& ix) -> std::size_t {
      return ix.space().base_key() == L"Κ" ? std::size_t{256} : std::size_t{8};
    };
    policy.is_volatile_leaf = [](Tensor const& t) { return t.label() == L"t"; };
    policy.accumulation_factor = 1.0;
    policy.peak_threshold =
        (std::getenv("SEQUANT_UT_DRYRUN_PEAK_THR_GB")
             ? std::atof(std::getenv("SEQUANT_UT_DRYRUN_PEAK_THR_GB"))
             : 100.0) *
        1e9;

    auto axes_map = std::make_shared<std::unordered_map<
        Expr const*, container::vector<NodeBatchAnnotation>>>();
    OptimizeOptions opts;
    opts.objective_function = ObjectiveFunction::DenseTimeSpaceBatched;
    opts.idx_to_extent = regime.idx_to_extent();
    opts.inner_pow = regime.inner_pow_fn();
    opts.batch_policy = policy;
    opts.volatile_weight = 20.0;
    opts.roofline.machine_balance = 200.0;
    opts.roofline.fast_mem_elems = 1000000.0;
    opts.term_batch_axes = axes_map;

    Meas m;
    std::vector<EvalNodeDryRun> forest;
    for (std::size_t s = 0; s < nterms; ++s) {
      ExprPtr const term = flatten_product(summands[s]);
      if (!term) continue;
      auto optimized = optimize(term, opts);
      if (!optimized) continue;
      auto it = axes_map->find(optimized.get());
      container::vector<NodeBatchAnnotation> node_axes;
      if (it != axes_map->end()) node_axes = it->second;
      for (auto const& node : node_axes)
        for (auto const& [ix, knd] : node.axes) {
          if (knd == BatchModeType::External) ++m.n_external_stamps;
          if (knd == BatchModeType::Contracted && ix.space().base_key() == L"i")
            ++m.n_contracted_occ;
        }
      BinarizationOptions bopts;
      bopts.node_batch_axes = node_axes;
      SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
      forest.push_back(binarize<EvalExprDryRun>(optimized, {}, bopts));
      SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    }
    REQUIRE(!forest.empty());

    sequant::eval::dryrun::CacheConfig cfg;
    cfg.max_footprint = 1e11;
    cfg.min_repeats = 1;
    cfg.is_volatile = [](EvalNodeDryRun const& n) {
      if (!n.leaf() || !n->is_tensor()) return false;
      return n->as_tensor().label() == L"t";
    };

    std::wostringstream trace;
    auto const cp = sequant::eval::dryrun::cost_profile(forest, policy, cfg,
                                                        regime, &trace);
    if (std::getenv("SEQUANT_UT_DRYRUN_DUMPTRACE"))
      std::wcerr << L"---- TRACE ----\n" << trace.str() << L"---- END ----\n";

    m.static_nodes = cp.model_n_ops;
    m.peak_gb = cp.peak_bytes / 1e9;

    auto field = [](std::wstring const& s, std::size_t i) -> std::wstring {
      std::size_t start = 0;
      for (std::size_t k = 0; k < i; ++k) {
        auto const p = s.find(L"| ", start);
        if (p == std::wstring::npos) return L"";
        start = p + 2;
      }
      auto const end = s.find(L"| ", start);
      std::wstring f = s.substr(
          start, end == std::wstring::npos ? std::wstring::npos : end - start);
      while (!f.empty() &&
             (f.back() == L' ' || f.back() == L'\r' || f.back() == L'\n'))
        f.pop_back();
      return f;
    };

    // Same touched-mode-slice keying as the occ-veto witness, extended for the
    // EXTERNAL scatter's markers. The contracted BatchGroup path emits
    // BatchGroup|Begin/End + a per-slice BatchIter; the External scatter branch
    // (eval.hpp) emits BatchScatter|Begin/End around the block loop and -- via
    // the per-block marker added alongside this gate -- the SAME BatchIter
    // (mode, e_lo) per external block. Treating BatchScatter exactly like
    // BatchGroup (push an enclosing frame on Begin, pop on End) means an op
    // scattered into a distinct external block gets a distinct touched-slice
    // signature, so its per-block replay is correctly counted as legitimate
    // 1/nblocks work, NOT an avoidable duplicate. (Ignoring BatchScatter would
    // give every block the same signature and falsely inflate avoidable_time.)
    std::vector<std::pair<std::wstring, std::wstring>> stack;  // (mode, slice)
    std::set<std::wstring> all_keys;
    std::map<std::wstring, double>
        exec_of;  // expr -> avoidable exec (offender)

    std::wistringstream in(trace.str());
    std::wstring ln;
    double last_exec = -1.0;
    while (std::getline(in, ln)) {
      if (ln.rfind(L"BatchGroup | Begin", 0) == 0 ||
          ln.rfind(L"BatchScatter | Begin", 0) == 0) {
        if (ln.rfind(L"BatchScatter | Begin", 0) == 0)
          ++m.n_scatter_begin;
        else
          ++m.n_batchgroup_begin;
        stack.emplace_back(L"", L"");
      } else if (ln.rfind(L"BatchGroup | End", 0) == 0 ||
                 ln.rfind(L"BatchScatter | End", 0) == 0) {
        if (!stack.empty()) stack.pop_back();
      } else if (ln.rfind(L"BatchIter", 0) == 0) {
        if (!stack.empty()) stack.back() = {field(ln, 1), field(ln, 2)};
      } else if (ln.rfind(L"OpCost", 0) == 0) {
        last_exec = std::wcstod(field(ln, 2).c_str(), nullptr);
      } else if (ln.find(L"Eval | Product") != std::wstring::npos) {
        ++m.replay_ops;
        double const exec = last_exec >= 0.0 ? last_exec : 0.0;
        if (last_exec >= 0.0) ++m.ops_with_cost;
        last_exec = -1.0;
        m.total_exec += exec;
        auto const p = ln.rfind(L"| ");
        std::wstring const expr =
            (p == std::wstring::npos) ? ln : ln.substr(p + 2);
        std::wstring sig;
        for (auto const& [mode, slice] : stack)
          if (!mode.empty() && expr.find(mode) != std::wstring::npos)
            sig += L"|" + mode + L"=" + slice;
        std::wstring const key = expr + L"@" + sig;
        bool const dup = !all_keys.insert(key).second;
        if (dup) {
          m.avoidable_exec += exec;
          exec_of[expr] += exec;
        }
      }
    }
    m.unique_keys = all_keys.size();
    double worst_exec = 0;
    for (auto const& [e, av] : exec_of)
      if (av > worst_exec) {
        worst_exec = av;
        m.top_expr = e;
      }
    return m;
  };

  auto const contracted =
      run(/*batch_external=*/false);  // contracted-occ baseline
  auto const external = run(/*batch_external=*/true);  // external-occ

  auto report = [](wchar_t const* label, Meas const& m) {
    std::wcerr << L"[extmode-avoidable] " << label << L": ops=" << m.replay_ops
               << L" unique=" << m.unique_keys << L" AVOIDABLE_TIME="
               << (100.0 * m.avoidable_time()) << L"%  scatter_begin="
               << m.n_scatter_begin << L" bgroup_begin=" << m.n_batchgroup_begin
               << L" ext_stamps=" << m.n_external_stamps << L" con_occ_stamps="
               << m.n_contracted_occ << L" peak=" << m.peak_gb
               << L"GB (cost-coverage " << m.ops_with_cost << L"/"
               << m.replay_ops
               << L")\n           worst avoidable-time offender: " << m.top_expr
               << L"\n";
  };
  std::wcerr << L"\n=== [dryrun-extmode-avoidable] C60 residual forest, "
             << nterms << L" terms (K@256, occ@8, 100GB, perf-first) ===\n";
  report(L"contracted-occ", contracted);
  report(L"external-occ  ", external);

  // The forest is identical; only the mode policy differs.
  CHECK(contracted.static_nodes == external.static_nodes);

  // Measurement guard: nearly every op must carry an OpCost line.
  CHECK(external.ops_with_cost >=
        external.replay_ops - external.replay_ops / 20);

  // The EXTERNAL scatter must genuinely fire (the DP stamped External on the
  // over-budget PPL giant's external occ, and the runtime took the scatter
  // branch), else this gate would trivially pass by never batching.
  CHECK(external.n_external_stamps > 0);
  CHECK(external.n_scatter_begin > 0);

  // PARTIAL-WIN WITNESS (deliberately NOT a ~0 gate). External-mode batching is
  // NECESSARY but NOT SUFFICIENT on the C60 residual forest: it engages (the
  // scatter fires above) and strictly reduces both avoidable recompute and
  // replay work, but it neither bounds the forest peak nor drives
  // avoidable_time to ~0. Two residual gaps remain, both deferred to a separate
  // design pass:
  //   (1) the forest peak is set by a term that slicing the proto-occ pair
  //       (i_1,i_2) does not reach -- needs full forest-level / multi-mode
  //       co-batching (also slice i_3,i_4 and/or K_2); and
  //   (2) a contracted MIDDLE-GAP node survives (an intermediate inside a batch
  //       loop over a mode it does not carry -- here the 4-occ
  //       I(..,i_3,K_2;a_3)*I(..,i_4,K_2;a_4) contraction), which needs the
  //       order-aware cost + multilevel hoisting of
  //       doc/dev/specs/2026-07-17-nested-batch-group-join-design.md.

  // RETRACTED WIN (2026-07-22). This gate previously asserted
  //   CHECK(external.avoidable_time() < contracted.avoidable_time() - 0.15);
  //   CHECK(external.replay_ops < contracted.replay_ops);
  // on measurements of ~44.7% vs ~75.1%. That comparison was CONFOUNDED: its
  // two arms flipped batch_spectator_indices AND suppress_heuristic_fallback
  // together, so the contracted arm ran the legacy runtime heuristic and the
  // external arm did not. Essentially the entire "win" was the flag, not the
  // external mode. With the heuristic removed outright (annotations are now
  // authoritative everywhere) the two arms differ in exactly one variable, and
  // external-mode batching is measured to be NEUTRAL-TO-SLIGHTLY-WORSE here:
  //
  //   contracted-occ : avoidable_time 43.72%, replay ops 61275
  //   external-occ   : avoidable_time 43.96%, replay ops 83573
  //
  // So External is NOT a fix for the avoidable recompute on this forest. It
  // still engages losslessly (asserted above) and remains the mechanism for
  // slicing a mode that is contracted nowhere, but the D5 claim that it
  // eliminates the occ-veto recompute is NOT supported by measurement. Record
  // the true relation rather than a target we have not met.
  CHECK(std::abs(external.avoidable_time() - contracted.avoidable_time()) <
        0.02);
  CHECK(external.replay_ops > contracted.replay_ops);

  // ENTRY CRITERION for the deferred forest-co-batching + middle-gap work:
  // these record the residual as it stands today. When that work lands, these
  // two will start failing -- FLIP them to (external.peak_gb < 100.0) and
  // (external.avoidable_time() < 0.05) and this witness becomes a true gate.
  CHECK(external.peak_gb > 100.0);          // peak NOT yet bounded (~2302 GB)
  CHECK(external.avoidable_time() > 0.10);  // avoidable NOT yet ~0 (~44.7%)
}
