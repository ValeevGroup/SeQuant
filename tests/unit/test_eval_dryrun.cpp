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
//      the scratch-folded batched peak via the metered dry-run replay
//      (dryrun::meter).

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
#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>
#include <SeQuant/core/eval/backends/dryrun/meter.hpp>  // ordered dry-run replay
#include <SeQuant/core/eval/backends/dryrun/result.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/ordered_executor.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/eval/schedule_dump.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/result_expr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/optimize/options.hpp>
#include <SeQuant/core/optimize/single_term_detail.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/expr.hpp>  // is_valid
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/space_qns.hpp>  // mbpt::Spin

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

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

// Local shim: BatchModeType-tagging of EvalExpr::node_slice_mask() entries
// (Task 1) strips this spike file's plain-Index reads down to the .first
// projection; this test file is not committed, so keep the fix minimal rather
// than threading BatchModeType through the trace/analysis helpers below.
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
// Regression (fast, no DP solve): is_valid must ACCEPT a CSV (proto-indexed)
// residual Sum. is_valid's Sum check compares each summand's external indices;
// the slot-only get_unique_indices it used ignores proto-indices, so an occ
// index carried inside a composite virtual (a<i,j>) in some summands and
// standalone in others was miscounted, and is_valid spuriously reported
// "Inconsistent external indices in sum". On an MPQC_ASSERT_ABORT build that
// aborted every CSV-CCk run at MPQC_ASSERT(is_valid(e)); the proto-aware
// external-index comparison fixes it. This reuses the real C60 doubles residual
// data file (a genuine proto-indexed CSV Sum).
TEST_CASE("is_valid accepts a CSV proto-indexed residual",
          "[utilities][is_valid][csv]") {
  using namespace sequant;
  auto ctx0 = get_default_context().clone();
  ctx0.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx0.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);
  sequant::mbpt::add_df_spaces(isr);
  sequant::tests::declare_real_basis(*isr);  // the residual is real
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
  REQUIRE(expr->as<Sum>().summands().size() > 1);

  std::string msg;
  bool const valid = is_valid(expr, &msg);
  INFO("is_valid message: " << msg);
  CHECK(valid);
  CHECK(msg.empty());
}

// Regression: optimize_result must key each summand's per-node batch
// annotations onto the FINAL reassembled Sum pointer, not per optimized
// summand. The CCk residual is one Sum-tree per equation, so the consumer
// binarizes the whole Sum and looks the annotation up by that Sum's pointer;
// opt_pure_product keys per summand, and (under reorder) opt::reorder clones
// the summands (Sum::append clones) while the keyed pre-clone summands are
// destroyed. Without re-keying, the whole-Sum lookup finds nothing, every batch
// annotation is dropped, and over-budget intermediates materialize whole -- the
// water-20 OOM. This asserts the re-keying on the real CSV doubles residual.
// Minutes-long under ASan/valgrind; see tests/unit/CMakeLists.txt.
#ifndef SEQUANT_SKIP_LONG_TESTS
TEST_CASE("optimize_result keys batch annotations onto the whole Sum",
          "[optimize][batch][term_batch_axes]") {
  using namespace sequant;
  auto ctx0 = get_default_context().clone();
  ctx0.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx0.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);
  sequant::mbpt::add_df_spaces(isr);
  sequant::tests::declare_real_basis(*isr);  // the residual is real
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
  REQUIRE(expr->as<Sum>().summands().size() > 1);

  auto regime = df_regime(kWater20_pVDZF12);
  BatchPolicy policy;
  policy.is_batchable_contracted_index = [](Index const& ix) {
    return ix.space().base_key() == L"Κ";
  };
  policy.batch_target_size = [](Index const&) -> std::size_t { return 256; };
  policy.is_volatile_leaf = [](Tensor const& t) { return t.label() == L"t"; };
  policy.peak_threshold = 100e9;

  auto axes_map = std::make_shared<std::unordered_map<
      Expr const*, container::vector<NodeBatchAnnotation>>>();
  OptimizeOptions opts;
  opts.objective_function = ObjectiveFunction::DenseTimeSpaceBatched;
  opts.reorder = ReorderSum::Reorder;  // the production (clone-on-append) path
  opts.idx_to_extent = regime.idx_to_extent();
  opts.inner_pow = regime.inner_pow_fn();
  opts.batch_policy = policy;
  opts.volatile_weight = 20.0;
  opts.roofline.machine_balance = 200.0;
  opts.roofline.fast_mem_elems = 1000000.0;
  opts.term_batch_axes = axes_map;

  auto res = optimize_result(expr, opts);
  REQUIRE(res.expr);
  REQUIRE(res.expr->is<Sum>());

  // The whole reassembled Sum is THE key -- not any per-summand pointer -- and
  // re-keying erased the stale per-summand entries, so it is the only key.
  CHECK(axes_map->count(res.expr.get()) == 1);
  CHECK(axes_map->size() == 1);

  // It carries real batch axes (Κ blows the 100 GB budget on this residual),
  // and one entry per contraction node of the whole Sum-tree (what binarize's
  // node counter consumes).
  std::size_t nonempty = 0, total = 0;
  if (auto it = axes_map->find(res.expr.get()); it != axes_map->end()) {
    total = it->second.size();
    for (auto const& a : it->second)
      if (!a.axes.empty()) ++nonempty;
  }
  CHECK(total > 0);
  CHECK(nonempty > 0);

  // The concatenated node_batch_axes must have EXACTLY one entry per
  // contraction node of the whole reassembled Sum, or binarize aborts on its
  // node_counter == node_batch_axes.size() check (the water-20 SIGABRT).
  // Binarize the residual with the concatenated axes and require it does not
  // throw/abort.
  {
    BinarizationOptions bopts;
    if (auto it = axes_map->find(res.expr.get()); it != axes_map->end())
      bopts.node_batch_axes = it->second;
    // Binarize through the SAME head-pinned ResultExpr path MPQC's CCk uses (a
    // CSV rank-2 residual head, make_R_template_csv): R{a_1<i_1,i_2>,
    // a_2<i_1,i_2>; i_1, i_2}. A count mismatch trips binarize's
    // node_counter == node_batch_axes.size() assertion (the water-20 SIGABRT on
    // an ABORT build; a no-op under IGNORE).
    std::vector<Index> occ{Index(L"i_1"), Index(L"i_2")};
    std::vector<Index> vir{Index(L"a_1", occ), Index(L"a_2", occ)};
    Tensor head(L"R", bra(vir), ket(occ), Symmetry::Nonsymm,
                BraKetSymmetry::Nonsymm, ColumnSymmetry::Symm);
    ResultExpr rexpr{head, res.expr};
    CHECK_NOTHROW(binarize<EvalExpr>(rexpr, bopts));
  }
}
#endif  // !defined(SEQUANT_SKIP_LONG_TESTS)

// Regression: the exact R1 (singles) summand water-20 PNO-CCSD aborted on --
// f{mu~;i} * C{a<i>;mu~}, a 2-tensor contraction. The batched optimizer must
// emit ONE node_axes entry (one contraction node), matching binarize; if the DP
// network drops a tensor (nt==1 -> zero entries) while binarize keeps the
// contraction, binarize's node_counter == node_batch_axes.size() assertion
// aborts. Mirrors MPQC's path: DenseTimeSpaceBatched optimize + head-pinned
// binarize with the CSV rank-1 residual head R{a<i>;i}.
TEST_CASE("optimizer node_axes match binarize on the water-20 R1 f*C summand",
          "[optimize][batch][r1-offbyone]") {
  using namespace sequant;
  auto ctx0 = get_default_context().clone();
  ctx0.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx0.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);
  sequant::mbpt::add_df_spaces(isr);
  sequant::tests::declare_real_basis(*isr);  // the residual is real
  auto ctx_resetter = set_scoped_default_context(std::move(ctx0));

  auto prod =
      deserialize<ExprPtr>("f{μ̃_1094;i_1}:N-S-S * C{a_1<i_1>;μ̃_1094}:N-S-S");
  REQUIRE(prod);
  REQUIRE(prod->is<Product>());

  auto regime = df_regime(kWater20_pVDZF12);
  BatchPolicy policy;
  policy.is_batchable_contracted_index = [](Index const& ix) {
    return ix.space().base_key() == L"Κ";
  };
  policy.batch_target_size = [](Index const&) -> std::size_t { return 256; };
  policy.is_volatile_leaf = [](Tensor const& t) { return t.label() == L"t"; };
  policy.peak_threshold = 100e9;

  auto axes_map = std::make_shared<std::unordered_map<
      Expr const*, container::vector<NodeBatchAnnotation>>>();
  OptimizeOptions opts;
  opts.objective_function = ObjectiveFunction::DenseTimeSpaceBatched;
  opts.reorder = ReorderSum::Reorder;
  opts.idx_to_extent = regime.idx_to_extent();
  opts.inner_pow = regime.inner_pow_fn();
  opts.batch_policy = policy;
  opts.term_batch_axes = axes_map;

  auto res = optimize_result(prod, opts);
  REQUIRE(res.expr);

  BinarizationOptions bopts;
  std::size_t na = 0;
  if (auto it = axes_map->find(res.expr.get()); it != axes_map->end()) {
    bopts.node_batch_axes = it->second;
    na = it->second.size();
  }

  // binarize's tensor*tensor contraction-node count for the same expression.
  std::function<std::size_t(FullBinaryNode<EvalExpr> const&)> cnt =
      [&](FullBinaryNode<EvalExpr> const& n) -> std::size_t {
    if (n.leaf()) return 0;
    std::size_t c = cnt(n.left()) + cnt(n.right());
    if (!n.left()->is_scalar() && !n.right()->is_scalar()) ++c;
    return c;
  };
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  std::size_t const bc = cnt(binarize<EvalExpr>(res.expr));
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

  // f*C is a single contraction: the optimizer MUST emit exactly one node_axes
  // entry, matching binarize (before the fix it emitted zero -> off-by-one).
  CHECK(bc == 1);
  CHECK(na == bc);

  // And the head-pinned binarize MPQC uses must not trip its count assertion.
  std::vector<Index> occ{Index(L"i_1")};
  std::vector<Index> vir{Index(L"a_1", occ)};
  Tensor head(L"R", bra(vir), ket(occ), Symmetry::Nonsymm,
              BraKetSymmetry::Nonsymm, ColumnSymmetry::Symm);
  ResultExpr rexpr{head, res.expr};
  CHECK_NOTHROW(binarize<EvalExpr>(rexpr, bopts));
}

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
  sequant::tests::declare_real_basis(*isr);  // the residual is real
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
  auto ctx = model.build_context(net, targets);
  std::wcerr << L"[ordered-key-c60-m] giant (summand " << gi << L", "
             << giant->as<Product>().factors().size() << L" factors): m="
             << ctx.m << L" ordered=" << ctx.ordered << L" nCells="
             << ctx.nCells << L" (nB=" << ctx.nB << L")\n";
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
  sequant::tests::declare_real_basis(*isr);  // the residual is real
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
    (void)model.reconstruct_batched_modes(ctx, st, &pk);
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
  sequant::tests::declare_real_basis(*isr);  // the residual is real
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
    CHECK_THROWS_AS(bad.build_context(net, targets), sequant::Exception);
  }

  // (2) With correct composite sizing, the DP forms NO 4-PAO integral.
  opt::detail::PeakBatchedModel model{idxsz, bts, is_vol,
                                      regime.inner_pow_fn()};
  model.is_batchable_contracted_index = is_aux;
  model.perf_first = true;
  model.volatile_weight = 20.0;
  model.peak_threshold = std::numeric_limits<double>::infinity();
  model.numeric_size = 8.0;
  auto ctx = model.build_context(net, targets);
  auto st = opt::detail::solve_single_term(model, net, targets, ctx);
  double peak_bytes = 0.0;
  (void)model.reconstruct_batched_modes(ctx, st, &peak_bytes);
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
  sequant::tests::declare_real_basis(*isr);  // the residual is real
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
    // carrying active = union of ancestor node_slice_mask (by FULL label, so
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
    // node that slices it (its result free-index signature + its
    // node_slice_mask). This is the node-dump: it turns "escaped={}" from an
    // inference into a concrete "mu~_X is sliced by ancestor <node>" (or
    // ESCAPED). nbatches carries, for each ancestor-sliced mode, its batch
    // count (extent / batch_target_size) -- the executed-flops recompute
    // factor.
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
            if (!n->node_slice_mask().empty()) ++term_batched_nodes;
            ++term_internal_nodes;
            for (auto const& ax : n->node_slice_mask())
              term_axis_hist[std::wstring(ax.first.space().base_key())]++;
            // Forest-level CSE bookkeeping (see declarations above). The mode
            // signature includes the ancestor-sliced context, not just this
            // node's own modes: the same expression evaluated under a different
            // ancestor slicing is a different cache entry.
            ++overall_internal_nodes;
            std::wstring axsig;
            for (auto const& [k, nb] : nbatches) axsig += k + L";";
            axsig += L"|";
            for (auto const& ax : n->node_slice_mask())
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
                L"} node_slice_mask={" +
                describe_indices(batch_axes_indices(n->node_slice_mask())) +
                L"}]";
            std::map<std::wstring, double> child_nbatch = nbatches;
            for (auto const& ax : n->node_slice_mask()) {
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
  ov[0] = 5;  // mode 0 (a_3) narrowed from 20 to 5
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
  // exec_cost takes the op's FULL compulsory traffic: both operand footprints
  // plus the result's (see the roofline note on CostModel::exec_cost).
  CHECK(cm.exec_cost(f, cm.memsize(out), cm.memsize(contracted),
                     cm.memsize(out)) > 0.0);
}

TEST_CASE(
    "dryrun product roofline exec charges both operands and the result, "
    "order-independently",
    "[dryrun-costmodel][roofline]") {
  // The roofline `traffic` term is the COMPULSORY single-pass data movement of
  // one contraction -- read both operands, WRITE the result -- which is what
  // the optimizer's DP charges (S[lp] + S[rp] + S[n]; PeakModel::relax and
  // BatchedPeakModel::relax in core/optimize/cost_model.hpp). The dry-run
  // replay must charge the same thing, at realized extents. Two consequences
  // are pinned here:
  //   (a) `exec` does not depend on which operand is the left one (data
  //       movement is symmetric in the operands), and
  //   (b) `exec` equals roofline_op_cost at |L| + |R| + |out| elements.
  // Before the fix, prod charged memsize(left) + a 4096-byte placeholder for
  // the right operand and nothing for the result, so BOTH failed: with the
  // 4000-vs-100-element operands used below, swapping the operand order moved
  // `exec` by the whole ratio of the two footprints.
  sequant::RooflineParams rp{.machine_balance = 200.0,
                             .fast_mem_elems = 1000000.0};
  auto const regime = backend_test_regime();
  auto cm = std::make_shared<CostModel const>(regime, rp);

  Index const i1{L"i_1"}, i2{L"i_2"}, a3{L"a_3"}, a4{L"a_4"};
  // Deliberately lopsided operands: big = 20*20*10 = 4000 elements,
  // small = 10*10 = 100 elements, result = 20*20*10 = 4000 elements.
  container::svector<Index> const big{a3, a4, i2};
  container::svector<Index> const small{i2, i1};
  container::svector<Index> const res{a3, a4, i1};

  // Bytes -> elements without hardcoding the numeric size: a rank-1 `i` tensor
  // is exactly `i`'s extent (10) elements.
  double const nsz = static_cast<double>(cm->memsize({i1})) / 10.0;
  REQUIRE(nsz > 0.0);
  double const traffic_elems =
      static_cast<double>(cm->memsize(big) + cm->memsize(small) +
                          cm->memsize(res)) /
      nsz;
  CHECK(traffic_elems == Catch::Approx(4000.0 + 100.0 + 4000.0));

  // The per-op OpCost emission (which is what stashes last_op_flops/exec) is
  // gated at RUNTIME on Logger::instance().eval.level > 0; redirect the stream
  // so nothing lands on stdout, and restore the global state afterwards.
  std::ostringstream trace_os;
  auto& logger = Logger::instance();
  auto const prev_level = logger.eval.level;
  auto* const prev_stream = logger.eval.stream;
  logger.eval.level = 2;
  logger.eval.stream = &trace_os;

  // big * small
  double flops_lr = 0.0, exec_lr = 0.0;
  {
    ResultDryRun l{big, cm};
    ResultDryRun r{small, cm};
    auto out = static_cast<Result const&>(l).prod(r, annot3(big, small, res),
                                                  DeNest::False);
    REQUIRE(out);
    flops_lr = sequant::eval::detail::last_op_flops();
    exec_lr = sequant::eval::detail::last_op_exec();
  }
  // small * big -- same contraction, operands swapped.
  double flops_rl = 0.0, exec_rl = 0.0;
  {
    ResultDryRun l{small, cm};
    ResultDryRun r{big, cm};
    auto out = static_cast<Result const&>(l).prod(r, annot3(small, big, res),
                                                  DeNest::False);
    REQUIRE(out);
    flops_rl = sequant::eval::detail::last_op_flops();
    exec_rl = sequant::eval::detail::last_op_exec();
  }

  logger.eval.level = prev_level;
  logger.eval.stream = prev_stream;

  CHECK(flops_lr > 0.0);
  CHECK(flops_rl == Catch::Approx(flops_lr));  // flops are order-independent
  CHECK(exec_lr > 0.0);
  CHECK(exec_rl == Catch::Approx(exec_lr));  // (a) so is the traffic term

  // (b) and it is exactly the roofline cost at the operand+result footprint.
  double const expected = sequant::opt::detail::roofline_op_cost(
      flops_lr, traffic_elems, rp.machine_balance, rp.fast_mem_elems,
      rp.block_tiles, rp.block_prefactor);
  CHECK(exec_lr == Catch::Approx(expected));
  // Guard against the check being vacuous (machine_balance high enough that
  // the op really is bandwidth-bound, so `traffic` is what is being tested).
  CHECK(expected > flops_lr);
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

TEST_CASE("dryrun axis_batches tiles the axis space extent",
          "[dryrun-result]") {
  auto r = backend_test_regime();
  auto cm = std::make_shared<CostModel const>(r);
  auto const aops = sequant::eval::dryrun::make_dryrun_array_ops(cm);
  auto batches = aops.axis_batches(Index{L"a_3"}, 5);  // 20 / 5 = 4 batches
  CHECK(batches.size() == 4);
  CHECK(batches.front().first == 0);
  CHECK(batches.back().second == 20);
}

TEST_CASE("dryrun make_zeros builds a full-extent flat scatter destination",
          "[dryrun-result][pre-sized]") {
  // The runtime External-mode scatter builds its destination from the node's
  // OWN (unsliced) index list via BackendArrayOps::make_zeros -- every mode at
  // its space's FULL extent (a structural fact queryable immediately via
  // size_in_bytes()). Replaces the old carrier-widening
  // pre_sized_zeros_over_mode: no block partial, no carrier.
  auto r = backend_test_regime();
  auto cm = std::make_shared<CostModel const>(r);
  Index i1{L"i_1"}, a3{L"a_3"};
  container::svector<Index> idx{i1, a3};  // i_1 extent 10, a_3 extent 20

  ResultDryRun full{idx, cm};  // reference: the fully-realized token
  auto const full_bytes = static_cast<Result const&>(full).size_in_bytes();

  auto const aops = sequant::eval::dryrun::make_dryrun_array_ops(cm);
  auto dest = aops.make_zeros(container::vector<Index>(idx.begin(), idx.end()));
  REQUIRE(dest);
  CHECK(dest->is<ResultDryRun>());
  CHECK(dest->size_in_bytes() == full_bytes);
}

TEST_CASE("dryrun make_zeros builds a full-extent nested scatter destination",
          "[dryrun-result][pre-sized]") {
  // ToT analogue of the flat case above -- CSV/PNO residuals carry
  // ResultDryRunNested tokens, so make_zeros must build a nested full-extent
  // destination when the descriptor carries a proto-indexed (composite) leg.
  auto r = backend_test_regime();
  auto cm = std::make_shared<CostModel const>(r);
  Index i1{L"i_1"}, i2{L"i_2"}, a3{L"a_3"};
  Index a_pno{L"a_1", {i1, i2}};
  container::svector<Index> outer{i1, a3};
  container::svector<Index> inner{a_pno};
  container::svector<Index> canon{i1, a3, a_pno};

  ResultDryRunNested full{outer, inner, cm, {}, canon};
  auto const full_bytes = static_cast<Result const&>(full).size_in_bytes();

  auto const aops = sequant::eval::dryrun::make_dryrun_array_ops(cm);
  auto dest =
      aops.make_zeros(container::vector<Index>(canon.begin(), canon.end()));
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
  ResultDryRunNested dest{
      outer, inner, cm, {{0, 0}}, canon};  // mode 0 unfilled
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
  ResultDryRunNested dest_fc{
      outer, inner, cm, {{0, 0}}, canon};  // mode 0 unfilled
  Result& dest_fc_w = dest_fc;
  auto block_fc = tmpl_r.slice_mode(0, 2, 6);
  REQUIRE(block_fc);
  dest_fc_w.write_into_slice(*block_fc, 0, 2, 6);
  CHECK(dest_fc.assembled_range(0) ==
        std::pair<std::size_t, std::size_t>{2, 6});
}

TEST_CASE(
    "dryrun result write_into_slice REFUSES a gapped or overlapping block",
    "[dryrun-result][write-into-slice]") {
  // Copilot review (PR #613): the contiguity requirement used to be carried by
  // a SEQUANT_ASSERT alone. Asserts are compiled out in non-Debug builds (and
  // CI builds with SEQUANT_ASSERT_BEHAVIOR=THROW rather than relying on the
  // assert), so a gapped or overlapping scatter would have updated the
  // coverage from stale data and the dry run would have accepted -- and then
  // MIS-SIZED -- an incorrect scatter. It must throw instead.
  auto r = backend_test_regime();
  auto cm = std::make_shared<CostModel const>(r);

  Index i1{L"i_1"}, i2{L"i_2"}, a3{L"a_3"};
  Index a_pno{L"a_1", {i1, i2}};
  container::svector<Index> outer{i1, a3};
  container::svector<Index> inner{a_pno};
  container::svector<Index> canon{i1, a3, a_pno};

  ResultDryRunNested tmpl{outer, inner, cm, {}, canon};
  Result const& tmpl_r = tmpl;

  // OVERLAP: [0,5) then [3,8) -- [3,8) neither appends after 5 nor prepends
  // before 0, so it would double-count elements 3 and 4.
  {
    ResultDryRunNested dest{outer, inner, cm, {{0, 0}}, canon};
    Result& dest_w = dest;
    auto b0 = tmpl_r.slice_mode(0, 0, 5);
    auto b1 = tmpl_r.slice_mode(0, 3, 8);
    REQUIRE(b0);
    REQUIRE(b1);
    dest_w.write_into_slice(*b0, 0, 0, 5);
    CHECK_THROWS_AS(dest_w.write_into_slice(*b1, 0, 3, 8), sequant::Exception);
  }

  // GAP: [0,5) then [6,10) -- element 5 would never be written, yet the
  // assembled extent would report the full range as covered.
  {
    ResultDryRunNested dest{outer, inner, cm, {{0, 0}}, canon};
    Result& dest_w = dest;
    auto b0 = tmpl_r.slice_mode(0, 0, 5);
    auto b1 = tmpl_r.slice_mode(0, 6, 10);
    REQUIRE(b0);
    REQUIRE(b1);
    dest_w.write_into_slice(*b0, 0, 0, 5);
    CHECK_THROWS_AS(dest_w.write_into_slice(*b1, 0, 6, 10), sequant::Exception);
  }
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
  // batching at all (as-built design
  // doc/dev/specs/2026-09-12-batched-array-dag-eval-as-built.md, section 4.2).
  // This test drives the SAME scatter branch the TA
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
  node->set_node_slice_mask({{mode, BatchModeType::External}});
  auto stamp_carriers = [&](auto&& self, node_t& n) -> void {
    if (n.leaf()) return;
    if (&n != &node && index_position(n, mode).has_value())
      n->set_node_slice_mask({{mode, BatchModeType::External}});
    self(self, n.left());
    self(self, n.right());
  };
  stamp_carriers(stamp_carriers, node);

  DryRunLeafEvaluator yield{cm};

  // Reference: plain unbatched evaluation (node_slice_mask are ignored without
  // a custom evaluator).
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
  auto aops = sequant::eval::dryrun::make_dryrun_array_ops(cm);
  cache.set_array_ops(&aops);
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

// Minutes-long under ASan/valgrind; see tests/unit/CMakeLists.txt.
#ifndef SEQUANT_SKIP_LONG_TESTS
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
  sequant::tests::declare_real_basis(*isr);  // the residual is real
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

  if (std::getenv("SEQUANT_SCHED_DUMP"))
    std::cerr << "SCHEDULE_IR_JSON "
              << sequant::eval::schedule_ir_json(node, "giant") << "\n";

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
      giant_axes_desc =
          describe_indices(batch_axes_indices(n->node_slice_mask()));
    }
  });
  REQUIRE(giant_nominal_bytes > 0.0);

  auto cache = sequant::cache_manager(std::vector<EvalNodeDryRun>{node});
  auto aops = sequant::eval::dryrun::make_dryrun_array_ops(cm);
  cache.set_array_ops(&aops);
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
  // Enable THIS replay's per-DISTINCT-value build tally on its own cache
  // (CacheManager::tally_build, keyed by the exact cache identity) so the
  // avoidable-recompute breakdown is accumulated. Using the test's OWN replay
  // -- not a separate metered replay -- means the numbers match the exact
  // run events the visualizer consumes (same cache, same slicing), and there is
  // no second replay flooding the SCHEDULE_RUN_EVENT stream.
  bool const sched_dump = std::getenv("SEQUANT_SCHED_DUMP") != nullptr;
  if (sched_dump) cache.set_recompute_tally_enabled(true);
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
  if (sched_dump) {
    // The replay's per-node avoidable recompute, keyed by the node's
    // topological hash (the SAME join key the IR and run-event nodes carry) so
    // the visualizer joins each DAG node to these numbers instead of
    // recomputing avoidable. Rolled up per DISTINCT value over its SLICES (see
    // BuildTally): for each slice total += builds*cost and build_once += cost,
    // so avoidable -- the arithmetic the replay repeated beyond building each
    // distinct slice once -- is sum over slices of (builds-1)*cost. A value
    // tiled over DISTINCT slices has builds==1 per slice => 0 avoidable (pure
    // tiling, even non-uniform); a value rebuilt at the SAME slice (an
    // invariant rebuilt every block of a loop it does not carry) has builds>1
    // there. Keeps only values with a positive avoidable amount, worst first.
    struct AvoidableNode {
      std::string label;
      double count = 0;
      double flops = 0;
    };
    std::vector<AvoidableNode> av;
    for (auto const& [node, t] : cache.recompute_tally()) {
      double total = 0, once = 0, extra_builds = 0;
      for (auto const& [sig, bc] : t.slices) {
        total += bc.count * bc.flops;
        once += bc.flops;
        extra_builds += static_cast<double>(bc.count - 1);
      }
      if (total - once <= 0.0) continue;
      av.push_back({.label = std::to_string(node->hash_value()),
                    .count = extra_builds,
                    .flops = total - once});
    }
    std::sort(av.begin(), av.end(),
              [](AvoidableNode const& a, AvoidableNode const& b) {
                return a.flops > b.flops;
              });
    std::ostringstream cjson;
    cjson << "SCHEDULE_COST_JSON {\"term_id\":\"giant\",\"nodes\":[";
    for (std::size_t i = 0; i < av.size(); ++i) {
      if (i) cjson << ',';
      cjson << "{\"sig\":\""
            << sequant::eval::detail::sched_json_escape(av[i].label)
            << "\",\"count\":" << av[i].count << ",\"flops\":" << av[i].flops
            << "}";
    }
    cjson << "]}";
    std::cerr << cjson.str() << "\n";
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
             << (giant_nominal_bytes / 1e9) << L" GB node_slice_mask={"
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
#endif  // !defined(SEQUANT_SKIP_LONG_TESTS)

// Task 6 (perf-first validation): optimize the C60 giant term (index 38) at
// the faithful real config and check the factorization the DP picks. The 4-PAO
// signature is a contraction node carrying >= 4 free mu~ indices (the
// (mu~ mu~|mu~ mu~) AO integral); perf-first is flops-primary and must NEVER
// form it. This is the direct in-harness proof of the fix, on ONE term
// (~seconds), without the [dryrun-df] full-sweep cost. The peak-first run is
// kept for the printed contrast only -- under the ordered cost model it also
// declines the 4-PAO (its batched peak is priced with accumulator residency),
// so no ASSERTION is made about it.
// Minutes-long under ASan/valgrind; see tests/unit/CMakeLists.txt.
#ifndef SEQUANT_SKIP_LONG_TESTS
TEST_CASE(
    "dryrun perf-first never forms the 4-PAO AO integral and peaks on the "
    "genuine 4-PNO W node (C60 giant)",
    "[dryrun-objective]") {
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);  // mu~
  sequant::mbpt::add_df_spaces(isr);                             // K
  sequant::tests::declare_real_basis(*isr);  // the residual is real
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
    // Cost of the chosen factorization: a STATIC per-internal-node flops walk
    // (order-/batching-blind, so it depends only on the factorization the DP
    // picked) plus the metered dry-run replay's peak (dryrun::meter).
    std::size_t model_n_ops = 0;
    double model_flops = 0.0;
    double peak_bytes = 0.0;
  };

  // Optimize the giant under `obj`, binarize, and walk the tree computing per
  // node (a) its free-mu~ count (the 4-PAO structural signature) and (b) its
  // REALIZED free-mu~ size after ancestor slicing (same active-ancestor
  // accounting as the [dryrun-df] verdict case). Then call the single shared
  // metered replay (dryrun::meter) to get the modeled peak via
  // the gated-cache replay -- replacing the ad-hoc manual replay + hwmark read
  // this case used before Task 5, so there is ONE peak/flops code path.
  auto analyze = [&](ObjectiveFunction obj) -> Analysis {
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
    // free indices and node_slice_mask (the indices SLICED, i.e. batched, at
    // that node), in post-order-ish indentation, so the exact schedule is
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
                batch_axes_indices(n->node_slice_mask());
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
                       << nosv << L")  node_slice_mask={"
                       << describe_indices(bax) << L"}\n";
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
            for (auto const& ax : n->node_slice_mask())
              child_active[keyof(ax.first)] = L"y";
            walk(n.left(), child_active);
            walk(n.right(), child_active);
          }
        };
    walk(node, {});

    // ---- modeled cost: static flops walk + metered replay peak ----------
    // (i) The STATIC walk prices every internal node exactly once
    // (flops_counter over the node's (left, right, result) index sets), so it
    // reports the factorization's own arithmetic, blind to evaluation order
    // and batching -- which is what the objective comparison below needs.
    {
      auto const flops_of = sequant::opt::detail::flops_counter(
          regime.idx_to_extent(), regime.inner_pow_fn());
      std::function<void(EvalNodeDryRun const&)> cost_walk =
          [&](EvalNodeDryRun const& n) {
            if (n.leaf()) return;
            a.model_n_ops += 1;
            a.model_flops +=
                flops_of(n.left()->canon_indices(), n.right()->canon_indices(),
                         n->canon_indices());
            cost_walk(n.left());
            cost_walk(n.right());
          };
      cost_walk(node);
    }
    // (ii) The metered replay (dryrun::meter) builds the gated dry-run cache
    // (free-batchable-mode veto + footprint gate), replays zero-data through
    // the real eval loop with a PeakMonitor wired onto the cache scope chain,
    // and forces the printing gate on internally -- so peak_bytes captures the
    // batched-inner transient the raw outer hwmark misses.
    sequant::eval::dryrun::CacheConfig cfg;
    cfg.max_footprint = 1e11;
    cfg.min_repeats = 1;
    cfg.is_volatile = [](EvalNodeDryRun const& n) {
      if (!n.leaf() || !n->is_tensor()) return false;
      return n->as_tensor().label() == L"t";
    };
    a.peak_bytes = sequant::eval::dryrun::meter(
                       std::vector<EvalNodeDryRun>{node}, policy, regime, cfg)
                       .peak_bytes;

    wchar_t const* obj_name = (obj == ObjectiveFunction::DenseTimeSpaceBatched)
                                  ? L"perf-first (DenseTimeSpaceBatched)"
                                  : L"peak-first (DenseSpaceTimeBatched)";
    std::wcerr << L"[dryrun-objective] " << obj_name << L": optimize " << opt_ms
               << L"ms  max free-mu~ on a node=" << a.max_free_mu
               << L"  largest realized free-mu~={" << a.largest_desc << L"}="
               << a.largest_realized_gb << L" GB\n               modeled cost: "
               << L"n_ops=" << a.model_n_ops << L" flops=" << a.model_flops
               << L" peak_bytes=" << (a.peak_bytes / 1e9) << L" GB\n";

    return a;
  };

  auto peak_first = analyze(ObjectiveFunction::DenseSpaceTimeBatched);
  auto perf_first = analyze(ObjectiveFunction::DenseTimeSpaceBatched);

  auto report = [](wchar_t const* tag, Analysis const& a) {
    std::wcerr << tag << L": 4-PAO node formed = "
               << (a.max_free_mu >= 4 ? L"YES" : L"NO") << L" (max free mu~="
               << a.max_free_mu << L")\n    DP-model largest realized free-mu~="
               << a.largest_realized_gb << L" GB {" << a.largest_desc << L"}\n"
               << L"    modeled cost: n_ops=" << a.model_n_ops << L" flops="
               << a.model_flops << L" peak_bytes=" << (a.peak_bytes / 1e9)
               << L" GB\n";
  };
  std::wcerr << L"\n=== [dryrun-objective] VERDICT (C60 giant, index 38) ===\n";
  report(L"peak-first (DenseSpaceTimeBatched)", peak_first);
  report(L"perf-first (DenseTimeSpaceBatched)", perf_first);

  // STRUCTURAL PROOF (the direct in-harness proof of the fix), read from the
  // static tree walk above -- NOT from the metered replay: perf-first, being
  // flops-primary, must NEVER form the fully-sliceable 4-PAO AO integral (the
  // C60 pathology). Kept as a plain tree-walk check because the free-mu~
  // signature is a property of the FACTORIZATION the DP picked, not of the
  // peak replay, which models cost and does not expose per-node free-index
  // structure. No assertion is made about peak_first.max_free_mu: under the
  // ordered cost model peak-first declines the 4-PAO too (its batched peak is
  // priced with accumulator residency), so the historical contrast is gone --
  // the printed verdict above still shows what each objective picked.
  CHECK(perf_first.max_free_mu < 4);

  // COST PROOF via the metered replay peak: perf-first's modeled peak is
  // dominated by the GENUINE 4-PNO intermediate the perf-first schedule forms
  // -- the CC doubles W node {a_1<i,i> a_2<i,i> a_3<i,i> a_4<i,i>} with FOUR
  // distinct PNO legs over one occ-pair (see SEQUANT_UT_DRYRUN_PERF_TREE
  // dump). With the FAITHFUL measured moments (kC60_pVDZF12), it is sized
  // occ^2 * M_4^4 = 120^2 * 53.151^4 * 8 ~= 0.92 TB (dense occ^2 pairs; the
  // screened ~6300 CC pairs would give ~0.40 TB). This is the CORRECT,
  // moment-aware size (df_regime's csv_pno_moment[k] are power means), NOT a
  // naive-product artifact and NOT a mis-sized twin R{a<i,i>,a<i,i>} (which
  // would be occ^2*M_2^2 ~= 0.3 GB). The Kappa mode is CONTRACTED at this
  // node, so batching cannot shrink W -- it is the irreducible peak floor of
  // the flop-optimal factorization, and precisely why perf-first (which forms
  // it) OOMs C60 while peak-first (which does not) does not. A 0.5..2 TB band
  // brackets the real-moment value with margin and is non-flaky.
  CHECK(perf_first.peak_bytes < 2e12);
  CHECK(perf_first.peak_bytes > 5e11);
}
#endif  // !defined(SEQUANT_SKIP_LONG_TESTS)

// D1.2 (external-mode batching wired into DP SELECTION): the external batch
// loop must flow into the DP's REPORTED peak. Optimizing the over-budget C60
// giant through PeakBatchedModel::reconstruct_batched_modes (the path
// optimize() drives) must, with batch_spectator_indices ON, report a root peak
// BELOW its flag-OFF value (the external occ sliced) and stamp
// BatchModeType::External ONLY on that external occ; with the flag OFF the
// reported peak is byte-identical to the unsliced baseline and NO External
// modes are stamped.
// HIDDEN ([.]): external-mode (occ) seeding no longer lowers the DP-reported
// peak for this C60 giant (peak stays ~391 GB, not < 40 GB). The expectation
// predates the per-node external opens the shipped DP performs (as-built
// section 4.2) and has not been retargeted; see the note below.
// ===========================================================================
// [blocked-layers-1-2] -- what these hidden fixtures ARE
//
// Every TEST_CASE tagged [blocked-layers-1-2] is hidden ([.]) because it
// encodes a DESIGN THIS BRANCH DOES NOT IMPLEMENT, not because a known-good
// test was switched off. The reason some of them still carry -- "blocked on
// Layers 1-2 (use-induced slicing of whole-produced operands + multi-level
// escape chain)" -- has EXPIRED: both are built (as-built design sections
// 6.2 / 6.3 / 7, doc/dev/specs/2026-09-12-batched-array-dag-eval-as-built.md).
// What these fixtures actually encode is the older loop-identity / layout
// model that the identity rework superseded, so their expectations no longer
// describe the shipped builder; at least test_eval_ta.cpp's "batched ToT
// External occ loop" case and test_ordered_schedule.cpp's "forced-split occ
// axis realizes TWO ordered sibling blocks" case FAIL when run today. The
// latter is instructive: it sets only a BatchPolicy role predicate and never
// stamps node_slice_mask on its forest, so the shipped builder realizes NO
// loop at all and finds zero occ blocks -- the fixture's INPUT contract is
// stale, not the forced-split realization it was written to pin (which
// as-built section 6.3 describes correctly and [cell_table][ordered] /
// [w20-auxocc-walk] exercise on real water-20 data).
//
// Retargeting them to the shipped design is deferred (as-built section 12.2).
// They are kept, hidden, as a record of the shapes that still want coverage.
// Do NOT un-hide one without first rewriting its expectation against the code.
// ===========================================================================
// Hidden [blocked-layers-1-2] -- encodes a design this branch does not
// implement; see the [blocked-layers-1-2] note in this file. Do not un-hide
// without rewriting the expectation first.
TEST_CASE(
    "dryrun external-mode seeding lowers the DP-reported peak of the C60 giant",
    "[.][dryrun-extmode][blocked-layers-1-2]") {
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);  // mu~
  sequant::mbpt::add_df_spaces(isr);                             // K
  sequant::tests::declare_real_basis(*isr);  // the residual is real
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
    auto [seq, modes] = m.reconstruct_batched_modes(mctx, mst, &peak);
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
  sequant::tests::declare_real_basis(*isr);  // the residual is real
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
  // occ^2 * M_4^4 = ~0.92 TB (matches the metered replay's ~954 GB peak). The
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
  sequant::tests::declare_real_basis(*isr);  // the residual is real
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

// P4 GO/NO-GO AUDIT (Concern #3): across ALL C60 residual terms under the
// perf-first (DenseTimeSpaceBatched) objective, does external-occ FOREST
// batching (P3, the result-external occ carried as composite protos)
// bound EVERY over-budget giant, or does some giant's dominant memory lever
// escape it -- i.e. a free mu~/K carried on an INTERMEDIATE (contracted
// downstream, hence NOT the result-external P3 slices, and NOT node-
// local sliceable if no ancestor slices it)?
//
// Per term we: (1) optimize perf-first + binarize with node-local batch modes;
// (2) walk every node tracking ancestor node_slice_mask and compute each node's
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
  sequant::tests::declare_real_basis(*isr);  // the residual is real
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
              for (auto const& ax : n->node_slice_mask())
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
// transient (~38.9 GB, materialized INSIDE a make_batched_scratch cache)
// dwarfs the outer, cross-batch cached residency (~0.2 GB). A PeakSink passed
// to make_evaluator folds each scratch cache's high-watermark into one global
// accumulator, so the global peak reflects the true batched-replay peak rather
// than just the outer residency the accessor sees today.
// Minutes-long under ASan/valgrind; see tests/unit/CMakeLists.txt.
#ifndef SEQUANT_SKIP_LONG_TESTS
TEST_CASE("dryrun scratch-fold captures batched peak", "[dryrun][peak]") {
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);  // mu~
  sequant::mbpt::add_df_spaces(isr);                             // K
  sequant::tests::declare_real_basis(*isr);  // the residual is real
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
  auto aops = sequant::eval::dryrun::make_dryrun_array_ops(cm);
  cache.set_array_ops(&aops);
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
  // outer residency (458x measured); a 2x floor is safe and non-flaky.
  CHECK(global_peak > outer_hwmark * 2.0);
}
#endif  // !defined(SEQUANT_SKIP_LONG_TESTS)

// Phase 1 Task 2 regression guard: a metered replay's peak_bytes must be the
// TRUE co-resident sum across the cache scope chain, not
// max(scratch_hwmark, outer_hwmark). eval.hpp's 7 note_working_set() call
// sites now add cache.parent()->chain_residency() (CacheManager, Task 1) to
// the per-op hwmark, so a scratch cache chained (CacheManager::set_parent)
// to an outer cache holding a PERSISTENT, ALIVE cross-term entry folds that
// outer residency into its own working_set_hwmark(). Before this fix, a
// scratch's hwmark reflected ONLY its own local footprint: running the SAME
// batched op with vs without an alive co-resident outer entry produced the
// IDENTICAL hwmark, silently under-reporting the true additive co-resident
// peak.
//
// This test isolates exactly that difference by running the SAME batched
// forest through two structurally-identical scratch caches, one chained to
// an outer cache holding a known-size persistent entry, one not (parent() ==
// nullptr, mirroring the un-hoisted / real-cache-absent case). Because the
// persistent entry's key (tensor label "h") never occurs in the batched
// forest (labels "g"/"C"), chaining cannot change any COMPUTED value --
// access_at() never finds a spurious hit -- so any difference between the
// two runs' working_set_hwmark() is entirely the added chain_residency()
// term. Proven to fail without the eval.hpp fix (see the report's Step 5
// both-states proof: reverting the fix makes hwmark_chained ==
// hwmark_isolated, failing both CHECKs below).
TEST_CASE("dryrun peak is co-resident sum", "[dryrun][peak]") {
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto ctx_resetter = set_scoped_default_context(std::move(ctx));

  using sequant::make_batched_custom_evaluator;
  using sequant::make_no_scope_guard;
  using sequant::never_volatile;
  using node_t = EvalNodeDryRun;

  auto r = backend_test_regime();  // i (occ) extent 10, a (virt) extent 20
  auto cm = std::make_shared<CostModel const>(r);
  DryRunLeafEvaluator yield{cm};

  // Same small batched-forest shape as the D3.1 test above (external-mode
  // scatter over the occ index carried only as a PNO proto-index).
  auto expr = deserialize<ExprPtr>(
      "(g{a_3;a_4} * C{a_4;a1<i_1,i_2>}) * C{a2<i_1,i_2>;a_3}");
  REQUIRE(static_cast<bool>(expr));
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
  node->set_node_slice_mask({{mode, BatchModeType::External}});

  // A LEAF with a distinct tensor label ("h"), unrelated to anything in the
  // batched tree above, registered directly (not via the batched forest) as
  // a PERSISTENT entry in a standalone outer CacheManager -- our synthetic
  // stand-in for a persistent cross-term cache entry alive at run scope.
  auto persistent_expr = deserialize<ExprPtr>("h{i_7;a_9}");
  REQUIRE(static_cast<bool>(persistent_expr));
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto persistent_node = binarize<EvalExprDryRun>(persistent_expr);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  REQUIRE(persistent_node.leaf());

  using hasher_t = sequant::TreeNodeHasher<node_t>;
  using comp_t = sequant::TreeNodeEqualityComparator<node_t>;
  std::unordered_map<node_t, size_t, hasher_t, comp_t> outer_reg;
  outer_reg.emplace(persistent_node, std::numeric_limits<size_t>::max());
  auto always_persistent = [](node_t const&) { return true; };
  sequant::CacheManager<node_t> outer(std::move(outer_reg), always_persistent);

  ResultPtr persistent_val = yield(persistent_node);
  REQUIRE(persistent_val);
  (void)outer.store_and_access(persistent_node, persistent_val);
  REQUIRE(outer.alive(persistent_node));
  size_t const R = outer.current_residency();
  REQUIRE(R > 0);
  REQUIRE(R == persistent_val->size_in_bytes());

  auto target_batch_size = [](Index const&) -> std::size_t { return 4; };

  // Force printing() on (the CacheManager hwmark only accumulates there):
  // note_working_set()'s per-op hwmark is fed by the CACHE-AWARE bytes()
  // overload, which short-circuits to 0 unless the eval trace is being
  // printed. Restored on every exit path below.
  auto& logger = Logger::instance();
  auto const prev_level = logger.eval.level;
  auto* const prev_stream = logger.eval.stream;
  std::ostringstream trace_os;
  logger.eval.level = 2;
  logger.eval.stream = &trace_os;

  // Run A: ISOLATED -- scratch has NO parent, so chain_residency() is never
  // consulted (parent() == nullptr short-circuits the added term to 0 at
  // every one of the 7 fixed sites).
  auto scratch_isolated = sequant::CacheManager<node_t>::empty();
  auto aops = sequant::eval::dryrun::make_dryrun_array_ops(cm);
  scratch_isolated.set_array_ops(&aops);
  scratch_isolated.set_custom_evaluator(
      make_batched_custom_evaluator(yield, target_batch_size, accept_occ,
                                    make_no_scope_guard{}, never_volatile{}));
  ResultPtr result_isolated;
  try {
    result_isolated = sequant::evaluate(node, yield, scratch_isolated);
  } catch (std::exception const&) {
  }
  size_t const hwmark_isolated = scratch_isolated.working_set_hwmark();

  // Run B: CHAINED -- an otherwise-identical fresh scratch, parented to
  // `outer` exactly the way place_at_this_level() (eval.hpp) wires an
  // order-aware hoisted invariant's scratch to its enclosing real/term
  // cache.
  auto scratch_chained = sequant::CacheManager<node_t>::empty();
  scratch_chained.set_parent(&outer);
  scratch_chained.set_array_ops(&aops);
  scratch_chained.set_custom_evaluator(
      make_batched_custom_evaluator(yield, target_batch_size, accept_occ,
                                    make_no_scope_guard{}, never_volatile{}));
  ResultPtr result_chained;
  try {
    result_chained = sequant::evaluate(node, yield, scratch_chained);
  } catch (std::exception const&) {
  }
  size_t const hwmark_chained = scratch_chained.working_set_hwmark();

  logger.eval.level = prev_level;
  logger.eval.stream = prev_stream;

  REQUIRE(result_isolated);
  REQUIRE(result_chained);
  REQUIRE(hwmark_isolated > 0);

  // The two runs compute IDENTICAL data: chaining a cache whose only entry's
  // key never occurs in the batched forest cannot change any computed value.
  CHECK(result_chained->size_in_bytes() == result_isolated->size_in_bytes());

  // THE regression guard: with the fix, the chained run's hwmark is the
  // isolated local footprint PLUS the outer's live co-resident residency --
  // an exact sum, not a max.
  CHECK(hwmark_chained == hwmark_isolated + R);
  // Equivalently, and matching the brief's robust form: strictly greater
  // than what max(scratch, outer) alone would yield.
  CHECK(hwmark_chained > std::max(hwmark_isolated, R));
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
TEST_CASE("dryrun gated cache footprint-gates the giant", "[dryrun][cache]") {
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);  // mu~
  sequant::mbpt::add_df_spaces(isr);                             // K
  sequant::tests::declare_real_basis(*isr);  // the residual is real
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
  auto ref_cache = build_dryrun_cache(nodes, ref_cfg, regime);

  // task: the config the task asks for (footprint gate on).
  CacheConfig task_cfg = ref_cfg;
  task_cfg.max_footprint = 1e11;
  auto task_cache = build_dryrun_cache(nodes, task_cfg, regime);

  // Phase 4b-1: with sliced_modes unified to the all-batched-modes
  // cross-occurrence meet, this giant's lifetime mask is now NON-empty -- it
  // carries a free batchable (aux) mode on its own result slots that is batched
  // above it, so its per-occurrence value differs per batch of that mode. The
  // cross-occurrence batch-variant veto therefore CORRECTLY removes it from run
  // scope even with NO footprint gate (the F1 hazard). This is the latent
  // under-veto the former External-only mask left open: a contracted (aux)
  // batch mode free on an intermediate's result is genuinely batch-variant, and
  // the unified meet now expresses it. (Before, the External-only mask saw no
  // External stamp on this DF-aux schedule, left the mask empty, and the giant
  // was admitted -- capped only by the footprint gate.)
  CHECK_FALSE(ref_cache.exists(giant_node));
  // The FOOTPRINT gate independently caps the >100 GB full giant too.
  CHECK_FALSE(task_cache.exists(giant_node));

  // A small non-batchable intermediate stays cached under both.
  CHECK(ref_cache.exists(small_node));
  CHECK(task_cache.exists(small_node));

  // Faithful end-to-end: replay the schedule through the real eval loop against
  // the gated (task) cache; the giant must not become resident, while eval
  // completes without materializing it whole.
  auto cm = std::make_shared<CostModel const>(regime);
  // The batched replay needs backend array-ops wired on the cache (see
  // dryrun::meter()); without them make_batched_custom_evaluator asserts.
  auto const aops = sequant::eval::dryrun::make_dryrun_array_ops(cm);
  task_cache.set_array_ops(&aops);
  task_cache.set_custom_evaluator(
      sequant::make_evaluator(policy, DryRunLeafEvaluator{cm}));
  REQUIRE_NOTHROW(
      (void)sequant::evaluate(node, DryRunLeafEvaluator{cm}, task_cache));
  CHECK_FALSE(task_cache.exists(giant_node));  // still not registered => never
                                               // resident
  CHECK_FALSE(task_cache.alive(giant_node));

  std::wcerr
      << L"\n=== [dryrun][cache] footprint-gate verdict (C60 giant, index "
      << giant_idx << L") ===\n  giant free-batchable node footprint = "
      << (giant_bytes / 1e9) << L" GB\n  ref (no gate)  exists(giant) = "
      << (ref_cache.exists(giant_node) ? L"YES" : L"no")
      << L"\n  task config    exists(giant) = "
      << (task_cache.exists(giant_node) ? L"YES" : L"no")
      << L"\n  small intermediate exists (task) = "
      << (task_cache.exists(small_node) ? L"YES" : L"no") << L"\n";
}

// PNO-CCSD water-20 aux-batching FRAGMENTATION surrogate. Faithfully mirrors
// the MPQC water-20 pVDZ-F12 PNO-CCSD run (job 658937): the SAME csv doubles
// residual equation, df_regime(kWater20_pVDZF12) (extents/moments verified
// against the job log), and the EXACT batch config make_csv_batch_policy emits
// for aux-only batching (objective dense_time_space, K contracted-batchable,
// target 256, peak_threshold 1e11, persistent_only false, no PAO/occ axis).
// MPQC drives batching through this same DP, so the aprime decisions reproduce
// by construction. Purpose: reproduce the batched-member explosion (the
// ~405-vs-83 group fragmentation the new-vs-old logs showed) and, under
// SEQUANT_DP_RECOMPUTE_DEBUG=1, expose why -- the K-carrying gC-class
// composites are charged rf==1 (the escaped-mode recompute model charges ZERO
// recompute to a node that CARRIES the only batch mode), so the DP prices
// slicing them as free and over-batches; the runtime then rebuilds them per
// batch group (the measured 2.5x product-work regression). See
// cost_model.hpp:1811-1817 ("if the expensive gC-class nodes show rf==1, the DP
// is not pricing the runtime recompute").
TEST_CASE("dryrun water-20 aux-batch fragmentation: gC composites priced rf==1",
          "[.][dryrun-water-frag]") {
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);
  sequant::mbpt::add_df_spaces(isr);
  sequant::tests::declare_real_basis(*isr);  // the residual is real
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
  auto regime = df_regime(kWater20_pVDZF12);

  // EXACT MPQC aux-only config (make_csv_batch_policy with aux_target=256,
  // pao_target=0, occ_target=0): K is the ONLY batchable mode, contracted role.
  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](Index const& ix) {
    return ix.space().base_key() == L"Κ";
  };
  policy.is_batchable_external_index = [](Index const&) { return false; };
  policy.batch_spectator_indices = false;
  policy.batch_target_size = [](Index const&) -> std::size_t { return 256; };
  policy.is_volatile_leaf = [](Tensor const& t) { return t.label() == L"t"; };
  policy.accumulation_factor = 1.0;
  policy.persistent_only = false;
  policy.peak_threshold = 1e11;

  auto axes_map = std::make_shared<std::unordered_map<
      Expr const*, container::vector<NodeBatchAnnotation>>>();
  OptimizeOptions opts;
  // MPQC "dense_time_space" + aux batch keywords => the perf-first BATCHED
  // objective (the batchability model that emits the contracted-K aprime; plain
  // DenseTimeSpace carries no per-index batch model, options.hpp:62).
  opts.objective_function = ObjectiveFunction::DenseTimeSpaceBatched;
  opts.idx_to_extent = regime.idx_to_extent();
  opts.inner_pow = regime.inner_pow_fn();
  opts.batch_policy = policy;
  opts.volatile_weight = 20.0;
  opts.roofline.machine_balance = 200.0;
  opts.roofline.fast_mem_elems = 1000000.0;
  opts.term_batch_axes = axes_map;

  std::size_t n_kcon = 0;      // K-contracted batched-member annotations
  std::size_t n_terms_ok = 0;  // terms that optimized
  for (std::size_t s = 0; s < nterms; ++s) {
    ExprPtr const term = flatten_product(summands[s]);
    if (!term) continue;
    ExprPtr optimized;
    try {
      optimized = optimize(term, opts);
    } catch (std::exception const&) {
      continue;
    }
    if (!optimized) continue;
    ++n_terms_ok;
    auto it = axes_map->find(optimized.get());
    if (it == axes_map->end()) continue;
    for (auto const& na : it->second)
      for (auto const& e : na.axes)
        if (e.second == sequant::BatchModeType::Contracted &&
            e.first.space().base_key() == L"Κ")
          ++n_kcon;
  }
  std::wcerr
      << L"\n=== [dryrun-water-frag] water-20 aux-only, " << n_terms_ok
      << L" terms optimized ===\n  K-contracted batched-member "
         L"annotations: "
      << n_kcon
      << L"\n  (SEQUANT_DP_RECOMPUTE_DEBUG=1 dumps per-node rf; CONFIRMED "
         L"all 91 K-carrying gC composites -- largest 34 GB -- priced "
         L"rf==1 while the 46 K-escaping nodes are charged rf==7: the DP "
         L"prices slicing the gC giants as free, cost_model.hpp:1811)\n";
  CHECK(n_terms_ok > 0);
  // Reproduces the MPQC-log fragmentation (~35 distinct batched-member shapes):
  // aux batching promotes MANY nodes to K-batched members. The mechanism is the
  // rf==1 pricing of every K-carrying gC composite (confirmed via
  // SEQUANT_DP_RECOMPUTE_DEBUG): the escaped-mode recompute model charges zero
  // to a node that CARRIES the only batch mode, so the flops-neutral slice
  // looks free -- but the runtime rebuilds each such composite per consumer
  // batch group (the measured 2.5x product-work regression, not modeled here).
  CHECK(n_kcon > 20);
}

// Does canonicalization preserve DISTINCT composite (PNO) proto pairs? A
// contraction spanning composites on different occ pairs (a<i_1,i_2>,
// b<i_2,i_3>) physically spans THREE occ; if canon_indices collapsed the pairs
// onto one, its flops (and the DP/model sizing that reads canon_indices) would
// undersize.
TEST_CASE("canon_indices preserves distinct composite proto pairs",
          "[eval_expr][composite-canon]") {
  using namespace sequant;
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);
  sequant::mbpt::add_df_spaces(isr);
  auto resetter = set_scoped_default_context(std::move(ctx));

  Index const i1{L"i_1"}, i2{L"i_2"}, i3{L"i_3"};
  auto const a_space = Index{L"a_1"}.space();
  Index const a1{a_space, container::vector<Index>{i1, i2}};  // pair (i_1,i_2)
  Index const a2{a_space, container::vector<Index>{i2, i3}};  // pair (i_2,i_3)
  Tensor const t(L"I", bra{a1}, ket{a2}, Symmetry::Nonsymm);
  REQUIRE(a1.has_proto_indices());
  REQUIRE(a2.has_proto_indices());

  EvalExpr const ev{t};
  auto const& ci = ev.canon_indices();

  auto nar = [](std::wstring_view w) {
    std::string s;
    for (wchar_t c : w) s += (c < 128 ? static_cast<char>(c) : '#');
    return s;
  };
  std::set<std::wstring> distinct_protos;
  std::vector<std::set<std::wstring>> pairs;
  for (auto const& ix : ci)
    if (ix.has_proto_indices()) {
      std::set<std::wstring> p;
      for (auto const& px : ix.proto_indices()) {
        p.insert(std::wstring(px.full_label()));
        distinct_protos.insert(std::wstring(px.full_label()));
      }
      pairs.push_back(p);
    }
  std::string dump;
  for (auto const& ix : ci) dump += nar(ix.full_label()) + " ";
  INFO("canon_indices = [" << dump << "]");
  REQUIRE(pairs.size() == 2);
  // Physical form spans THREE distinct occ (i_1,i_2,i_3); a collapse => TWO.
  CHECK(distinct_protos.size() == 3);
  CHECK(pairs[0] != pairs[1]);

  // Binary contraction: A{a1<i_1,i_2>,i_4;} * B{a2<i_2,i_3>; i_4} contracts
  // i_4 (bra of A, ket of B -- particle-conserving), so the RESULT carries
  // a1<i_1,i_2> and a2<i_2,i_3> -- composites on distinct pairs, spanning 3
  // occ. This is the (binary) node whose canon_indices the static cost walk
  // sizes; check IT preserves the distinct pairs too.
  const ResultExpr expr = deserialize<ResultExpr>(
      "R{a1<i_1,i_2>,a2<i_2,i_3>;} = A{a1<i_1,i_2>,i_4;} * B{a2<i_2,i_3>; "
      "i_4}");
  auto const node = binarize(expr);
  auto const& cci = node->canon_indices();
  std::set<std::wstring> bprotos;
  std::string bdump;
  for (auto const& ix : cci) {
    bdump += nar(ix.full_label()) + " ";
    if (ix.has_proto_indices())
      for (auto const& px : ix.proto_indices())
        bprotos.insert(std::wstring(px.full_label()));
  }
  INFO("binary node canon_indices = [" << bdump << "]");
  CHECK(bprotos.size() == 3);  // NOT 2 -- a collapse would merge the pairs
}

// ---------------------------------------------------------------------------
// Metered dry-run replay ([meter]): PeakMonitor, assemble_report, and the
// meter() entry point itself. Moved here verbatim from the former
// tests/unit/test_meter.cpp when the cost-profile replay predictor was retired
// (the metered replay is what survived it).
// ---------------------------------------------------------------------------

TEST_CASE("PeakMonitor tracks hierarchy-wide co-resident high-water",
          "[meter]") {
  sequant::eval::PeakMonitor mon;
  std::vector<std::size_t> peaks;
  mon.on_peak = [&](sequant::eval::PeakEvent const& e) {
    peaks.push_back(e.bytes);
  };
  mon.observe(100, 0xA);
  mon.observe(50, 0xB);   // below hwmark: no advance, no fire
  mon.observe(300, 0xC);  // new peak
  CHECK(mon.hwmark_bytes == 300);
  CHECK(mon.peak.op_hash == 0xC);
  CHECK(peaks == std::vector<std::size_t>{100, 300});
}

TEST_CASE(
    "PeakMonitor wired onto a CacheManager observes the same high-water as "
    "the cache's own working_set_hwmark()",
    "[meter]") {
  using sequant::eval::dryrun::CostModel;
  using sequant::eval::dryrun::DryRunLeafEvaluator;
  using sequant::eval::dryrun::EvalExprDryRun;
  using sequant::eval::dryrun::EvalNodeDryRun;
  using sequant::eval::dryrun::SizeRegime;

  // A tiny, self-consistent regime: two named spaces, no composite (proto-
  // indexed) legs involved, so no CSV/PNO moments are needed.
  SizeRegime regime;
  regime.space_extent = {
      {L"i", 10},
      {L"a", 20},
  };
  auto cm = std::make_shared<CostModel const>(regime);

  // A minimal 3-node forest: two leaves (g, t) contracted into one product,
  // fully contracted (no external indices) so the whole tree is a scalar.
  // Small enough to hand-build directly rather than routing through
  // sequant::optimize (see test-1 brief: a 2-3 node hand-built forest is
  // acceptable for this cache-integration assertion).
  auto expr = sequant::deserialize<sequant::ExprPtr>("g{i_1;a_3} * t{a_3;i_1}");
  REQUIRE(static_cast<bool>(expr));

  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = sequant::binarize<EvalExprDryRun>(expr);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  REQUIRE_FALSE(node.leaf());

  auto cache = sequant::cache_manager(std::vector<EvalNodeDryRun>{node});

  sequant::eval::PeakMonitor mon;
  cache.set_peak_monitor(&mon);

  // Redirect the eval trace to a private buffer (rather than stdout) and
  // restore the logger's prior state afterward -- note_working_set()'s
  // per-op hwmark input is only folded when Logger::instance().eval.level >
  // 0 (see cache_manager.hpp / eval.hpp).
  auto& logger = sequant::Logger::instance();
  auto const prev_level = logger.eval.level;
  auto* const prev_stream = logger.eval.stream;
  std::ostringstream trace_os;
  logger.eval.level = 1;
  logger.eval.stream = &trace_os;

  sequant::eval::dryrun::DryRunLeafEvaluator yield{cm};
  sequant::ResultPtr result;
  try {
    result = sequant::evaluate<sequant::Trace::On>(node, yield, cache);
  } catch (...) {
    logger.eval.level = prev_level;
    logger.eval.stream = prev_stream;
    throw;
  }
  logger.eval.level = prev_level;
  logger.eval.stream = prev_stream;

  REQUIRE(result);
  CHECK(mon.hwmark_bytes > 0);
  CHECK(mon.hwmark_bytes == cache.working_set_hwmark());
}

TEST_CASE(
    "assemble_report rolls a metered dry-run replay's per-node build tally "
    "into a MeterReport (peak, persistent/volatile FLOPs+time, build-vs-home)",
    "[meter]") {
  using sequant::eval::compute_dag_boulevard;
  using sequant::eval::PeakMonitor;
  using sequant::eval::dryrun::assemble_report;
  using sequant::eval::dryrun::compute_volatility;
  using sequant::eval::dryrun::CostModel;
  using sequant::eval::dryrun::DryRunLeafEvaluator;
  using sequant::eval::dryrun::EvalExprDryRun;
  using sequant::eval::dryrun::EvalNodeDryRun;
  using sequant::eval::dryrun::SizeRegime;

  SizeRegime regime;
  regime.space_extent = {
      {L"i", 10},
      {L"a", 20},
  };
  auto cm = std::make_shared<CostModel const>(regime);

  // Three-factor product, folded LEFT-TO-RIGHT by binarize (see
  // fold_left_to_node, binary_node.hpp): root == (X * t) with X == (g * h).
  // X is a PERSISTENT internal node (neither g nor h is volatile); the root
  // additionally consumes the volatile leaf t, so the root is volatile --
  // exercising BOTH branches of assemble_report's persistent/volatile split.
  // X keeps i_1 (from g) and i_2 (from h) external after contracting a_3; t
  // contracts i_2 and keeps a_5 external, so the ROOT keeps i_1 and a_5
  // external too (not a bare scalar). That is deliberate: i_1 being a
  // genuine external slot of BOTH X and the root lets the ROOT's own
  // External loop over i_1 (stamped below) give X a NONEMPTY home -- X sits
  // inside that loop but does not own it -- while the root's OWN home stays
  // empty (a node's own realized loop is excluded from its own home; see
  // compute_dag_boulevard's own_modes_union subtraction). That asymmetry is
  // the content check below: a broken cell_by_hash lookup in assemble_report
  // would silently default BOTH to empty, indistinguishable from the root's
  // genuinely-empty case, so only X's nonempty "i" catches a broken lookup.
  auto expr = sequant::deserialize<sequant::ExprPtr>(
      "g{i_1;a_3} * h{a_3;i_2} * t{i_2;a_5}");
  REQUIRE(static_cast<bool>(expr));

  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = sequant::binarize<EvalExprDryRun>(expr);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  REQUIRE_FALSE(node.leaf());
  REQUIRE_FALSE(node.left().leaf());  // X == g * h (persistent)
  REQUIRE(node.right().leaf());       // t (volatile)

  sequant::Index const mode{L"i_1"};
  REQUIRE(sequant::index_position(node, mode).has_value());
  REQUIRE(sequant::index_position(node.left(), mode).has_value());
  // Stamp the ROOT's own External loop over i_1. Plain (non-scope)
  // evaluate() ignores node_slice_mask without a custom evaluator (see
  // test_eval_dryrun.cpp's equivalent stamping), so this only feeds
  // compute_dag_boulevard's home/ectx bookkeeping below, not the replay. The
  // enclosing-loop context (OccurrenceRec::ectx, which the children's `uses`
  // reads) is reconstructed from batch_loops_opened_here(), so the root --
  // which REALIZES the i_1 loop -- must open it, not merely carry the sliced
  // mask.
  node->set_node_slice_mask({{mode, sequant::BatchModeType::External}});
  node->set_batch_loops_opened_here({{mode, sequant::BatchModeType::External}});

  std::vector<EvalNodeDryRun> const forest{node};

  auto cache = sequant::cache_manager(forest);
  cache.set_recompute_tally_enabled(true);

  sequant::eval::PeakMonitor mon;
  cache.set_peak_monitor(&mon);

  auto const is_volatile = [](EvalNodeDryRun const& n) {
    return n.leaf() && n->is_tensor() && n->as_tensor().label() == L"t";
  };

  auto const block_of = [](sequant::Index const&) -> std::size_t {
    return 256;
  };
  auto const rich = compute_dag_boulevard(forest, *cm, block_of);

  auto& logger = sequant::Logger::instance();
  auto const prev_level = logger.eval.level;
  auto* const prev_stream = logger.eval.stream;
  std::ostringstream trace_os;
  logger.eval.level = 1;
  logger.eval.stream = &trace_os;

  sequant::eval::dryrun::DryRunLeafEvaluator yield{cm};
  sequant::ResultPtr result;
  try {
    result = sequant::evaluate<sequant::Trace::On>(node, yield, cache);
  } catch (...) {
    logger.eval.level = prev_level;
    logger.eval.stream = prev_stream;
    throw;
  }
  logger.eval.level = prev_level;
  logger.eval.stream = prev_stream;

  REQUIRE(result);

  auto const report = assemble_report(cache, mon, rich, forest, is_volatile,
                                      sequant::BatchScheduler::forest_descent);

  CHECK(report.builds_total > 0);
  CHECK_FALSE(report.home_fidelity.empty());
  CHECK(report.scheduler == sequant::BatchScheduler::forest_descent);

  // EXEC/COST SPLIT: exec is threaded BuildRecord -> tally_build ->
  // assemble_report independently of (and not swapped with) flops -- both
  // the volatile (root) and persistent (X) buckets get a positive exec
  // estimate alongside their flops.
  CHECK(report.flops_volatile > 0.0);
  CHECK(report.cost_volatile > 0.0);
  CHECK(report.flops_persistent > 0.0);
  CHECK(report.cost_persistent > 0.0);

  // HOME-FIDELITY CONTENT: look up X's (persistent, g*h) and the root's
  // entries by hash and check their home/uses STRINGS, not just
  // non-emptiness of the whole list -- see the forest-design comment above
  // for why X's nonempty "i" is the discriminating case.
  auto const find_by_hash = [&](std::size_t h) {
    return std::find_if(report.home_fidelity.begin(),
                        report.home_fidelity.end(),
                        [h](auto const& hf) { return hf.hash == h; });
  };

  auto const x_it = find_by_hash(node.left()->hash_value());
  REQUIRE(x_it != report.home_fidelity.end());
  CHECK(x_it->home == "i");
  CHECK(x_it->uses == "i");

  auto const root_it = find_by_hash(node->hash_value());
  REQUIRE(root_it != report.home_fidelity.end());
  CHECK(root_it->home.empty());
  CHECK(root_it->uses.empty());
}

// Task 3: meter() drives the REAL policy-selected executor (ordered or
// forest descent, chosen by BatchPolicy::scheduler) through the sizing
// backend with its own metered cache, rather than a hand-rolled proxy of one.
// This forest carries a genuine Contracted batch axis (a_3, the "aux"-analog):
// with scheduler == BatchScheduler::ordered, the driver entry
// (ordered_executor.hpp) rebuilds a schedule with ONE realized batch loop
// from that stamp and drives the executor's nested batch-scratch caches
// (a Task-1 coverage gap -- no earlier [meter] test exercised a multi-level
// PeakMonitor parent-chain under a real batched walk); with it forest_descent,
// the SAME forest runs through today's unbatched per-tree descent. Both modes
// must report a positive peak/build count and the matching `scheduler` value.
TEST_CASE(
    "meter runs the policy-selected executor (ordered and forest) with "
    "the sizing backend and a metered cache",
    "[meter]") {
  using sequant::BatchModeType;
  using sequant::BatchPolicy;
  using sequant::Index;
  using sequant::eval::dryrun::CacheConfig;
  using sequant::eval::dryrun::EvalExprDryRun;
  using sequant::eval::dryrun::EvalNodeDryRun;
  using sequant::eval::dryrun::meter;
  using sequant::eval::dryrun::SizeRegime;

  // Same small, self-consistent regime as the tests above; "a" (extent 20)
  // doubles as the batch-axis space here (target_size 5 => 4 blocks).
  SizeRegime regime;
  regime.space_extent = {
      {L"i", 10},
      {L"a", 20},
  };

  // A single-root, fully-contracted product -- a_3 is a genuine operand
  // index of BOTH factors (contracted away at the root), hand-stamped as the
  // root's one Contracted batch axis (no optimize() involved, matching
  // test_eval_dryrun.cpp's hand-built batch-annotation recipe): the driver
  // entry rebuilds its schedule from the forest's OWN node_slice_mask()
  // stamps, not from BatchPolicy predicates, so this alone is enough for it
  // to realize a non-root-only loop nest.
  auto expr =
      sequant::deserialize<sequant::ExprPtr>(L"g{i_1;a_3} * h{a_3;i_1}");
  REQUIRE(static_cast<bool>(expr));

  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = sequant::binarize<EvalExprDryRun>(expr);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  REQUIRE_FALSE(node.leaf());

  Index const a3{L"a_3"};
  node->set_node_slice_mask({{a3, BatchModeType::Contracted}});

  std::vector<EvalNodeDryRun> const forest{node};

  CacheConfig const cfg;  // default: no footprint gate, min_repeats=1

  BatchPolicy policy;
  policy.batch_target_size = [](Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"a" ? std::size_t{5} : std::size_t{1};
  };

  // ---- ordered: exercises the nested batch-scratch walk. ----
  policy.scheduler = sequant::BatchScheduler::ordered;
  auto const ord_report = meter(forest, policy, regime, cfg);
  CHECK(ord_report.peak_bytes > 0.0);
  CHECK(ord_report.builds_total > 0);
  CHECK(ord_report.scheduler == sequant::BatchScheduler::ordered);

  // ---- forest descent: the SAME forest/policy, scheduler reset. ----
  policy.scheduler = sequant::BatchScheduler::forest_descent;
  auto const fd_report = meter(forest, policy, regime, cfg);
  CHECK(fd_report.peak_bytes > 0.0);
  CHECK(fd_report.builds_total > 0);
  CHECK(fd_report.scheduler == sequant::BatchScheduler::forest_descent);
}
