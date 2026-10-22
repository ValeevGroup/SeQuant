//
// Named CSV (PNO/OSV) problem sizes for the dry-run cost model, shared by the
// eval/dry-run/ordered-executor/legality/schedule-document unit tests, plus
// the occupancy adjustment MPQC applies to the moment tables it hands SeQuant.
//
// Lives in a named namespace (not an anonymous one) so it can be included by
// several test files under a unity build without ODR collisions; the per-file
// aliases keep each file's historical identifiers.
//

#ifndef SEQUANT_TESTS_UNIT_CSV_SIZE_FIXTURES_HPP
#define SEQUANT_TESTS_UNIT_CSV_SIZE_FIXTURES_HPP

#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>

#include <array>
#include <cmath>
#include <cstddef>

namespace sequant::tests::csv_sizes {

/// One named (molecule, basis, parameter-set) problem size for the dry-run cost
/// model: base-space extents, the measured heavy-tailed CSV power means, and
/// the per-rank cluster-grid occupancy.
///
/// `pno_M[k]` / `osv_M[k]` (k in [1, kMaxCompositeGroup]; index 0 unused,
/// kept 1) are the
/// per-NONNULL-cluster power means `M_k = (mean over nonnull clusters of
/// d^k)^(1/k)` of the per-pair PNO / per-orbital OSV domain size d, exactly as
/// mpqc's PaoPnoRMP2 prints them ("PNO domain power means M_1..M_4 per pair",
/// "OSV domain power means M_1..M_4 per orbital"). `pno_occupancy` /
/// `osv_occupancy` are the fractions of the dense `nocc^2` / `nocc` cluster
/// grids whose clusters are genuinely nonnull (element level: a nested TA tile
/// stores nothing for a null inner cell).
struct CsvProblemSize {
  std::size_t mu_tilde;             // PAO extent (the DP's extent of mu~ and a)
  std::size_t aux;                  // DF aux (K) extent
  std::size_t i_occ;                // active occupied extent
  eval::dryrun::MomentTable pno_M;  // per-nonnull-pair PNO power means M_0..M_K
  eval::dryrun::MomentTable osv_M;  // per-nonnull-orbital OSV power means
  double pno_occupancy;             // nonnull ordered pairs / i_occ^2
  double osv_occupancy;  // orbitals with a nonnull OSV domain / i_occ
};

/// M_5..M_8 were not printed by the runs these fixtures were taken from (the
/// producer measured M_1..M_4 at the time). They are filled FLAT with M_4 --
/// a lower bound for a heavy-tailed domain (M_k is non-decreasing in k), i.e.
/// exactly what the former column-4 clamp priced -- so a group of five to
/// eight composites on one cluster is sized as the clamp used to size it
/// rather than rejected. Replace with measured values when a run prints
/// them; a test that depends on a k >= 5 group should say so.
inline constexpr double kUnmeasured_w20_pno_M4 = 30.03848680550367;
inline constexpr double kUnmeasured_w20_osv_M4 = 59.872014818179686;
inline constexpr double kUnmeasured_c60_pno_M4 = 53.151291880343109;
inline constexpr double kUnmeasured_c60_osv_M4 = 166.85553430303926;

/// (H2O)20 / cc-pVDZ-F12 / aug-cc-pVDZ-RI, frozen core (80 active occupied),
/// PNO-CCSD with tIJ = 3e-5. Three provenances, kept apart on purpose:
/// - `i_occ`: the active occupied count (exact).
/// - `mu_tilde`, `aux`: EFFECTIVE extents backed out of realized leaf sizes in
///   Owl job 649160 -- ext(K) = 1682 from g(i,i,K) = 86.1 MB, ext(mu~) = 896
///   from g(mu~,mu~,K) = 10.8 GB. They are NOT the basis-function counts (960
///   PAOs, 2360 aux functions for this system): the difference is the block
///   sparsity of those two leaves, folded into a dense extent so the dense
///   cost model reproduces the realized leaf sizes.
/// - `pno_M`, `osv_M`, `pno_occupancy`: the z820 water-20 runs at tIJ = 3e-5
///   (2026-09-14/15, `w20-task15-occ{8,16}-tTA1e-10` logs), which print these
///   exact moments and "Number of {all,unique} CSV pairs: {1354,717}" ->
///   1354 ordered pairs of 80^2 = 0.2115625. Every active orbital has OSVs.
inline constexpr CsvProblemSize kWater20_pVDZF12{
    /*mu_tilde=*/896u,
    /*aux=*/1682u,
    /*i_occ=*/80u,
    /*pno_M=*/
    {1.0, 23.175775480059084, 25.865548281212597, 28.171416142614103,
     30.03848680550367, kUnmeasured_w20_pno_M4, kUnmeasured_w20_pno_M4,
     kUnmeasured_w20_pno_M4, kUnmeasured_w20_pno_M4},
    /*osv_M=*/
    {1.0, 58.987499999999997, 59.289227520688783, 59.584437469011633,
     59.872014818179686, kUnmeasured_w20_osv_M4, kUnmeasured_w20_osv_M4,
     kUnmeasured_w20_osv_M4, kUnmeasured_w20_osv_M4},
    /*pno_occupancy=*/1354.0 / 6400.0,
    /*osv_occupancy=*/1.0};

/// C60 / cc-pVDZ-F12, TCUTPNO = 1e-8 / TCUTOSV = 1e-9, 120 active occupied.
/// - extents: the Owl job log (1800 PAOs: no PAO-domain truncation on the
///   delocalized cage; 4320 aux functions).
/// - `pno_M`, `osv_M`: measured by mpqc PaoPnoRMP2 in Owl job 617809 (the
///   true pVDZ-F12 values; the earlier 48.60/64.35 PNO and 206/234 OSV
///   numbers came from job 617653, mis-configured with an aug-cc-pVTZ OBS).
/// - `pno_occupancy`: ~6300 significant element pairs of 120^2 = 14400 as
///   logged by Owl job 606308 and recorded in the C60 PNO-CCSD OOM analysis;
///   a logged count, rounded to the hundreds, so treat 0.4375 as +-1%.
inline constexpr CsvProblemSize kC60_pVDZF12{
    /*mu_tilde=*/1800u,
    /*aux=*/4320u,
    /*i_occ=*/120u,
    /*pno_M=*/
    {1.0, 42.029069767441861, 46.039206412923569, 49.766252354482994,
     53.151291880343109, kUnmeasured_c60_pno_M4, kUnmeasured_c60_pno_M4,
     kUnmeasured_c60_pno_M4, kUnmeasured_c60_pno_M4},
    /*osv_M=*/
    {1.0, 148.25, 155.04434849422921, 161.33527408797721, 166.85553430303926,
     kUnmeasured_c60_osv_M4, kUnmeasured_c60_osv_M4, kUnmeasured_c60_osv_M4,
     kUnmeasured_c60_osv_M4},
    /*pno_occupancy=*/6300.0 / 14400.0,
    /*osv_occupancy=*/1.0};

/// The occupancy adjustment mpqc applies before handing the moment tables to
/// SeQuant (mpqc `CSVSizeMoments::extent(rank, k)`): each member of a k-group
/// of composites sharing one cluster is sized by `p^(1/k) * M_k`, so the
/// group's per-member product (`inner_aware_volume`) comes out as
/// `p * M_k^k`, and `nocc^rank * p * M_k^k = sum over nonnull clusters of d^k`
/// is the block-sparse volume. This is the contract the fixtures share with
/// mpqc; `[dryrun][sizing]` checks the group-volume identity.
///
/// \param M per-nonnull-cluster power means, index 0 unused
/// \param p the cluster grid's occupancy, in (0, 1]
/// \return the adjusted table, index 0 untouched
inline eval::dryrun::MomentTable occupancy_adjusted_moments(
    eval::dryrun::MomentTable M, double p) {
  for (std::size_t k = 1; k <= eval::dryrun::kMaxCompositeGroup; ++k)
    M[k] *= std::pow(p, 1.0 / static_cast<double>(k));
  return M;
}

/// The dry-run regime for the post-transform (mu~ + K) residual of \p p:
/// extents for i, mu~, K and a (the latter sized as mu~, the PAO basis a CSV
/// composite is expanded in), and the occupancy-adjusted moment tables.
inline eval::dryrun::SizeRegime df_size_regime(CsvProblemSize const& p) {
  eval::dryrun::SizeRegime r;
  r.space_extent = {
      {L"i", p.i_occ},
      {L"μ̃", p.mu_tilde},
      {L"Κ", p.aux},
      {L"a", p.mu_tilde},
  };
  r.csv_pno_moment = occupancy_adjusted_moments(p.pno_M, p.pno_occupancy);
  r.csv_osv_moment = occupancy_adjusted_moments(p.osv_M, p.osv_occupancy);
  return r;
}

}  // namespace sequant::tests::csv_sizes

#endif  // SEQUANT_TESTS_UNIT_CSV_SIZE_FIXTURES_HPP
