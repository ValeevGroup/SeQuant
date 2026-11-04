//
// Integral projection (PNO domain projection and integral-basis stamps): a
// post-Wick rewrite of R2 terms; run after spintrace, before density fitting.
// Riplinger & Neese, JCP 138, 034106 (2013); Jiang et al., JCP 161, 082502.
//

#ifndef SEQUANT_DOMAIN_MBPT_RULES_PNO_HPP
#define SEQUANT_DOMAIN_MBPT_RULES_PNO_HPP

#include <SeQuant/core/basis.hpp>
#include <SeQuant/core/expr.hpp>

#include <cstdint>
#include <functional>
#include <optional>
#include <string>

namespace sequant::mbpt {

/// @return true if @p idx sits in a pure-occupied space of the default
/// context's index-basis registry
[[nodiscard]] bool is_occupied(Index const& idx);

/// (ov|ov): each Mulliken column (bra[k], ket[k]) holds exactly one occupied
/// (exchange, K)
[[nodiscard]] bool is_ovov_integral(AbstractTensor const& tnsr);

/// (oo|vv): homogeneous Mulliken columns that differ (Coulomb, J)
[[nodiscard]] bool is_oovv_integral(AbstractTensor const& tnsr);

/// the R2 terms the projection fires on; the four cells partition eligible R2
/// terms by (Mulliken class, order in the amplitudes), where linear means <= 1
/// in the amplitudes
enum class ProjectionTerms : std::uint32_t {
  None = 0u,
  /// `(ov|ov)` integrals in R2 terms of order <= 1 in the amplitudes
  ExchangeLinear = 1u << 0,
  /// `(ov|ov)` integrals in R2 terms of order >= 2 in the amplitudes
  ExchangeNonlinear = 1u << 1,
  /// `(oo|vv)` integrals in R2 terms of order <= 1 in the amplitudes
  CoulombLinear = 1u << 2,
  /// `(oo|vv)` integrals in R2 terms of order >= 2 in the amplitudes
  CoulombNonlinear = 1u << 3,

  Exchange = ExchangeLinear | ExchangeNonlinear,
  Coulomb = CoulombLinear | CoulombNonlinear,
  Linear = ExchangeLinear | CoulombLinear,
  Nonlinear = ExchangeNonlinear | CoulombNonlinear,
  All = Exchange | Coulomb
};

[[nodiscard]] constexpr ProjectionTerms operator|(ProjectionTerms a,
                                                  ProjectionTerms b) noexcept {
  return static_cast<ProjectionTerms>(static_cast<std::uint32_t>(a) |
                                      static_cast<std::uint32_t>(b));
}
[[nodiscard]] constexpr ProjectionTerms operator&(ProjectionTerms a,
                                                  ProjectionTerms b) noexcept {
  return static_cast<ProjectionTerms>(static_cast<std::uint32_t>(a) &
                                      static_cast<std::uint32_t>(b));
}
/// @return true if @p a selects at least one cell
[[nodiscard]] constexpr bool any(ProjectionTerms a) noexcept {
  return static_cast<std::uint32_t>(a) != 0u;
}
/// @return true if every cell of @p sub is selected by @p set
[[nodiscard]] constexpr bool contains(ProjectionTerms set,
                                      ProjectionTerms sub) noexcept {
  return (set & sub) == sub;
}

/// where a moved leg lands: its partner's pair (a basis stamp) or the
/// integral's own pair, its two distinct occupied indices (PNO)
enum class ProjectionDomain { PartnerPair, OwnPair };

/// the basis instance a moved leg takes: its partner's (the amplitude family)
/// or the one ProjectionOptions::cell_instance names for the integral's cell
enum class ProjectionBasis { Amplitude, Integral };

struct ProjectionOptions {
  ProjectionTerms terms = ProjectionTerms::None;
  ProjectionDomain domain = ProjectionDomain::OwnPair;
  ProjectionBasis basis = ProjectionBasis::Integral;
  /// the instance of a selected cell's integral legs; required for `Integral`
  std::function<IndexBasis::instance_type(ProjectionTerms cell)> cell_instance =
      {};
  std::wstring integral_label = L"g";
  /// also move a virtual leg that is a residual external (the overlap then
  /// carries the external): every leg of a selected integral lands on the
  /// target pair, so a density-fitted `(vv|Κ)` factor of an own-pair
  /// integral can fold into that pair's leaf; an approximation beyond the
  /// CSV truncation itself, exact at complete pair spaces. Off by default.
  bool move_externals = false;
  /// under `Amplitude`, the basis every moved leg is minted in, instead of
  /// the leg's own space and instance: a Kramers-restricted kit passes its
  /// spin-free union expansion basis, so a flavoured leg moves into the
  /// pair's whole PNS space (the pair vectors have support on both Kramers
  /// flavours; the flavoured half alone is not a projector). Ignored under
  /// `Integral`.
  std::optional<IndexBasis> leg_basis = {};
};

/// In each R2 term (residual rank read off the `Ŝ`/`Â` symmetrizer; other
/// terms are returned as is; a term without a symmetrizer, a residual block
/// whose `Â` the caller expanded or stripped as the Kramers-restricted trace
/// does, takes its externals from its own external indices and its rank from
/// the number of unoccupied ones) every `integral_label` tensor of a selected
/// cell gets each virtual leg `x` that is not a residual external replaced by a
/// fresh `x'` on the target pair and basis, times the overlap `s{x';x}`
/// (`s{x;x'}` for a bra leg); a leg already there stays. An approximation,
/// exact only at complete domains, and opt-in: with the default
/// ProjectionOptions::terms (`None`) nothing changes. Run it after spin
/// tracing and before density fitting.
/// @return @p expr itself if nothing changed
/// @note the result holds fresh temporary indices: canonicalize it before
/// serializing (deserialization rejects ordinals >= Index::min_tmp_index())
/// A leg shared by two integrals of one term moves once per integral, each
/// to a fresh index, so the overlaps chain through the leg: `x' s{x';x}`
/// ... `s{x;x''} x''`.
/// @throw Exception if an integral of an R2 term carries auxiliary indices,
/// `Integral` lacks `cell_instance`, or `PartnerPair` is combined with
/// `Amplitude`
[[nodiscard]] ExprPtr project_integral_domains(ExprPtr const& expr,
                                               ProjectionOptions const& opts);

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_RULES_PNO_HPP
