//
// Created by Eduard Valeyev on 2019-04-01.
//

#ifndef SEQUANT_CONVENTION_HPP
#define SEQUANT_CONVENTION_HPP

namespace sequant {
namespace mbpt {

enum class Convention { QCiFS };

// clang-format off
/// @brief Loads defaults for Convention @c conv

/// This registers IndexSpace objects standard for the chosen convention,
/// updates default context's IndexRegistry object, and
/// updates TensorCanonicalizer cardinal labels
/// @warning should be only called once
/// @param conv convention to load; the only supported convention at the moment
///        is QCiFS
/// @note The QCiFS convention introduced the following definitions:
///
/// | Label | Space                               | Comment                                                                                                                                                      |
/// |-------|-------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------|
/// | `p`   | IndexSpace::all                     |  p,q,r... for OBS spstates introduced in [DOI 10.1063/1.444231 (QCiFS I)](https://dx.doi.org/10.1063/1.444231)                                               |
/// | `i`   | IndexSpace::active_occupied         |  i,j,k... for active occupied spstates introduced in [DOI 10.1063/1.446736 (QCiFS III)](https://dx.doi.org/10.1063/1.446736)                                 |
/// | `a`   | IndexSpace::active_unoccupied       |  a,b,c... for active unoccupied spstates introduced in [DOI 10.1063/1.446736 (QCiFS III)](https://dx.doi.org/10.1063/1.446736)                               |
/// | `α`   | IndexSpace::complete_unoccupied     |  α,β... for complete unoccupied spstates introduced in [DOI 10.1063/1.459921 (MP2-R12 I)](https://dx.doi.org/10.1063/1.459921)                               |
/// | `κ`   | IndexSpace::complete                |  κ,𝛌... for complete spstates introduced in [DOI 10.1063/1.459921 (MP2-R12 I)](https://dx.doi.org/10.1063/1.459921)                                          |
/// | \c α' | IndexSpace::other_unoccupied        |  α',β'... for orthogonal complement to OBS (CABS) introduced in [DOI 10.1016/j.cplett.2004.07.061 (CABS)](https://dx.doi.org/10.1016/j.cplett.2004.07.061)   |
/// | `m`   | IndexSpace::occupied                |  m,n.. for all occupied (including inactive/frozen orbitals) de facto introduced in [DOI 10.1016/j.cplett.2004.07.061 (CABS)](https://dx.doi.org/10.1016/j.cplett.2004.07.061), though formally not explicitly defined so |
/// | `e`   | IndexSpace::unoccupied              |  e,f... for all unoccupied (including inactive/frozen orbitals) used internally in MPQC LCAOWavefunction |
/// | `g`   | IndexSpace::OBS_unfrozen            |  unfrozen orbitals in the OBS
/// | `x`   | IndexSpace::all_active              |  used internally in MPQC for GF, CT-F12, and other ad hoc uses |
/// | `γ`   | IndexSpace::complete_inactive_unoccupied | useful in CT-F12 with active space geminals  |
/// | `z`   | IndexSpace::complete_unfrozen       | union of all unfrozen orbitals including the CABS orbitals. may be useful towards universal projector. |
/// | `u`   | IndexSpace::active                  |  origin unknown, for recent use see [DOI 10.1063/5.0067511](https://dx.doi.org/10.1063/5.0067511) |
/// | `I`   | IndexSpace::active_maybe_occupied   |  origin unknown, for recent use see [DOI 10.1063/5.0067511](https://dx.doi.org/10.1063/5.0067511); N.B. although QCiFS uses capital letters for spin-free indices, since there is usually no need to distinguish spin-orbital from spin-free indices (the type is deduced from the context), this use of capital letters seems preferable |
/// | `A`   | IndexSpace::active_maybe_unoccupied |  origin unknown, for recent use see [DOI 10.1063/5.0067511](https://dx.doi.org/10.1063/5.0067511); N.B. although QCiFS uses capital letters for spin-free indices, since there is usually no need to distinguish spin-orbital from spin-free indices (the type is deduced from the context), this use of capital letters seems preferable |
/// | `M`   | IndexSpace::maybe_occupied          |  combination of of `I` for `m` |
/// | `E`   | IndexSpace::maybe_unoccupied        |  combination of of `A` for `e` |
/// | `Δ`   | IndexSpace::complete_maybe_unoccupied | combination of `α` and `A`; since capital Alpha looks indistinguishable from `A` use `Δ`                          |
// clang-format on
void set_default_convention(Convention conv = Convention::QCiFS);

}  // namespace mbpt
}  // namespace sequant

#endif  // SEQUANT_CONVENTION_HPP
