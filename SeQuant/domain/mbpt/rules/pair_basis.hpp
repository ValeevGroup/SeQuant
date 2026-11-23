#ifndef SEQUANT_DOMAIN_MBPT_RULES_PAIR_BASIS_HPP
#define SEQUANT_DOMAIN_MBPT_RULES_PAIR_BASIS_HPP

#include <SeQuant/core/expr.hpp>

#include <cstddef>
#include <string>

namespace sequant::mbpt {

struct PairBasisFoldOptions {
  /// the label of the density-fitted integrals, `g{a;b;Κ}`
  std::wstring integral_label = L"g";
  /// the label of the CSV coefficients the CSV transform introduced
  std::wstring coeff_label = L"C";
};

/// what fold_pair_basis_integrals() did to an expression
struct PairBasisFoldReport {
  /// integrals replaced, with their two coefficients, by a pair-basis leaf
  std::size_t folded = 0;
  /// integrals left alone because the conjugate coefficient sits on the ket
  /// leg and the stored-orientation one on the bra leg: `Σ_ab C*{a;x} g{b;a;Κ}
  /// C{b;y}` is not the leaf and gets no spelling of its own here
  std::size_t mirrored = 0;
};

/// Folds a density-fitted integral and the two CSV coefficients of one pair
/// that project its virtual legs into the integral in that pair's basis.
///
/// In every product of @p expr (nested sums and products are visited), an
/// integral `g{a;b;Κ}` (one bra, one ket and one auxiliary index, both legs
/// unoccupied and free of proto indices, no conjugation mark) whose bra leg
/// `a` is the ket of a coefficient `C{x<P>;a}` -- the proto index in the bra:
/// the conjugate orientation -- and whose ket leg `b` is the bra of a
/// coefficient `C{b;y<P>}` -- the proto index in the ket: the stored
/// orientation -- with the same proto list `P` on both, is replaced together
/// with the two coefficients by the leaf. The conjugation of a coefficient is
/// its slot orientation, as the CSV transform spells it; a coefficient
/// carrying the adjoint mark with the proto index in the bra is accepted as
/// the same conjugate.
///
///     g{x<P>;y<P>;Κ} = Σ_ab conj(C{a;x<P>}) g{a;b;Κ} C{b;y<P>}
///
/// The leaf keeps the integral's label, conjugation parity and Kramers
/// symmetry; its permutational and bra-ket symmetries are `Nonsymm` (the two
/// legs run over different PNS sets). Nothing else in the product changes; the
/// scalar stays. The leaf's value depends on the coefficients and the
/// integrals only, so a consumer can build it once and keep it across
/// iterations.
///
/// Run it after the density fitting and the CSV transform.
/// @return @p expr itself if nothing changed
[[nodiscard]] ExprPtr fold_pair_basis_integrals(
    ExprPtr const& expr, PairBasisFoldOptions const& opts = {},
    PairBasisFoldReport* report = nullptr);

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_RULES_PAIR_BASIS_HPP
