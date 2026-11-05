#ifndef SEQUANT_CORE_IO_PRINT_OPTIONS_HPP
#define SEQUANT_CORE_IO_PRINT_OPTIONS_HPP

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/options.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <optional>

namespace sequant {

/// options that control how expressions are printed as LaTeX
/// @note these are presentation settings only: no library logic depends on
/// them
/// @sa Context::print_options(), io::latex::to_string()
struct PrintOptions {
  SEQUANT_DESIGNATED_INIT_ONLY;

  /// whether tensor indices of ket (covariant, primal) modes are typeset as
  /// superscript or as subscript
  BraKetTypesetting braket_typesetting = BraKetTypesetting::ContraSub;
  /// whether tensor indices are typeset using the `tensor` LaTeX package
  /// (preferred for alignment of superscript with subscript slots) or
  /// natively
  BraKetSlotTypesetting braket_slot_typesetting =
      BraKetSlotTypesetting::TensorPackage;
  /// if set, the top-level printers (io::latex::to_string(),
  /// to_latex_align()) print a copy of the expression whose summands and
  /// commuting factors are sorted in this order (see sort_expr()); for a
  /// complete normalization of the form use prettify()
  std::optional<LexicographicCompareOptions> ordering = std::nullopt;

  bool operator==(const PrintOptions&) const = default;
};

}  // namespace sequant

#endif  // SEQUANT_CORE_IO_PRINT_OPTIONS_HPP
