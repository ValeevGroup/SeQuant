#ifndef SEQUANT_CORE_EXPORT_MARKED_NAME_HPP
#define SEQUANT_CORE_EXPORT_MARKED_NAME_HPP

#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <string>

namespace sequant {

/// @return Tensor::label() followed by one ASCII suffix per set core state,
///         `_adj` for the adjointed `⁺` and `_conj` for the K-conjugated `꙳`.
///         A marked tensor is a different array than the one label() names, so
///         the marks belong in the array name; spelled in ASCII, that name is
///         an identifier in every export target, unlike
///         Tensor::decorated_label().
inline std::wstring export_label(const Tensor &tensor) {
  std::wstring result(tensor.label());
  if (tensor.adjointed()) result += L"_adj";
  if (tensor.kconjugated()) result += L"_conj";
  return result;
}

/// @return export_label() as UTF-8
inline std::string export_name(const Tensor &tensor) {
  return toUtf8(export_label(tensor));
}

/// @brief renames @p tensor to export_label() and clears its core states, so
///        that a marked tensor is an array of its own wherever the export
///        machinery keys on a label: the block-keyed maps (declarations,
///        reference counts, load strategies, import names) and the backends'
///        own label matching, such as the ITF two-electron integral remap.
///        A tensor with neither state set is left untouched.
/// @note clearing both states consumes no sign, so the tensor still denotes
///       the value it denoted
inline void fold_marks_into_label(Tensor &tensor) {
  if (!tensor.adjointed() && !tensor.kconjugated()) return;

  tensor.set_label(export_label(tensor));
  [[maybe_unused]] const auto sign = tensor.set_states(false, false);
  SEQUANT_ASSERT(sign == 1);
}

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_MARKED_NAME_HPP
