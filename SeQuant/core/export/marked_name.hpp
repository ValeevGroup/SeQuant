#ifndef SEQUANT_CORE_EXPORT_MARKED_NAME_HPP
#define SEQUANT_CORE_EXPORT_MARKED_NAME_HPP

#include <SeQuant/core/expressions/tensor.hpp>
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

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_MARKED_NAME_HPP
