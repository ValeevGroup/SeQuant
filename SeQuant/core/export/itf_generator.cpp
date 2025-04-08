#include <SeQuant/core/export/itf_generator.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <cctype>
#include <limits>
#include <locale>
#include <optional>
#include <string>

namespace sequant {

std::string ItfGeneratorContext::index_name(const IndexSpace &space,
                                            std::size_t ordinal) const {
  std::string base_key = toUtf8(space.base_key());

  char limit = [&]() -> char {
    auto it = m_index_label_limits.find(space);
    if (it == m_index_label_limits.end()) {
      // Note that in between capital and lowercase letters there are symbols
      // [,\,],^,_ and `. Hence, this default is not quite safe but should work
      // in most cases anyway
      return 'z';
    }

    return it->second;
  }();

  if (base_key.empty() || base_key.size() > 1 ||
      base_key[0] + ordinal > limit) {
    return base_key + std::to_string(ordinal);
  }

  // Merge base key and ordinal into a single letter. That is
  // i1 -> i
  // i2 -> j
  // i3 -> k
  // etc.
  char target = base_key[0] + ordinal;

  return std::string(1, target);
}

std::string ItfGeneratorContext::get_name(const IndexSpace &space) const {
  auto it = m_space_names.find(space);

  if (it == m_space_names.end()) {
    throw std::runtime_error("No name known for index space '" +
                             toUtf8(space.base_key()) + "'");
  }

  return it->second;
}

std::string ItfGeneratorContext::get_tag(const IndexSpace &space) const {
  auto it = m_tags.find(space);

  if (it == m_tags.end()) {
    throw std::runtime_error("No tag known for index space '" +
                             toUtf8(space.base_key()) + "'");
  }

  return it->second;
}

std::optional<std::string> ItfGeneratorContext::import_name(
    const Tensor &tensor) const {
  auto it = m_tensor_imports.find(tensor);

  if (it == m_tensor_imports.end()) {
    return {};
  }

  return it->second;
}

std::optional<std::string> ItfGeneratorContext::import_name(
    const Variable &variable) const {
  auto it = m_variable_imports.find(variable);

  if (it == m_variable_imports.end()) {
    return {};
  }

  return it->second;
}

void ItfGeneratorContext::set_name(const IndexSpace &space, std::string name) {
  m_space_names[space] = std::move(name);
}

void ItfGeneratorContext::set_tag(const IndexSpace &space, std::string tag) {
  m_tags[space] = std::move(tag);
}

void ItfGeneratorContext::set_import_name(const Tensor &tensor,
                                          std::string name) {
  m_tensor_imports[tensor] = std::move(name);
}

void ItfGeneratorContext::set_import_name(const Variable &variable,
                                          std::string name) {
  m_variable_imports[variable] = std::move(name);
}

bool ItfGeneratorContext::rewrite(Tensor &tensor) const {
  // TODO: implement integral remapping
  // Maybe also do this generally, based on the specified tensor symmetries
  return false;
}

}  // namespace sequant
