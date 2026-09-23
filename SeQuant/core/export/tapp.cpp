#include <SeQuant/core/export/tapp.hpp>
#include <SeQuant/core/export/utils.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/space.hpp>

#include <string>
#include <utility>

namespace sequant {

std::string TAPPGeneratorContext::get_tag(const IndexSpace &space) const {
  auto it = m_index_tags.find(space);

  if (it == m_index_tags.end()) {
    return detail::sanitize_identifier(space.base_key(),
                                       detail::NonAsciiPolicy::Escape);
  }

  return it->second;
}

std::string TAPPGeneratorContext::get_dim(const IndexSpace &space) const {
  auto it = m_index_dims.find(space);

  if (it == m_index_dims.end()) {
    return "dim_" + get_tag(space);
  }

  return it->second;
}

std::string TAPPGeneratorContext::get_tags(const Tensor &tensor) const {
  std::string tags;

  for (const Index &idx : tensor.const_indices()) {
    tags += get_tag(idx.space());
  }

  return tags;
}

void TAPPGeneratorContext::set_tag(const IndexSpace &space, std::string tag) {
  m_index_tags[space] = std::move(tag);
}

void TAPPGeneratorContext::set_dim(const IndexSpace &space, std::string dim) {
  m_index_dims[space] = std::move(dim);
}

const std::string &TAPPGeneratorContext::prefix() const { return m_prefix; }

void TAPPGeneratorContext::set_prefix(std::string prefix) {
  m_prefix = std::move(prefix);
}

TAPPScalarType TAPPGeneratorContext::scalar_type() const {
  return m_scalar_type;
}

void TAPPGeneratorContext::set_scalar_type(TAPPScalarType type) {
  m_scalar_type = type;
}

}  // namespace sequant
