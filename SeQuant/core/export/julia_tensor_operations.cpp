#include <SeQuant/core/export/julia_tensor_operations.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/utility/string.hpp>

namespace sequant {

JuliaTensorOperationsGeneratorContext::JuliaTensorOperationsGeneratorContext(
    TagMap index_tags, DimMap index_dims)
    : m_index_tags(std::move(index_tags)),
      m_index_dims(std::move(index_dims)) {}

std::string JuliaTensorOperationsGeneratorContext::get_tag(
    const IndexSpace &space) const {
  auto it = m_index_tags.find(space);

  if (it == m_index_tags.end()) {
    throw std::runtime_error("No known tags for indices of space \"" +
                             toUtf8(space.base_key()) + "\"");
  }

  return it->second;
}

std::string JuliaTensorOperationsGeneratorContext::get_dim(
    const IndexSpace &space) const {
  auto it = m_index_dims.find(space);

  if (it == m_index_dims.end()) {
    // Auto-generate a variable name for the space's dimension
    return "dim_" + get_tag(space);
  }

  return it->second;
}

std::string JuliaTensorOperationsGeneratorContext::get_tags(
    const Tensor &tensor) const {
  std::string tags;

  for (const Index &idx : tensor.const_braket()) {
    tags += get_tag(idx.space());
  }

  return tags;
}

}  // namespace sequant
