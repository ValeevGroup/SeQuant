#include <SeQuant/core/export/julia_tensoroperationsgen.hpp>

namespace sequant {

JuliaTensorOperationsGenContext::JuliaTensorOperationsGenContext(
    std::map<IndexSpace, std::string> index_tags,
    std::map<IndexSpace, std::string> index_dims,
    bool print_intermediate_comments)
    : m_index_tags(std::move(index_tags)),
      m_index_dims(std::move(index_dims)),
      m_print_intermediate_comments(print_intermediate_comments) {}

std::string JuliaTensorOperationsGenContext::get_tag(
    const IndexSpace &indexspace) const {
  auto it = m_index_tags.find(indexspace);
  if (it != m_index_tags.end()) {
    return it->second;
  } else {
    std::string bad_index_string = toUtf8(indexspace.base_key());
    throw std::runtime_error("Tag of index \"" + bad_index_string +
                             "\" not found");
  }
}

std::string JuliaTensorOperationsGenContext::get_dim(
    const IndexSpace &indexspace) const {
  auto it = m_index_dims.find(indexspace);
  if (it != m_index_dims.end()) {
    return it->second;
  } else {
    return "dim_" + get_tag(indexspace);
  }
}

std::string JuliaTensorOperationsGenContext::get_tags(
    const Tensor &tensor) const {
  std::string tagstring;
  const auto &indices = tensor.const_braket();

  for (const auto &index : indices) {
    std::string tag = get_tag(index.space());
    tagstring += tag;
  }

  return tagstring;
}

bool JuliaTensorOperationsGenContext::print_intermediate_comments() const {
  return m_print_intermediate_comments;
}
}  // namespace sequant