#include "utils.hpp"

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/index_basis_registry.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <system_error>

using namespace sequant;

std::size_t IndexSpaceMeta::getSize(const IndexSpace &space) const {
  return space.dimension();
}

std::size_t IndexSpaceMeta::getSize(const Index &index) const {
  return index.basis().extent();
}

std::string IndexSpaceMeta::getName(const IndexSpace &space) const {
  auto iter = m_entries.find(space);
  if (iter == m_entries.end()) {
    throw Exception("No known name for index space " +
                    toUtf8(space.base_key()));
  }

  return iter->second.name;
}

std::string IndexSpaceMeta::getTag(const IndexSpace &space) const {
  auto iter = m_entries.find(space);
  if (iter == m_entries.end()) {
    throw Exception("No known tag for index space " + toUtf8(space.base_key()));
  }

  return iter->second.tag;
}

void IndexSpaceMeta::registerSpace(IndexSpace space, Entry entry) {
  m_entries.insert({std::move(space), std::move(entry)});
}

std::set<std::wstring> registered_space_labels() {
  std::set<std::wstring> labels;
  for (const IndexSpace &space :
       get_default_context().index_basis_registry()->spaces()) {
    labels.insert(std::wstring(space.base_key()));
  }
  return labels;
}

void create_parent_directories(const std::filesystem::path &path,
                               std::string_view step_kind) {
  if (!path.has_parent_path()) {
    return;
  }
  std::error_code ec;
  std::filesystem::create_directories(path.parent_path(), ec);
  if (ec) {
    throw Exception("Failed to create output directory '" +
                    path.parent_path().string() + "' for " +
                    std::string(step_kind) + ": " + ec.message());
  }
}
