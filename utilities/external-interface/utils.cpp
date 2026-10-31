#include "utils.hpp"

#include <SeQuant/core/index.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/string.hpp>

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
