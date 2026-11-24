#ifndef SEQUANT_EXTERNAL_INTERFACE_UTILS_HPP
#define SEQUANT_EXTERNAL_INTERFACE_UTILS_HPP

#include <SeQuant/core/index.hpp>
#include <SeQuant/core/space.hpp>

#include <filesystem>
#include <map>
#include <set>
#include <string>
#include <string_view>

class IndexSpaceMeta {
 public:
  struct Entry {
    std::string tag;
    std::string name;
  };

  IndexSpaceMeta() = default;

  std::size_t getSize(const sequant::IndexSpace &space) const;

  std::size_t getSize(const sequant::Index &index) const;

  std::string getName(const sequant::IndexSpace &space) const;

  std::string getTag(const sequant::IndexSpace &space) const;

  void registerSpace(sequant::IndexSpace space, Entry entry);

 private:
  std::map<sequant::IndexSpace, Entry> m_entries;
};

/// @return the base_key() of every index space in the default context's
/// registry
std::set<std::wstring> registered_space_labels();

/// Creates the parent directories of @p path, if it names any
/// @throw sequant::Exception mentioning @p step_kind if that fails
void create_parent_directories(const std::filesystem::path &path,
                               std::string_view step_kind);

#endif
