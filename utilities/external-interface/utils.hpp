#ifndef SEQUANT_EXTERNAL_INTERFACE_UTILS_HPP
#define SEQUANT_EXTERNAL_INTERFACE_UTILS_HPP

#include <SeQuant/core/index.hpp>
#include <SeQuant/core/space.hpp>

#include <map>
#include <string>

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

#endif
