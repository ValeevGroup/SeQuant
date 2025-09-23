#include <SeQuant/core/context.hpp>
#include <SeQuant/core/index_space_registry.hpp>
#include <SeQuant/core/space.hpp>

#include <string_view>

namespace sequant {

// Note: we can't add a templated ctor for IndexSpace as that would have to be
// defined in the header, requiring an include to index_space_registry,
// resulting in a circular include issue.
template <typename StrView>
IndexSpace retrieve(StrView label) {
  auto registry_ptr = get_default_context().index_space_registry();

  if (!registry_ptr) {
    throw std::runtime_error(
        "Can't use IndexSpace(string) without an active index space registry");
  }

  return registry_ptr->retrieve(label);
}

IndexSpace::IndexSpace(std::string_view label) : IndexSpace(retrieve(label)) {}

IndexSpace::IndexSpace(std::wstring_view label) : IndexSpace(retrieve(label)) {}

std::string to_string(IndexSpace::Type type) {
  std::ostringstream oss;
  oss << "0x" << std::hex << type.to_int32();
  return oss.str();
}

std::string to_string(IndexSpace::QuantumNumbers qns) {
  std::ostringstream oss;
  oss << "0x" << std::hex << qns.to_int32();
  return oss.str();
}

std::string to_string(IndexSpace::Attr attr) {
  std::ostringstream oss;
  oss << "{type=" << to_string(attr.type()) << ",qns=" << to_string(attr.qns())
      << "}";
  return oss.str();
}

std::string to_string(const IndexSpace& space) {
  std::ostringstream oss;
  oss << "{attr=" << to_string(space.attr());
  if (space.base_key().empty() == false) {
    oss << ",base_key=" << to_string(space.base_key());
  }
  oss << ",approximate_size=" << std::to_string(space.approximate_size())
      << "}";
  return oss.str();
}

}  // namespace sequant
