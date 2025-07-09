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

}  // namespace sequant
