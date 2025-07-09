#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/export/export.hpp>
#include <SeQuant/core/expr.hpp>

#include <cwctype>
#include <string>

namespace sequant::detail {

auto get_with_default(const auto &map, const auto &key, auto def_value) {
  auto it = map.find(key);

  return it == map.end() ? def_value : it->second;
}

bool rename(Tensor &tensor, PreprocessResult &result) {
  assert(result.tensorReferences.at(tensor) > 0);

  std::wstring basename = [&tensor]() {
    std::size_t size = tensor.label().size();

    while (size > 0 && std::iswdigit(tensor.label().at(size - 1))) {
      size--;
    }
    assert(size > 0);

    return std::wstring(tensor.label().substr(0, size));
  }();

  tensor.set_label(basename);

  std::size_t counter = 2;
  while (get_with_default(result.tensorReferences, tensor, 0) > 0) {
    tensor.set_label(basename + std::to_wstring(counter));
    counter++;
  }

  return result.tensorReferences.find(tensor) != result.tensorReferences.end();
}

bool rename(Variable &variable, PreprocessResult &result) {
  assert(result.variableReferences.at(variable) > 0);

  std::wstring basename = [&variable]() {
    std::size_t size = variable.label().size();

    while (size > 0 && std::iswdigit(variable.label().at(size - 1))) {
      size--;
    }
    assert(size > 0);

    return std::wstring(variable.label().substr(0, size));
  }();

  variable.set_label(basename);

  std::size_t counter = 2;
  while (get_with_default(result.variableReferences, variable, 0) > 0) {
    variable.set_label(basename + std::to_wstring(counter));
    counter++;
  }

  return result.variableReferences.find(variable) !=
         result.variableReferences.end();
}

}  // namespace sequant::detail
