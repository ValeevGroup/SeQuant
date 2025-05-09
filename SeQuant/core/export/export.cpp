#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/export/export.hpp>
#include <SeQuant/core/tensor.hpp>

#include <string>

namespace sequant::details {

void rename(Tensor &tensor, std::size_t counter) {
  tensor.set_label(std::wstring(tensor.label()) + std::to_wstring(counter));
}

void rename(Variable &variable, std::size_t counter) {
  variable.set_label(std::wstring(variable.label()) + std::to_wstring(counter));
}

}  // namespace sequant::details
