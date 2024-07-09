#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/export/export.hpp>
#include <SeQuant/core/tensor.hpp>

#include <string>

namespace sequant::details {

void rename(Tensor &tensor, std::size_t counter) {
  tensor = Tensor(std::wstring(tensor.label()) + std::to_wstring(counter),
                  tensor.bra(), tensor.ket(), tensor.symmetry(),
                  tensor.braket_symmetry(), tensor.particle_symmetry());
}

void rename(Variable &variable, std::size_t counter) {
  variable =
      Variable(std::wstring(variable.label()) + std::to_wstring(counter));
}

}  // namespace sequant::details
