#include <SeQuant/core/export/itf_generator.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/tensor.hpp>

#include <optional>
#include <string>

namespace sequant {

std::string ItfGeneratorContext::index_name(const IndexSpace &space,
                                            std::size_t ordinal) const {
  // TODO
  return "Steve";
}

std::string ItfGeneratorContext::get_name(const IndexSpace &space) const {
  // TODO
  return "MySpace";
}

std::string ItfGeneratorContext::get_tag(const IndexSpace &space) const {
  // TODO
  return "T";
}

std::optional<std::string> ItfGeneratorContext::import_name(
    const Tensor &tensor) const {
  return {};
}

std::optional<std::string> ItfGeneratorContext::import_name(
    const Variable &variable) const {
  return {};
}

}  // namespace sequant
