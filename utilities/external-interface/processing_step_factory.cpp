#include "processing_step_factory.hpp"
#include "processing_step.hpp"

#include <SeQuant/core/utility/exception.hpp>

#include <map>
#include <memory>
#include <string>
#include <string_view>

namespace sequant::util::extint {

ProcessingStepFactory &ProcessingStepFactory::instance() {
  static ProcessingStepFactory factory;

  return factory;
}

std::unique_ptr<ProcessingStep> ProcessingStepFactory::instantiate(
    std::string_view type) const {
  auto it = instatiators_.find(type);

  if (it == instatiators_.end()) {
    throw Exception(
        "Attempted to instantiated unknown processing step of kind '" +
        std::string(type) + "'");
  }

  return it->second();
}

bool ProcessingStepFactory::register_class(std::string type,
                                           InstatiateFunc func) {
  auto [it, inserted] = instatiators_.emplace(std::move(type), std::move(func));

  if (!inserted) {
    throw Exception("Duplicate processing step type name '" + it->first + "'");
  }

  return inserted;
}

}  // namespace sequant::util::extint
