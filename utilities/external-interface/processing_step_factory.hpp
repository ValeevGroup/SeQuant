#ifndef SEQUANT_EXTERNAL_INTERFACE_PROCESSING_STEP_FACTORY_HPP
#define SEQUANT_EXTERNAL_INTERFACE_PROCESSING_STEP_FACTORY_HPP

#include "processing_step.hpp"

#include <SeQuant/core/utility/macros.hpp>

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <string_view>

namespace sequant::util::extint {

class ProcessingStepFactory {
 public:
  using InstatiateFunc = std::function<std::unique_ptr<ProcessingStep>()>;

  static ProcessingStepFactory &instance();

  std::unique_ptr<ProcessingStep> instantiate(std::string_view type) const;

  bool register_class(std::string type, InstatiateFunc func);

 private:
  std::map<std::string, InstatiateFunc, std::less<>> instatiators_;

  ProcessingStepFactory() = default;
};

#define SEQUANT_EXTINT_REGISTER_STEP_TYPE(class_name, type_string)         \
  namespace {                                                              \
  static const bool SEQUANT_CONCAT(SEQUANT_CONCAT(registered_step_type_,   \
                                                  class_name),             \
                                   __LINE__) =                             \
      ProcessingStepFactory::instance().register_class(type_string, []() { \
        return std::make_unique<class_name>();                             \
      });                                                                  \
  }

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_PROCESSING_STEP_FACTORY_HPP
