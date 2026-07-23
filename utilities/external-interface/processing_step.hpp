#ifndef SEQUANT_EXTERNAL_INTERFACE_PROCESSINGSTEP_HPP
#define SEQUANT_EXTERNAL_INTERFACE_PROCESSINGSTEP_HPP

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/utility/exception.hpp>

#include <nlohmann/json_fwd.hpp>

#include <concepts>
#include <string>
#include <string_view>
#include <type_traits>
#include <variant>
#include <vector>

namespace sequant::util::extint {

struct VoidInput {
  static constexpr std::string_view name{"VoidInput"};
};

struct ExpressionInput {
  static constexpr std::string_view name{"ExpressionInput"};

  std::vector<ResultExpr> expressions;
};

struct ExportTreeInput {
  static constexpr std::string_view name{"ExportTreeInput"};
};

using ProcessingInput =
    std::variant<VoidInput, ExpressionInput, ExportTreeInput>;

template <typename Target, typename Input>
  requires(std::is_assignable_v<std::remove_cvref_t<Input>, Target>)
auto convert_input(Input &&input) {
  if constexpr (std::same_as<Target, VoidInput>) {
    return VoidInput{};
  }

  using RetType = decltype(meta::forward_like<Input>(std::declval<Target>()));

  auto res = std::visit(
      [](auto &&inp) -> RetType {
        using Current = std::remove_cvref_t<decltype(inp)>;
        if constexpr (std::is_same_v<Target, Current>) {
          return std::forward<decltype(inp)>(inp);
        }

        throw Exception("Can't convert " + std::string(Current::name) + " to " +
                        std::string(Target::name));
      },
      std::forward<Input>(input));

  return res;
}

struct VoidOutput {};

struct ExpressionOutput {
  std::vector<ResultExpr> expressions;
};

struct ExportTreeOutput {};

using ProcessingOutput =
    std::variant<VoidOutput, ExpressionOutput, ExportTreeOutput>;

template <typename Target, typename Output>
  requires(std::is_assignable_v<std::remove_cvref_t<Output>, Target>)
auto convert_output(Output &&input) {
  if constexpr (std::same_as<Target, VoidOutput>) {
    return VoidOutput{};
  }

  struct InputConverter {
    auto operator()(VoidOutput) const {
      throw Exception("Can't convert VoidOutput to " +
                      std::string(Target::name));
    }

    auto operator()(ExpressionOutput &&inp) const {
      if constexpr (std::is_same_v<Target, ExpressionOutput>) {
        return std::get<Target>(std::forward<ExpressionOutput>(inp));
      }

      throw Exception("Can't convert ExpressionOutput to " +
                      std::string(Target::name));
    }

    auto operator()(ExportTreeOutput &&inp) const {
      if constexpr (std::is_same_v<Target, ExportTreeOutput>) {
        return std::get<Target>(std::forward<ExportTreeOutput>(inp));
      }

      throw Exception("Can't convert ExportTreeOutput to " +
                      std::string(Target::name));
    }
  };

  return std::visit(InputConverter{}, std::forward<Output>(input));
}

class ProcessingStep {
 public:
  ProcessingStep() = default;
  virtual ~ProcessingStep();

  virtual std::string kind() const = 0;

  virtual bool accepts_options() const = 0;
  virtual bool requires_options() const = 0;
  virtual void set_options(const nlohmann::json &options) = 0;

  virtual ProcessingOutput run(const ProcessingInput &input) = 0;
};

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_PROCESSINGSTEP_HPP
