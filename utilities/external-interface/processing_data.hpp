#ifndef SEQUANT_EXTERNAL_INTERFACE_PROCESSINGDATA_HPP
#define SEQUANT_EXTERNAL_INTERFACE_PROCESSINGDATA_HPP

#include <SeQuant/core/export/export_node.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/utility/exception.hpp>

#include <optional>
#include <string_view>
#include <type_traits>
#include <variant>
#include <vector>

namespace sequant::util::extint {

struct ExpressionData {
  static constexpr std::string_view name{"ExpressionData"};

  std::vector<ResultExpr> expressions;
};

struct ExportTreeData {
  static constexpr std::string_view name{"ExportTreeData"};

  struct Entry {
    ExportNode<> tree;
    /// The result the current tree contributes to _after symmetrizing_ over the
    /// external indices. This implies that if this is set, symmetrization is
    /// required.
    std::optional<Tensor> symm_contribution_target;
  };

  std::vector<Entry> entries;
};

using ProcessingData = std::variant<ExpressionData, ExportTreeData>;

template <typename Target, typename Input>
  requires(std::is_assignable_v<std::remove_cvref_t<Input>, Target>)
decltype(auto) convert_data(Input &&input) {
  using RetType = decltype(meta::forward_like<Input>(std::declval<Target>()));

  return std::visit(
      [](auto &&inp) -> RetType {
        using Current = std::remove_cvref_t<decltype(inp)>;
        if constexpr (std::is_same_v<Target, Current>) {
          return std::forward<decltype(inp)>(inp);
        }

        throw Exception("Can't convert " + std::string(Current::name) + " to " +
                        std::string(Target::name));
      },
      std::forward<Input>(input));
}

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_PROCESSINGDATA_HPP
