#ifndef SEQUANT_CORE_EXPORT_CONTEXT_HPP
#define SEQUANT_CORE_EXPORT_CONTEXT_HPP

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/utility/tensor.hpp>

#include <limits>
#include <map>
#include <optional>
#include <string>
#include <type_traits>

namespace sequant {

/// Different options for what to do when loading a tensor in the sense of
/// bringing it into memory to make it available for subsequent operations
enum class LoadStrategy {
  Load,
  Create,
};

/// Different options for when to (re)set a tensor's value to zero
enum class ZeroStrategy {
  NeverZero = 0,
  ZeroOnLoad = 1,
  ZeroOnCreate = 1 << 1,
  ZeroOnReuse = 1 << 2,
  AlwaysZero = ZeroOnLoad | ZeroOnCreate | ZeroOnReuse,
};

ZeroStrategy operator|(ZeroStrategy lhs, ZeroStrategy rhs);
ZeroStrategy operator&(ZeroStrategy lhs, ZeroStrategy rhs);
bool operator!(ZeroStrategy strategy);

/// A combination of a Load and ZeroStrategy
/// @sa LoadStrategy
/// @sa ZeroStrategy
struct StrategyPair {
  StrategyPair() = default;
  StrategyPair(LoadStrategy load, ZeroStrategy zero);
  StrategyPair(LoadStrategy load);
  StrategyPair(ZeroStrategy zero);

  LoadStrategy load = LoadStrategy::Create;
  ZeroStrategy zero = ZeroStrategy::ZeroOnCreate;
};

/// Different ways a given tensor/variable can be used in a given code.
enum class Usage {
  None = 0,
  Result = 0b001,
  Intermediate = 0b010,
  Terminal = 0b100,
};

/// A set of Usage scenarios for a given tensor/variable. Contrary to Usage,
/// this can hold a combination of usage scenarios.
/// @sa Usage
class UsageSet {
 public:
  using Value = std::underlying_type_t<Usage>;
  UsageSet() = default;

  UsageSet &operator=(Usage usage);

  UsageSet &operator|=(UsageSet usage);
  UsageSet &operator&=(UsageSet usage);
  UsageSet &operator|=(Usage usage);
  UsageSet &operator&=(Usage usage);

  friend bool operator&(UsageSet set, Usage usage);
  friend bool operator&(Usage usage, UsageSet set);

  friend bool operator==(UsageSet set, Usage usage);
  friend bool operator==(Usage usage, UsageSet set);

  friend bool operator!=(UsageSet set, Usage usage);
  friend bool operator!=(Usage usage, UsageSet set);

 private:
  Value m_val = static_cast<Value>(Usage::None);
};

/// The ExportContext is intended to provide context during export operations.
/// This task is achieved by two core tasks:
/// 1. Acting as a meta-information storage. Before export starts, certain meta
///    information are collected and stored on the context object.
/// 2. Providing a customization point for users to influence any export
///    behavior that is not universal for a given backend. In this sense, the
///    context object can be seen as a configuration storage system.
/// The context object has to be user-provided when initiating an export and
/// is passed down to every single function during the export context. Hence,
/// it acts as a communication channel down into the export process.
class ExportContext {
 public:
  using TensorStrategyMap =
      std::map<Tensor, StrategyPair, TensorBlockLessThanComparator>;
  using VariableStrategyMap = std::map<Variable, StrategyPair>;

  ExportContext() = default;
  ExportContext(TensorStrategyMap tensorMap,
                VariableStrategyMap variableMap = {});
  ExportContext(VariableStrategyMap map);
  virtual ~ExportContext();

  /// @returns The LoadStrategy for the given tensor
  /// @sa LoadStrategy
  [[nodiscard]] virtual LoadStrategy loadStrategy(const Tensor &tensor) const;
  /// @returns The LoadStrategy for the given variable
  /// @sa LoadStrategy
  [[nodiscard]] virtual LoadStrategy loadStrategy(
      const Variable &variable) const;
  /// Sets the LoadStrategy for the given tensor
  /// @param tensor The Tensor to set the strategy for
  /// @param strategy The LoadStrategy to set
  /// @param expression_id The expression for which to set the strategy for. If
  ///        left unspecified, the strategy will become global. The expression
  ///        ID is the ID of the root of the associated ExportNode
  /// @sa LoadStrategy
  /// @sa ExportNode
  virtual void setLoadStrategy(
      const Tensor &tensor, LoadStrategy strategy,
      const std::optional<std::size_t> &expression_id = {});
  /// Sets the LoadStrategy for the given variable
  /// @param variable The Variable to set the strategy for
  /// @param strategy The LoadStrategy to set
  /// @param expression_id The expression for which to set the strategy for. If
  ///        left unspecified, the strategy will become global. The expression
  ///        ID is the ID of the root of the associated ExportNode
  /// @sa LoadStrategy
  /// @sa ExportNode
  virtual void setLoadStrategy(
      const Variable &variable, LoadStrategy strategy,
      const std::optional<std::size_t> &expression_id = {});

  /// @returns The ZeroStrategy for the given tensor
  /// @sa ZeroStrategy
  [[nodiscard]] virtual ZeroStrategy zeroStrategy(const Tensor &tensor) const;
  /// @returns The ZeroStrategy for the given variable
  /// @sa ZeroStrategy
  [[nodiscard]] virtual ZeroStrategy zeroStrategy(
      const Variable &variable) const;
  /// Sets the ZeroStrategy for the given tensor
  /// @param tensor The Tensor to set the strategy for
  /// @param strategy The ZeroStrategy to set
  /// @param expression_id The expression for which to set the strategy for. If
  ///        left unspecified, the strategy will become global. The expression
  ///        ID is the ID of the root of the associated ExportNode
  /// @sa ZeroStrategy
  /// @sa ExportNode
  virtual void setZeroStrategy(
      const Tensor &tensor, ZeroStrategy strategy,
      const std::optional<std::size_t> &expression_id = {});
  /// Sets the ZeroStrategy for the given variable
  /// @param variable The Variable to set the strategy for
  /// @param strategy The ZeroStrategy to set
  /// @param expression_id The expression for which to set the strategy for. If
  ///        left unspecified, the strategy will become global. The expression
  ///        ID is the ID of the root of the associated ExportNode
  /// @sa ZeroStrategy
  /// @sa ExportNode
  virtual void setZeroStrategy(
      const Variable &variable, ZeroStrategy strategy,
      const std::optional<std::size_t> &expression_id = {});

  /// Performs necessary (backend-specific) modifications (if any) on the
  /// given Tensor in-place
  /// @param tensor The Tensor to be modified
  /// @returns Whether any modification has been performed
  virtual bool rewrite(Tensor &tensor) const;
  /// Performs necessary (backend-specific) modifications (if any) on the
  /// given Variable in-place
  /// @param variable The Variable to be modified
  /// @returns Whether any modification has been performed
  virtual bool rewrite(Variable &variable) const;

  /// @returns Whether the current export workflow is generating code inside a
  /// named section
  [[nodiscard]] virtual bool inside_named_section() const;
  /// @returns The name of the current section. If not inside a named section,
  /// this function will throw.
  [[nodiscard]] virtual const std::string &current_section_name() const;
  /// Sets the name of the section that subsequent export shall be for.
  /// Implicitly notifies the context that from now on we're in a named section.
  /// @param name The name of the current section
  virtual void set_current_section_name(std::string name);
  /// Resets the current section name. Implicitly notifies the context that
  /// we're no longer inside a named section.
  virtual void clear_current_section_name();

  /// @returns Whether or not the ID of the current expression is available
  [[nodiscard]] virtual bool has_current_expression_id() const;
  /// @returns The ID of the current expression. If none is available, this
  /// function will throw.
  [[nodiscard]] virtual std::size_t current_expression_id() const;
  /// Sets the current expression ID
  /// @param id The ID of the current expression
  virtual void set_current_expression_id(std::size_t id);
  /// Resets the ID of the current expression
  virtual void clear_current_expression_id();

 private:
  static constexpr std::size_t GLOBAL = std::numeric_limits<std::size_t>::max();

  std::map<std::size_t, TensorStrategyMap> m_tensorStrategies;
  std::map<std::size_t, VariableStrategyMap> m_variableStrategies;
  std::optional<std::string> m_currentSection;
  std::optional<std::size_t> m_currentExpressionID;
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_CONTEXT_HPP
