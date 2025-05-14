#ifndef SEQUANT_CORE_EXPORT_CONTEXT_HPP
#define SEQUANT_CORE_EXPORT_CONTEXT_HPP

#include <SeQuant/core/export/utils.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/tensor.hpp>

#include <limits>
#include <map>
#include <optional>
#include <string>
#include <type_traits>

namespace sequant {

enum class LoadStrategy {
  Load,
  Create,
};

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

struct StrategyPair {
  StrategyPair() = default;
  StrategyPair(LoadStrategy load, ZeroStrategy zero);
  StrategyPair(LoadStrategy load);
  StrategyPair(ZeroStrategy zero);

  LoadStrategy load = LoadStrategy::Create;
  ZeroStrategy zero = ZeroStrategy::ZeroOnCreate;
};

enum class Usage {
  None = 0,
  Result = 0b001,
  Intermediate = 0b010,
  Terminal = 0b100,
};

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

class ExportContext {
 public:
  using TensorStrategyMap = std::map<Tensor, StrategyPair, TensorBlockCompare>;
  using VariableStrategyMap = std::map<Variable, StrategyPair>;

  ExportContext() = default;
  ExportContext(TensorStrategyMap tensorMap,
                VariableStrategyMap variableMap = {});
  ExportContext(VariableStrategyMap map);
  virtual ~ExportContext();

  virtual LoadStrategy loadStrategy(const Tensor &tensor) const;
  virtual LoadStrategy loadStrategy(const Variable &variable) const;
  virtual void setLoadStrategy(
      const Tensor &tensor, LoadStrategy strategy,
      const std::optional<std::size_t> &expression_id = {});
  virtual void setLoadStrategy(
      const Variable &variable, LoadStrategy strategy,
      const std::optional<std::size_t> &expression_id = {});

  virtual ZeroStrategy zeroStrategy(const Tensor &tensor) const;
  virtual ZeroStrategy zeroStrategy(const Variable &variable) const;
  virtual void setZeroStrategy(
      const Tensor &tensor, ZeroStrategy strategy,
      const std::optional<std::size_t> &expression_id = {});
  virtual void setZeroStrategy(
      const Variable &variable, ZeroStrategy strategy,
      const std::optional<std::size_t> &expression_id = {});

  virtual bool rewrite(Tensor &tensor) const;
  virtual bool rewrite(Variable &variable) const;

  virtual bool inside_named_section() const;
  virtual const std::string &current_section_name() const;
  virtual void set_current_section_name(std::string name);
  virtual void clear_current_section_name();

  virtual bool has_current_expression_id() const;
  virtual std::size_t current_expression_id() const;
  virtual void set_current_expression_id(std::size_t id);
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
