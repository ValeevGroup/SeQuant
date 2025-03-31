#ifndef SEQUANT_CORE_EXPORT_CONTEXT_HPP
#define SEQUANT_CORE_EXPORT_CONTEXT_HPP

#include <SeQuant/core/export/utils.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/tensor.hpp>

#include <map>
#include <optional>
#include <string>

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
  virtual void setLoadStrategy(const Tensor &tensor, LoadStrategy strategy);
  virtual void setLoadStrategy(const Variable &variable, LoadStrategy strategy);

  virtual ZeroStrategy zeroStrategy(const Tensor &tensor) const;
  virtual ZeroStrategy zeroStrategy(const Variable &variable) const;
  virtual void setZeroStrategy(const Tensor &tensor, ZeroStrategy strategy);
  virtual void setZeroStrategy(const Variable &variable, ZeroStrategy strategy);

  virtual bool rewrite(Tensor &tensor) const;
  virtual bool rewrite(Variable &variable) const;

  virtual bool inside_named_section() const;
  virtual const std::string &current_section_name() const;
  virtual void set_current_section_name(std::string name);
  virtual void clear_current_section_name();

 private:
  TensorStrategyMap m_tensorStrategies;
  VariableStrategyMap m_variableStrategies;
  std::optional<std::string> m_currentSection;
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_CONTEXT_HPP
