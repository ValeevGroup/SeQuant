#ifndef SEQUANT_CORE_EXPORT_CONTEXT_HPP
#define SEQUANT_CORE_EXPORT_CONTEXT_HPP

#include <SeQuant/core/export/utils.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/tensor.hpp>

#include <map>

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
  using StrategyMap = std::map<Tensor, StrategyPair, TensorBlockCompare>;

  ExportContext() = default;
  ExportContext(StrategyMap map);
  virtual ~ExportContext();

  virtual LoadStrategy loadStrategy(const Tensor &tensor) const;
  virtual void setLoadStrategy(const Tensor &tensor, LoadStrategy strategy);

  virtual ZeroStrategy zeroStrategy(const Tensor &tensor) const;
  virtual void setZeroStrategy(const Tensor &tensor, ZeroStrategy strategy);

  virtual bool rewrite(Tensor &tensor) const;
  virtual bool rewrite(Variable &variable) const;

 private:
  StrategyMap m_tensorStrategies;
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_CONTEXT_HPP
