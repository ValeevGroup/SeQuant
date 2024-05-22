#include <SeQuant/core/export/context.hpp>

#include <type_traits>

namespace sequant {

ZeroStrategy operator|(ZeroStrategy lhs, ZeroStrategy rhs) {
  using underlying_type = std::underlying_type_t<ZeroStrategy>;
  return static_cast<ZeroStrategy>(static_cast<underlying_type>(lhs) |
                                   static_cast<underlying_type>(rhs));
}

ZeroStrategy operator&(ZeroStrategy lhs, ZeroStrategy rhs) {
  using underlying_type = std::underlying_type_t<ZeroStrategy>;
  return static_cast<ZeroStrategy>(static_cast<underlying_type>(lhs) &
                                   static_cast<underlying_type>(rhs));
}

bool operator!(ZeroStrategy strategy) {
  return !static_cast<std::underlying_type_t<ZeroStrategy>>(strategy);
}

StrategyPair::StrategyPair(LoadStrategy load, ZeroStrategy zero)
    : load(load), zero(zero) {}

StrategyPair::StrategyPair(LoadStrategy load) : StrategyPair() {
  this->load = load;
}

StrategyPair::StrategyPair(ZeroStrategy zero) : StrategyPair() {
  this->zero = zero;
}

ExportContext::ExportContext(StrategyMap map)
    : m_tensorStrategies(std::move(map)) {}

ExportContext::~ExportContext() = default;

LoadStrategy ExportContext::loadStrategy(const Tensor &tensor) const {
  auto iter = m_tensorStrategies.find(tensor);
  return iter == m_tensorStrategies.end() ? LoadStrategy::Create
                                          : iter->second.load;
}

void ExportContext::setLoadStrategy(const Tensor &tensor,
                                    LoadStrategy strategy) {
  auto iter = m_tensorStrategies.find(tensor);

  if (iter == m_tensorStrategies.end()) {
    m_tensorStrategies[tensor] = StrategyPair(strategy);
  } else {
    iter->second.load = strategy;
  }
}

ZeroStrategy ExportContext::zeroStrategy(const Tensor &tensor) const {
  auto iter = m_tensorStrategies.find(tensor);
  return iter == m_tensorStrategies.end() ? ZeroStrategy::ZeroOnCreate
                                          : iter->second.zero;
}

void ExportContext::setZeroStrategy(const Tensor &tensor,
                                    ZeroStrategy strategy) {
  auto iter = m_tensorStrategies.find(tensor);

  if (iter == m_tensorStrategies.end()) {
    m_tensorStrategies[tensor] = StrategyPair(strategy);
  } else {
    iter->second.zero = strategy;
  }
}

bool ExportContext::rewrite(Tensor &tensor) const { return false; }

bool ExportContext::rewrite(Variable &variable) const { return false; }

}  // namespace sequant
