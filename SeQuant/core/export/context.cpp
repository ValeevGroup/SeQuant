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

ExportContext::ExportContext(TensorStrategyMap tensorMap,
                             VariableStrategyMap variableMap)
    : m_tensorStrategies(std::move(tensorMap)),
      m_variableStrategies(std::move(variableMap)) {}

ExportContext::ExportContext(VariableStrategyMap map)
    : m_variableStrategies(std::move(map)) {}

ExportContext::~ExportContext() = default;

LoadStrategy ExportContext::loadStrategy(const Tensor &tensor) const {
  auto iter = m_tensorStrategies.find(tensor);
  return iter == m_tensorStrategies.end() ? LoadStrategy::Create
                                          : iter->second.load;
}

LoadStrategy ExportContext::loadStrategy(const Variable &variable) const {
  auto iter = m_variableStrategies.find(variable);
  return iter == m_variableStrategies.end() ? LoadStrategy::Create
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

void ExportContext::setLoadStrategy(const Variable &variable,
                                    LoadStrategy strategy) {
  auto iter = m_variableStrategies.find(variable);

  if (iter == m_variableStrategies.end()) {
    m_variableStrategies[variable] = StrategyPair(strategy);
  } else {
    iter->second.load = strategy;
  }
}

ZeroStrategy ExportContext::zeroStrategy(const Tensor &tensor) const {
  auto iter = m_tensorStrategies.find(tensor);
  return iter == m_tensorStrategies.end() ? ZeroStrategy::ZeroOnCreate
                                          : iter->second.zero;
}

ZeroStrategy ExportContext::zeroStrategy(const Variable &variable) const {
  auto iter = m_variableStrategies.find(variable);
  return iter == m_variableStrategies.end() ? ZeroStrategy::ZeroOnCreate
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

void ExportContext::setZeroStrategy(const Variable &variable,
                                    ZeroStrategy strategy) {
  auto iter = m_variableStrategies.find(variable);

  if (iter == m_variableStrategies.end()) {
    m_variableStrategies[variable] = StrategyPair(strategy);
  } else {
    iter->second.zero = strategy;
  }
}

bool ExportContext::rewrite(Tensor &tensor) const { return false; }

bool ExportContext::rewrite(Variable &variable) const { return false; }

bool ExportContext::inside_named_section() const {
  return m_currentSection.has_value();
}

const std::string &ExportContext::current_section_name() const {
  return m_currentSection.value();
}

void ExportContext::set_current_section_name(std::string name) {
  m_currentSection = std::move(name);
}

void ExportContext::clear_current_section_name() { m_currentSection.reset(); }

}  // namespace sequant
