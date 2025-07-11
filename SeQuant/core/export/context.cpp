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

UsageSet &UsageSet::operator=(Usage usage) {
  m_val = static_cast<Value>(usage);
  return *this;
}

UsageSet &UsageSet::operator|=(UsageSet usage) {
  m_val |= usage.m_val;
  return *this;
}

UsageSet &UsageSet::operator&=(UsageSet usage) {
  m_val &= usage.m_val;
  return *this;
}

UsageSet &UsageSet::operator|=(Usage usage) {
  m_val |= static_cast<Value>(usage);
  return *this;
}

UsageSet &UsageSet::operator&=(Usage usage) {
  m_val &= static_cast<Value>(usage);
  return *this;
}

bool operator&(UsageSet set, Usage usage) {
  return set.m_val & static_cast<UsageSet::Value>(usage);
}

bool operator&(Usage usage, UsageSet set) { return set & usage; }

bool operator==(UsageSet set, Usage usage) {
  return set.m_val == static_cast<UsageSet::Value>(usage);
}

bool operator==(Usage usage, UsageSet set) { return set == usage; }

bool operator!=(UsageSet set, Usage usage) { return !(set == usage); }

bool operator!=(Usage usage, UsageSet set) { return set != usage; }

ExportContext::ExportContext(TensorStrategyMap tensorMap,
                             VariableStrategyMap variableMap)
    : m_tensorStrategies({{GLOBAL, std::move(tensorMap)}}),
      m_variableStrategies({{GLOBAL, std::move(variableMap)}}) {}

ExportContext::ExportContext(VariableStrategyMap map)
    : m_variableStrategies({{GLOBAL, std::move(map)}}) {}

ExportContext::~ExportContext() = default;

LoadStrategy ExportContext::loadStrategy(const Tensor &tensor) const {
  if (auto map_iter = m_tensorStrategies.find(GLOBAL);
      map_iter != m_tensorStrategies.end()) {
    auto iter = map_iter->second.find(tensor);
    if (iter != map_iter->second.end()) {
      return iter->second.load;
    }
  }

  if (has_current_expression_id()) {
    if (auto map_iter = m_tensorStrategies.find(current_expression_id());
        map_iter != m_tensorStrategies.end()) {
      auto iter = map_iter->second.find(tensor);
      if (iter != map_iter->second.end()) {
        return iter->second.load;
      }
    }
  }

  return LoadStrategy::Create;
}

LoadStrategy ExportContext::loadStrategy(const Variable &variable) const {
  if (auto map_iter = m_variableStrategies.find(GLOBAL);
      map_iter != m_variableStrategies.end()) {
    auto iter = map_iter->second.find(variable);
    if (iter != map_iter->second.end()) {
      return iter->second.load;
    }
  }

  if (has_current_expression_id()) {
    if (auto map_iter = m_variableStrategies.find(current_expression_id());
        map_iter != m_variableStrategies.end()) {
      auto iter = map_iter->second.find(variable);
      if (iter != map_iter->second.end()) {
        return iter->second.load;
      }
    }
  }

  return LoadStrategy::Create;
}

void ExportContext::setLoadStrategy(
    const Tensor &tensor, LoadStrategy strategy,
    const std::optional<std::size_t> &expression_id) {
  std::size_t id = expression_id.value_or(
      has_current_expression_id() ? current_expression_id() : GLOBAL);

  auto iter = m_tensorStrategies[id].find(tensor);

  if (iter == m_tensorStrategies[id].end()) {
    m_tensorStrategies[id][tensor] = StrategyPair(strategy);
  } else {
    iter->second.load = strategy;
  }
}

void ExportContext::setLoadStrategy(
    const Variable &variable, LoadStrategy strategy,
    const std::optional<std::size_t> &expression_id) {
  std::size_t id = expression_id.value_or(
      has_current_expression_id() ? current_expression_id() : GLOBAL);

  auto iter = m_variableStrategies[id].find(variable);

  if (iter == m_variableStrategies[id].end()) {
    m_variableStrategies[id][variable] = StrategyPair(strategy);
  } else {
    iter->second.load = strategy;
  }
}

ZeroStrategy ExportContext::zeroStrategy(const Tensor &tensor) const {
  if (auto map_iter = m_tensorStrategies.find(GLOBAL);
      map_iter != m_tensorStrategies.end()) {
    auto iter = map_iter->second.find(tensor);
    if (iter != map_iter->second.end()) {
      return iter->second.zero;
    }
  }

  if (has_current_expression_id()) {
    if (auto map_iter = m_tensorStrategies.find(current_expression_id());
        map_iter != m_tensorStrategies.end()) {
      auto iter = map_iter->second.find(tensor);
      if (iter != map_iter->second.end()) {
        return iter->second.zero;
      }
    }
  }

  return ZeroStrategy::ZeroOnCreate;
}

ZeroStrategy ExportContext::zeroStrategy(const Variable &variable) const {
  if (auto map_iter = m_variableStrategies.find(GLOBAL);
      map_iter != m_variableStrategies.end()) {
    auto iter = map_iter->second.find(variable);
    if (iter != map_iter->second.end()) {
      return iter->second.zero;
    }
  }

  if (has_current_expression_id()) {
    if (auto map_iter = m_variableStrategies.find(current_expression_id());
        map_iter != m_variableStrategies.end()) {
      auto iter = map_iter->second.find(variable);
      if (iter != map_iter->second.end()) {
        return iter->second.zero;
      }
    }
  }

  return ZeroStrategy::ZeroOnCreate;
}

void ExportContext::setZeroStrategy(
    const Tensor &tensor, ZeroStrategy strategy,
    const std::optional<std::size_t> &expression_id) {
  std::size_t id = expression_id.value_or(
      has_current_expression_id() ? current_expression_id() : GLOBAL);

  auto iter = m_tensorStrategies[id].find(tensor);

  if (iter == m_tensorStrategies[id].end()) {
    m_tensorStrategies[id][tensor] = StrategyPair(strategy);
  } else {
    iter->second.zero = strategy;
  }
}

void ExportContext::setZeroStrategy(
    const Variable &variable, ZeroStrategy strategy,
    const std::optional<std::size_t> &expression_id) {
  std::size_t id = expression_id.value_or(
      has_current_expression_id() ? current_expression_id() : GLOBAL);

  auto iter = m_variableStrategies[id].find(variable);

  if (iter == m_variableStrategies[id].end()) {
    m_variableStrategies[id][variable] = StrategyPair(strategy);
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

bool ExportContext::has_current_expression_id() const {
  return m_currentExpressionID.has_value();
}

std::size_t ExportContext::current_expression_id() const {
  return m_currentExpressionID.value();
}

void ExportContext::set_current_expression_id(std::size_t id) {
  m_currentExpressionID = id;
}

void ExportContext::clear_current_expression_id() {
  m_currentExpressionID.reset();
}

}  // namespace sequant
