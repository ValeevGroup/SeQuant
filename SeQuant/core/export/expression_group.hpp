#ifndef SEQUANT_CORE_EXPORT_EXPRESSION_GROUP_HPP
#define SEQUANT_CORE_EXPORT_EXPRESSION_GROUP_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval_node.hpp>

#include <optional>
#include <string>

namespace sequant {

/// Groups multiple expressions together that are supposed to be exported as a
/// single (named) group (e.g. in a function).
template <typename T>
class ExpressionGroup {
 private:
  using container_type = container::svector<EvalNode<T>>;

 public:
  using iterator = container_type::iterator;
  using const_iterator = container_type::const_iterator;

  explicit ExpressionGroup(EvalNode<T> expression)
      : ExpressionGroup(container_type{std::move(expression)}, {}) {}
  ExpressionGroup(EvalNode<T> expression, std::optional<std::string> name)
      : ExpressionGroup(container_type{std::move(expression)},
                        std::move(name)) {}
  explicit ExpressionGroup(container_type expressions)
      : ExpressionGroup(std::move(expressions), {}) {}
  ExpressionGroup(container_type expressions, std::optional<std::string> name)
      : m_expressions(std::move(expressions)), m_name(std::move(name)) {
    assert(!m_name.has_value() || !m_name.value().empty());
  }

  bool is_named() const { return m_name.has_value(); }

  const std::string name() const { return m_name.value(); }

  iterator begin() { return m_expressions.begin(); }
  iterator end() { return m_expressions.end(); }
  const_iterator begin() const { return m_expressions.begin(); }
  const_iterator end() const { return m_expressions.end(); }
  const_iterator cbegin() const { return m_expressions.cbegin(); }
  const_iterator cend() const { return m_expressions.cend(); }

 private:
  container::svector<EvalNode<T>> m_expressions;
  std::optional<std::string> m_name;
};

}  // namespace sequant

#endif
