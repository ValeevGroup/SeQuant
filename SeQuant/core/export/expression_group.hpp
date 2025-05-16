#ifndef SEQUANT_CORE_EXPORT_EXPRESSION_GROUP_HPP
#define SEQUANT_CORE_EXPORT_EXPRESSION_GROUP_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/export/export_expr.hpp>
#include <SeQuant/core/export/export_node.hpp>

#include <optional>
#include <string>

namespace sequant {

/// Groups multiple expressions together that are supposed to be exported as a
/// single (named) group (e.g. in a function).
template <typename T = ExportExpr>
class ExpressionGroup {
 private:
  using container_type = container::svector<ExportNode<T>>;
  using node_data_type = T;

 public:
  using iterator = container_type::iterator;
  using const_iterator = container_type::const_iterator;

  explicit ExpressionGroup(std::optional<std::string> name = {})
      : ExpressionGroup(container_type{}, std::move(name)) {}
  explicit ExpressionGroup(ExportNode<T> expression)
      : ExpressionGroup(container_type{std::move(expression)}, {}) {}
  ExpressionGroup(ExportNode<T> expression, std::optional<std::string> name)
      : ExpressionGroup(container_type{std::move(expression)},
                        std::move(name)) {}
  explicit ExpressionGroup(container_type expressions)
      : ExpressionGroup(std::move(expressions), {}) {}
  ExpressionGroup(container_type expressions, std::optional<std::string> name)
      : m_expressions(std::move(expressions)), m_name(std::move(name)) {
    assert(!m_name.has_value() || !m_name.value().empty());
  }

  template <typename Range>
    requires std::ranges::range<Range> &&
             std::is_same_v<std::ranges::range_value_t<Range>,
                            typename container_type::value_type>
  explicit ExpressionGroup(const Range &rng,
                           std::optional<std::string> name = {})
      : ExpressionGroup(
            container_type(std::ranges::begin(rng), std::ranges::end(rng)),
            std::move(name)) {}

  bool is_named() const { return m_name.has_value(); }

  const std::string name() const { return m_name.value(); }

  iterator begin() { return m_expressions.begin(); }
  iterator end() { return m_expressions.end(); }
  const_iterator begin() const { return m_expressions.begin(); }
  const_iterator end() const { return m_expressions.end(); }
  const_iterator cbegin() const { return m_expressions.cbegin(); }
  const_iterator cend() const { return m_expressions.cend(); }

  std::size_t size() const { return m_expressions.size(); }

  void add(ExportNode<T> expr) { m_expressions.emplace_back(std::move(expr)); }

 private:
  container::svector<ExportNode<T>> m_expressions;
  std::optional<std::string> m_name;
};

}  // namespace sequant

#endif
