#ifndef SEQUANT_UTILS_BINARY_EXPR_HPP
#define SEQUANT_UTILS_BINARY_EXPR_HPP

#include <functional>
#include <iostream>
#include <memory>
#include <string>

namespace sequant::utils {

template <typename T>
class binary_expr {
 private:
  T data_;

 protected:
  explicit binary_expr(const T& data) : data_{data} {}

  explicit binary_expr(T&& data) : data_{std::move(data)} {}

 public:
  using type = T;

  using node_ptr = std::unique_ptr<binary_expr<T>>;

  virtual ~binary_expr() = default;

  [[nodiscard]] virtual bool leaf() const = 0;

  virtual const node_ptr& left() const = 0;

  virtual const node_ptr& right() const = 0;

  const T& data() const { return data_; }

  T& data() { return data_; }
};

namespace detail {

template <typename T>
class binary_expr_internal final : public binary_expr<T> {
 private:
  using typename binary_expr<T>::node_ptr;

  node_ptr left_;

  node_ptr right_;

 public:
  binary_expr_internal(T&& data, node_ptr&& left, node_ptr&& right)
      : binary_expr<T>(std::forward<T>(data)),
        left_{std::move(left)},
        right_{std::move(right)} {}

  binary_expr_internal(binary_expr_internal<T> const&) = delete;

  binary_expr_internal<T>& operator=(binary_expr_internal<T> const&) = delete;

  [[nodiscard]] bool leaf() const override { return false; }

  const node_ptr& left() const override { return left_; }

  const node_ptr& right() const override { return right_; }
};

template <typename T>
class binary_expr_leaf final : public binary_expr<T> {
 private:
  using typename binary_expr<T>::node_ptr;

 public:
  explicit binary_expr_leaf(T&& data) : binary_expr<T>(std::forward<T>(data)) {}

  binary_expr_leaf(binary_expr_leaf<T> const&) = delete;

  binary_expr_leaf<T>& operator=(binary_expr_leaf<T> const&) = delete;

  [[nodiscard]] bool leaf() const override { return true; }

  const node_ptr& left() const override {
    throw std::logic_error("left() called on leaf node");
  }

  const node_ptr& right() const override {
    throw std::logic_error("right() called on leaf node");
  }
};

}  // namespace detail

template <typename T, typename V>
bool operator==(const binary_expr<T>& lhs, const binary_expr<V>& rhs) {
  return lhs.data() == rhs.data() &&
         ((lhs.leaf() && rhs.leaf()) ||
          (*lhs.left() == *rhs.left()) && (*lhs.right() == *rhs.right()));
}

template <typename T>
typename binary_expr<T>::node_ptr make_binary_expr(T&& data) {
  return std::make_unique<detail::binary_expr_leaf<T>>(std::forward<T>(data));
}

template <typename T>
typename binary_expr<T>::node_ptr make_binary_expr(
    T&& data, typename binary_expr<T>::node_ptr&& left,
    typename binary_expr<T>::node_ptr&& right) {
  return std::make_unique<detail::binary_expr_internal<T>>(
      std::forward<T>(data), std::move(left), std::move(right));
}

template <typename T>
typename binary_expr<T>::node_ptr make_binary_expr(T&& p, T&& l, T&& r) {
  return make_binary_expr<T>(std::forward<T>(p),
                             make_binary_expr<T>(std::forward<T>(l)),
                             make_binary_expr<T>(std::forward<T>(r)));
}

template <typename T, typename F>
void visit_inorder_binary_expr(typename binary_expr<T>::node_ptr const& node,
                               F&& visitor) {
  static_assert(
      std::is_invocable_v<F, typename binary_expr<T>::node_ptr const&>,
      "visitor signature not matched");
  if (!node) return;  // if nullptr passed explicitly

  if (node->leaf()) {
    visitor(node);
    return;
  }

  // visit left node
  visit_inorder_binary_expr<T, F>(node->left(), std::forward<F>(visitor));

  // visit this node
  visitor(node);

  // visit right node
  visit_inorder_binary_expr<T, F>(node->right(), std::forward<F>(visitor));
}

template <typename T, typename F>
void visit_preorder_binary_expr(typename binary_expr<T>::node_ptr const& node,
                                F&& visitor) {
  static_assert(
      std::is_invocable_v<F, typename binary_expr<T>::node_ptr const&>,
      "visitor signature not matched");
  if (!node) return;  // if nullptr passed explicitly

  if (node->leaf()) {
    visitor(node);
    return;
  }

  // visit this node
  visitor(node);

  // visit left node
  visit_preorder_binary_expr<T, F>(node->left(), std::forward<F>(visitor));

  // visit right node
  visit_preorder_binary_expr<T, F>(node->right(), std::forward<F>(visitor));
}

template <typename T, typename F>
void visit_postorder_binary_expr(typename binary_expr<T>::node_ptr const& node,
                                 F&& visitor) {
  static_assert(
      std::is_invocable_v<F, typename binary_expr<T>::node_ptr const&>,
      "visitor signature not matched");

  if (node->leaf()) {
    visitor(node);
    return;
  }
  // visit left node
  visit_postorder_binary_expr<T, F>(node->left(), std::forward<F>(visitor));

  // visit right node
  visit_postorder_binary_expr<T, F>(node->right(), std::forward<F>(visitor));

  // visit this node
  visitor(node);
}

template <typename T, typename R, typename F>
R evaluate_binary_expr(typename binary_expr<T>::node_ptr const& node,
                       F&& evaluator) {
  static_assert(
      std::is_convertible_v<
          F, std::function<R(typename binary_expr<T>::node_ptr const&)>>,
      "evaluator signature does not match");
  static_assert(std::is_convertible_v<
                    F, std::function<R(typename binary_expr<T>::node_ptr const&,
                                       const R&, const R&)>>,
                "evaluator signature does not match");

  if (!node) throw std::logic_error("Called evaluate on nullptr");

  if (node->leaf()) return evaluator(node);

  return evaluator(
      node,
      evaluate_binary_expr<T, R, F>(node->left(), std::forward<F>(evaluator)),
      evaluate_binary_expr<T, R, F>(node->right(), std::forward<F>(evaluator)));
}

template <typename T, typename R, typename F>
typename binary_expr<R>::node_ptr transform_binary_expr(
    typename binary_expr<T>::node_ptr const& node, F&& transformer) {
  static_assert(std::is_convertible_v<F, std::function<R(const T&)>>,
                "transformer signature does not match");
  if (!node) return nullptr;

  auto parent_node = make_binary_expr<R>(transformer(node->data()));

  if (node->leaf()) return parent_node;

  return make_binary_expr<R>(transformer(node->data()),
                             transform_binary_expr<T, R, F>(
                                 node->left(), std::forward<F>(transformer)),
                             transform_binary_expr<T, R, F>(
                                 node->right(), std::forward<F>(transformer)));
}

template <typename T, typename Os, typename F>
void node_connect(Os& out, typename binary_expr<T>::node_ptr const& node,
                  size_t& node_count, F&& pred) {
  auto label = "node" + std::to_string(node_count);
  out << label.data() << "[label=" << pred(node) << "];\n";

  if (node->leaf()) return;

  out << label.data() << " -> "
      << "node" << ++node_count << ";\n";

  node_connect<T, Os>(out, node->left(), node_count, std::forward<F>(pred));

  out << label.data() << " -> "
      << "node" << ++node_count << ";\n";

  node_connect<T, Os>(out, node->right(), node_count, std::forward<F>(pred));
}

template <typename T, typename Os, typename F>
void digraph_binary_expr(
    Os& out, typename binary_expr<T>::node_ptr const& node,
    F&& label_gen = [](typename binary_expr<T>::node_ptr const&) {
      return "";
    }) {
  static_assert(
      std::is_invocable_v<F, typename binary_expr<T>::node_ptr const&>,
      "label generator signature not matched");

  out << "digraph binary_expr {\n";
  size_t count = 0;
  node_connect<T>(out, node, count, std::forward<F>(label_gen));
  out << "}\n";
}

}  // namespace sequant::utils

#endif  // SEQUANT_UTILS_BINARY_EXPR_HPP
