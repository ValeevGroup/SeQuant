#ifndef SEQUANT_UTILS_BINARY_EXPR_HPP
#define SEQUANT_UTILS_BINARY_EXPR_HPP

#include <range/v3/numeric/accumulate.hpp>

#include <memory>
#include <string>

namespace sequant::utils {

template <typename T>
class binary_expr {
 private:
  T data_;

 protected:
  explicit binary_expr(T data) : data_{std::move(data)} {}

 public:
  using type = T;

  using node_ptr = std::unique_ptr<binary_expr<T>>;

  virtual ~binary_expr() = default;

  [[nodiscard]] virtual bool leaf() const = 0;

  virtual const node_ptr& left() const = 0;

  virtual const node_ptr& right() const = 0;

  virtual node_ptr clone() const = 0;

  const T& data() const { return data_; }

  T& data() { return data_; }
};

namespace detail {

template <typename T>
class binary_expr_leaf final : public binary_expr<T> {
 private:
  using typename binary_expr<T>::node_ptr;

 public:
  explicit binary_expr_leaf(T data) : binary_expr<T>(std::move(data)) {}

  binary_expr_leaf(binary_expr_leaf<T> const& other)
      : binary_expr<T>{other.data()} {
    static_assert(std::is_copy_constructible<T>::value,
                  "Cannot copy construct binary_expr<T> "
                  "for non-copy-constructible T.");
  }

  binary_expr_leaf<T>& operator=(binary_expr_leaf<T> const& other) {
    static_assert(std::is_copy_assignable<T>::value,
                  "Cannot copy assign binary_expr<T> "
                  "for non-copy-assignable T.");
    return *other.clone();
  }

  binary_expr_leaf(binary_expr_leaf<T>&&) = default;

  binary_expr_leaf<T>& operator=(binary_expr_leaf<T>&&) = default;

  [[nodiscard]] bool leaf() const override { return true; }

  const node_ptr& left() const override {
    throw std::logic_error("left() called on leaf node");
  }

  const node_ptr& right() const override {
    throw std::logic_error("right() called on leaf node");
  }

  node_ptr clone() const override {
    return std::make_unique<binary_expr_leaf<T>>(
        binary_expr_leaf<T>{this->data()});
  }
};

template <typename T>
class binary_expr_internal final : public binary_expr<T> {
 private:
  using typename binary_expr<T>::node_ptr;

  node_ptr left_;

  node_ptr right_;

 public:
  binary_expr_internal(T data, node_ptr&& left, node_ptr&& right)
      : binary_expr<T>(std::move(data)),
        left_{std::move(left)},
        right_{std::move(right)} {}

  binary_expr_internal(binary_expr_internal<T> const& other)
      : binary_expr<T>{other.data()} {
    static_assert(std::is_copy_constructible<T>::value,
                  "Cannot copy construct binary_expr<T> "
                  "for non-copy-constructible T.");
    left_ = std::move(other.left()->clone());
    right_ = std::move(other.right()->clone());
  }

  binary_expr_internal<T>& operator=(binary_expr_internal<T> const& other) {
    static_assert(std::is_copy_assignable<T>::value,
                  "Cannot copy assign binary_expr<T> "
                  "for non-copy-assignable T.");
    return *other.clone();
  }

  binary_expr_internal(binary_expr_internal<T>&&) = default;

  binary_expr_internal<T>& operator=(binary_expr_internal<T>&&) = default;

  [[nodiscard]] bool leaf() const override { return false; }

  const node_ptr& left() const override { return left_; }

  const node_ptr& right() const override { return right_; }

  node_ptr clone() const override {
    return std::make_unique<binary_expr_internal>(
        binary_expr_internal{this->data(), left()->clone(), right()->clone()});
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

template <typename Container, typename F>
auto make_binary_expr(Container const& container, F&& function) {
  static_assert(std::is_invocable_v<F, typename Container::value_type>,
                "Binarizer to handle terminal nodes missing");

  using return_data_t = std::invoke_result_t<F, typename Container::value_type>;

  static_assert(!std::is_same_v<return_data_t, void>,
                "binary_expr<void> not allowed");

  static_assert(std::is_invocable_v<F, return_data_t, return_data_t>,
                "Binarizer to handle non-terminal nodes missing");

  static_assert(
      std::is_same_v<return_data_t,
                     std::invoke_result_t<F, return_data_t, return_data_t>>,
      "function(...) and function(..., ...) have different return types");

  using node_ptr = typename binary_expr<return_data_t>::node_ptr;
  using value_type = typename Container::value_type;

  return ranges::accumulate(
      ranges::begin(container) + 1, ranges::end(container),
      make_binary_expr<return_data_t>(function(*(ranges::begin(container)))),
      [&function](auto&& acc, auto const& val) {
        auto right_node = make_binary_expr<return_data_t>(function(val));
        auto data = function(acc->data(), right_node->data());
        return make_binary_expr<return_data_t>(std::move(data), std::move(acc),
                                               std::move(right_node));
      });
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
auto evaluate_binary_expr(typename binary_expr<T>::node_ptr const& root,
                          F&& evaluator) {
  using node_ptr = typename binary_expr<T>::node_ptr;

  static_assert(std::is_invocable_v<F, node_ptr>,
                "terminal node evaluator not found");

  using ret_type = std::invoke_result_t<F, node_ptr>;

  static_assert(std::is_invocable_v<F, node_ptr, ret_type, ret_type>,
                "non-terminal node evaluator not found");
  static_assert(
      std::is_same_v<ret_type,
                     std::invoke_result_t<F, node_ptr, ret_type, ret_type>>,
      "evaluator(node_ptr) and evaluator(node_ptr, l_result, r_result) have "
      "different return types");

  // if nullptr is passed explicitly
  if (!root) throw std::logic_error("called evaluator on nullptr");

  if (root->leaf()) return evaluator(root);

  return evaluator(
      root,
      evaluate_binary_expr<T, F>(root->left(), std::forward<F>(evaluator)),
      evaluate_binary_expr<T, F>(root->right(), std::forward<F>(evaluator)));
}

namespace detail {

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

}  // namespace detail

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
  detail::node_connect<T>(out, node, count, std::forward<F>(label_gen));
  out << "}\n";
}

}  // namespace sequant::utils

#endif  // SEQUANT_UTILS_BINARY_EXPR_HPP
