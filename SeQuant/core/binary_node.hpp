#ifndef SEQUANT_BINARY_NODE_HPP
#define SEQUANT_BINARY_NODE_HPP

#include <iostream>
#include <memory>
#include <range/v3/numeric/accumulate.hpp>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace sequant {

template <typename>
class FullBinaryNode;

/// Full-binary node visit orders.

/// Visit children nodes followed by parent node.
/// Left will be visited before right.
struct PreOrder {};

/// Visit parent node followed by children nodes.
/// Left will be visited before right.
struct PostOrder {};

/// Visit left node then parent node followed by right node.
struct InOrder {};

namespace {

/// Visit leaf nodes only.
struct VisitLeaf {};

/// Visit internal nodes only.
struct VisitInternal {};

/// Visit both leaf and internal nodes.
struct VisitAll {};

///
/// \brief Visit a full binary node.
///
/// \tparam V Visitor type. Must be invocable with a FullBinaryNode<T> argument.
/// \tparam Order Visit order.
/// \param node Node to visit.
/// \param f Visitor.
/// \param lf If true, only leaf nodes will be visited.
/// \param in If true, parent node will be visited between left and right.
///
template <
    typename T, typename V, typename Order, typename NodeType,
    typename = std::enable_if_t<std::is_invocable_v<V, FullBinaryNode<T>>>>
void visit(FullBinaryNode<T> const& node, V const& f, Order, NodeType) {
  static_assert(std::is_same_v<Order, PreOrder> ||
                std::is_same_v<Order, InOrder> ||
                std::is_same_v<Order, PostOrder> && "Unsupported visit order");
  static_assert(std::is_same_v<NodeType, VisitLeaf> ||
                std::is_same_v<NodeType, VisitInternal> ||
                std::is_same_v<NodeType, VisitAll> &&
                    "Not sure which nodes to visit");
  if (node.leaf()) {
    if constexpr (!std::is_same_v<NodeType, VisitInternal>) f(node);
  } else {
    if constexpr (std::is_same_v<Order, PreOrder> &&
                  !std::is_same_v<NodeType, VisitLeaf>)
      f(node);

    visit(node.left(), f, Order{}, NodeType{});

    if constexpr (std::is_same_v<Order, InOrder> &&
                  !std::is_same_v<NodeType, VisitLeaf>)
      f(node);

    visit(node.right(), f, Order{}, NodeType{});

    if constexpr (std::is_same_v<Order, PostOrder> &&
                  !std::is_same_v<NodeType, VisitLeaf>)
      f(node);
  }
}

}  // namespace

///
/// @brief Represents a node with data of @c T type in a full-binary tree.
///
/// A full binary tree is a binary tree in which each node has two children
/// or no children.
///
template <typename T>
class FullBinaryNode {
 public:
  using node_ptr = std::unique_ptr<FullBinaryNode<T>>;
  using value_type = T;

 private:
  T data_;

  node_ptr left_{nullptr};

  node_ptr right_{nullptr};

  node_ptr deep_copy() const {
    return leaf() ? std::make_unique<FullBinaryNode<T>>(data_)
                  : std::make_unique<FullBinaryNode<T>>(
                        data_, left_->deep_copy(), right_->deep_copy());
  }

  static node_ptr const& checked_ptr_access(node_ptr const& n) {
    if (n)
      return n;
    else
      throw std::runtime_error(
          "Dereferenced nullptr: use leaf() method to check leaf node.");
  }

 public:
  ///
  /// Construct an internal node with emtpy left and right nodes.
  ///
  /// \param d Data in the internal node.
  FullBinaryNode(T d) : data_{std::move(d)} {}

  ///
  /// Construct an internal node with left and right node data.
  ///
  /// \param d Data in the internal node.
  /// \param l Data in the left node.
  /// \param r Data in the right node.
  FullBinaryNode(T d, T l, T r)
      : data_{std::move(d)},
        left_{std::make_unique<FullBinaryNode>(std::move(l))},
        right_{std::make_unique<FullBinaryNode>(std::move(r))} {}

  ///
  /// Constructs an internal node with left and right nodes.
  ///
  /// \param d Data in the internal node.
  /// \param l Left node.
  /// \param r Right node
  FullBinaryNode(T d, FullBinaryNode<T> l, FullBinaryNode<T> r)
      : data_{std::move(d)},
        left_{std::make_unique<FullBinaryNode<T>>(std::move(l))},
        right_{std::make_unique<FullBinaryNode<T>>(std::move(r))} {}

  ///
  /// Constructs an internal node with left and right node pointers.
  ///
  /// \param d Data in the internal node.
  /// \param l Left node pointer.
  /// \param r Right node pointer.
  FullBinaryNode(T d, node_ptr&& l, node_ptr&& r)
      : data_{std::move(d)}, left_{std::move(l)}, right_{std::move(r)} {}

  FullBinaryNode(FullBinaryNode<T> const& other)
      : data_{other.data_},
        left_{other.left_ ? other.left_->deep_copy() : nullptr},
        right_{other.right_ ? other.right_->deep_copy() : nullptr} {}

  FullBinaryNode& operator=(FullBinaryNode<T> const& other) {
    auto temp = other.deep_copy();
    data_ = std::move(temp->data_);
    left_ = std::move(temp->left_);
    right_ = std::move(temp->right_);
    return *this;
  }

  FullBinaryNode(FullBinaryNode<T>&&) = default;

  FullBinaryNode& operator=(FullBinaryNode<T>&&) = default;

  FullBinaryNode const& left() const { return *checked_ptr_access(left_); }

  FullBinaryNode const& right() const { return *checked_ptr_access(right_); }

  bool leaf() const { return !(left_ || right_); }

  ///
  /// \return Returns the data stored by the node.
  T const& operator*() const { return data_; }

  ///
  /// \return Returns the data stored by the node.
  T& operator*() { return data_; }

  ///
  /// \return Returns the pointer to the data stored by the node.
  T const* operator->() const { return &data_; }

  ///
  /// \return Returns the pointer to the data stored by the node.
  T* operator->() { return &data_; }

  ///
  /// Left-fold a container to make a full-binary node.
  ///
  /// \param container To be binarized.
  /// \param binarize Fold function.
  ///        \c binarize needs to support:
  ///          - unary function call with a return value (say of type R)
  ///             to the element type of the container (say of type V)
  ///          - binary function call of kind f(R,V) that returns R type
  ///
  template <typename Cont, typename F>
  FullBinaryNode(Cont const& container, F&& binarize) {
    using value_type = decltype(*ranges::begin(container));
    static_assert(std::is_invocable_v<F, value_type const&>,
                  "Binarizer to handle terminal nodes missing");

    using return_data_t = std::invoke_result_t<F, value_type const&>;

    static_assert(
        std::is_invocable_v<F, return_data_t const&, return_data_t const&>,
        "Binarizer to handle non-terminal nodes missing");

    static_assert(
        std::is_same_v<return_data_t,
                       std::invoke_result_t<F, return_data_t const&,
                                            return_data_t const&>>,
        "function(...) and function(..., ...) have different return types");

    using ranges::accumulate;
    using ranges::begin;
    using ranges::end;

    auto node =
        accumulate(begin(container) + 1, end(container),         // range
                   FullBinaryNode{binarize(*begin(container))},  // init
                   [&binarize](auto&& acc, const auto& val) {    // predicate
                     auto rnode = FullBinaryNode{binarize(val)};
                     return FullBinaryNode{binarize(*acc, *rnode),
                                           std::move(acc), std::move(rnode)};
                   });

    *this = std::move(node);
  }

  ///
  /// \brief Visit the tree in the order specified by the Order argument.
  /// \tparam F Type of the visitor.
  /// \tparam Order Type of the order. Can be PreOrder, InOrder or PostOrder.
  ///               By default it is PostOrder. That means parent node will
  ///               be visited after visiting left and right children in that
  ///               order.
  /// \param visitor Visitor to be invoked on each node.
  ///
  template <
      typename F, typename Order = PostOrder,
      std::enable_if_t<
          std::is_void_v<std::invoke_result_t<F, FullBinaryNode<T> const&>>,
          bool> = true>
  void visit(F const& visitor, Order = {}) const {
    sequant::visit(*this,    //
                   visitor,  //
                   Order{},  //
                   VisitAll{});
  }

  ///
  /// \tparam F Type of the visitor.
  /// \param visitor Visitor to be invoked on each node.
  /// \brief Visit the children nodes only if the visitor returns true upon
  ///        visiting the parent node.
  /// \details This is a pre-order traversal with short-circuit behavior.
  ///
  template <
      typename F,
      std::enable_if_t<std::is_invocable_r_v<bool, F, FullBinaryNode<T> const&>,
                       bool> = true>
  void visit(F const& visitor) const {
    if (visitor(*this) && !leaf()) {
      left().visit(visitor);
      right().visit(visitor);
    }
  }

  ///
  /// \brief Visit the internal nodes of the tree in the order specified by the
  ///        Order argument.
  /// \tparam F Type of the visitor.
  /// \tparam Order Type of the order. Can be PreOrder, InOrder or PostOrder.
  ///               By default it is PostOrder. That means parent node will
  ///               be visited after visiting left and right children in that
  ///               order.
  /// \param visitor Visitor to be invoked on each node.
  ///
  template <typename F, typename Order = PostOrder,
            std::enable_if_t<std::is_invocable_v<F, FullBinaryNode<T> const&>,
                             bool> = true>
  void visit_internal(F const& visitor, Order = {}) const {
    sequant::visit(*this,    //
                   visitor,  //
                   Order{},  //
                   VisitInternal{});
  }

  ///
  /// \brief Visit the leaf nodes of the tree in the order specified by the
  ///        Order argument.
  /// \tparam F Type of the visitor.
  /// \tparam Order Type of the order. Can be PreOrder, InOrder or PostOrder.
  ///               By default it is PostOrder. That means parent node will
  ///               be visited after visiting left and right children in that
  ///               order.
  /// \param visitor Visitor to be invoked on each node.
  ///
  template <typename F, typename Order = PostOrder,
            std::enable_if_t<std::is_invocable_v<F, FullBinaryNode<T> const&>,
                             bool> = true>
  void visit_leaf(F const& visitor, Order = {}) const {
    sequant::visit(*this,    //
                   visitor,  //
                   Order{},  //
                   VisitLeaf{});
  }

  template <
      typename F,
      std::enable_if_t<
          std::is_invocable_v<F, FullBinaryNode<T> const&> &&
              std::is_invocable_v<
                  F, FullBinaryNode<T> const&,
                  std::invoke_result_t<F, FullBinaryNode<T> const&> const&,
                  std::invoke_result_t<F, FullBinaryNode<T> const&> const&>,
          bool> = true>
  auto evaluate(F const& func) const {
    if (leaf()) return func(*this);
    return func(*this, left().evaluate(func), right().evaluate(func));
  }

 private:
  template <typename Ostream, typename F>
  [[maybe_unused]] int digraph(Ostream& os, F const& label_gen,
                               int count = 0) const {
    os << "node" << count << "[label=" << label_gen(*this) << "];\n";

    if (this->leaf()) return count;

    auto lcount = left().digraph(os, label_gen, count + 1);
    auto rcount = right().digraph(os, label_gen, lcount + 1);
    os << "node" << count << " -> "
       << "node" << count + 1 << ";\n";
    os << "node" << count << " -> "
       << "node" << lcount + 1 << ";\n";

    return rcount;
  }

  template <typename Ostream, typename F, typename G>
  void tikz(Ostream& os, F const& label_gen, G const& spec_gen,
            size_t indent = 2) const {
    auto pad = [](Ostream& o, size_t i) {
      for (auto j = 0; j < i; ++j) o << " ";
    };

    // pad(os, indent);

    os << "node [" << spec_gen(*this) << "]"
       << " {" << label_gen(*this) << "}";
    if (leaf()) return;
    os << "\n";

    pad(os, indent);
    os << "child {";
    left().tikz(os, label_gen, spec_gen, indent + 2);
    os << "}";
    os << "\n";

    pad(os, indent);
    os << "child {";
    right().tikz(os, label_gen, spec_gen, indent + 2);
    os << "}";
  }

 public:
  template <typename string_t, typename F>
  string_t digraph(F const& label_gen, string_t const& graph_name = {}) const {
    static_assert(std::is_invocable_r_v<string_t, F, FullBinaryNode<T> const&>,
                  "node label generator F(FullBinaryNode<T> const &) should "
                  "return string_t");

    auto oss = std::basic_ostringstream{string_t{}};

    oss << "digraph " << graph_name << "{\n";
    this->digraph(oss, label_gen, 0);
    oss << "}";
    oss.flush();

    return oss.str();
  }

  template <typename string_t>
  string_t tikz(
      std::function<string_t(FullBinaryNode<T> const&)> label_gen,
      std::function<string_t(FullBinaryNode<T> const&)> spec_gen) const {
    auto oss = std::basic_ostringstream{string_t{}};
    oss << "\\tikz{\n\\";
    tikz(oss, label_gen, spec_gen);
    oss << "\n}";
    oss.flush();
    return oss.str();
  }

};  // FullBinaryNode<T>

template <typename T, typename U>
bool operator==(FullBinaryNode<T> const& lhs, FullBinaryNode<U> const& rhs) {
  return ((*lhs == *rhs) &&
          ((lhs.leaf() && rhs.leaf()) ||
           (lhs.left() == rhs.left() && lhs.right() == rhs.right())));
}

}  // namespace sequant

#endif  // SEQUANT_BINARY_NODE_HPP
