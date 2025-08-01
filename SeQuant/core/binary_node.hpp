#ifndef SEQUANT_BINARY_NODE_HPP
#define SEQUANT_BINARY_NODE_HPP

#include <iostream>
#include <memory>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/view.hpp>

namespace sequant {

template <typename>
class FullBinaryNode;

namespace meta {

template <typename>
constexpr bool is_full_binary_node{false};

template <typename T>
constexpr bool is_full_binary_node<FullBinaryNode<T>>{true};

}  // namespace meta

/// Full-binary node visit orders.

/// Visit parent node followed by children nodes.
/// Left will be visited before right.
struct PreOrder {};

/// Visit children nodes followed by parent node.
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
/// \tparam NodeType The type of nodes to be visited.
///                  Can be VisitInternal, VisitLeaf, or VisitAll.
/// \param node Node to visit.
/// \param f Visitor.
///
template <
    typename T, typename V, typename Order, typename NodeType,
    typename = std::enable_if_t<std::is_invocable_v<V, FullBinaryNode<T>>>>
void visit(FullBinaryNode<T> const& node, V f, Order, NodeType) {
  static_assert(std::is_same_v<Order, PreOrder> ||
                    std::is_same_v<Order, InOrder> ||
                    std::is_same_v<Order, PostOrder>,
                "Unsupported visit order");
  static_assert(std::is_same_v<NodeType, VisitLeaf> ||
                    std::is_same_v<NodeType, VisitInternal> ||
                    std::is_same_v<NodeType, VisitAll>,
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
  explicit FullBinaryNode(T d) : data_{std::move(d)} {}

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

  ///
  /// \return Left node if this is an internal node, throws otherwise.
  ///
  FullBinaryNode const& left() const { return *checked_ptr_access(left_); }

  ///
  /// \return Right node if this is an internal node, throws otherwise.
  ///
  FullBinaryNode const& right() const { return *checked_ptr_access(right_); }

  ///
  /// \brief Check if the object is a leaf node.
  ///
  /// \return True if this object represents a terminal binary node.
  /// \note Left and right children are nullptr like
  ///
  [[nodiscard]] bool leaf() const { return !(left_ || right_); }

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
  void visit(F visitor, Order = {}) const {
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
  void visit(F visitor) const {
    if (visitor(*this) && !leaf()) {
      left().visit(visitor);
      right().visit(visitor);
    }
  }

  ///
  /// \tparam F Type of the visitor.
  /// \param visitor Visitor to be invoked on each node.
  /// \brief Visit the children nodes only if the visitor returns true upon
  ///        visiting the parent node, and current node is not a leaf node.
  /// \details This is a pre-order traversal with short-circuit behavior.
  ///
  template <
      typename F,
      std::enable_if_t<std::is_invocable_r_v<bool, F, FullBinaryNode<T> const&>,
                       bool> = true>
  void visit_internal(F visitor) const {
    if (leaf()) return;
    if (visitor(*this)) {
      left().visit_internal(visitor);
      right().visit_internal(visitor);
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
  template <
      typename F, typename Order = PostOrder,
      std::enable_if_t<
          std::is_void_v<std::invoke_result_t<F, FullBinaryNode<T> const&>>,
          bool> = true>
  void visit_internal(F visitor, Order = {}) const {
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
  template <
      typename F, typename Order = PostOrder,
      std::enable_if_t<
          std::is_void_v<std::invoke_result_t<F, FullBinaryNode<T> const&>>,
          bool> = true>
  void visit_leaf(F visitor, Order = {}) const {
    sequant::visit(*this,    //
                   visitor,  //
                   Order{},  //
                   VisitLeaf{});
  }

 private:
  template <typename F, typename Os>
  [[maybe_unused]] int digraph(F label_gen, Os& os, int count) const {
    os << "node" << count << "[label=" << label_gen(*this) << "];\n";
    if (this->leaf()) return count;
    auto lcount = left().digraph(label_gen, os, count + 1);
    auto rcount = right().digraph(label_gen, os, lcount + 1);
    os << "node" << count << " -> "
       << "node" << count + 1 << ";\n";
    os << "node" << count << " -> "
       << "node" << lcount + 1 << ";\n";

    return rcount;
  }

  template <typename Ostream, typename F, typename G>
  void tikz(Ostream& os, F label_gen, G spec_gen, size_t indent = 2) const {
    auto pad = [](Ostream& o, size_t i) {
      for (size_t j = 0; j < i; ++j) o << " ";
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
  template <typename F,
            typename String = std::invoke_result_t<F, FullBinaryNode>>
  String digraph(F label_gen,
                 std::basic_string_view<typename String::value_type>
                     graph_name = {}) const {
    auto oss = std::basic_ostringstream{String{}};

    oss << "digraph " << graph_name << "{\n";
    this->digraph(label_gen, oss, 0);
    oss << "}";
    oss.flush();

    return oss.str();
  }

  ///
  /// \param label_gen Generates the label for tikz nodes.
  /// \param spec_gen Generates the node spec that goes into the square
  ///                 bracket of tikz node statement.
  /// \note Make sure the following are present in the preamble.
  ///
  ///                  @c \usepackage{tikz}
  ///                  @c \usetikzlibrary{graphs,graphdrawing}
  ///                  @c \usegdlibrary{trees}
  ///
  /// \return TikZ graph.

  template <typename L,
            typename String = std::invoke_result_t<L, FullBinaryNode>,
            typename S = std::function<String(FullBinaryNode)>>
  String tikz(
      L label_gen, S spec_gen = [](auto&&) { return String{}; }) const {
    auto oss = std::basic_ostringstream{String{}};
    oss << "\\tikz[binary tree layout]{\n\\";
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

///
/// \return A new binary node where each node value is the result of applying @c
///         fun on the argument node value.
///
template <typename T, typename F,
          typename = std::enable_if_t<std::is_invocable_v<F, T>>>
auto transform_node(FullBinaryNode<T> const& node, F fun) {
  if (node.leaf())
    return FullBinaryNode(fun(*node));
  else {
    return FullBinaryNode(fun(*node), transform_node(node.left(), fun),
                          transform_node(node.right(), fun));
  }
}

///
/// \brief Accumulates the given range of binary nodes into a single binary node
///        using a given operation to generate the internal node values.
///
/// \tparam Node The FullBinaryNode type.
/// \param rng A range of Node objects.
/// \param op A binary function that returns Node::value_type.
///           The signature could be `op(Node, Node) -> Node::value_type` or
///           `op(Node::value_type, Node::value_type) -> Node::value_type`.
/// \return A binary node with subtrees built by left-folding given range.
///
template <typename Rng,                                //
          typename F,                                  //
          typename Node = ranges::range_value_t<Rng>,  //
          typename = std::enable_if_t<meta::is_full_binary_node<Node>>>
Node fold_left_to_node(Rng rng, F op) {
  using Value = typename Node::value_type;

  constexpr bool invoke_on_value =  //
      std::is_invocable_r_v<Value, F, Value, Value>;

  constexpr bool invoke_on_node =  //
      std::is_invocable_r_v<Value, F, Node, Node>;

  static_assert(invoke_on_value || invoke_on_node);

  using ranges::views::move;
  using ranges::views::tail;
  return ranges::accumulate(
      rng | tail | move, std::move(ranges::front(rng)),
      [&op](auto&& l, auto&& r) {
        if constexpr (invoke_on_node) {
          auto&& val = op(l, r);
          return FullBinaryNode(std::move(val), std::move(l), std::move(r));

        } else {
          auto&& val = op(*l, *r);
          return FullBinaryNode(std::move(val), std::move(l), std::move(r));
        }
      });
}

}  // namespace sequant

#endif  // SEQUANT_BINARY_NODE_HPP
