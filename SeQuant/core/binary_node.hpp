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
enum class TreeTraversal {
  None = 0,
  PreOrder = 0b001,
  PostOrder = 0b010,
  InOrder = 0b100,

  PreAndPostOrder = PreOrder | PostOrder,
  PreAndInOrder = PreOrder | InOrder,
  PostAndInOrder = PostOrder | InOrder,

  Any = PreOrder | PostOrder | InOrder,
};

constexpr TreeTraversal operator|(TreeTraversal lhs, TreeTraversal rhs) {
  return static_cast<TreeTraversal>(
      static_cast<std::underlying_type_t<TreeTraversal>>(lhs) |
      static_cast<std::underlying_type_t<TreeTraversal>>(rhs));
}

constexpr TreeTraversal operator&(TreeTraversal lhs, TreeTraversal rhs) {
  return static_cast<TreeTraversal>(
      static_cast<std::underlying_type_t<TreeTraversal>>(lhs) &
      static_cast<std::underlying_type_t<TreeTraversal>>(rhs));
}

#define TRAVERSAL_TO_TEMPLATE_ARG(order, functionName, functionArgs) \
  switch (order) {                                                   \
    case TreeTraversal::None:                                        \
      break;                                                         \
    case TreeTraversal::PreOrder:                                    \
      functionName<TreeTraversal::PreOrder> functionArgs;            \
      break;                                                         \
    case TreeTraversal::PostOrder:                                   \
      functionName<TreeTraversal::PostOrder> functionArgs;           \
      break;                                                         \
    case TreeTraversal::InOrder:                                     \
      functionName<TreeTraversal::InOrder> functionArgs;             \
      break;                                                         \
    case TreeTraversal::PreAndPostOrder:                             \
      functionName<TreeTraversal::PreAndPostOrder> functionArgs;     \
      break;                                                         \
    case TreeTraversal::PreAndInOrder:                               \
      functionName<TreeTraversal::PreAndInOrder> functionArgs;       \
      break;                                                         \
    case TreeTraversal::PostAndInOrder:                              \
      functionName<TreeTraversal::PostAndInOrder> functionArgs;      \
      break;                                                         \
    case TreeTraversal::Any:                                         \
      functionName<TreeTraversal::Any> functionArgs;                 \
      break;                                                         \
  }

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
/// \tparam Order The kinds of tree traversals the visitor shall be notified
/// about \tparam V Visitor type. Must be invocable with a FullBinaryNode<T>
/// argument. Optionally, can take a second argument of type TreeTraversal which
/// indicates the context of the visitor invocation. Can optionally return a
/// value convertible to bool to indicate whether subtree exploration for the
/// current node shall proceed. \tparam NodeType The type of nodes to be
/// visited.
///                  Can be VisitInternal, VisitLeaf, or VisitAll.
/// \param node Node to visit.
/// \param f Visitor.
///
template <TreeTraversal order, typename T, typename V, typename NodeType>
void visit(FullBinaryNode<T> const& node, V const& f, NodeType) {
  using Node = FullBinaryNode<T>;
  static_assert(std::is_same_v<NodeType, VisitLeaf> ||
                    std::is_same_v<NodeType, VisitInternal> ||
                    std::is_same_v<NodeType, VisitAll>,
                "Not sure which nodes to visit");
  const auto invoke = [](const V& f, const Node& node,
                         TreeTraversal context) -> bool {
    if constexpr (std::is_invocable_v<decltype(f), decltype(node),
                                      decltype(context)>) {
      using result_type =
          std::invoke_result_t<decltype(f), decltype(node), decltype(context)>;
      if constexpr (std::is_same_v<result_type, void>) {
        f(node, context);
        return true;
      } else {
        return static_cast<bool>(f(node, context));
      }
    } else {
      static_assert(
          std::is_invocable_v<decltype(f), decltype(node)>,
          "Visitor must be a (const) callable that takes a FullBinaryNode<T> "
          "and optionally a TreeTraversal argument");
      using result_type = std::invoke_result_t<decltype(f), decltype(node)>;
      if constexpr (std::is_same_v<result_type, void>) {
        f(node);
        return true;
      } else {
        return f(node);
      }
    }
  };

  if (node.leaf()) {
    if constexpr (!std::is_same_v<NodeType, VisitInternal>) {
      invoke(f, node, TreeTraversal::Any);
    }
  } else {
    if constexpr ((order & TreeTraversal::PreOrder) ==
                      TreeTraversal::PreOrder &&
                  !std::is_same_v<NodeType, VisitLeaf>) {
      if (!invoke(f, node, TreeTraversal::PreOrder)) return;
    }

    visit<order>(node.left(), f, NodeType{});

    if constexpr ((order & TreeTraversal::InOrder) == TreeTraversal::InOrder &&
                  !std::is_same_v<NodeType, VisitLeaf>) {
      if (!invoke(f, node, TreeTraversal::InOrder)) return;
    }

    visit<order>(node.right(), f, NodeType{});

    if constexpr ((order & TreeTraversal::PostOrder) ==
                      TreeTraversal::PostOrder &&
                  !std::is_same_v<NodeType, VisitLeaf>) {
      if (!invoke(f, node, TreeTraversal::PostOrder)) return;
    }
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
  FullBinaryNode& left() { return *checked_ptr_access(left_); }

  ///
  /// \return Right node if this is an internal node, throws otherwise.
  ///
  FullBinaryNode const& right() const { return *checked_ptr_access(right_); }
  FullBinaryNode& right() { return *checked_ptr_access(right_); }

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
  /// \brief Visit the tree in the order specified by the order argument.
  /// \tparam F Type of the visitor.
  /// \param visitor Visitor to be invoked on each node. The visitor can
  /// optionally return a value convertible to bool to indicate whether the
  /// subtree of the current node shall be explored.
  /// \param order Tree traversal order(s) to invoke the visitor for.
  ///
  template <typename F>
  void visit(F const& visitor,
             TreeTraversal order = TreeTraversal::PreOrder) const {
    TRAVERSAL_TO_TEMPLATE_ARG(order, sequant::visit,
                              (*this, visitor, VisitAll{}));
  }

  ///
  /// \brief Visit the internal nodes of the tree in the order specified by the
  /// order argument.
  /// \tparam F Type of the visitor.
  /// \param visitor Visitor to
  /// be invoked on each node. The visitor can optionally return a value
  /// convertible to bool to indicate whether the subtree of the current node
  /// shall be explored.
  /// \param order Tree traversal order(s) to invoke the visitor for.
  ///
  template <typename F>
  void visit_internal(F const& visitor,
                      TreeTraversal order = TreeTraversal::PreOrder) const {
    TRAVERSAL_TO_TEMPLATE_ARG(order, sequant::visit,
                              (*this, visitor, VisitInternal{}));
  }

  ///
  /// \brief Visit the leaf nodes of the tree in the order specified by the
  /// order argument.
  /// \tparam F Type of the visitor.
  /// \param visitor Visitor to
  /// be invoked on each node. The visitor can optionally return a value
  /// convertible to bool to indicate whether the subtree of the current node
  /// shall be explored.
  /// \param order Tree traversal order(s) to invoke the visitor for.
  ///
  template <typename F>
  void visit_leaf(F const& visitor,
                  TreeTraversal order = TreeTraversal::PreOrder) const {
    TRAVERSAL_TO_TEMPLATE_ARG(order, sequant::visit,
                              (*this, visitor, VisitLeaf{}));
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

#undef TRAVERSAL_TO_TEMPLATE_ARG

#endif  // SEQUANT_BINARY_NODE_HPP
