#ifndef SEQUANT_BINARY_NODE_HPP
#define SEQUANT_BINARY_NODE_HPP

#include <cassert>
#include <memory>
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
void visit(FullBinaryNode<T> const& node, V f, NodeType) {
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

  FullBinaryNode<T>* parent_{nullptr};

  node_ptr deep_copy() const {
    return leaf() ? std::make_unique<FullBinaryNode<T>>(data_)
                  : std::make_unique<FullBinaryNode<T>>(
                        data_, left_->deep_copy(), right_->deep_copy());
  }

  template <typename Ptr>
  static Ptr const& checked_ptr_access(Ptr const& n) {
    if (n)
      return n;
    else
      throw std::runtime_error(
          "Dereferenced nullptr: use leaf() or root() methods to check for "
          "leaf and root nodes");
  }

 public:
  ///
  /// Construct an internal node with empty left and right nodes.
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
      : FullBinaryNode(std::move(d),
                       std::make_unique<FullBinaryNode>(std::move(l)),
                       std::make_unique<FullBinaryNode>(std::move(r))) {}

  ///
  /// Constructs an internal node with left and right nodes.
  ///
  /// \param d Data in the internal node.
  /// \param l Left node.
  /// \param r Right node
  FullBinaryNode(T d, FullBinaryNode<T> l, FullBinaryNode<T> r)
      : FullBinaryNode(std::move(d),
                       std::make_unique<FullBinaryNode<T>>(std::move(l)),
                       std::make_unique<FullBinaryNode<T>>(std::move(r))) {}

  ///
  /// Constructs an internal node with left and right node pointers.
  ///
  /// \param d Data in the internal node.
  /// \param l Left node pointer.
  /// \param r Right node pointer.
  FullBinaryNode(T d, node_ptr&& l, node_ptr&& r)
      : data_{std::move(d)}, left_{std::move(l)}, right_{std::move(r)} {
    if (left_) {
      left_->parent_ = this;
    }
    if (right_) {
      right_->parent_ = this;
    }
  }

  FullBinaryNode(FullBinaryNode<T> const& other)
      : data_{other.data_},
        left_{other.left_ ? other.left_->deep_copy() : nullptr},
        right_{other.right_ ? other.right_->deep_copy() : nullptr},
        parent_(nullptr) {
    if (left_) {
      left_->parent_ = this;
    }
    if (right_) {
      right_->parent_ = this;
    }
  }

  FullBinaryNode& operator=(FullBinaryNode<T> const& other) {
    auto temp = other.deep_copy();
    data_ = std::move(temp->data_);
    left_ = std::move(temp->left_);
    right_ = std::move(temp->right_);
    if (left_) {
      left_->parent_ = this;
    }
    if (right_) {
      right_->parent_ = this;
    }
    // parent_ remains unchanged
    return *this;
  }

  FullBinaryNode(FullBinaryNode<T>&& other) : data_(other.data_) {
    // Note: we explicitly have to initialize data_ in the member initializer
    // list as not doing that would impose default-constructibility on T. The
    // assumption is that all halfway decent compilers will optimize this extra
    // copy away.
    *this = std::move(other);
  }

  FullBinaryNode& operator=(FullBinaryNode<T>&& node) {
    data_ = std::move(node.data_);

    // We have to save a temporary copy of these, in case the node we're moving
    // from is pointed to (and thus owned) by either left_.
    // If we don't do this, overwriting of the owning pointer leads to deleting
    // node, in which case subsequent accesses to it are invalid.
    auto left_tmp = std::move(left_);

    left_ = std::move(node.left_);
    right_ = std::move(node.right_);

    if (left_) {
      left_->parent_ = this;
    }
    if (right_) {
      right_->parent_ = this;
    }

    // parent_ remains unchanged

    return *this;
  }

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
  /// \brief Check if the object is a root node
  ///
  /// \return Whether this node represents the root in its respective tree
  [[nodiscard]] bool root() const { return !parent_; }

  ///
  /// \return The parent node
  /// \note This assumes that this object is not a root node
  ///
  [[nodiscard]] const FullBinaryNode<T>& parent() const {
    return *checked_ptr_access(parent_);
  }

  ///
  /// \return The parent node
  /// \note This assumes that this object is not a root node
  ///
  [[nodiscard]] FullBinaryNode<T>& parent() {
    return *checked_ptr_access(parent_);
  }

  ///
  /// \return Size of the tree rooted at this node
  ///
  [[nodiscard]] std::size_t size() const {
    if (leaf()) {
      return 1;
    }

    return left().size() + right().size() + 1;
  }

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
  void visit(F visitor, TreeTraversal order = TreeTraversal::PreOrder) const {
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
  void visit_internal(F visitor,
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
  void visit_leaf(F visitor,
                  TreeTraversal order = TreeTraversal::PreOrder) const {
    TRAVERSAL_TO_TEMPLATE_ARG(order, sequant::visit,
                              (*this, visitor, VisitLeaf{}));
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

  using ranges::size;
  assert(size(rng) > 0);

  using ranges::views::move;
  using ranges::views::tail;
  return ranges::accumulate(
      rng | tail | move, std::move(ranges::front(rng)),
      [&op, invoke_on_node](auto&& l, auto&& r) {
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

#undef TRAVERSAL_TO_TEMPLATE_ARG

#endif  // SEQUANT_BINARY_NODE_HPP
