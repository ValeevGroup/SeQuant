#ifndef SEQUANT_BINARY_NODE_HPP
#define SEQUANT_BINARY_NODE_HPP

#include <memory>
#include <range/v3/numeric/accumulate.hpp>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace sequant {

template <typename T>
class BinaryNode {
  template <typename U>
  struct DataNode {
    virtual ~DataNode() = default;
    virtual U const& data() const = 0;
    virtual U& data() = 0;
    [[nodiscard]] virtual bool leaf() const = 0;
    virtual BinaryNode<U> const& left() const = 0;
    virtual BinaryNode<U> const& right() const = 0;
    virtual BinaryNode<U> clone() const = 0;
  };  // DataNode<U>

  template <typename U>
  class data_node_internal final : public DataNode<U> {
    U data_;
    BinaryNode<U> left_;
    BinaryNode<U> right_;

   public:
    data_node_internal(U d, BinaryNode<U>&& l, BinaryNode<U>&& r)
        : data_{std::move(d)}, left_{std::move(l)}, right_{std::move(r)} {}

    U const& data() const override { return data_; }

    U& data() override { return data_; }

    [[nodiscard]] bool leaf() const override { return false; }

    BinaryNode<U> const& left() const override { return left_; }

    BinaryNode<U> const& right() const override { return right_; }

    BinaryNode<U> clone() const override {
      return BinaryNode<U>{data_, left_.clone(), right_.clone()};
    }
  };  // data_node_internal<U>

  template <typename U>
  class data_node_leaf final : public DataNode<U> {
    U data_;

   public:
    explicit data_node_leaf(U d) : data_{std::move(d)} {}

    U const& data() const override { return data_; }

    U& data() override { return data_; }

    [[nodiscard]] bool leaf() const override { return true; }

    BinaryNode<U> const& left() const override {
      throw std::logic_error("left() called on leaf node");
    }

    BinaryNode<U> const& right() const override {
      throw std::logic_error("right() called on leaf node");
    }

    BinaryNode<U> clone() const override { return BinaryNode<U>{data_}; }
  };  // data_node_leaf<U>

 private:
  std::unique_ptr<DataNode<T>> dnode;

 public:
  T const& operator*() const { return dnode->data(); }

  T& operator*() { return dnode->data(); }

  T* operator->() const { return &dnode->data(); }

  [[nodiscard]] bool leaf() const { return dnode->leaf(); }

  BinaryNode<T> const& left() const { return dnode->left(); }

  BinaryNode<T> const& right() const { return dnode->right(); }

  BinaryNode<T> clone() const { return dnode->clone(); }

  BinaryNode(BinaryNode<T> const& other) noexcept {
    *this = other.clone();
  }

  BinaryNode& operator=(BinaryNode<T> const& other) noexcept {
    *this = other.clone();
    return *this;
  }

  BinaryNode(BinaryNode<T>&& other) noexcept = default;

  BinaryNode& operator=(BinaryNode<T>&& other) noexcept = default;

  explicit BinaryNode(T d)
      : dnode{std::make_unique<data_node_leaf<T>>(std::move(d))} {}

  BinaryNode(T d, BinaryNode<T>&& ln, BinaryNode<T>&& rn)
      : dnode{std::make_unique<data_node_internal<T>>(std::move(d),   //
                                                      std::move(ln),  //
                                                      std::move(rn))} {}

  BinaryNode(T d, T ld, T rd)
      : BinaryNode<T>{std::move(d), BinaryNode<T>{ld}, BinaryNode<T>{rd}} {}

  template <typename Cont, typename F>
  BinaryNode(Cont const& container, F&& binarize) {
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
        accumulate(begin(container) + 1, end(container),       // range
                   BinaryNode{binarize(*begin(container))},    // init
                   [&binarize](auto&& acc, const auto& val) {  // predicate
                     auto rnode = BinaryNode{binarize(val)};
                     return BinaryNode{binarize(*acc, *rnode), std::move(acc),
                                       std::move(rnode)};
                   });

    *this = std::move(node);
  }

  template <typename F,
            std::enable_if_t<std::is_invocable_v<F, BinaryNode<T> const&>,
                             bool> = true>
  void visit(F&& pred) const {
    pred(*this);
    if (leaf()) return;

    left().visit(std::forward<F>(pred));
    right().visit(std::forward<F>(pred));
  }

  template <typename F,
            std::enable_if_t<std::is_invocable_v<F, BinaryNode<T> const&>,
                             bool> = true>
  void visit_internal(F&& pred) const {
    if (!leaf()) {
      left().visit_internal(std::forward<F>(pred));
      pred(*this);
      right().visit_internal(std::forward<F>(pred));
    }
  }

  template <typename F,
            std::enable_if_t<std::is_invocable_v<F, BinaryNode<T> const&>,
                             bool> = true>
  void visit_leaf(F&& pred) const {
    if (leaf()) {
      pred(*this);
    } else {
      left().visit_leaf(std::forward<F>(pred));
      right().visit_leaf(std::forward<F>(pred));
    }
  }

  template <typename F,
            std::enable_if_t<
                std::is_invocable_v<F, BinaryNode<T> const&> &&
                    std::is_invocable_v<
                        F, BinaryNode<T> const&,
                        std::invoke_result_t<F, BinaryNode<T> const&> const&,
                        std::invoke_result_t<F, BinaryNode<T> const&> const&>,
                bool> = true>
  auto evaluate(F&& evaluator) const {
    if (leaf()) return evaluator(*this);
    return evaluator(*this, left().evaluate(std::forward<F>(evaluator)),
                     right().evaluate(std::forward<F>(evaluator)));
  }

 private:
  template <typename Ostream, typename F>
  [[maybe_unused]] int digraph(Ostream& os, F&& label_gen,
                               int count = 0) const {
    os << "node" << count << "[label=" << label_gen(*this) << "];\n";

    if (this->leaf()) return count;

    auto lcount = left().digraph(os, std::forward<F>(label_gen), count + 1);
    auto rcount = right().digraph(os, std::forward<F>(label_gen), lcount + 1);
    os << "node" << count << " -> "
       << "node" << count + 1 << ";\n";
    os << "node" << count << " -> "
       << "node" << lcount + 1 << ";\n";

    return rcount;
  }

  template <typename Ostream, typename F, typename G>
  void tikz(Ostream& os, F&& label_gen, G&& spec_gen, size_t indent = 2) const {
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
    left().tikz(os, std::forward<F>(label_gen), std::forward<G>(spec_gen),
                indent + 2);
    os << "}";
    os << "\n";

    pad(os, indent);
    os << "child {";
    right().tikz(os, std::forward<F>(label_gen), std::forward<G>(spec_gen),
                 indent + 2);
    os << "}";
  }

 public:
  template <typename string_t, typename F>
  string_t digraph(F&& label_gen, string_t const& graph_name = {}) const {
    static_assert(std::is_invocable_r_v<string_t, F, BinaryNode<T> const&>,
                  "node label generator F(BinaryNode<T> const &) should "
                  "return string_t");

    auto oss = std::basic_ostringstream{string_t{}};

    oss << "digraph " << graph_name << "{\n";
    this->digraph(oss, std::forward<F>(label_gen), 0);
    oss << "}";
    oss.flush();

    return oss.str();
  }

  template <typename string_t>
  string_t tikz(std::function<string_t(BinaryNode<T> const&)> label_gen,
                std::function<string_t(BinaryNode<T> const&)> spec_gen) const {
    auto oss = std::basic_ostringstream{string_t{}};
    oss << "\\tikz{\n\\";
    tikz(oss, label_gen, spec_gen);
    oss << "\n}";
    oss.flush();
    return oss.str();
  }

};  // BinaryNode<T>

template <typename T, typename U>
bool operator==(BinaryNode<T> const& lhs, BinaryNode<U> const& rhs) {
  return ((*lhs == *rhs) &&
          ((lhs.leaf() && rhs.leaf()) ||
           (lhs.left() == rhs.left() && lhs.right() == rhs.right())));
}

}  // namespace sequant

#endif  // SEQUANT_BINARY_NODE_HPP
