#ifndef SEQUANT_UTILS_BINARY_NODE_HPP
#define SEQUANT_UTILS_BINARY_NODE_HPP

#include <memory>
#include <range/v3/numeric/accumulate.hpp>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace sequant::utils {

template <typename T>
class binary_node {
  template <typename U>
  struct data_node {
    virtual ~data_node() = default;
    virtual U const& data() const = 0;
    virtual U& data() = 0;
    [[nodiscard]] virtual bool leaf() const = 0;
    virtual binary_node<U> const& left() const = 0;
    virtual binary_node<U> const& right() const = 0;
    virtual binary_node<U> clone() const = 0;
  };  // data_node<U>

  template <typename U>
  class data_node_internal final : public data_node<U> {
    U data_;
    binary_node<U> left_;
    binary_node<U> right_;

   public:
    data_node_internal(U d, binary_node<U>&& l, binary_node<U>&& r)
        : data_{std::move(d)}, left_{std::move(l)}, right_{std::move(r)} {}

    U const& data() const override { return data_; }

    U& data() override { return data_; }

    [[nodiscard]] bool leaf() const override { return false; }

    binary_node<U> const& left() const override { return left_; }

    binary_node<U> const& right() const override { return right_; }

    binary_node<U> clone() const override {
      return binary_node<U>{data_, *left_, *right_};
    }
  };  // data_node_internal<U>

  template <typename U>
  class data_node_leaf final : public data_node<U> {
    U data_;

   public:
    explicit data_node_leaf(U d) : data_{std::move(d)} {}

    U const& data() const override { return data_; }

    U& data() override { return data_; }

    [[nodiscard]] bool leaf() const override { return true; }

    binary_node<U> const& left() const override {
      throw std::logic_error("left() called on leaf node");
    }

    binary_node<U> const& right() const override {
      throw std::logic_error("right() called on leaf node");
    }

    binary_node<U> clone() const override { return binary_node<U>{data_}; }
  };  // data_node_leaf<U>

 private:
  std::unique_ptr<data_node<T>> dnode;

 public:
  T const& operator*() const { return dnode->data(); }

  T& operator*() { return dnode->data(); }

  T* operator->() const { return &dnode->data(); }

  [[nodiscard]] bool leaf() const { return dnode->leaf(); }

  binary_node<T> const& left() const { return dnode->left(); }

  binary_node<T> const& right() const { return dnode->right(); }

  binary_node<T> clone() const {
    return leaf() ? binary_node<T>{**this}
                  : binary_node<T>{**this, *left(), *right()};
  }

  explicit binary_node(T d)
      : dnode{std::make_unique<data_node_leaf<T>>(std::move(d))} {}

  binary_node(T d, binary_node<T>&& ln, binary_node<T>&& rn)
      : dnode{std::make_unique<data_node_internal<T>>(std::move(d),   //
                                                      std::move(ln),  //
                                                      std::move(rn))} {}

  binary_node(T d, T ld, T rd)
      : binary_node<T>{std::move(d), binary_node<T>{ld}, binary_node<T>{rd}} {}

  template <typename Cont, typename F>
  binary_node(Cont const& container, F&& binarize) {
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
                   binary_node{binarize(*begin(container))},   // init
                   [&binarize](auto&& acc, const auto& val) {  // predicate
                     auto rnode = binary_node{binarize(val)};
                     return binary_node{binarize(*acc, *rnode), std::move(acc),
                                        std::move(rnode)};
                   });

    *this = std::move(node);
  }

  template <typename F,
            std::enable_if_t<std::is_invocable_v<F, T const&>, bool> = true>
  void visit(F&& pred) const {
    if (leaf()) {
      pred(**this);
      return;
    }

    left().visit(std::forward<F>(pred));
    pred(**this);
    right().visit(std::forward<F>(pred));
  }

  template <typename F,
            std::enable_if_t<std::is_invocable_v<F, T const&>, bool> = true>
  void visit_internal(F&& pred) const {
    if (leaf()) {
      return;
    }

    left().visit_internal(std::forward<F>(pred));
    pred(**this);
    right().visit_internal(std::forward<F>(pred));
  }

  template <typename F,
            std::enable_if_t<std::is_invocable_v<F, T const&>, bool> = true>
  void visit_leaf(F&& pred) const {
    if (leaf()) {
      pred(**this);
    }

    left().visit_leaf(std::forward<F>(pred));
    right().visit_leaf(std::forward<F>(pred));
  }

  template <typename F,
            std::enable_if_t<
                std::is_invocable_v<F, T const&> &&
                    std::is_invocable_v<
                        F, T const&, std::invoke_result_t<F, T const&> const&,
                        std::invoke_result_t<F, T const&> const&>,
                bool> = true>
  auto evaluate(F&& evaluator) const {
    using ret_type = std::invoke_result_t<F, T const&>;

    // static_assert(
    //     std::is_same_v<ret_type,
    //                 std::invoke_result_t<F, ret_type const&,
    //                                         ret_type const&>>,
    //     "evaluator(T const&)"
    //     " and evaluator(T const&, ret_type const&, ret_type const&"
    //     " have different return types");

    if (leaf()) return evaluator(**this);
    return evaluator(**this, left().evaluate(std::forward<F>(evaluator)),
                     right().evaluate(std::forward<F>(evaluator)));
  }

  template <typename F,
            std::enable_if_t<
                std::is_invocable_v<F, binary_node<T> const&> &&
                    std::is_invocable_v<
                        F, binary_node<T> const&,
                        std::invoke_result_t<F, binary_node<T> const&> const&,
                        std::invoke_result_t<F, binary_node<T> const&> const&>,
                bool> = true>
  auto evaluate(F&& evaluator) const {
    using ret_type = std::invoke_result_t<F, binary_node<T> const&>;

    // static_assert(
    //     std::is_same_v<ret_type,
    //                    std::invoke_result_t<F, binary_node<T> const&,
    //                                         ret_type const&, ret_type
    //                                         const&>>,
    //     "evaluator(binary_node<T> const&)"
    //     " and evaluator(binary_node<T> const&, ret_type const&, ret_type
    //     const&" " have different return types");

    if (leaf()) return evaluator(*this);
    return evaluator(*this, left().evaluate(std::forward<F>(evaluator)),
                     right().evaluate(std::forward<F>(evaluator)));
  }

 private:
  template <typename Ostream, typename F>
  [[maybe_unused]] int digraph(Ostream& os, F&& label_gen,
                               int count = 0) const {
    os << "node" << count << "[label=" << label_gen(**this) << "];\n";

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
    static_assert(std::is_invocable_r_v<string_t, F, T const&>,
                  "node labels generator F(const T&) should return string_t");

    auto oss = std::basic_ostringstream{string_t{}};

    oss << "digraph " << graph_name << "{\n";
    this->digraph(oss, std::forward<F>(label_gen), 0);
    oss << "}";
    oss.flush();

    return oss.str();
  }

  template <typename string_t>
  string_t tikz(std::function<string_t(binary_node<T> const&)> label_gen,
                std::function<string_t(binary_node<T> const&)> spec_gen) const {
    auto oss = std::basic_ostringstream{string_t{}};
    oss << "\\tikz{\n\\";
    tikz(oss, label_gen, spec_gen);
    oss << "\n}";
    oss.flush();
    return oss.str();
  }

};  // binary_node<T>

template <typename T, typename U>
bool operator==(binary_node<T> const& lhs, binary_node<U> const& rhs) {
  return  // (&lhs == &rhs) ||
      ((*lhs == *rhs) &&
       ((lhs.leaf() && rhs.leaf()) ||
        (lhs.left() == rhs.left() && lhs.right() == rhs.right())));
}

}  // namespace sequant::utils

#endif  // SEQUANT_UTILS_BINARY_NODE_HPP
