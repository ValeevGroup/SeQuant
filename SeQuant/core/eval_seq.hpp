#ifndef SEQUANT_UTILS_EVAL_SEQ_HPP
#define SEQUANT_UTILS_EVAL_SEQ_HPP

#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view.hpp>
#include <vector>

#include <SeQuant/core/binary_node.hpp>

namespace sequant {

/**
 * General ordered tree.
 *
 * @author Bimal Gaudel
 *
 * @version 07 Oct 2020
 */
template <typename T>
class EvalSeq {
  /** Label of the node. */
  T label_;

  /** Children of the node. */
  std::vector<EvalSeq> nodes_;

 public:
  /** Construct EvalSeq.
   *
   * @param l Root label.
   */
  explicit EvalSeq(T l) : label_{std::move(l)} {}

  /** Construct EvalSeq.
   *
   * @param l Root label.
   * @param children Children EvalSeq objects.
   */
  EvalSeq(T l, std::vector<EvalSeq<T>> children) : EvalSeq(std::move(l)) {
    nodes_ = std::move(children);
  }

  /**
   * Construct EvalSeq.
   *
   * @param l Root label.
   * @param labels Children labels.
   */
  EvalSeq(T l, std::initializer_list<T> &&labels) : EvalSeq(std::move(l)) {
    nodes_.reserve(labels.size());
    for (auto &&lbl : labels) nodes_.emplace_back(std::move(lbl));
  }

  const T &label() const { return label_; }

  const std::vector<EvalSeq<T>> &nodes() const { return nodes_; }

  /**
   * Append a node to the sequence of evaluation.
   */
  void seque(T n) { nodes_.emplace_back(std::move(n)); }

  /**
   * Append a node to the sequence.
   */
  void seque(EvalSeq<T> n) { nodes_.emplace_back(std::move(n)); }

  /**
   * Check if the node is a terminal in the evaluation sequence.
   */
  [[nodiscard]] bool terminal() const { return nodes_.empty(); }

  template <typename F>
  auto evaluate(F &&fun) const {
    static_assert(std::is_invocable_v<F, T const &>,
                  "F can't be called on T const&");

    using return_data_t = std::invoke_result_t<F, T const &>;

    static_assert(
        std::is_invocable_v<F, return_data_t const &, return_data_t const &>,
        "non-terminal node evaluator missing");

    static_assert(std::is_same_v<return_data_t,
                                 std::invoke_result_t<F, return_data_t const &,
                                                      return_data_t const &>>,
                  "fun(T) and fun(T, T) have different return types");

    auto parent_result = fun(label());

    if (terminal()) return std::move(parent_result);

    return ranges::accumulate(
        nodes().begin(), nodes().end(), std::move(parent_result),
        [&fun](auto &&leval, auto const &rseq) {
          //
          auto reval = rseq.evaluate(std::forward<F>(fun));
          return fun(std::move(leval), std::move(reval));
        });
  }

  template <typename F>
  auto binarize(F &&binarizer) const {
    static_assert(std::is_invocable_v<F, T const &>,
                  "terminal node binarizer missing");

    using return_data_t = std::invoke_result_t<F, T const &>;

    static_assert(
        std::is_invocable_v<F, return_data_t const &, return_data_t const &>,
        "non-terminal node binarizer missing");
    static_assert(
        std::is_same_v<return_data_t,
                       std::invoke_result_t<F, return_data_t const &,
                                            return_data_t const &>>,
        "binarizer(T) and binarizer(T, T) have different return types");

    //
    // todo
    //
    // struct {
    //   auto operator()(T const &node) const {
    //     return BinaryNode<return_data_t>{binarizer(node)};
    //   }
    //
    //   auto operator()(BinaryNode<return_data_t> &&lnode,
    //                   BinaryNode<return_data_t> &&rnode) const {
    //     auto pres = BinaryNode<return_data_t>{binarizer(*lnode, *rnode)};
    //     return BinaryNode<return_data_t>{std::move(pres), std::move(lnode),
    //                                       std::move(rnode)};
    //   }
    // } evaluator;
    //
    // return evaluate(std::forward<decltype(evaluator)>(evaluator));

    auto parent_result = BinaryNode<return_data_t>{binarizer(label())};

    if (terminal()) return std::move(parent_result);

    return ranges::accumulate(
        nodes().begin(), nodes().end(), std::move(parent_result),
        [&binarizer](auto &&lexpr, const auto &rseq) {
          //
          auto rnode = rseq.binarize(std::forward<F>(binarizer));

          auto bin_res = binarizer(*lexpr, *rnode);

          return BinaryNode<return_data_t>{std::move(bin_res),
                                            std::move(lexpr), std::move(rnode)};
        });
  }

  template <typename F>
  auto transform(F &&pred) const {
    static_assert(std::is_invocable_v<F, T const &>,
                  "predicate to transform missing");

    using return_data_t = std::invoke_result_t<F, T const &>;

    using ranges::views::transform;

    auto parent_result = pred(label());
    if (terminal()) return EvalSeq<return_data_t>{parent_result};

    auto children = nodes() | transform([&pred](const auto &x) {
                      return x.transform(std::forward<F>(pred));
                    }) |
                    ranges::to<std::vector<EvalSeq<return_data_t>>>;

    return EvalSeq<return_data_t>{parent_result, std::move(children)};
  }

};  // EvalSeq<T>

template <typename T, typename V>
bool operator==(const EvalSeq<T> &lhs, const EvalSeq<V> &rhs) {
  return  // (&lhs == &rhs) ||
      (lhs.label() == rhs.label() && lhs.nodes() == rhs.nodes());
}

/**
 * Enumerates all possible rooted trees with the elements of @c paths as the
 * leaves and call @c callback on each such enumerations.
 *
 * @param nodes Leaves of the trees being enumerated. All enumerations have all
 * of these leaves.
 *
 * @param callback Function to be called on each enumeration. Arity one function
 * that takes const EvalSeq<T>& argument.
 *
 *
 * There are (2k - 1)!! number of ways to evaluate a product of k factors.
 * The notation x!! implies 'double factorial' of x which is defined as the
 * product of all integers from 1 upto x whose parity matches with that of x.
 *
 * Eg. 5!! =     5 * 3 * 1.
 *     7!! = 7 * 5 * 3 * 1.
 *     8!! = 8 * 6 * 4 * 2.
 *     4!! =         4 * 2.
 */
template <typename T, typename Comp = std::less<T>, typename F>
void enumerate_eval_seq(
    const std::vector<EvalSeq<T>> &nodes,
    F &&callback = [](const EvalSeq<T> &) -> void {}) {
  static_assert(std::is_invocable_v<F, const EvalSeq<T> &>,
                "callback function signature doesn't match");
  if (nodes.size() == 1) callback(*nodes.begin());
  for (auto i = 0; i < nodes.size(); ++i) {
    for (auto j = i + 1; j < nodes.size(); ++j) {
      std::vector<EvalSeq<T>> new_args{nodes[i]};
      new_args.begin()->seque(nodes[j]);

      bool skip_recursive_call = false;
      for (auto k = 0; k < nodes.size(); ++k)
        if (k != i && k != j) {
          // remove redundancy by lexicographic comparison
          if ((!nodes[k].terminal()) &&
              (Comp{}(nodes[k].label(), nodes[i].label()))) {
            skip_recursive_call = true;
            new_args.clear();
            break;
          }
          new_args.emplace_back(nodes[k]);
        }

      if (!skip_recursive_call)
        enumerate_eval_seq(new_args, std::forward<F>(callback));
    }  // for j
  }    // for i
}  // enumerate_eval_seq

}  // namespace sequant

#endif  // SEQUANT_EVAL_SEQ_HPP
