#ifndef SEQUANT_UTILS_EVAL_SEQUENCE_HPP
#define SEQUANT_UTILS_EVAL_SEQUENCE_HPP

#include <functional>
#include <numeric>
#include <ostream>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view.hpp>
#include <vector>

#include "binary_expr.hpp"

namespace sequant::utils {

/**
 * General ordered tree.
 *
 * @author Bimal Gaudel
 *
 * @version 07 Oct 2020
 */
template <typename T>
class eval_sequence {
  /** Label of the node. */
  T label_;

  /** Children of the node. */
  std::vector<eval_sequence> nodes_;

 public:
  /** Construct eval_sequence.
   *
   * @param l Root label.
   */
  explicit eval_sequence(T l) : label_{std::move(l)} {}

  /** Construct eval_sequence.
   *
   * @param l Root label.
   * @param children Children eval_sequence objects.
   */
  eval_sequence(T l, std::vector<eval_sequence<T>> children)
      : eval_sequence(std::move(l)) {
    nodes_ = std::move(children);
  }

  /**
   * Construct eval_sequence.
   *
   * @param l Root label.
   * @param labels Children labels.
   */
  eval_sequence(T l, std::initializer_list<T> &&labels)
      : eval_sequence(std::move(l)) {
    nodes_.reserve(labels.size());
    for (auto &&lbl : labels) nodes_.emplace_back(std::move(lbl));
  }

  const T &label() const { return label_; }

  const std::vector<eval_sequence<T>> &nodes() const { return nodes_; }

  /**
   * Append a node to the sequence of evaluation.
   */
  void seque(T n) { nodes_.emplace_back(std::move(n)); }

  /**
   * Append a node to the sequence.
   */
  void seque(eval_sequence<T> n) { nodes_.emplace_back(std::move(n)); }

  /**
   * Check if the node is a terminal in the evaluation sequence.
   */
  [[nodiscard]] bool terminal() const { return nodes_.empty(); }
};

template <typename T, typename V>
bool operator==(const eval_sequence<T> &lhs, const eval_sequence<V> &rhs) {
  return ((&lhs == &rhs) ||
          (lhs.label() == rhs.label() && lhs.nodes() == rhs.nodes()));
}

template <typename T, typename Os>
Os &operator<<(Os &os, const eval_sequence<T> &seq) {
  if (seq.terminal()) {
    os << seq.label();
    return os;
  }
  os << "(" << seq.label() << " ";

  for (auto ii = 0; ii < seq.nodes().size() - 1; ++ii)
    os << seq.nodes()[ii] << " ";
  os << *(seq.nodes().end() - 1);

  os << ")";
  return os;
}

/**
 * Enumerates all possible rooted trees with the elements of @c paths as the
 * leaves and call @c callback on each such enumerations.
 *
 * @param nodes Leaves of the trees being enumerated. All enumerations have all
 * of these leaves.
 *
 * @param callback Function to be called on each enumeration. Arity one function
 * that takes const eval_sequence<T>& argument.
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
void enumerate_eval_sequence(
    const std::vector<eval_sequence<T>> &nodes,
    F &&callback = [](const eval_sequence<T> &) -> void {}) {
  static_assert(std::is_invocable_v<F, const eval_sequence<T> &>,
                "callback function signature doesn't match");
  if (nodes.size() == 1) callback(*nodes.begin());
  for (auto i = 0; i < nodes.size(); ++i) {
    for (auto j = i + 1; j < nodes.size(); ++j) {
      std::vector<eval_sequence<T>> new_args{nodes[i]};
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
        enumerate_eval_sequence(new_args, std::forward<F>(callback));
    }  // for j
  }    // for i
}

template <typename T, typename F>
auto binarize_eval_sequence(const eval_sequence<T> &seq, F &&binarizer) {
  static_assert(std::is_invocable_v<F, T>, "terminal node binarizer missing");

  using return_data_t = std::invoke_result_t<F, T>;

  static_assert(std::is_invocable_v<F, return_data_t, return_data_t>,
                "non-terminal node binarizer missing");
  static_assert(
      std::is_same_v<return_data_t,
                     std::invoke_result_t<F, return_data_t, return_data_t>>,
      "binarizer(T) and binarizer(T, T) have different return types");

  auto parent_result = make_binary_expr<return_data_t>(binarizer(seq.label()));

  if (seq.terminal()) return std::move(parent_result);

  return ranges::accumulate(
      seq.nodes().begin(), seq.nodes().end(), std::move(parent_result),
      [&binarizer](auto &&lexpr, const auto &rseq) {
        auto bin_res = binarizer(
            lexpr->data(),
            binarize_eval_sequence<T, F>(rseq, std::forward<F>(binarizer))
                ->data());

        return make_binary_expr<return_data_t>(
            std::move(bin_res), std::move(lexpr),
            binarize_eval_sequence<T, F>(rseq, std::forward<F>(binarizer)));
      });
}

template <typename T, typename F>
auto transform_eval_sequence(const eval_sequence<T> &seq, F &&pred) {
  static_assert(std::is_invocable_v<F, T>, "predicate to transform missing");

  using return_data_t = std::invoke_result_t<F, T>;

  using ranges::views::transform;

  auto parent_result = pred(seq.label());
  if (seq.terminal()) return eval_sequence<return_data_t>{parent_result};

  auto children =
      seq.nodes() | transform([&pred](const auto &x) {
        return transform_eval_sequence<T, F>(x, std::forward<F>(pred));
      }) |
      ranges::to<std::vector<eval_sequence<return_data_t>>>;

  return eval_sequence<return_data_t>{parent_result, std::move(children)};
}

}  // namespace sequant::utils

#endif  // SEQUANT_UTILS_EVAL_SEQUENCE_HPP
