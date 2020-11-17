#ifndef SEQUANT_FACTORIZE_EVAL_SEQUENCE_HPP
#define SEQUANT_FACTORIZE_EVAL_SEQUENCE_HPP

#include <functional>
#include <ostream>
#include <vector>

namespace sequant::factorize {

/** Almost like a rooted binary tree.
 *
 * @author Bimal Gaudel
 *
 * @version 07 Oct 2020
 */
struct eval_sequence {
  /** Label of the node. */
  size_t label;

  /** Children of the node. */
  std::vector<eval_sequence> children;

  /** Construct eval_sequence.
   *
   * @param l Root label.
   */
  eval_sequence(size_t l) : label{l} {}  // not explicit

  /** Construct eval_sequence.
   *
   * @param l Root label.
   * @param children Children eval_sequence objects.
   */
  eval_sequence(size_t l, std::vector<eval_sequence> &&children)
      : label{l}, children{std::move(children)} {}

  /**
   * Construct eval_sequence.
   * @param l Root label.
   * @param labels Children labels.
   */
  eval_sequence(size_t l, std::initializer_list<size_t> &&labels);

  /** Default ctor. */
  eval_sequence() = default;

  /** Default dtor. */
  ~eval_sequence() = default;

  /** Default copy ctor. */
  eval_sequence(const eval_sequence &) = default;

  /** Default copy assign. */
  eval_sequence &operator=(const eval_sequence &) = default;

  /** Default move ctor. */
  eval_sequence(eval_sequence &&) = default;

  /** Default move assign. */
  eval_sequence &operator=(eval_sequence &&) = default;
};

/** Compare two eval_sequence objects. */
bool operator==(const eval_sequence &, const eval_sequence &);

/** Stream out an eval_sequence. */
std::wostream &operator<<(std::wostream &, const eval_sequence &);

/**
 * Enumerates all possible evaluation sequences with the elements of @c leaves as the
 * leaves and call @c callback on each such enumerations.
 *
 * @param leaves Leaves of the trees being enumerated. All enumerations have all
 * of these leaves.
 *
 * @param callback Function to be called on each enumeration.
 *
 *
 * There are (2k - 1)!! number of ways to evaluate a product of k factors.
 * The notation x!! implies 'double factorial' of x which is defined as the
 * product of all integers from 1 upto x whose parity matches that of x.
 *
 * Eg. 5!! =     5 * 3 * 1.
 *     7!! = 7 * 5 * 3 * 1.
 *     8!! = 8 * 6 * 4 * 2.
 *     4!! =         4 * 2.
 */
void enumerate_eval_sequence(
    const std::vector<eval_sequence> &leaves,
    const std::function<void(const eval_sequence &)> &callback = {});

}  // namespace sequant::factorize

#endif  // SEQUANT_FACTORIZE_EVAL_SEQUENCE_HPP
