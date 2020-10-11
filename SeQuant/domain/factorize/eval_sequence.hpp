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
struct rooted_tree {
  /** Label of the node. */
  size_t label;

  /** Children of the node. */
  std::vector<rooted_tree> children;

  /** Construct rooted_tree.
   *
   * @param l Root label.
   */
  rooted_tree(size_t l) : label{l} {}  // not explicit

  /** Construct rooted_tree.
   *
   * @param l Root label.
   * @param children Children rooted_tree objects.
   */
  rooted_tree(size_t l, std::vector<rooted_tree> &&children)
      : label{l}, children{std::move(children)} {}

  /**
   * Construct rooted_tree.
   * @param l Root label.
   * @param labels Children labels.
   */
  rooted_tree(size_t l, std::initializer_list<size_t> &&labels);

  /** Default ctor. */
  rooted_tree() = default;

  /** Default dtor. */
  ~rooted_tree() = default;

  /** Default copy ctor. */
  rooted_tree(const rooted_tree &) = default;

  /** Default copy assign. */
  rooted_tree &operator=(const rooted_tree &) = default;

  /** Default move ctor. */
  rooted_tree(rooted_tree &&) = default;

  /** Default move assign. */
  rooted_tree &operator=(rooted_tree &&) = default;
};

/** Compare two rooted_tree objects. */
bool operator==(const rooted_tree &, const rooted_tree &);

/** Stream out a rooted tree. */
std::wostream &operator<<(std::wostream &, const rooted_tree &);

/**
 * Enumerates all possible rooted trees with the elements of @c leaves as the
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
 * product of all integers from 1 upto x whose parity matches with that of x.
 *
 * Eg. 5!! =     5 * 3 * 1.
 *     7!! = 7 * 5 * 3 * 1.
 *     8!! = 8 * 6 * 4 * 2.
 *     4!! =         4 * 2.
 */
void enumerate_eval_sequence(
    const std::vector<rooted_tree> &leaves,
    const std::function<void(const rooted_tree &)> &callback = {});

}  // namespace sequant::factorize

#endif  // SEQUANT_FACTORIZE_EVAL_SEQUENCE_HPP
