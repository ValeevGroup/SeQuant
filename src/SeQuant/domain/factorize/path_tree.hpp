//
// Created by Bimal Gaudel on 12/26/19.
//
// A tree data structure of positions of SeQuant Expr objects
// has to be implemented. Choosing positions to be ordinals makes it
// easier for eliminating redundant path computations.
//
// Such a tree for a product of four tensors will have paths that look like
// this:
// (0, 1, 2, 3)
// (0, 1, 3, 2)
// (0, 2, 1, 3)
// (0, 2, 3, 1)
// (0, 3, 1, 2)
// (0, 3, 2, 1)
// (1, 2, 0, 3)
// (1, 2, 3, 0)
// ((0, 3), (1, 2))
// (1, 3, 0, 2)
// (1, 3, 2, 0)
// ((0, 2), (1, 3))
// (2, 3, 0, 1)
// (2, 3, 1, 0)
// ((0, 1), (2, 3))
//

#ifndef SEQUANT_FACTORIZE_PATH_TREE_HPP
#define SEQUANT_FACTORIZE_PATH_TREE_HPP

#include <SeQuant/core/container.hpp>

#include <cstddef>
#include <memory>

namespace sequant {
namespace factorize {
namespace detail {

class PathTree;

class PathTree {
 public:
  PathTree() = default;

  explicit PathTree(size_t);

  ~PathTree() = default;

  size_t get_label() const;

  const sequant::container::svector<std::shared_ptr<PathTree>>& get_children()
      const;

  sequant::container::svector<std::shared_ptr<PathTree>>& get_children();

  void add_child(std::shared_ptr<PathTree>&);

  void pop_last_child();

  bool is_leaf() const;

  std::string print_tree() const;

 private:
  size_t label_;

  sequant::container::svector<std::shared_ptr<PathTree>> children_;
};

}  // namespace detail
}  // namespace factorize
}  // namespace sequant

#endif  // SEQUANT_FACTORIZE_PATH_TREE_HPP
