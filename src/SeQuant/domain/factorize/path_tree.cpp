//
// Created by Bimal Gaudel on 12/26/19.
//

#include "path_tree.hpp"

#include <SeQuant/core/container.hpp>

#include <cstddef>
#include <memory>
#include <string>

namespace sequant {
namespace factorize {
namespace detail {

// constructor
PathTree::PathTree(size_t x) : label_{x} { children_.clear(); }

size_t PathTree::get_label() const { return label_; }

void PathTree::add_child(std::shared_ptr<PathTree>& ptr) {
  children_.push_back(ptr);
}

const sequant::container::svector<std::shared_ptr<PathTree>>&
PathTree::get_children() const {
  return children_;
}
sequant::container::svector<std::shared_ptr<PathTree>>&
PathTree::get_children() {
  return children_;
}

void PathTree::pop_last_child() { children_.pop_back(); }

bool PathTree::is_leaf() const { return children_.empty(); }

std::string PathTree::print_tree() const {
  if (is_leaf()) {
    return " " + std::to_string(get_label());
  } else {
    std::string result = " (" + std::to_string(get_label());
    for (const auto& ch : get_children()) {
      result += ch->print_tree();
    }
    result += ")";
    return result;
  }
}

}  // namespace detail
}  // namespace factorize
}  // namespace sequant
