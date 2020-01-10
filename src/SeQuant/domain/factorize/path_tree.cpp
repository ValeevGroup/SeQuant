//
// Created by Bimal Gaudel on 12/26/19.
//

#include "path_tree.hpp"

#include <SeQuant/core/container.hpp>

#include <cstddef>
#include <memory>
#include <string>

namespace sequant::factorize::detail {

// constructor
PathTree::PathTree(size_t x) : label_{x} { children_.clear(); }

// copy constructor
// PathTree::PathTree(const std::shared_ptr<PathTree>& rhs) {
//   label_ = rhs->label_;
//   if (rhs->is_leaf())
//     return;
//   for(const auto& rch: rhs->children_)
//     children_.push_back(std::make_shared<PathTree>(PathTree(rch)));
// }

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

std::wstring PathTree::print_tree() const {
  if (is_leaf()) {
    return L" " + std::to_wstring(get_label());
  } else {
    std::wstring result = L" (" + std::to_wstring(get_label());
    for (const auto& ch : get_children()) {
      result += ch->print_tree();
    }
    result += L")";
    return result;
  }
}

}  // namespace sequant
