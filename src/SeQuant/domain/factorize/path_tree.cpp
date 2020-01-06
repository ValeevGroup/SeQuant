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
PathTree::PathTree(const PathTree& rhs) {
  label_ = rhs.get_label();
  auto chsize = rhs.get_children().size();
  children_.resize(chsize);
  for (auto i=0; i<chsize; ++i){
    children_[i] = std::make_shared<PathTree>(*rhs.get_children()[i]);
  }
}

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
