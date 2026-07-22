#ifndef SEQUANT_EXTERNAL_INTERFACE_PROCESSINGTREE_HPP
#define SEQUANT_EXTERNAL_INTERFACE_PROCESSINGTREE_HPP

#include <cassert>
#include <cstddef>
#include <vector>

namespace sequant::util::extint {

class ProcessingTree {
 public:
  ProcessingTree() = default;

  bool is_root() const { return parent_; }
  bool is_leaf() const { return children_.empty(); }

  std::size_t num_children() const { return children_.size(); }

  ProcessingTree &parent() {
    assert(parent_);
    return *parent_;
  }

  std::size_t step_id() const { return step_id_; }

  auto begin() { return children_.begin(); }
  auto end() { return children_.end(); }
  auto begin() const { return children_.begin(); }
  auto end() const { return children_.end(); }
  auto cbegin() const { return children_.begin(); }
  auto cend() const { return children_.end(); }

 private:
  ProcessingTree *parent_ = nullptr;
  std::vector<std::size_t> children_;
  std::size_t step_id_;
};

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_PROCESSINGTREE_HPP
