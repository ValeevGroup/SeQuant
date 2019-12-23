//
// Created by Bimal Gaudel on 12/22/19.
//

// #include "../sequant_setup.hpp"
// #include <SeQuant/domain/factorize/factorizer.hpp>
#include <iostream>
#include <string>
#include <cstddef>
#include <vector>

class Path;

using PathPtr    = Path *;
using PathPtrVec = std::vector<PathPtr>;

class Path {
 public:
  Path() = default;

  explicit Path(size_t x): label_{x} {
    parent_ = nullptr;
    children_.empty();
  }

  auto get_label() const { return label_;}

  void set_parent(PathPtr ptr) { parent_ = ptr; }

  auto get_parent() const { return parent_; }

  void add_child(PathPtr ptr) {
    children_.push_back(ptr);
    ptr->set_parent(this);
  }

  void add_sibling(PathPtr ptr) { parent_->add_child(ptr); }

  bool is_leaf() const { return children_.size() == 0; }

  bool is_root() const { return parent_ == nullptr; }

 private:

  size_t label_;

  PathPtr parent_;

  PathPtrVec children_;
};

std::string path_is_leaf(const Path& path) {
  if (path.is_leaf()) return "True";
  return "False";
}

std::string path_is_root(const Path& path) {
  if (path.is_root()) return "True";
  return "False";
}

int main() {
  // global sequant setup...
  // std::setlocale(LC_ALL, "en_US.UTF-8");
  // std::wcout.precision(std::numeric_limits<double>::max_digits10);
  // std::wcerr.precision(std::numeric_limits<double>::max_digits10);
  // std::wcout.sync_with_stdio(true);
  // std::wcerr.sync_with_stdio(true);
  // std::wcout.imbue(std::locale("en_US.UTF-8"));
  // std::wcerr.imbue(std::locale("en_US.UTF-8"));
  // std::wcout.sync_with_stdio(true);
  // std::wcerr.sync_with_stdio(true);
  // sequant::detail::OpIdRegistrar op_id_registrar;

  // sequant::mbpt::set_default_convention();

  // TensorCanonicalizer::register_instance(
  //     std::make_shared<DefaultTensorCanonicalizer>());
  // Logger::get_instance().wick_stats = false;
  // auto cc_r = cceqvec{ 2, 2 }(true, true, true, true);

  using std::cout;

  auto p1 = Path{0};
  auto p2 = Path{1};
  auto p3 = Path{2};

  cout << "All must be leaves..\n";
  cout << "p1.is_leaf() :" << path_is_leaf(p1) << "\n";
  cout << "p2.is_leaf() :" << path_is_leaf(p2) << "\n";
  cout << "p3.is_leaf() :" << path_is_leaf(p3) << "\n";
  
  cout << "\nAll must be roots..\n";
  cout << "p1.is_root() :" << path_is_root(p1) << "\n";
  cout << "p2.is_root() :" << path_is_root(p2) << "\n";
  cout << "p3.is_root() :" << path_is_root(p3) << "\n";

  cout << "\np1 -> p2 -> p3\n";
  p1.add_child(&p2);
  p2.add_child(&p3);

  cout << "\nOnly p3 shall be leaf..\n";
  cout << "p1.is_leaf() :" << path_is_leaf(p1) << "\n";
  cout << "p2.is_leaf() :" << path_is_leaf(p2) << "\n";
  cout << "p3.is_leaf() :" << path_is_leaf(p3) << "\n";
  
  cout << "\nOnly p1 shall be root..\n";
  cout << "p1.is_root() :" << path_is_root(p1) << "\n";
  cout << "p2.is_root() :" << path_is_root(p2) << "\n";
  cout << "p3.is_root() :" << path_is_root(p3) << "\n";

  return  0;
}
