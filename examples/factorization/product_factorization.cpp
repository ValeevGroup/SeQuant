//
// Created by Bimal Gaudel on 12/22/19.
//

// #include "../sequant_setup.hpp"
#include <SeQuant/domain/factorize/factorizer.hpp>
#include <SeQuant/domain/factorize/path_tree.hpp>
#include <SeQuant/core/container.hpp>

#include <memory>
#include <iostream>
#include <string>

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
using sequant::factorize::detail::PathTree;
using PathPtr      = std::shared_ptr<PathTree>;
using vec_path_ptr = sequant::container::svector<PathPtr>;

  PathPtr tree0, tree1, tree2, tree3, tree4, tree5;
  tree0 = std::make_shared<PathTree>(PathTree{0});
  tree1 = std::make_shared<PathTree>(PathTree{1});
  tree2 = std::make_shared<PathTree>(PathTree{2});
  tree3 = std::make_shared<PathTree>(PathTree{3});
  tree4 = std::make_shared<PathTree>(PathTree{4});
  tree5 = std::make_shared<PathTree>(PathTree{5});

  cout << "Created trees..\n";
  cout << "tree0: " << tree0->print_tree() << "\n";
  cout << "tree1: " << tree1->print_tree() << "\n";
  cout << "tree2: " << tree2->print_tree() << "\n";
  cout << "tree3: " << tree3->print_tree() << "\n";
  cout << "tree4: " << tree4->print_tree() << "\n";
  cout << "tree5: " << tree5->print_tree() << "\n";

  cout << "\nComputing contraction paths..\n";
  vec_path_ptr args{tree0, tree1, tree2, tree3, tree4};
  sequant::factorize::detail::path_scanner(args);

  return  0;
}