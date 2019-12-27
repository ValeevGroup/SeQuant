//
// Created by Bimal Gaudel on 12/20/19.
//

#include "factorizer.hpp"
#include "path_tree.hpp"

#include <SeQuant/core/container.hpp>

#include <iostream>

namespace sequant {
  namespace factorize {
    namespace detail {

using vec_path_ptr = sequant::container::svector<std::shared_ptr<sequant::factorize::detail::PathTree>>;

void path_scanner(vec_path_ptr& paths) {
  if (paths.size() == 1)
    std::cout << paths[0]->print_tree() << "\n";

  for (auto i=0; i < paths.size(); ++i){
    for (auto j=i+1; j < paths.size(); ++j) {

      // form a product of ith and jth tensors
      paths[i]->add_child(paths[j]);

      // set up new arguments for the recursive call
      vec_path_ptr new_args;
      // newly formed prod. is the first element
      new_args.push_back(paths[i]);
      
      // all tensors from paths argument should be
      // in the new_args except the ith and the jth 
      // since they are absorbed in the product formed above
      bool skip_recursive_call = false; // unlike in Python there's no 'for else' loop in C++
      for (auto k = 0; k < paths.size(); ++k){
        if ((k != i) && (k !=j)){

          // remove redundancy by lexicographic comparison
          if ((! paths[k]->is_leaf()) &&
              (paths[k]->get_label() < paths[i]->get_label())){
            skip_recursive_call = true;
            new_args.clear();
            break;
          }
          new_args.push_back(paths[k]);
        }
      }
      if (! skip_recursive_call) path_scanner(new_args);

      paths[i]->pop_last_child();
    } // for j
  } // for i
}

} // namespace detail
} // namespace factorize
} // namespace sequant
