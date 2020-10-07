#include "eval_sequence.hpp"

#include <ostream>

namespace sequant::factorize {

std::ostream &operator<<(std::ostream &os, const rooted_tree &tree) {
  if (tree.children.empty()) {
    os << tree.label;
    return os;
  }
  os << "(" << tree.label << " ";

  for (auto ii = 0; ii < tree.children.size() - 1; ++ii)
    os << tree.children.at(ii) << " ";
  os << *(tree.children.end() - 1);

  os << ")";
  return os;
}

void enumerate_eval_sequence(
    const std::vector<rooted_tree> &paths,
    const std::function<void(const rooted_tree &)> &callback) {
  if (paths.size() == 1) callback(*paths.begin());
  for (auto i = 0; i < paths.size(); ++i) {
    for (auto j = i + 1; j < paths.size(); ++j) {
      std::vector<rooted_tree> new_args{paths[i]};
      new_args.begin()->children.push_back(paths[j]);

      bool skip_recursive_call = false;
      for (auto k = 0; k < paths.size(); ++k)
        if (k != i && k != j) {
          // remove redundancy by lexicographic comparison
          if ((!paths[k].children.empty()) &&
              (paths[k].label < paths[i].label)) {
            skip_recursive_call = true;
            new_args.clear();
            break;
          }
          new_args.push_back(paths[k]);
        }

      if (!skip_recursive_call)
        enumerate_eval_sequence(new_args, std::move(callback));
    }  // for j
  }    // for i
}

}  // namespace sequant::factorize
