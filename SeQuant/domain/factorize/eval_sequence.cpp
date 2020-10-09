#include "eval_sequence.hpp"

#include <ostream>

namespace sequant::factorize {

std::wostream &operator<<(std::wostream &os, const rooted_tree &tree) {
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
    const std::vector<rooted_tree> &leaves,
    const std::function<void(const rooted_tree &)> &callback) {
  if (leaves.size() == 1) callback(*leaves.begin());
  for (auto i = 0; i < leaves.size(); ++i) {
    for (auto j = i + 1; j < leaves.size(); ++j) {
      std::vector<rooted_tree> new_args{leaves[i]};
      new_args.begin()->children.push_back(leaves[j]);

      bool skip_recursive_call = false;
      for (auto k = 0; k < leaves.size(); ++k)
        if (k != i && k != j) {
          // remove redundancy by lexicographic comparison
          if ((!leaves[k].children.empty()) &&
              (leaves[k].label < leaves[i].label)) {
            skip_recursive_call = true;
            new_args.clear();
            break;
          }
          new_args.push_back(leaves[k]);
        }

      if (!skip_recursive_call) enumerate_eval_sequence(new_args, callback);
    }  // for j
  }    // for i
}

}  // namespace sequant::factorize
