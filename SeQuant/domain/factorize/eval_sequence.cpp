#include "eval_sequence.hpp"

#include <ostream>

namespace sequant::factorize {

eval_sequence::eval_sequence(size_t l, std::initializer_list<size_t> &&labels)
    : label{l} {
  children.reserve(labels.size());
  for (const auto &l : labels) children.emplace_back(l);
}

bool operator==(const eval_sequence &lhs, const eval_sequence &rhs) {
  return lhs.label == rhs.label && lhs.children == rhs.children;
}

std::wostream &operator<<(std::wostream &os, const eval_sequence &tree) {
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
    const std::vector<eval_sequence> &leaves,
    const std::function<void(const eval_sequence &)> &callback) {
  if (leaves.size() == 1) callback(*leaves.begin());
  for (auto i = 0; i < leaves.size(); ++i) {
    for (auto j = i + 1; j < leaves.size(); ++j) {
      std::vector<eval_sequence> new_args{leaves[i]};
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
