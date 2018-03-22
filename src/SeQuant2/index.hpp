//
// Created by Eduard Valeyev on 3/20/18.
//

#ifndef SEQUANT2_INDEX_H
#define SEQUANT2_INDEX_H

#include "space.hpp"
#include <set>

namespace sequant2 {

/// @brief Index = label + IndexSpace
/// @note Unlike SeQuant1's ParticleIndex, this Index supports dependencies
/// between indices to be able to express
///       e.g. hiearchical partitioning of index spaces or hiearchical nesting
///       of spaces
class Index {
 public:
  Index() = default;

  /// @param label the label, does not need to be unique
  /// @param space (a const ref to) the IndexSpace object that specifies to this
  /// space this object belongs
  /// @param proto_index_labels labels of proto indices (all must be unique,
  /// i.e. duplicates are not allowed)
  /// @note Index does not manage lifetime of index spaces, so it is user's
  /// responsibility to ensure that @c space will
  ///       outlive this object
  Index(std::wstring_view label, const IndexSpace &space,
        std::initializer_list<Index> proto_indices)
      : label_(label), space_(&space), proto_indices_(proto_indices) {
    check_for_duplicate_proto_indices();
  }

  /// @param label the index label, does not need to be unique, but must be
  /// convertible into an IndexSpace (@sa IndexSpace::instance )
  explicit Index(std::wstring_view label)
      : Index(label, IndexSpace::instance(label), {}) {}

  /// @param label the label, does not need to be unique
  /// @param space (a const ref to) the IndexSpace object that specifies to this
  /// space this object belongs
  /// @note Index does not manage lifetime of index spaces, so it is user's
  /// responsibility to ensure that @c space will
  ///       outlive this object
  Index(std::wstring_view label, const IndexSpace &space) noexcept
      : label_(label), space_(&space), proto_indices_() {}

  /// @param label the label, does not need to be unique
  /// @param proto_index_labels labels of proto indices (all must be unique,
  /// i.e. duplicates are not allowed)
  Index(std::wstring_view label,
        std::initializer_list<Index> proto_indices)
      : Index(label, IndexSpace::instance(label), proto_indices) {
    check_for_duplicate_proto_indices();
  }

  /// @param label the label, does not need to be unique
  /// @param proto_index_labels labels of proto indices (all must be unique,
  /// i.e. duplicates are not allowed)
  Index(std::wstring_view label,
        std::initializer_list<std::wstring_view> proto_index_labels)
      : label_(label), space_(&IndexSpace::instance(label)) {
    if (proto_index_labels.size() != 0) {
      proto_indices_.reserve(proto_index_labels.size());
      for (const auto &plabel : proto_index_labels)
        proto_indices_.push_back(Index(plabel));
    }
    check_for_duplicate_proto_indices();
  }

  /// @return the label
  std::wstring_view label() const { return label_; }
  /// @return the IndexSpace object
  const IndexSpace &space() const {
    assert(space_ != nullptr);
    return *space_;
  }

  /// @return true if this index has proto indices
  bool has_proto_indices() const { return !proto_indices_.empty(); }
  /// @return the list of proto indices of this index
  const std::vector<Index> &proto_indices() const { return proto_indices_; }

 private:
  std::wstring label_{};
  const IndexSpace *space_{};          // pointer to allow default initialization
  std::vector<Index> proto_indices_{}; // an unordered set of unique indices on
  // which this index depends on

  /// throws std::invalid_argument if have duplicates in proto_indices_
  inline void check_for_duplicate_proto_indices();
};

/// @return true if @c index1 is identical to @c index2 , i.e. they belong to
/// the same space, they have the same label, and the same proto-indices (if
/// any)
inline bool operator==(const Index &i1, const Index &i2) {
  return i1.space() == i2.space() && i1.label() == i2.label() &&
      i1.proto_indices() == i2.proto_indices();
}

/// @return false if @c index1 is identical to @c index2 , i.e. they belong to
/// different spaces or they have different labels or they have different
/// proto-indices (if any)
inline bool operator!=(const Index &i1, const Index &i2) { return !(i1 == i2); }

/// @brief The ordering operator

/// @return true if @c i1 preceeds @c i2 in the canonical order; Index objects
/// are ordered lexicographically,
///         first by space, then by label, then by protoindices (if any)
inline bool operator<(const Index &i1, const Index &i2) {
  if (i1.space() < i2.space()) {
    return true;
  } else if (i1.space() == i2.space()) {
    if (i1.label() < i2.label()) {
      return true;
    } else return i1.label() == i2.label() && i1.proto_indices() < i2.proto_indices();
  } else { // i1.space > i2.space
    return false;
  }
}

void Index::check_for_duplicate_proto_indices() {
  std::vector<Index const *> vp;
  vp.reserve(proto_indices_.size());
  for (size_t i = 0; i < proto_indices_.size(); ++i)
    vp.push_back(&proto_indices_[i]);
  std::sort(vp.begin(), vp.end(),
            [](Index const *l, Index const *r) { return *l < *r; });
  if (std::adjacent_find(vp.begin(), vp.end(),
                         [](Index const *l, Index const *r) {
                           return *l == *r;
                         }) != vp.end()) {
    throw std::invalid_argument("Index ctor: duplicate proto indices detected");
  }
}

inline std::wstring to_latex(const Index &index) {
  std::wstring result;
  result = L"{";
  result += index.label();
  if (index.has_proto_indices()) {
    result += L"^{";
    for (const auto &pi: index.proto_indices()) {
      result += to_latex(pi);
    }
    result += L"}";
  }
  result += L"}";
  return result;
}

} // namespace sequant2

#endif // SEQUANT2_INDEX_H
