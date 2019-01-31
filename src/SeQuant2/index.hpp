//
// Created by Eduard Valeyev on 3/20/18.
//

#ifndef SEQUANT2_INDEX_H
#define SEQUANT2_INDEX_H

#include <atomic>
#include <optional>
#include <mutex>
#include <set>
#include <string>

#include <boost/container_hash/hash.hpp>

#include "space.hpp"
#include "vector.hpp"
#include "tag.hpp"

namespace sequant2 {

class Index;
using WstrList = std::initializer_list<std::wstring_view>;
using IndexList = std::initializer_list<Index>;

/// @brief Index = label + IndexSpace
/// @note Unlike SeQuant1's ParticleIndex, this Index supports dependencies
/// between indices to be able to express
///       e.g. hierarchical partitioning of index spaces or hiearchical nesting
///       of spaces
/// @note label has format "label_index" where "label" is a string of characters excluding '_', and "index"
///       is an integer less than the value returned by min_tmp_label() .
class Index : public Taggable {
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
        IndexList proto_indices)
      : label_(label), space_(&space), proto_indices_(proto_indices) {
    check_for_duplicate_proto_indices();
    check_nontmp_label();
  }

  /// @param label the index label, does not need to be unique, but must be
  /// convertible into an IndexSpace (@sa IndexSpace::instance )
  explicit Index(std::wstring_view label)
      : Index(label, IndexSpace::instance(label), {}) {
    check_nontmp_label();
  }

  /// @param label the label, does not need to be unique
  /// @param space (a const ref to) the IndexSpace object that specifies to this
  /// space this object belongs
  /// @note Index does not manage lifetime of index spaces, so it is user's
  /// responsibility to ensure that @c space will
  ///       outlive this object
  Index(std::wstring_view label, const IndexSpace &space) noexcept
      : label_(label), space_(&space), proto_indices_() {
    check_nontmp_label();
  }

  /// @param label the label, does not need to be unique
  /// @param proto_index_labels labels of proto indices (all must be unique,
  /// i.e. duplicates are not allowed)
  Index(std::wstring_view label,
        IndexList proto_indices)
      : Index(label, IndexSpace::instance(label), proto_indices) {
    check_for_duplicate_proto_indices();
    check_nontmp_label();
  }

  /// @param label the label, does not need to be unique
  /// @param proto_index_labels labels of proto indices (all must be unique,
  /// i.e. duplicates are not allowed)
  Index(std::wstring_view label,
        WstrList proto_index_labels)
      : label_(label), space_(&IndexSpace::instance(label)) {
    if (proto_index_labels.size() != 0) {
      proto_indices_.reserve(proto_index_labels.size());
      for (const auto &plabel : proto_index_labels)
        proto_indices_.push_back(Index(plabel));
    }
    check_for_duplicate_proto_indices();
    check_nontmp_label();
  }

  /// creates a globaly-unique temporary index in space @c space . The label of the resulting index =
  /// @c IndexSpace::base_key(space) + '_' + temporary counter.
  /// Each call increments the current tmp counter (see next_tmp_index() ) . To make neater temporary indices
  /// unique in a given scope (e.g. a single term in an expression) use IndexFactory.
  /// @param space an IndexSpace object
  /// @return a unique temporary index in space @c space
  static Index make_tmp_index(const IndexSpace& space) {
    return Index(IndexSpace::base_key(space) + L'_' + std::to_wstring(Index::next_tmp_index()), &space);
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
  const auto& proto_indices() const { return proto_indices_; }

  std::wstring to_latex() const {
    auto protect_subscript = [](const std::wstring_view str) {
      auto subsc_pos = str.find(L'_');
      if (subsc_pos == std::wstring_view::npos)
        return std::wstring(str);
      else {
        assert(subsc_pos + 1 < str.size());
        if (subsc_pos + 2 == str.size())  // don't protect single character
          return std::wstring(str);
        std::wstring
            result = std::wstring(str.substr(0, subsc_pos + 1)) + L"{" + std::wstring(str.substr(subsc_pos + 1)) + L"}";
        return result;
      }
    };

    std::wstring result;
    result = L"{";
    result += protect_subscript(this->label());
    if (this->has_proto_indices()) {
      result += L"^{";
      for (const auto &pi: this->proto_indices()) {
        result += pi.to_latex();
      }
      result += L"}";
    }
    result += L"}";
    return result;
  }

  /// @return the smallest index of a generated index
  static constexpr std::size_t min_tmp_index() {
    return 100;
  }

  /// @return a unique temporary index, its value is equal to or greater than that
  static std::size_t next_tmp_index() {
    static std::atomic<std::size_t> index = min_tmp_index() - 1;
    return ++index;
  }

 private:
  std::wstring label_{};
  const IndexSpace *space_{};          // pointer to allow default initialization
  container::vector<Index> proto_indices_{}; // an unordered set of unique indices on
  // which this index depends on

  Index(std::wstring_view label, const IndexSpace *space) : label_(label), space_(space) {}

  /// throws std::invalid_argument if have duplicates in proto_indices_
  inline void check_for_duplicate_proto_indices();

  /// throws std::invalid_argument if label_ is in reserved
  void check_nontmp_label() {
    const auto index = label_index(label_);
    if (index && index > min_tmp_index()) {
      throw std::invalid_argument("Index ctor: label index must be less than the value returned by min_tmp_index()");
    }
  }

  static std::optional<std::size_t> label_index(std::wstring_view label) {
    const auto underscore_position = label.find(L'_');
    if (underscore_position != std::wstring::npos) {
      assert(underscore_position+1 < label.size());  // check that there is at least one char past the underscore
      return std::wcstol(label.substr(underscore_position+1, std::wstring::npos).data(), NULL, 10);
    }
    else
      return {};
  }

  friend class IndexFactory;
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
  container::vector<Index const *> vp;
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
  return index.to_latex();
}

class IndexSwapper {
 public:
  IndexSwapper() : even_num_of_swaps_(true) {}
  static IndexSwapper &thread_instance() {
    static thread_local IndexSwapper instance_{};
    return instance_;
  }

  bool even_num_of_swaps() const { return even_num_of_swaps_; }
  void reset() { even_num_of_swaps_ = true; }

 private:
  std::atomic<bool> even_num_of_swaps_;
  void toggle() { even_num_of_swaps_ = !even_num_of_swaps_; }

  friend inline void swap(Index &, Index &);
};

/// swap operator helps tracking # of swaps
inline void swap(Index &first, Index &second) {
  std::swap(first, second);
  IndexSwapper::thread_instance().toggle();
}

/// Generates temporary indices
class IndexFactory {
 public:

  /// creates a temporary index in space @c space . The label of the resulting index =
  /// @c IndexSpace::base_key(space) + '_' + temporary counter.
  /// Each call increments the current tmp counter (see next_tmp_index() ) . To make cleaner temporary indices
  /// with one counter per space use IndexFactory
  /// @param space an IndexSpace object
  /// @return a unique temporary index in space @c space
  Index make(const IndexSpace &space) {
    auto counter_it = counters_.begin();
    {  // if don't have a counter for this space
      std::scoped_lock lock(mutex_);
      if ((counter_it = counters_.find(space)) == counters_.end()) {
        counters_[space] = Index::min_tmp_index() - 1;
        counter_it = counters_.find(space);
      }
    }
    return Index(IndexSpace::base_key(space) + L'_' + std::to_wstring(++(counter_it->second)), &space);
  }

 private:
  std::mutex mutex_;
  std::map<IndexSpace, std::atomic<std::size_t>> counters_;
};

inline auto hash_value(const Index &idx) {
  return boost::hash_value(idx.label());
}

} // namespace sequant2

#endif // SEQUANT2_INDEX_H
