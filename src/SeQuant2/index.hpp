//
// Created by Eduard Valeyev on 3/20/18.
//

#ifndef SEQUANT2_INDEX_H
#define SEQUANT2_INDEX_H

#include <atomic>
#include <mutex>
#include <optional>
#include <set>
#include <string>

#include <range/v3/all.hpp>

#include <boost/container_hash/hash.hpp>

#include "attr.hpp"
#include "container.hpp"
#include "space.hpp"
#include "tag.hpp"

// change to 1 to make thread-safe
#define SEQUANT2_INDEX_THREADSAFE 1

namespace sequant2 {

class Index;
using WstrList = std::initializer_list<std::wstring_view>;
using IndexList = std::initializer_list<Index>;

/// @brief Index = label + IndexSpace
/// @note Unlike SeQuant1's ParticleIndex, this Index supports dependencies
/// between indices to be able to express
///       e.g. hierarchical partitioning of index spaces or hiearchical nesting
///       of spaces
/// @note label has format "label_index" where "label" is a string of characters
/// excluding '_', and "index"
///       is an integer less than the value returned by min_tmp_label() .
class Index : public Taggable {
 public:
  Index() = default;

  /// @param label the label, does not need to be unique
  /// @param space (a const ref to) the IndexSpace object that specifies to this
  /// space this object belongs
  /// @param proto_index labels of proto indices (all must be unique,
  /// i.e. duplicates are not allowed)
  /// @param symmetric_proto_indices if true, proto_indices can be permuted at
  /// will and will always be sorted
  Index(std::wstring_view label, const IndexSpace &space,
        IndexList proto_indices, bool symmetric_proto_indices = true)
      : label_(label),
        space_(space),
        proto_indices_(proto_indices),
        symmetric_proto_indices_(symmetric_proto_indices) {
    canonicalize_proto_indices();
    check_for_duplicate_proto_indices();
    check_nontmp_label();
  }

  /// @param label the index label, does not need to be unique, but must be
  /// convertible into an IndexSpace (@sa IndexSpace::instance )
  Index(const std::wstring_view label)
      : Index(label, IndexSpace::instance(label), {}) {
    check_nontmp_label();
  }

  /// @param label the index label, does not need to be unique, but must be
  /// convertible into an IndexSpace (@sa IndexSpace::instance )
  template <size_t N>
  Index(const wchar_t (&label)[N])
      : Index(std::wstring_view(&label[0]), IndexSpace::instance(&label[0]),
              {}) {
    check_nontmp_label();
  }

  /// @param label the index label, does not need to be unique, but must be
  /// convertible into an IndexSpace (@sa IndexSpace::instance )
  Index(const wchar_t *label)
      : Index(std::wstring_view(label), IndexSpace::instance(label), {}) {
    check_nontmp_label();
  }

  /// @param label the label, does not need to be unique
  /// @param space (a const ref to) the IndexSpace object that specifies to this
  /// space this object belongs
  /// @note Index does not manage lifetime of index spaces, so it is user's
  /// responsibility to ensure that @c space will
  ///       outlive this object
  Index(std::wstring_view label, const IndexSpace &space) noexcept
      : label_(label), space_(space), proto_indices_() {
    check_nontmp_label();
  }

  /// @param label the label, does not need to be unique
  /// @param proto_indices list of proto indices, or their labels (all must be
  /// unique, i.e. duplicates are not allowed)
  /// @param symmetric_proto_indices if true, proto_indices can be permuted at
  /// will and will always be sorted
  template <typename I1, typename I2>
  Index(I1 &&index, std::initializer_list<I2> proto_indices,
        bool symmetric_proto_indices = true)
      : symmetric_proto_indices_(symmetric_proto_indices) {
    if constexpr (!std::is_same_v<std::decay_t<I1>, Index>) {
      label_ = index;
      space_ = IndexSpace::instance(label_);
    } else {
      label_ = index.label();
      space_ = index.space();
    }
    if constexpr (!std::is_same_v<std::decay_t<I2>, Index>) {
      if (proto_indices.size() != 0) {
        proto_indices_.reserve(proto_indices.size());
        for (const auto &plabel : proto_indices)
          proto_indices_.push_back(Index(plabel));
      }
    } else
      proto_indices_ = proto_indices;
    canonicalize_proto_indices();
    check_for_duplicate_proto_indices();
    check_nontmp_label();
  }

  /// @param label the label, does not need to be unique
  /// @param proto_indices list of proto indices (all must be unique,
  /// i.e. duplicates are not allowed)
  /// @param symmetric_proto_indices if true, proto_indices can be permuted at
  ///  will and will always be sorted
  template <
      typename I1, typename IndexContainer,
      typename = std::enable_if_t<std::is_convertible_v<
          std::remove_reference_t<IndexContainer>, container::vector<Index>>>>
  Index(I1 &&index, IndexContainer &&proto_indices,
        bool symmetric_proto_indices = true)
      : proto_indices_(std::forward<IndexContainer>(proto_indices)),
        symmetric_proto_indices_(symmetric_proto_indices) {
    if constexpr (!std::is_same_v<std::decay_t<I1>, Index>) {
      label_ = index;
      check_nontmp_label();
      space_ = IndexSpace::instance(label_);
    } else {
      label_ = index.label();
      space_ = index.space();
    }
    canonicalize_proto_indices();
    check_for_duplicate_proto_indices();
  }

  /// creates a globally-unique temporary index in space @c space . The label of
  /// the resulting index =
  /// @c IndexSpace::base_key(space) + '_' + temporary counter.
  /// Each call increments the current tmp counter (see next_tmp_index() ) . To
  /// make neater temporary indices unique in a given scope (e.g. a single term
  /// in an expression) use IndexFactory.
  /// @param space an IndexSpace object
  /// @return a unique temporary index in space @c space
  static Index make_tmp_index(const IndexSpace &space) {
    Index result;
    result.label_ = IndexSpace::base_key(space) + L'_' +
                    std::to_wstring(Index::next_tmp_index());
    result.space_ = space;
    return result;
  }

  /// creates a globaly-unique temporary index in space @c space . The label of
  /// the resulting index =
  /// @c IndexSpace::base_key(space) + '_' + temporary counter.
  /// Each call increments the current tmp counter (see next_tmp_index() ) . To
  /// make neater temporary indices unique in a given scope (e.g. a single term
  /// in an expression) use IndexFactory.
  /// @param space an IndexSpace object
  /// @param proto_indices list of proto indices (all must be unique,
  /// i.e. duplicates are not allowed)
  /// @param symmetric_proto_indices if true, proto_indices can be permuted at
  ///  will and will always be sorted
  /// @return a unique temporary index in space @c space
  template <
      typename IndexContainer,
      typename = std::enable_if_t<std::is_convertible_v<
          std::remove_reference_t<IndexContainer>, container::vector<Index>>>>
  static Index make_tmp_index(const IndexSpace &space,
                              IndexContainer &&proto_indices,
                              bool symmetric_proto_indices = true) {
    Index result;
    result.label_ = IndexSpace::base_key(space) + L'_' +
                    std::to_wstring(Index::next_tmp_index());
    result.space_ = space;
    result.proto_indices_ = std::forward<IndexContainer>(proto_indices);
    result.symmetric_proto_indices_ = symmetric_proto_indices;
    result.canonicalize_proto_indices();
    result.check_for_duplicate_proto_indices();
    return result;
  }

  /// @return the label
  /// @warning this does not include the proto index labels, use
  /// Index::full_label() instead
  std::wstring_view label() const { return label_; }
  /// @return the full label
  /// @warning this includes the proto index labels (if any), use
  /// Index::label() instead if only want the label
  std::wstring_view full_label() const {
    if (!has_proto_indices()) return label();
    if (full_label_) return *full_label_;
    std::wstring result = label_;
    ranges::for_each(proto_indices_, [&result](const Index &idx) {
      result += idx.full_label();
    });
    full_label_ = result;
    return *full_label_;
  }
  /// @return the IndexSpace object
  const IndexSpace &space() const {
    assert(space_.attr().is_valid());
    return space_;
  }

  /// @return true if this index has proto indices
  bool has_proto_indices() const { return !proto_indices_.empty(); }
  /// @return the list of proto indices of this index
  const auto &proto_indices() const { return proto_indices_; }

  std::wstring to_latex() const {
    auto protect_subscript = [](const std::wstring_view str) {
      auto subsc_pos = str.find(L'_');
      if (subsc_pos == std::wstring_view::npos)
        return std::wstring(str);
      else {
        assert(subsc_pos + 1 < str.size());
        if (subsc_pos + 2 == str.size())  // don't protect single character
          return std::wstring(str);
        std::wstring result = std::wstring(str.substr(0, subsc_pos + 1)) +
                              L"{" + std::wstring(str.substr(subsc_pos + 1)) +
                              L"}";
        return result;
      }
    };

    std::wstring result;
    result = L"{";
    result += protect_subscript(this->label());
    if (this->has_proto_indices()) {
      result += L"^{";
      for (const auto &pi : this->proto_indices()) {
        result += pi.to_latex();
      }
      result += L"}";
    }
    result += L"}";
    return result;
  }

  template <typename ... Attrs> std::wstring to_wolfram(Attrs && ... attrs) const {
    auto protect_subscript = [](const std::wstring_view str) {
      auto subsc_pos = str.find(L'_');
      if (subsc_pos == std::wstring_view::npos)
        return std::wstring(str);
      else{
        assert(subsc_pos + 1 < str.size());
        if (subsc_pos + 2 == str.size())  // don't protect single character
          return std::wstring(str);
        std::wstring result = L"\\!\\(\\*SubscriptBox[\\(";
        result += std::wstring(str.substr(0, subsc_pos));
        result += L"\\), \\(";
        result += std::wstring(str.substr(subsc_pos + 1));
        result += L"\\)]\\)";
        return result;
      }
    };

    using namespace std::literals;
    std::wstring result =
        L"particleIndex[\""s + protect_subscript(this->label()) + L"\"";
    if (this->has_proto_indices()) {
      assert(false && "not yet supported");
    }
    using namespace std::literals;
    result += L","s + ::sequant2::to_wolfram(space());
    ((result += ((L","s + ::sequant2::to_wolfram(std::forward<Attrs>(attrs))))),
     ...);
    result += L"]";
    return result;
  }

  /// @return the smallest index of a generated index
  static constexpr std::size_t min_tmp_index() { return 100; }

  /// @return a unique temporary index, its value is equal to or greater than
  /// that
  static std::size_t next_tmp_index() {
    static std::atomic<std::size_t> index = min_tmp_index() - 1;
    return ++index;
  }

  /// @brief index replacement
  /// replaces this object with its image in the Index map.
  /// If this object was not found in the map, tries replacing its subindices.
  /// @param index_map maps Index to Index
  /// @return false if no replacements were made
  template <template <typename, typename, typename... Args> class Map,
            typename... Args>
  bool transform(const Map<Index, Index, Args...> &index_map) {
    // first try replacing this first; if not found, try replacing protoindices
    // first, then try replacing this again
    auto it = index_map.find(*this);
    bool mutated = false;
    if (it != index_map.end()) {
      *this = it->second;
      mutated = true;
    } else {
      bool proto_indices_transformed = false;
      for (auto &&subidx : proto_indices_) {
        if (subidx.transform(index_map)) proto_indices_transformed = true;
      }
      if (proto_indices_transformed) {
        mutated = true;
        canonicalize_proto_indices();
      }
      auto it = index_map.find(*this);
      if (it != index_map.end()) {
        *this = it->second;
        mutated = true;
      }
    }
    if (mutated) full_label_.reset();
    return mutated;
  }

  /// compares Index objects using labels only
  struct LabelCompare {
    using is_transparent = void;
    bool operator()(const Index &first, const Index &second) const {
      return first.label() < second.label();
    }
    bool operator()(const Index &first, const std::wstring_view &second) const {
      return first.label() < second;
    }
    bool operator()(const std::wstring_view &first, const Index &second) const {
      return first < second.label();
    }
  };

  /// compares Index objects using type only (but since type is defined by the
  /// *values* of proto indices those are not ignored)
  struct TypeCompare {
    bool operator()(const Index &first, const Index &second) const {
      bool result;
      if (first.space() == second.space()) {
        result = first.proto_indices() < second.proto_indices();
      } else
        result = first.space() < second.space();
      return result;
    }
  };

  /// tests equality of Index objects using type only (but since type is defined
  /// by the *values* of proto indices those are not ignored)
  struct TypeEquality {
    bool operator()(const Index &first, const Index &second) const {
      bool result = (first.space() == second.space()) &&
                    (first.proto_indices() == second.proto_indices());
      return result;
    }
  };

 private:
  std::wstring label_{};
  IndexSpace space_{};
  container::vector<Index>
      proto_indices_{};  // an unordered set of unique indices on which this
                         // index depends on
  // whether proto_indices_ is symmetric w.r.t. permutations; if true,
  // proto_indices_ will be ordered
  bool symmetric_proto_indices_ = true;

  mutable std::optional<std::wstring> full_label_;

  /// sorts proto_indices_ if symmetric_proto_indices_
  inline void canonicalize_proto_indices();

  /// @warning disabled if NDEBUG is defined
  /// @throw std::invalid_argument  have duplicates in proto_indices_
  inline void check_for_duplicate_proto_indices();

  /// throws std::invalid_argument if label_ is in reserved
  void check_nontmp_label() {
    const auto index = label_index(label_);
    if (index && index > min_tmp_index()) {
      throw std::invalid_argument(
          "Index ctor: label index must be less than the value returned by "
          "min_tmp_index()");
    }
  }

  static std::optional<std::size_t> label_index(std::wstring_view label) {
    const auto underscore_position = label.find(L'_');
    if (underscore_position != std::wstring::npos) {
      assert(underscore_position + 1 <
             label.size());  // check that there is at least one char past the
                             // underscore
      return std::wcstol(
          label.substr(underscore_position + 1, std::wstring::npos).data(),
          NULL, 10);
    } else
      return {};
  }

  friend class IndexFactory;

  // this ctor is only used by make_tmp_index and IndexFactory and bypasses
  // check for nontmp index
  Index(std::wstring_view label, const IndexSpace *space) noexcept
      : label_(label), space_(*space), proto_indices_() {}
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
/// are ordered lexicographically, first by tags (if defined for both), then
/// by space, then by label, then by protoindices (if any)
inline bool operator<(const Index &i1, const Index &i2) {
  // compare tags first
  const bool have_tags = i1.has_tag() && i2.has_tag();
  if (!have_tags || i1.tag_equal(i2)) {
    if (i1.space() == i2.space()) {
      if (i1.label() == i2.label()) {
        return i1.proto_indices() < i2.proto_indices();
      } else {
        return i1.label() < i2.label();
      }
    } else {
      return i1.space() < i2.space();
    }
  } else {
    return i1.tag_less_than(i2);
  }
}

void Index::check_for_duplicate_proto_indices() {
#ifndef NDEBUG
  if (!symmetric_proto_indices_) {  // if proto indices not symmetric, sort via
                                    // ptrs
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
      throw std::invalid_argument(
          "Index ctor: duplicate proto indices detected");
    }
  } else {  // else search directly
    if (std::adjacent_find(begin(proto_indices_), end(proto_indices_)) !=
        proto_indices_.end()) {
      throw std::invalid_argument(
          "Index ctor: duplicate proto indices detected");
    }
  }
#endif
}

void Index::canonicalize_proto_indices() {
  if (symmetric_proto_indices_)
    std::stable_sort(begin(proto_indices_), end(proto_indices_));
}

inline std::wstring to_latex(const Index &index) { return index.to_latex(); }

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
  IndexFactory() = default;
  /// @tparam IndexValidator IndexValidator(const Index&) -> bool is valid and
  /// returns true generated index is valid
  /// @param min_index start indexing indices for each space with this value;
  /// must be greater than 0; the default is to use Index::min_tmp_index()
  template <typename IndexValidator>
  explicit IndexFactory(IndexValidator validator,
                        size_t min_index = Index::min_tmp_index())
      : min_index_(min_index), validator_(validator) {
    assert(min_index_ > 0);
  }

  /// creates a temporary index in space @c space . The label of the resulting
  /// index =
  /// @c IndexSpace::base_key(space) + '_' + temporary counter.
  /// Each call increments the current tmp counter (see next_tmp_index() ) .
  /// @param space an IndexSpace object
  /// @return a unique temporary index in space @c space
  Index make(const IndexSpace &space) {
    Index result;
    bool valid = false;
    do {
      auto counter_it = counters_.begin();
      {  // if don't have a counter for this space
#if SEQUANT2_INDEX_THREADSAFE
        std::scoped_lock lock(mutex_);
#endif
        if ((counter_it = counters_.find(space)) == counters_.end()) {
          counters_[space] = min_index_ - 1;
          counter_it = counters_.find(space);
        }
      }
      result = Index(IndexSpace::base_key(space) + L'_' +
                         std::to_wstring(++(counter_it->second)),
                     &space);
      valid = validator_ ? validator_(result) : true;
    } while (!valid);
    return result;
  }

  /// creates a temporary index that inherits the space and protoindices of @c
  /// idx . The label of the resulting index =
  /// @c IndexSpace::base_key(space) + '_' + temporary counter.
  /// Each call increments the current tmp counter (see next_tmp_index() ) .
  /// @param idx an Index object
  /// @return a unique temporary index in space @c space with same protoindices
  /// as @c idx
  Index make(const Index &idx) {
    const auto &space = idx.space();
    Index result;
    bool valid = false;
    do {
      auto counter_it = counters_.begin();
      {  // if don't have a counter for this space
#if SEQUANT2_INDEX_THREADSAFE
        std::scoped_lock lock(mutex_);
#endif
        if ((counter_it = counters_.find(space)) == counters_.end()) {
          counters_[space] = min_index_ - 1;
          counter_it = counters_.find(space);
        }
      }
      result = Index(Index(IndexSpace::base_key(space) + L'_' +
                               std::to_wstring(++(counter_it->second)),
                           &space),
                     idx.proto_indices());
      valid = validator_ ? validator_(result) : true;
    } while (!valid);
    return result;
  }

 private:
  std::size_t min_index_ = Index::min_tmp_index();
  std::function<bool(const Index &)> validator_ = {};
#if SEQUANT2_INDEX_THREADSAFE
  std::mutex mutex_;
  // boost::container::flat_map needs copyable value, which std::atomic is not,
  // so must use std::map
  std::map<IndexSpace, std::atomic<std::size_t>> counters_;
#else
  // until multithreaded skip atomic
  container::map<IndexSpace, std::size_t> counters_;
#endif
};

inline auto hash_value(const Index &idx) {
  const auto &proto_indices = idx.proto_indices();
  using std::begin;
  using std::end;
  auto val = boost::hash_range(begin(proto_indices), end(proto_indices));
  boost::hash_combine(val, idx.label());
  return val;
}

template <typename Container>
auto make_indices(WstrList index_labels = {}) {
  Container result;
  for (const auto &label : index_labels) {
    result.push_back(Index{label});
  }
  return result;
}

}  // namespace sequant2

#endif  // SEQUANT2_INDEX_H
