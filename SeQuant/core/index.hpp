//
// Created by Eduard Valeyev on 3/20/18.
//

#ifndef SEQUANT_INDEX_H
#define SEQUANT_INDEX_H

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/hash.hpp>
#include <iostream>
// #include <SeQuant/core/space.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/tag.hpp>
#include <SeQuant/core/utility/string.hpp>
// Only needed due to a (likely) compiler bug in Apple Clang
// #include <SeQuant/core/attr.hpp>

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <cwchar>
#include <functional>
#include <initializer_list>
#include <iterator>
#include <map>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <range/v3/all.hpp>

// change to 1 to make thread-safe
#define SEQUANT_INDEX_THREADSAFE 1

namespace sequant {

class Index;
using WstrList = std::initializer_list<std::wstring_view>;
using IndexList = std::initializer_list<Index>;

/// @brief Index = label + IndexSpace
/// @note Unlike SeQuant1's ParticleIndex, this Index supports dependencies
/// between indices to be able to express
///       e.g. hierarchical partitioning of index spaces or hierarchical nesting
///       of spaces
/// @note Index label can be plain (`label`) or composite (`label_ordinal`)
/// where `label` is a string of characters excluding '_', and `ordinal` is
/// an integer less than the value returned by min_tmp_label() .
class Index : public Taggable {
  static auto &tmp_index_accessor() {
    // initialized so that the first call to next_tmp_index will return
    // min_tmp_index()
    static std::atomic<std::size_t> index = min_tmp_index() - 1;
    return index;
  }

 public:
  using index_vector = container::vector<Index>;

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

  /// @param label the label, does not need to be unique
  /// @param space (a const ref to) the IndexSpace object that specifies to this
  /// space this object belongs
  /// @param proto_index labels of proto indices (all must be unique,
  /// i.e. duplicates are not allowed)
  /// @param symmetric_proto_indices if true, proto_indices can be permuted at
  /// will and will always be sorted
  Index(std::wstring_view label, const IndexSpace &space,
        container::vector<Index> proto_indices,
        bool symmetric_proto_indices = true)
      : label_(label),
        space_(space),
        proto_indices_(std::move(proto_indices)),
        symmetric_proto_indices_(symmetric_proto_indices) {
    canonicalize_proto_indices();
    check_for_duplicate_proto_indices();
    check_nontmp_label();
  }

  /// @param label the index label, does not need to be unique, but must be
  /// convertible into an IndexSpace (@sa IndexSpace::instance )
  Index(const std::wstring_view label)
      : Index(label,
              get_default_context().index_space_registry()->retrieve(label),
              {}) {
    check_nontmp_label();
  }

  /// @param label the index label, does not need to be unique, but must be
  /// convertible into an IndexSpace (@sa IndexSpace::instance )
  template <size_t N>
  Index(const wchar_t (&label)[N])
      : Index(std::wstring_view(&label[0]),
              get_default_context().index_space_registry()->retrieve(&label[0]),
              {}) {
    check_nontmp_label();
  }

  /// @param label the index label, does not need to be unique, but must be
  /// convertible into an IndexSpace (@sa IndexSpace::instance )
  Index(const wchar_t *label)
      : Index(std::wstring_view(label),
              get_default_context().index_space_registry()->retrieve(label),
              {}) {
    check_nontmp_label();
  }

  /// @brief constructs an Index using an existing Index's label and space and a
  /// list of proto indices

  /// @tparam IndexOrIndexLabel either Index or std::wstring or
  /// std::wstring_view
  /// @tparam I either Index or a type that can be converted to Index
  /// @param label the label, does not need to be unique
  /// @param proto_indices list of proto indices, or their labels (all must be
  /// unique, i.e. duplicates are not allowed)
  /// @param symmetric_proto_indices if true, proto_indices can be permuted at
  /// will and will always be sorted
  template <
      typename IndexOrIndexLabel, typename I,
      typename = std::enable_if_t<
          (std::is_same_v<std::decay_t<IndexOrIndexLabel>, Index> ||
           meta::is_wstring_convertible_v<std::decay_t<IndexOrIndexLabel>>)>>
  Index(IndexOrIndexLabel &&index, std::initializer_list<I> proto_indices,
        bool symmetric_proto_indices = true)
      : symmetric_proto_indices_(symmetric_proto_indices) {
    if constexpr (!std::is_same_v<std::decay_t<IndexOrIndexLabel>, Index>) {
      label_ = index;
      space_ = get_default_context().index_space_registry()->retrieve(label_);
    } else {
      label_ = index.label();
      space_ = index.space();
    }
    if constexpr (!std::is_same_v<std::decay_t<I>, Index>) {
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

  /// @brief constructs an Index using an existing Index's label and space and a
  /// list of proto indices

  /// @tparam IndexOrIndexLabel either Index or std::wstring or
  /// std::wstring_view
  /// @param label the label, does not need to be unique
  /// @param proto_indices list of proto indices (all must be unique,
  /// i.e. duplicates are not allowed)
  /// @param symmetric_proto_indices if true, proto_indices can be permuted at
  ///  will and will always be sorted
  template <
      typename IndexOrIndexLabel, typename IndexContainer,
      typename = std::enable_if_t<
          std::is_convertible_v<std::remove_reference_t<IndexContainer>,
                                container::vector<Index>> &&
          (std::is_same_v<std::decay_t<IndexOrIndexLabel>, Index> ||
           meta::is_wstring_convertible_v<std::decay_t<IndexOrIndexLabel>>)>>
  Index(IndexOrIndexLabel &&index, IndexContainer &&proto_indices,
        bool symmetric_proto_indices = true)
      : proto_indices_(std::forward<IndexContainer>(proto_indices)),
        symmetric_proto_indices_(symmetric_proto_indices) {
    if constexpr (!std::is_same_v<std::decay_t<IndexOrIndexLabel>, Index>) {
      label_ = index;
      check_nontmp_label();
      space_ = get_default_context().index_space_registry()->retrieve(label_);
    } else {
      label_ = index.label();
      space_ = index.space();
    }
    canonicalize_proto_indices();
    check_for_duplicate_proto_indices();
  }

  /// @brief constructs an Index using an existing Index's label and proto
  /// indices (if any) and an IndexSpace

  /// @tparam IndexOrIndexLabel either Index or std::wstring or
  /// std::wstring_view
  /// @param[in] index an Index object
  /// @param space an IndexSpace object
  /// @param label the label, does not need to be unique
  /// @param space (a const ref to) the IndexSpace object that specifies to this
  /// space this object belongs
  template <typename IndexOrIndexLabel>
  Index(IndexOrIndexLabel &&index_or_index_label, IndexSpace space) {
    if constexpr (std::is_same_v<IndexOrIndexLabel, Index>) {
      label_ = std::move(index_or_index_label.label());
      space_ = std::move(space);
      proto_indices_ = std::move(index_or_index_label.proto_indices_);
      symmetric_proto_indices_ = index_or_index_label.symmetric_proto_indices_;
      canonicalize_proto_indices();
      check_for_duplicate_proto_indices();
    } else if constexpr (std::is_same_v<std::decay_t<IndexOrIndexLabel>,
                                        Index>) {
      label_ = index_or_index_label.label();
      space_ = std::move(space);
      proto_indices_ = index_or_index_label.proto_indices_;
      symmetric_proto_indices_ = index_or_index_label.symmetric_proto_indices_;
      canonicalize_proto_indices();
      check_for_duplicate_proto_indices();
    } else {
      label_ = index_or_index_label;
      space_ = std::move(space);
    }
    check_nontmp_label();
  }

  /// @brief constructs an Index using this object's label and proto indices (if
  /// any) and a new IndexSpace
  /// @param space an IndexSpace object
  [[nodiscard]] Index replace_space(IndexSpace space) const {
    return Index(*this, std::move(space));
  }

  /// @brief constructs an Index using this object's label and proto indices (if
  /// any), its IndexSpaceType, and a new set of QuantumNumbers
  [[nodiscard]] Index replace_qns(QuantumNumbersAttr qns) const {
    return Index(*this, IndexSpace(this->space().get_base_key(),
                                   this->space().attr(), std::move(qns)));
  }

  /// @return this cast to Taggable&
  Taggable &tag() { return static_cast<Taggable &>(*this); }
  /// @return this cast to const Taggable&
  const Taggable &tag() const { return static_cast<const Taggable &>(*this); }
  /// resets tag of this and its protoindices (if any)
  /// @note do @c this->tag().reset() if you only want to reset tag on this (not
  /// its protoindices)
  void reset_tag() const {
    this->tag().reset();
    ranges::for_each(proto_indices_,
                     [](const Index &idx) { idx.tag().reset(); });
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
    result.label_ =
        space.get_base_key() + L'_' + std::to_wstring(Index::next_tmp_index());
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
    result.label_ =
        space.get_base_key() + L'_' + std::to_wstring(Index::next_tmp_index());
    result.space_ = space;
    result.proto_indices_ = std::forward<IndexContainer>(proto_indices);
    result.symmetric_proto_indices_ = symmetric_proto_indices;
    result.canonicalize_proto_indices();
    result.check_for_duplicate_proto_indices();
    return result;
  }

  /// @return the label as a UTF-8 encoded wide-character string
  /// @warning this does not include the proto index labels, use
  /// Index::full_label() instead
  std::wstring_view label() const { return label_; }

  /// @return A string label representable in ASCII encoding
  /// @warning not to be used with proto indices
  /// @brief Replaces non-ascii wstring characters with human-readable analogs,
  ///        each such UTF-8 character will be encoded by one or more chars.
  /// @note Maps: `↑` -> `a`, `↓` -> `b`, and all greek characters to their
  ///       english language equivalents (e.g. `α` -> `alpha`, `Ξ` -> `XI`,
  ///       etc.)
  [[deprecated(
      "use to_string to produce TiledArray-compatible index label "
      "representation")]] std::string
  ascii_label() const;

  /// @return A UTF-8 encoded narrow-character string label
  /// @warning not to be used with proto indices
  /// @note equivalent to `sequant::to_string(this->label())`
  std::string to_string() const;

  /// @return the full label as a UTF-8 encoded wide-character string
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

  /// @brief makes a new label by appending a suffix to the label

  /// Appends @p suffix to the label itself (if plain) or to its core (if
  /// composite)
  /// @param suffix a string to append to the label
  /// @return `this->label()` with @p suffix appended
  template <typename WS, typename = std::enable_if_t<(
                             meta::is_wstring_convertible_v<std::decay_t<WS>>)>>
  [[nodiscard]] std::wstring make_label_plus_suffix(WS &&suffix) const {
    return Index::make_label_plus_suffix(this->label(),
                                         std::forward<WS>(suffix));
  }

  /// @brief makes a new label by appending a suffix to the label

  /// Appends @p suffix to the label itself (if plain) or to its core (if
  /// composite)
  /// @param suffix a string to append to the label
  /// @return `this->label()` with @p suffix appended
  template <typename WS1, typename WS2,
            typename = std::enable_if_t<
                (meta::is_wstring_or_view_v<std::decay_t<WS1>> &&
                 meta::is_wstring_convertible_v<std::decay_t<WS2>>)>>
  [[nodiscard]] static std::wstring make_label_plus_suffix(WS1 &&label,
                                                           WS2 &&suffix) {
    auto underscore_position = label.find(L'_');
    std::wstring result;
    if (underscore_position == std::wstring::npos) {
      result = label;
      result += suffix;
    } else {
      result = label.substr(0, underscore_position);
      result += suffix;
      result += label.substr(underscore_position);
    }
    return result;
  }

  /// @brief makes a new label by removing a substring from the label

  /// Removes @p substr from the label itself (if plain) or from its core (if
  /// composite)
  /// @param substr a string to remove from the label
  /// @return `this->label()` with @p substr removed
  template <typename WS, typename = std::enable_if_t<(
                             meta::is_wstring_convertible_v<std::decay_t<WS>>)>>
  [[nodiscard]] std::wstring make_label_minus_substring(WS &&substr) const {
    return Index::make_label_minus_substring(this->label(),
                                             std::forward<WS>(substr));
  }

  /// @brief makes a new label by removing a substring from the label

  /// Removes @p substr from the label itself (if plain) or from its core (if
  /// composite)
  /// @param substr a string to remove from the label
  /// @return `this->label()` with @p substr removed
  template <typename WS1, typename WS2,
            typename = std::enable_if_t<
                (meta::is_wstring_or_view_v<std::decay_t<WS1>> &&
                 meta::is_wstring_convertible_v<std::decay_t<WS2>>)>>
  [[nodiscard]] static std::wstring make_label_minus_substring(WS1 &&label,
                                                               WS2 &&substr) {
    auto underscore_position = label.find(L'_');
    std::wstring result;

    auto erase = [](auto &result, const auto &substr) {
      auto pos = result.find(substr);
      if (pos != std::wstring::npos) {
        if constexpr (std::is_same_v<std::decay_t<WS2>, std::wstring> ||
                      std::is_same_v<std::decay_t<WS2>, std::wstring_view>) {
          result.erase(pos, substr.size());
        } else if constexpr (std::is_same_v<std::decay_t<WS2>,
                                            const wchar_t[]> ||
                             std::is_same_v<std::decay_t<WS2>, wchar_t[]> ||
                             std::is_same_v<std::decay_t<WS2>,
                                            const wchar_t *> ||
                             std::is_same_v<std::decay_t<WS2>, wchar_t *>) {
          result.erase(pos, std::strlen(substr));
        } else {
          result.erase(pos, 1);
        }
      }
    };

    if (underscore_position == std::wstring::npos) {
      result = label;
      erase(result, substr);
    } else {
      result = label.substr(0, underscore_position);
      erase(result, substr);
      result += label.substr(underscore_position);
    }
    return result;
  }

  /// @return the IndexSpace object
  const IndexSpace &space() const {
    assert(space_.attr().is_valid());
    return space_;
  }

  /// @return true if this index has proto indices
  bool has_proto_indices() const { return !proto_indices_.empty(); }
  /// @return the list of proto indices of this index
  const index_vector &proto_indices() const { return proto_indices_; }
  /// @return true if the index is symmetric with respect to the permutation of
  /// protoindices
  bool symmetric_proto_indices() const { return symmetric_proto_indices_; }
  /// drops the proto indices from this Index
  /// @return a copy of this Index without proto indices
  Index drop_proto_indices() const {
    return Index(this->label(), this->space());
  }

  std::wstring to_latex() const;

  /*template <typename... Attrs>
  std::wstring to_wolfram(Attrs &&...attrs) const {
    auto protect_subscript = [](const std::wstring_view str) {
      auto subsc_pos = str.rfind(L'_');
      if (subsc_pos == std::wstring_view::npos)
        return std::wstring(str);
      else {
        assert(subsc_pos + 1 < str.size());
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
    result += L","s + ::sequant::to_wolfram(space());
    ((result += ((L","s + ::sequant::to_wolfram(std::forward<Attrs>(attrs))))),
     ...);
    result += L"]";
    return result;
  }*/

  /// @param protoindex_range a range of Index objects
  /// @return the color of the protoindices
  /// @sa Index::color()
  template <typename Range>
  static auto proto_indices_color(const Range &protoindex_range) {
    auto space_attr_view =
        protoindex_range | ranges::views::transform([](const Index &idx) {
          return int64_t(idx.space().attr());
        });
    return hash::range(ranges::begin(space_attr_view),
                       ranges::end(space_attr_view));
  }

  /// @return the color of the protoindices
  /// @sa Index::color()
  auto proto_indices_color() const {
    return proto_indices_color(proto_indices_);
  }

  /// Color of an Index = hashed IndexSpace + IndexSpace objects of the
  /// protoindices
  /// @return the color of this object
  auto color() const {
    if (has_proto_indices()) {
      auto result = proto_indices_color();
      hash::combine(result, int64_t(space().attr()));
      return result;
    } else {
      auto result = hash::value(int64_t(space().attr()));
      return result;
    }
  }

  /// @return the smallest index of a generated index
  static const std::size_t min_tmp_index();

  /// @return a unique temporary index, its value is equal to or greater than
  /// that returned by min_tmp_index()
  static std::size_t next_tmp_index() { return ++tmp_index_accessor(); }

  /// resets the temporary index counter so that the next call to
  /// next_tmp_index() will return the value returned by min_tmp_index()
  /// @warning should only to be used when reproducibility matters (e.g. unit
  /// testing)
  static void reset_tmp_index();

  /// @brief index replacement
  /// replaces this object with its image in the Index map.
  /// If this object was not found in the map, tries replacing its subindices.
  /// @param index_map maps Index to Index
  /// @return false if no replacements were made
  /// @pre  \code this->tag().has_value() == false || (this->tag().has_value()
  /// == true && this->tag().value<int>() == 0) \endcode
  /// @post if return value is true: \code this->tag().has_value() == true &&
  /// this->tag().value<int>() == 0 \endcode
  template <template <typename, typename, typename... Args> class Map,
            typename... Args>
  bool transform(const Map<Index, Index, Args...> &index_map) {
    bool mutated = false;

    // outline:
    // - try replacing this first
    //   - if this is replaced by an index with protoindices, the protoindices
    //   should not be tagged since they are original and may need to be
    //   replaced also
    // - if not found, try replacing protoindices
    // - if protoindices mutated, try replacing this again

    // is this tagged already? if yes, can't skip, the protoindices may need to
    // be transformed also
    const auto this_is_tagged = this->tag().has_value();
    // sanity check that tag = 0
    if (this_is_tagged) {
      assert(this->tag().value<int>() == 0);
    } else {  // only try replacing this if not already tagged
      for (auto it = index_map.begin(); it != index_map.end(); ++it) {
        LabelCompare comp;
        bool test_comp_ab = comp(it->first, it->second);
        bool test_comp_ba = comp(it->second, it->first);
        std::wcout << "Key: " << it->first.label() << std::endl;
        std::wcout << "Value: " << it->second.label() << std::endl;
      }
      auto it = index_map.find(*this);
      if (it != index_map.end()) {
        *this = it->second;
        this->tag().assign(0);
        mutated = true;
      }
    }

    if (!mutated) {
      bool proto_indices_transformed = false;
      for (auto &&subidx : proto_indices_) {
        if (subidx.transform(index_map)) proto_indices_transformed = true;
      }
      if (proto_indices_transformed) {
        mutated = true;
        canonicalize_proto_indices();
        if (!this_is_tagged) {  // if protoindices were mutated, try again, but
                                // only if no tag yet
          auto it = index_map.find(*this);
          if (it != index_map.end()) {
            *this = it->second;
            this->tag().assign(0);
            mutated = true;
          }
        }
      }
    }
    if (mutated) {
      full_label_.reset();
    }
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
  // an unordered set of unique indices on which this index depends on
  // whether proto_indices_ is symmetric w.r.t. permutations; if true,
  // proto_indices_ will be ordered
  index_vector proto_indices_{};
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
    const auto underscore_position = label.rfind(L'_');
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

  /// @return true if @c index1 is identical to @c index2 , i.e. they belong to
  /// the same space, they have the same label, and the same proto-indices (if
  /// any)
  friend bool operator==(const Index &i1, const Index &i2) {
    return i1.space() == i2.space() && i1.label() == i2.label() &&
           i1.proto_indices() == i2.proto_indices();
  }

  /// @return false if @c index1 is identical to @c index2 , i.e. they belong to
  /// different spaces or they have different labels or they have different
  /// proto-indices (if any)
  friend bool operator!=(const Index &i1, const Index &i2) {
    return !(i1 == i2);
  }

  /// @brief The ordering operator

  /// @return true if @c i1 preceeds @c i2 in the canonical order; Index objects
  /// are ordered lexicographically, first by qns, followed by tags (if defined
  /// for both), then by space, then by label, then by protoindices (if any)
  friend bool operator<(const Index &i1, const Index &i2) {
    // compare qns, tags and spaces in that sequence
    assert(i1.space().attr().is_valid());
    assert(i2.space().attr().is_valid());

    auto i1_Q = i1.space().qns();
    auto i2_Q = i2.space().qns();
    const bool have_qns =
        i1_Q != IndexSpace::nullqns || i2_Q != IndexSpace::nullqns;

    auto compare_space = [&i1, &i2]() {
      if (i1.space() == i2.space()) {
        if (i1.label() == i2.label()) {
          return i1.proto_indices() < i2.proto_indices();
        } else {
          return i1.label() < i2.label();
        }
      } else {
        return i1.space() < i2.space();
      }
    };

    if (have_qns || (i1_Q != i2_Q)) {
      return compare_space();
    }
    const bool have_tags = i1.tag().has_value() && i2.tag().has_value();

    if (!have_tags || i1.tag() == i2.tag()) {
      return compare_space();
    } else {
      return i1.tag() < i2.tag();
    }
  }

};  // class Index

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
#if SEQUANT_INDEX_THREADSAFE
        std::scoped_lock lock(mutex_);
#endif
        if ((counter_it = counters_.find(space)) == counters_.end()) {
          counters_[space] = min_index_ - 1;
          counter_it = counters_.find(space);
        }
      }
      result = Index(
          space.get_base_key() + L'_' + std::to_wstring(++(counter_it->second)),
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
#if SEQUANT_INDEX_THREADSAFE
        std::scoped_lock lock(mutex_);
#endif
        if ((counter_it = counters_.find(space)) == counters_.end()) {
          counters_[space] = min_index_ - 1;
          counter_it = counters_.find(space);
        }
      }
      result = Index(Index(space.get_base_key() + L'_' +
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
#if SEQUANT_INDEX_THREADSAFE
  std::mutex mutex_;
  // boost::container::flat_map needs copyable value, which std::atomic is not,
  // so must use std::map
  std::map<IndexSpace, std::atomic<std::size_t>> counters_;
#else
  // until multithreaded skip atomic
  container::map<IndexSpace, std::size_t> counters_;
#endif
};

/// @brief hashing function

/// @paramp[in] idx a const reference to an Index object
/// @return the hash value of the object referred to by idx
inline auto hash_value(const Index &idx) {
  const auto &proto_indices = idx.proto_indices();
  using std::begin;
  using std::end;
  auto val = hash::range(begin(proto_indices), end(proto_indices));
  hash::combine(val, idx.label());
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

class IndexRegistry {
 public:
  using Record =
      std::tuple<std::function<long(const Index &)>>;  // index record = {sizer}

  IndexRegistry() = default;

  /// updates an existing entry, or creates a new one if it does not exist
  template <typename... Args>
  void update(const Index &idx, Args &&...args) {
    auto it = registry_.find(idx);
    if (it != registry_.end()) {
      registry_.erase(it);
    }
    auto insertion_result =
        registry_.try_emplace(idx, std::forward<Args>(args)...);
  }
  /// creates a new entry
  template <typename... Args>
  void make(const Index &idx, Args &&...args) {
    auto insertion_result =
        registry_.try_emplace(idx, std::forward<Args>(args)...);
    assert(insertion_result.second);
  }

  /// retrieves the pointer to the Record object for Index @idx , or nullptr if
  /// not found
  const Record *retrieve(const Index &idx) const {
    auto result = registry_.find(idx);
    if (result != registry_.end())
      return &(result->second);
    else
      return nullptr;
  }

 private:
  container::map<Index, Record> registry_;
};

}  // namespace sequant

#endif  // SEQUANT_INDEX_H
