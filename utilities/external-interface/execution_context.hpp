#ifndef SEQUANT_EXTERNAL_INTERFACE_EXECUTIONCONTEXT_HPP
#define SEQUANT_EXTERNAL_INTERFACE_EXECUTIONCONTEXT_HPP

#include "processing_data.hpp"

#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <deque>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <ranges>
#include <set>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace sequant::util::extint {

/// A data object that a step is free to modify, obtained from
/// ExecutionContext::working_copy
template <typename DataT>
struct WorkingCopy {
  DataT data;
  /// The data object @c data was copied from, which remains untouched in the
  /// ExecutionContext. If this is null, @c data has been moved out of the
  /// context instead and the original no longer exists.
  const DataT *original = nullptr;

  [[nodiscard]] bool is_copy() const { return original != nullptr; }

  /// @returns Whether the original exists and compares equal to @c data. Hence,
  /// this is always false for a working copy that was moved out of the context.
  [[nodiscard]] bool unchanged() const {
    return original != nullptr && data == *original;
  }
};

class ExecutionContext {
 public:
  template <typename DataT>
  struct Data {
    /// The actual data object (reference)
    std::reference_wrapper<DataT> data;
    /// List of IDs that this data object is assigned to
    std::vector<std::string_view> associated_ids = {};
    /// List of IDs that are associated with a group of data objects which the
    /// present one is a part of
    std::vector<std::string_view> associated_group_ids = {};
    /// Position of the data object in the context's store
    std::size_t index = std::numeric_limits<std::size_t>::max();

    operator DataT &() { return data.get(); }
    operator std::add_const_t<DataT> &() const { return data.get(); }

    bool operator==(const Data &other) const {
      return data.get() == other.data.get();
    }
  };

  ExecutionContext() = default;

  void set_data(std::string_view prefix, std::size_t counter,
                ProcessingData data);

  template <std::ranges::range CounterValues, std::ranges::range Data>
    requires(std::integral<std::ranges::range_value_t<CounterValues>> &&
             std::same_as<std::ranges::range_value_t<Data>, ProcessingData>)
  void set_data(std::string_view prefix, CounterValues &&counters,
                Data &&data) {
    if (std::ranges::empty(counters)) {
      throw Exception("Attempted to add data without specifying any counter");
    }
    if (std::ranges::empty(data)) {
      throw Exception("Attempted to register empty dataset");
    }

    std::vector<std::string> ids;
    ids.reserve(std::ranges::size(counters));

    for (auto &&current : counters) {
      ids.emplace_back(std::string(prefix) + "." + std::to_string(current));

      if (id_to_data_indices_.find(ids.back()) != id_to_data_indices_.end()) {
        throw Exception("Duplicate data ID '" + ids.back() + "'");
      }
      if (!is_valid_id(ids.back(), false)) {
        throw Exception("Invalid id '" + ids.back() + "'");
      }
    }

    std::vector<std::size_t> data_indices(std::ranges::size(data));
    std::iota(data_indices.begin(), data_indices.end(), data_.size());

    for (std::string &current : ids) {
      set_id_idx_assoc(std::move(current), data_indices);
    }

    for (auto &&current : data) {
      data_.emplace_back(std::move(current));
    }
  }

  void add_data_alias(std::string_view id, std::string alias);

  template <std::ranges::range IDs>
    requires(std::is_convertible_v<std::ranges::range_value_t<IDs>,
                                   std::string_view>)
  void add_data_alias(IDs &&ids, std::string alias) {
    if (id_to_data_indices_.find(alias) != id_to_data_indices_.end()) {
      throw Exception("Alias '" + alias + "' already exists as a data ID");
    }

    for (const auto &current_id : ids) {
      for (const std::string &expanded : expand_id(current_id)) {
        auto it = id_to_data_indices_.find(expanded);

        if (it == id_to_data_indices_.end()) {
          throw Exception("Attempted to alias non-existent ID '" + expanded +
                          "'");
        }

        SEQUANT_ASSERT(!it->second.empty());
        set_id_idx_assoc(alias, it->second);
      }
    }
  }

  bool has_data(std::string_view id) const;

  bool ids_are_equivalent(std::string_view lhs, std::string_view rhs) const;

  template <std::ranges::range LhsIDs, std::ranges::range RhsIDs>
    requires(std::convertible_to<std::ranges::range_value_t<LhsIDs>,
                                 std::string_view> &&
             std::convertible_to<std::ranges::range_value_t<RhsIDs>,
                                 std::string_view>)
  bool ids_are_equivalent(LhsIDs &&lhs, RhsIDs &&rhs) {
    return get_data_indices(std::forward<LhsIDs>(lhs)) ==
           get_data_indices(std::forward<RhsIDs>(rhs));
  }

  template <std::ranges::range IDs>
    requires(
        std::convertible_to<std::ranges::range_value_t<IDs>, std::string_view>)
  bool ids_are_equivalent(IDs &&lhs, std::string_view rhs) {
    return ids_are_equivalent(std::forward<IDs>(lhs),
                              std::ranges::single_view(rhs));
  }

  template <std::ranges::range IDs>
    requires(
        std::convertible_to<std::ranges::range_value_t<IDs>, std::string_view>)
  bool ids_are_equivalent(std::string_view lhs, IDs &&rhs) {
    return ids_are_equivalent(std::ranges::single_view(lhs),
                              std::forward<IDs>(rhs));
  }

  std::size_t dataset_size(std::string_view id) const;

  std::vector<Data<const ProcessingData>> get_data(std::string_view id) const;
  std::vector<Data<ProcessingData>> get_data(std::string_view id);

  /// Installs a predicate that tells whether the data registered under a given
  /// ID is still needed once the currently running step has finished. This
  /// allows steps to find out whether they may consume their inputs (see
  /// working_copy).
  void set_needed_later_predicate(std::function<bool(std::string_view)> pred);

  /// Removes the predicate installed via set_needed_later_predicate
  void clear_needed_later_predicate();

  /// @returns Whether @p data may still be accessed after the currently running
  /// step. This is the case if the predicate set via
  /// set_needed_later_predicate holds for any of the IDs associated with the
  /// data. Without a predicate (e.g. when a step is run outside of an
  /// Executor), the data is conservatively assumed to be needed.
  template <typename DataT>
  bool is_needed_later(const Data<DataT> &data) const {
    return is_needed_later(data.index);
  }

  /// Provides the calling step with its own copy of @p data that it can modify
  /// at will. If the data is not needed later (see is_needed_later), it is
  /// moved out of the context instead of being copied. In this case, the data
  /// is no longer available from the context afterwards, so that any further
  /// attempt to retrieve it (e.g. by a second request for a working copy)
  /// throws.
  /// @tparam T The type of data that @p data is expected to hold
  /// @param keep_original If true, the data is always copied, so that the
  /// original remains available for the step to compare against
  template <typename T>
  WorkingCopy<T> working_copy(const Data<ProcessingData> &data,
                              bool keep_original = false) {
    T &payload = convert_data<T>(live_data(data.index));

    if (keep_original || is_needed_later(data.index)) {
      return {.data = clone(std::as_const(payload)), .original = &payload};
    }

    WorkingCopy<T> working{.data = std::move(payload)};
    consume(data.index);
    return working;
  }

  /// Frees all data objects for which @p is_needed returns false for every ID
  /// associated with them. The IDs themselves remain registered (keeping
  /// alias and group bookkeeping intact) but retrieving the released data via
  /// get_data throws.
  /// @returns The number of released data objects
  template <std::predicate<std::string_view> Predicate>
  std::size_t release_data(Predicate &&is_needed) {
    std::size_t n_released = 0;

    for (const auto &[idx, ids] : data_idx_to_ids_) {
      std::optional<ProcessingData> &slot = data_.at(idx);

      if (!slot.has_value() || std::ranges::any_of(ids, is_needed)) {
        continue;
      }

      slot.reset();
      ++n_released;
    }

    return n_released;
  }

  static bool is_valid_id(std::string_view id, bool allow_selectors);

  static std::vector<std::string> expand_id(std::string_view id);

 private:
  std::deque<std::optional<ProcessingData>> data_;
  std::map<std::string, std::vector<std::size_t>, std::less<>>
      id_to_data_indices_;
  std::map<std::size_t, std::deque<std::string>> data_idx_to_ids_;
  std::function<bool(std::string_view)> needed_later_;
  /// Positions of the data objects that were moved out by working_copy
  std::set<std::size_t> consumed_;

  bool is_needed_later(std::size_t index) const;

  /// @returns The data at the given position, throwing if it is no longer
  /// available
  ProcessingData &live_data(std::size_t index);

  /// Frees the data at the given position, recording that it was consumed
  void consume(std::size_t index);

  template <typename Indices>
  void set_id_idx_assoc(std::string id, Indices &&indices) {
    if constexpr (!std::ranges::range<Indices>) {
      set_id_idx_assoc(std::move(id),
                       std::ranges::single_view(std::move(indices)));
    } else {
      for (std::size_t idx : indices) {
        id_to_data_indices_[id].emplace_back(idx);
        data_idx_to_ids_[idx].emplace_back(id);
      }
    }
  }

  template <std::ranges::range IDs>
    requires(
        std::convertible_to<std::ranges::range_value_t<IDs>, std::string_view>)
  std::vector<std::size_t> get_data_indices(
      IDs &&ids, bool *contained_unknown_id = nullptr) const {
    std::vector<std::size_t> indices;

    if (contained_unknown_id) {
      *contained_unknown_id = false;
    }

    for (std::string_view id : ids) {
      if (!ExecutionContext::is_valid_id(id, true)) {
        throw Exception("Invalid id '" + std::string(id) + "'");
      }

      for (const std::string_view current : ExecutionContext::expand_id(id)) {
        auto it = id_to_data_indices_.find(current);

        if (it == std::ranges::end(id_to_data_indices_)) {
          if (contained_unknown_id) {
            *contained_unknown_id = true;
          } else {
            throw Exception("No data available for ID '" +
                            std::string(current) + "'");
          }

          continue;
        }

        std::ranges::copy_if(it->second, std::back_inserter(indices),
                             [&](std::size_t idx) {
                               return std::ranges::find(indices, idx) ==
                                      std::ranges::end(indices);
                             });
      }
    }

    std::ranges::sort(indices);

    return indices;
  }
};

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_EXECUTIONCONTEXT_HPP
