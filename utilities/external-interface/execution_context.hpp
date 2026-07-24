#ifndef SEQUANT_EXTERNAL_INTERFACE_EXECUTIONCONTEXT_HPP
#define SEQUANT_EXTERNAL_INTERFACE_EXECUTIONCONTEXT_HPP

#include "processing_data.hpp"

#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <functional>
#include <map>
#include <ranges>
#include <string>
#include <type_traits>
#include <vector>

namespace sequant::util::extint {

class ExecutionContext {
 public:
  ExecutionContext() = default;

  void set_data(std::string id, ProcessingData data);

  template <std::ranges::range IDs, std::ranges::range Data>
    requires(
        std::is_convertible_v<std::ranges::range_value_t<IDs>, std::string> &&
        std::same_as<std::ranges::range_value_t<Data>, ProcessingData>)
  void set_data(IDs &&ids, Data &&data) {
    if (std::ranges::empty(ids)) {
      throw Exception("Attempted to add data without specifying any id");
    }
    if (std::ranges::empty(data)) {
      throw Exception("Attempted to register empty dataset");
    }

    for (auto &&current : ids) {
      if (data_ids_.find(current) != data_ids_.end()) {
        throw Exception("Duplicate data ID '" + std::string(current) + "'");
      }
      if (!is_valid_id(current, false)) {
        throw Exception("Invalid id '" + std::string(current) + "'");
      }
    }

    std::vector<std::size_t> data_indices(std::ranges::size(data));
    std::iota(data_indices.begin(), data_indices.end(), data_.size());

    for (auto &&current : ids) {
      data_ids_.emplace(std::string(std::move(current)), data_indices);
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
    if (data_ids_.find(alias) != data_ids_.end()) {
      throw Exception("Alias '" + alias + "' already exists as a data ID");
    }

    std::vector<std::size_t> indices;

    for (std::string_view current_id : ids) {
      for (const std::string &expanded : expand_id(current_id)) {
        auto it = data_ids_.find(expanded);

        if (it == data_ids_.end()) {
          throw Exception("Attempted to alias non-existent ID '" + expanded +
                          "'");
        }

        indices.insert(indices.end(), it->second.begin(), it->second.end());
      }
    }

    SEQUANT_ASSERT(!indices.empty());

    data_ids_.emplace(std::move(alias), std::move(indices));
  }

  bool has_data(std::string_view id) const;

  std::size_t dataset_size(std::string_view id) const;

  std::vector<std::reference_wrapper<const ProcessingData>> get_data(
      std::string_view id) const;
  std::vector<std::reference_wrapper<ProcessingData>> get_data(
      std::string_view id);

  static bool is_valid_id(std::string_view id, bool allow_selectors);

  static std::vector<std::string> expand_id(std::string_view id);

 private:
  std::vector<ProcessingData> data_;
  std::map<std::string, std::vector<std::size_t>, std::less<>> data_ids_;
};

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_EXECUTIONCONTEXT_HPP
