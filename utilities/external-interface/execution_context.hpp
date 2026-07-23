#ifndef SEQUANT_EXTERNAL_INTERFACE_EXECUTIONCONTEXT_HPP
#define SEQUANT_EXTERNAL_INTERFACE_EXECUTIONCONTEXT_HPP

#include "processing_data.hpp"

#include <SeQuant/core/utility/exception.hpp>

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

  template <std::ranges::range Range>
    requires(
        std::is_convertible_v<std::ranges::range_value_t<Range>, std::string>)
  void set_data(Range &&ids, ProcessingData data) {
    if (std::ranges::empty(ids)) {
      throw Exception("Attempted to add data without specifying any id");
    }

    for (auto &&current : ids) {
      if (data_ids_.find(current) != data_ids_.end()) {
        throw Exception("Duplicate data ID '" + std::string(current) + "'");
      }
      if (!is_valid_id(current, false)) {
        throw Exception("Invalid id '" + std::string(current) + "'");
      }
    }

    for (auto &&current : ids) {
      data_ids_.emplace(std::string(std::move(current)), data_.size());
    }

    data_.emplace_back(std::move(data));
  }

  void add_data_alias(std::string_view id, std::string alias);

  bool has_data(std::string_view id) const;

  const ProcessingData &get_data(std::string_view id) const;
  ProcessingData &get_data(std::string_view id);

  static bool is_valid_id(std::string_view id, bool allow_selectors);

  static std::vector<std::string> expand_selectors(std::string_view id);

 private:
  std::vector<ProcessingData> data_;
  std::map<std::string, std::size_t, std::less<>> data_ids_;
};

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_EXECUTIONCONTEXT_HPP
