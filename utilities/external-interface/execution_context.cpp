#include "execution_context.hpp"
#include "processing_data.hpp"

#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/utility/conversion.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <algorithm>
#include <cctype>
#include <limits>
#include <ranges>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace sequant::util::extint {

void ExecutionContext::set_data(std::string_view prefix, std::size_t counter,
                                ProcessingData data) {
  set_data(prefix, std::ranges::single_view{counter},
           std::ranges::single_view{std::move(data)});
}

void ExecutionContext::add_data_alias(std::string_view id, std::string alias) {
  add_data_alias(std::ranges::single_view(std::move(id)), std::move(alias));
}

bool ExecutionContext::has_data(std::string_view id) const {
  return id_to_data_indices_.find(id) != id_to_data_indices_.end();
}

std::size_t ExecutionContext::dataset_size(std::string_view id) const {
  if (!is_valid_id(id, false)) {
    throw Exception("Invalid ID '" + std::string(id) + "'");
  }

  auto it = id_to_data_indices_.find(id);

  return it == id_to_data_indices_.end() ? 0 : it->second.size();
}

template <std::ranges::random_access_range DataVec,
          std::ranges::range ID2DataIdxMap, std::ranges::range DataIdx2IDMap>
std::vector<
    ExecutionContext::Data<meta::mimic_constness_t<DataVec, ProcessingData>>>
get_data_impl(DataVec &&data, ID2DataIdxMap &&id2dataidx,
              DataIdx2IDMap &&dataidx2id, std::string_view id) {
  if (!ExecutionContext::is_valid_id(id, true)) {
    throw Exception("Invalid id '" + std::string(id) + "'");
  }

  using RefT = meta::mimic_constness_t<DataVec, ProcessingData>;

  std::vector<ExecutionContext::Data<RefT>> selected_data;

  for (const std::string_view current : ExecutionContext::expand_id(id)) {
    auto it = id2dataidx.find(current);

    if (it == id2dataidx.end()) {
      throw Exception("No data available for ID '" + std::string(current) +
                      "'");
    }

    for (std::size_t idx : it->second) {
      ExecutionContext::Data<RefT> ret_data{.data = data.at(idx)};

      for (std::string_view alias : dataidx2id.at(idx)) {
#if defined(__cpp_lib_associative_heterogeneous_insertion) && \
    __cpp_lib_associative_heterogeneous_insertion >= 202311L
        if (id2dataidx.at(alias).size() > 1) {
#else
        // No heterogenous lookup overload for at() until C++26
        if (id2dataidx.at(std::string(alias)).size() > 1) {
#endif
          ret_data.associated_group_ids.emplace_back(std::move(alias));
        } else {
          ret_data.associated_ids.emplace_back(std::move(alias));
        }
      }

      // Ensure that we order the IDs such that the generic index-IDs come after
      // potentially manually assigned (that will have more meaningful names)
      auto not_ends_with_num = [](std::string_view id) -> bool {
        SEQUANT_ASSERT(!id.empty());
        if (id.empty()) {
          return false;
        }

        auto it = id.rfind('.');
        if (it == std::string_view::npos) {
          // The generic number-IDs always contain a period
          return true;
        }

        std::string_view suffix = id.substr(it + 1);
        return suffix.find_first_not_of("0123456789") != std::string_view::npos;
      };

      std::ranges::stable_partition(ret_data.associated_ids, not_ends_with_num);
      std::ranges::stable_partition(ret_data.associated_group_ids,
                                    not_ends_with_num);

      selected_data.emplace_back(std::move(ret_data));
    }
  }

  return selected_data;
}

std::vector<ExecutionContext::Data<const ProcessingData>>
ExecutionContext::get_data(std::string_view id) const {
  return get_data_impl(data_, id_to_data_indices_, data_idx_to_ids_, id);
}

std::vector<ExecutionContext::Data<ProcessingData>> ExecutionContext::get_data(
    std::string_view id) {
  return get_data_impl(data_, id_to_data_indices_, data_idx_to_ids_, id);
}

bool ExecutionContext::is_valid_id(std::string_view id, bool allow_selectors) {
  auto validate_non_selector = [](std::string_view part) -> bool {
    for (char c : part) {
      if (!std::isalnum(c) & c != '_' && c != '.') {
        return false;
      }
    }

    return true;
  };

  auto validate_range_component = [](std::string_view comp) -> bool {
    for (char c : comp) {
      if (!std::isdigit(c)) {
        return false;
      }
    }

    return true;
  };

  auto validate_selector = [&validate_non_selector, &validate_range_component](
                               std::string_view selector) -> bool {
    for (auto &&comp : selector | std::ranges::views::split(',')) {
      std::string_view part(comp.begin(), comp.end());

      if (validate_non_selector(part)) {
        continue;
      }

      auto dash_pos = part.find('-');
      if (dash_pos == std::string_view::npos ||
          part.find('-', dash_pos + 1) != std::string_view::npos) {
        // Either no dash or more than one dash in single selector component
        return false;
      }

      std::string_view prefix = part.substr(0, dash_pos);
      std::string_view suffix = part.substr(dash_pos + 1);

      if (prefix.empty() || suffix.empty()) {
        return false;
      }

      if (!validate_range_component(prefix) ||
          !validate_range_component(suffix)) {
        return false;
      }
    }

    return true;
  };

  if (id.empty()) {
    return false;
  }

  std::string_view::size_type bracket_begin = 0;
  std::string_view::size_type prev_pos = 0;
  do {
    bracket_begin = id.find('[', prev_pos);

    if (!validate_non_selector(id.substr(prev_pos, bracket_begin - prev_pos))) {
      return false;
    }

    if (bracket_begin == std::string_view::npos) {
      continue;
    }

    if (!allow_selectors) {
      return false;
    }

    std::string_view::size_type bracket_end = id.find(']', bracket_begin);

    if (bracket_end == std::string_view::npos) {
      return false;
    }
    if (bracket_begin + 1 == bracket_end) {
      return false;
    }

    if (!validate_selector(
            id.substr(bracket_begin + 1, bracket_end - bracket_begin - 1))) {
      return false;
    }

    prev_pos = bracket_end + 1;
  } while (bracket_begin != std::string_view::npos && prev_pos < id.size());

  return true;
}

std::vector<std::string> expand_selector(std::string_view selector) {
  std::vector<std::string> processed;

  for (auto &&comp : selector | std::ranges::views::split(',')) {
    std::string_view current(comp.begin(), comp.end());
    // trim whitespace
    while (!current.empty() && current.front() == ' ') {
      current.remove_prefix(1);
    }
    while (!current.empty() && current.back() == ' ') {
      current.remove_suffix(1);
    }

    if (current.empty()) {
      throw Exception("Empty selector component in '[" + std::string(selector) +
                      "']");
    }

    if (auto dash_pos = current.find('-'); dash_pos != std::string::npos) {
      // Numeric ranges such as "1-3"
      std::size_t from = string_to<std::size_t>(current.substr(0, dash_pos));
      std::size_t to = string_to<std::size_t>(current.substr(dash_pos + 1));

      if (from > to) {
        std::swap(from, to);
      }

      for (std::size_t val : std::ranges::views::iota(from, to + 1)) {
        processed.emplace_back(std::to_string(val));
      }
    } else {
      processed.emplace_back(std::move(current));
    }
  }

  SEQUANT_ASSERT(std::none_of(processed.begin(), processed.end(),
                              [](const auto &p) { return p.empty(); }));

  return processed;
}

std::vector<std::vector<std::string>> create_id_partitions(
    std::string_view id) {
  std::vector<std::vector<std::string>> partitions;

  std::size_t begin = -1;
  std::size_t prev_pos = 0;
  do {
    begin = id.find('[', begin + 1);

    std::vector<std::string> part = {std::string(id.substr(prev_pos, begin))};
    SEQUANT_ASSERT(!part.back().empty());
    partitions.emplace_back(std::move(part));

    if (begin == std::string_view::npos) {
      continue;
    }

    auto end = id.find(']', begin);

    if (end == std::string_view::npos) {
      throw Exception("Unbalanced brackets in selector '" + std::string(id) +
                      "'");
    }
    if (begin + 1 == end) {
      throw Exception("Empty selector in '" + std::string(id) + "'");
    }

    std::string_view selector = id.substr(begin + 1, end - begin - 1);
    partitions.emplace_back(expand_selector(selector));

    prev_pos = end + 1;
  } while (begin != std::string_view::npos && prev_pos < id.size());

  SEQUANT_ASSERT(std::none_of(partitions.begin(), partitions.end(),
                              [](const auto &p) { return p.empty(); }));

  return partitions;
}

std::vector<std::string> ExecutionContext::expand_id(std::string_view id) {
  if (!is_valid_id(id, true)) {
    throw Exception("Invalid id '" + std::string(id) + "'");
  }

  std::vector<std::vector<std::string>> partitions = create_id_partitions(id);
  std::vector<std::size_t> indices(partitions.size(), 0);

  auto has_more = [&partitions, &indices]() {
    for (std::size_t i = 0; i < indices.size(); ++i) {
      SEQUANT_ASSERT(partitions[i].size() > 0);
      if (indices[i] >= partitions[i].size()) {
        return false;
      }
    }

    return true;
  };

  auto increment = [&partitions, &indices]() mutable {
    for (std::size_t i = 0; i < indices.size(); ++i) {
      if (indices[i] < partitions[i].size() - 1) {
        ++indices[i];
        return;
      }
      indices[i] = 0;
    }

    // indicate end has been reached
    std::fill(indices.begin(), indices.end(),
              std::numeric_limits<std::size_t>::max());
  };

  SEQUANT_ASSERT(has_more());

  std::vector<std::string> expanded;

  while (has_more()) {
    std::stringstream stream;
    for (std::size_t i = 0; i < indices.size(); ++i) {
      stream << partitions[i].at(indices[i]);

      if (i + 1 < indices.size()) {
        stream << ".";
      }
    }

    expanded.emplace_back(stream.str());

    increment();
  }

  SEQUANT_ASSERT(!expanded.empty());

  return expanded;
}

}  // namespace sequant::util::extint
