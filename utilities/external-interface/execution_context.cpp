#include "execution_context.hpp"
#include "processing_data.hpp"

#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/utility/conversion.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <algorithm>
#include <limits>
#include <ranges>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace sequant::util::extint {

void ExecutionContext::set_data(std::string id, ProcessingData data) {
  set_data(std::ranges::single_view{std::move(id)},
           std::ranges::single_view{std::move(data)});
}

void ExecutionContext::add_data_alias(std::string_view id, std::string alias) {
  add_data_alias(std::ranges::single_view(std::move(id)), std::move(alias));
}

bool ExecutionContext::has_data(std::string_view id) const {
  return data_ids_.find(id) != data_ids_.end();
}

std::size_t ExecutionContext::dataset_size(std::string_view id) const {
  auto it = data_ids_.find(id);

  return it == data_ids_.end() ? 0 : it->second.size();
}

template <std::ranges::random_access_range DataVec,
          std::ranges::range DataIDMap>
std::vector<
    std::reference_wrapper<meta::mimic_constness_t<DataVec, ProcessingData>>>
get_data_impl(DataVec &&data, DataIDMap &&data_ids, std::string_view id) {
  if (!ExecutionContext::is_valid_id(id, true)) {
    throw Exception("Invalid id '" + std::string(id) + "'");
  }

  using RefWrapT =
      std::reference_wrapper<meta::mimic_constness_t<DataVec, ProcessingData>>;

  std::vector<RefWrapT> selected_data;

  for (const std::string_view current : ExecutionContext::expand_id(id)) {
    auto it = data_ids.find(current);

    if (it == data_ids.end()) {
      throw Exception("No data available for ID '" + std::string(current) +
                      "'");
    }

    for (std::size_t idx : it->second) {
      selected_data.emplace_back(data.at(idx));
    }
  }

  return selected_data;
}

std::vector<std::reference_wrapper<const ProcessingData>>
ExecutionContext::get_data(std::string_view id) const {
  return get_data_impl(data_, data_ids_, id);
}

std::vector<std::reference_wrapper<ProcessingData>> ExecutionContext::get_data(
    std::string_view id) {
  return get_data_impl(data_, data_ids_, id);
}

bool ExecutionContext::is_valid_id(std::string_view id, bool allow_selectors) {
  // TODO: implement
  (void)id;
  (void)allow_selectors;
  // Rules:
  // - No dash outside of selectors
  // - No newlines or tabs
  // - square brackets may only represent selectors
  // - No asterix outside of selector
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
