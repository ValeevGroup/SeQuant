#include "execution_context.hpp"
#include "processing_data.hpp"

#include <SeQuant/core/utility/conversion.hpp>

#include <ranges>
#include <string>
#include <string_view>
#include <vector>

namespace sequant::util::extint {

void ExecutionContext::set_data(std::string id, ProcessingData data) {
  set_data(std::ranges::single_view{std::move(id)}, std::move(data));
}

void ExecutionContext::add_data_alias(std::string_view id, std::string alias) {
  auto it = data_ids_.find(id);

  if (it == data_ids_.end()) {
    throw Exception("Attempted to alias non-existent ID '" + std::string(id) +
                    "'");
  }

  if (data_ids_.find(alias) != data_ids_.end()) {
    throw Exception("Alias '" + alias + "' already exists as a data ID");
  }

  data_ids_.emplace(std::move(alias), it->second);
}

bool ExecutionContext::has_data(std::string_view id) const {
  return data_ids_.find(id) != data_ids_.end();
}

const ProcessingData &ExecutionContext::get_data(std::string_view id) const {
  if (!is_valid_id(id, false)) {
    throw Exception("Invalid id '" + std::string(id) + "'");
  }

  auto it = data_ids_.find(id);

  if (it == data_ids_.end()) {
    throw Exception("No data available for ID '" + std::string(id) + "'");
  }

  return data_.at(it->second);
}

ProcessingData &ExecutionContext::get_data(std::string_view id) {
  // const-cast is safe as we're calling from a non-const function (implying a
  // non-const this)
  return const_cast<ProcessingData &>(std::as_const(*this).get_data(id));
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

std::vector<std::string> ExecutionContext::expand_selectors(
    std::string_view id) {
  if (!is_valid_id(id, true)) {
    throw Exception("Invalid id '" + std::string(id) + "'");
  }

  std::vector<std::string> expanded;

  if (auto begin = id.find('['); begin != std::string_view::npos) {
    auto end = id.find(']', begin);

    if (end == std::string::npos) {
      throw Exception("Unbalanced brackets in selector '" + std::string(id) +
                      "'");
    }
    if (begin + 1 == end) {
      throw Exception("Empty selector in '" + std::string(id) + "'");
    }

    std::string_view selector = id.substr(begin + 1, end - begin - 1);

    std::vector<std::string_view> parts;
    for (const auto &current : selector | std::ranges::views::split(',')) {
      std::string_view current_view(current.begin(), current.end());
      // trim whitespace
      while (!current_view.empty() && current_view.front() == ' ') {
        current_view.remove_prefix(1);
      }
      while (!current_view.empty() && current_view.back() == ' ') {
        current_view.remove_suffix(1);
      }

      if (current_view.empty()) {
        throw Exception("Empty selector component in '" + std::string(id) +
                        "'");
      }

      if (auto dash_pos = current_view.find('-');
          dash_pos != std::string::npos) {
        // Numeric ranges such as "1-3"
        std::size_t from =
            string_to<std::size_t>(current_view.substr(0, dash_pos));
        std::size_t to =
            string_to<std::size_t>(current_view.substr(dash_pos + 1));

        if (from > to) {
          std::swap(from, to);
        }

        for (std::size_t val : std::ranges::views::iota(from, to + 1)) {
          parts.emplace_back(std::to_string(val));
        }
      } else {
        parts.emplace_back(std::move(current_view));
      }
    }

    std::string_view prefix = id.substr(0, begin);
    std::string_view suffix = id.substr(end + 1);
    for (std::string_view current : parts) {
      expanded.emplace_back(std::string(prefix) + "." + std::string(current) +
                            std::string(suffix));

      if (expanded.back().find('[') != std::string::npos) {
        throw Exception("Multiple selectors in single ID not yet supported");
      }
    }
  } else {
    expanded.emplace_back(std::string(id));
  }

  return expanded;
}

}  // namespace sequant::util::extint
