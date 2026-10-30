#include "json_include.hpp"

#include <SeQuant/core/utility/exception.hpp>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

namespace sequant::util::extint {

namespace {

using nlohmann::json;
namespace fs = std::filesystem;

constexpr const char *include_key = "$include";

json load(const fs::path &path, std::vector<fs::path> &stack);

/// Loads the file named by the directive @p directive, found in a file located
/// in @p dir
json load_directive(const json &directive, const fs::path &dir,
                    std::vector<fs::path> &stack) {
  if (!directive.is_string()) {
    throw Exception(std::string("The value of '") + include_key +
                    "' must be a string");
  }
  return load(dir / directive.get<std::string>(), stack);
}

json resolve(const json &value, const fs::path &dir,
             std::vector<fs::path> &stack) {
  if (value.is_array()) {
    json result = json::array();
    for (const json &element : value) {
      json resolved = resolve(element, dir, stack);
      const bool is_include = element.is_object() && element.size() == 1 &&
                              element.contains(include_key);
      if (is_include && resolved.is_array()) {
        for (json &sub : resolved) {
          result.push_back(std::move(sub));
        }
      } else {
        result.push_back(std::move(resolved));
      }
    }
    return result;
  }

  if (!value.is_object()) {
    return value;
  }

  json result = json::object();
  for (const auto &[key, member] : value.items()) {
    if (key != include_key) {
      result[key] = resolve(member, dir, stack);
    }
  }

  if (!value.contains(include_key)) {
    return result;
  }

  json included = load_directive(value.at(include_key), dir, stack);
  if (result.empty()) {
    return included;
  }

  if (!included.is_object()) {
    throw Exception(std::string("'") + include_key +
                    "' alongside other members requires the included file "
                    "to contain an object");
  }
  for (auto &[key, member] : included.items()) {
    if (!result.emplace(key, std::move(member)).second) {
      throw Exception("Key '" + key + "' is defined both locally and in '" +
                      value.at(include_key).get<std::string>() + "'");
    }
  }
  return result;
}

json load(const fs::path &path, std::vector<fs::path> &stack) {
  std::error_code ec;
  const bool regular = fs::is_regular_file(path, ec);
  if (ec || !regular) {
    throw Exception("JSON file '" + path.string() +
                    "' does not exist or is not a file");
  }

  const fs::path canonical = fs::canonical(path, ec);
  if (ec) {
    throw Exception("Unable to resolve JSON file '" + path.string() +
                    "': " + ec.message());
  }
  if (std::ranges::find(stack, canonical) != stack.end()) {
    throw Exception("Cyclic " + std::string(include_key) + " of '" +
                    path.string() + "'");
  }

  std::ifstream in(path);
  if (!in) {
    throw Exception("Unable to open JSON file '" + path.string() + "'");
  }
  json parsed;
  try {
    parsed = json::parse(in, /*callback*/ nullptr, /*allow_exceptions*/ true,
                         /*skip_comments*/ true);
  } catch (const json::parse_error &e) {
    throw Exception("Failed to parse '" + path.string() + "': " + e.what());
  }

  // Resolve against the (symlink-preserving) location of the file itself
  stack.push_back(canonical);
  json result = resolve(parsed, fs::absolute(path).parent_path(), stack);
  stack.pop_back();
  return result;
}

}  // namespace

json load_json_with_includes(const fs::path &path) {
  std::vector<fs::path> stack;
  return load(path, stack);
}

}  // namespace sequant::util::extint
