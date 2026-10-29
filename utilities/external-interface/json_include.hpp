#ifndef SEQUANT_EXTERNAL_INTERFACE_JSON_INCLUDE_HPP
#define SEQUANT_EXTERNAL_INTERFACE_JSON_INCLUDE_HPP

#include <nlohmann/json_fwd.hpp>

#include <filesystem>

namespace sequant::util::extint {

/// Parses the JSON file at @p path (C-style comments are allowed) and resolves
/// all @c $include directives contained in it, recursively.
///
/// An object of the form <tt>{"$include": "file.json"}</tt> is replaced by the
/// content of the referenced file. If the object has further members, the
/// included file must be an object whose members are merged with them (a key
/// present on both sides is an error). Inside an array, an include of an array
/// splices the included elements in place. Relative include paths are resolved
/// against the directory of the file containing the directive.
///
/// @throws Exception on I/O errors, malformed directives, merge conflicts and
///         include cycles.
nlohmann::json load_json_with_includes(const std::filesystem::path &path);

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_JSON_INCLUDE_HPP
