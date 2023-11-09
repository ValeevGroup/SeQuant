//
// Created by Robert Adam on 2023-10-10
//

#ifndef SEQUANT_CORE_UTILITY_STRING_HPP
#define SEQUANT_CORE_UTILITY_STRING_HPP

#include <string>
#include <string_view>

namespace sequant {

/// Converts the given wide-string to a UTF-8 encoded narrow string
std::string toUtf8(std::wstring_view str);

/// Converts the given UTF-8 encoded narrow-string to a UTF-16 encoded
/// wide-string
std::wstring toUtf16(std::string_view str);

}  // namespace sequant

#endif  // SEQUANT_CORE_UTILITY_STRING_HPP
