#ifndef SEQUANT_EXTERNAL_FORMAT_SUPPORT_HPP
#define SEQUANT_EXTERNAL_FORMAT_SUPPORT_HPP

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/format.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <range/v3/core.hpp>
#include <range/v3/view/join.hpp>
#include <range/v3/view/transform.hpp>

#include <format>
#include <ranges>
#include <string>
#include <string_view>

// Index
template <>
struct std::formatter<sequant::Index> : std::formatter<std::string_view> {
  template <typename FormatContext>
  auto format(const sequant::Index &idx, FormatContext &ctx) const
      -> decltype(ctx.out()) {
    if (idx.has_proto_indices()) {
      return std::format_to(ctx.out(), "{}", sequant::toUtf8(idx.full_label()));
    }

    return std::format_to(ctx.out(), "{}", sequant::toUtf8(idx.label()));
  }
};

// ResultExpr
template <>
struct std::formatter<sequant::ResultExpr> : std::formatter<std::string_view> {
  template <typename FormatContext>
  auto format(const sequant::ResultExpr &result, FormatContext &ctx) const
      -> decltype(ctx.out()) {
    std::string label =
        result.has_label() ? sequant::toUtf8(result.label()) : "?";

    if (result.bra().empty() && result.ket().empty() && result.aux().empty()) {
      return std::format_to(ctx.out(), "{} =\n{:s}", label,
                            result.expression());
    }

    auto idx_to_string = [](const sequant::Index &idx) {
      return std::format("{}", idx);
    };

    return std::format_to(
        ctx.out(), "{}[{};{};{}] =\n {:s}", label,
        result.bra() | ::ranges::views::transform(idx_to_string) |
            ::ranges::views::join(", "sv) | ::ranges::to<std::string>(),
        result.ket() | ::ranges::views::transform(idx_to_string) |
            ::ranges::views::join(", "sv) | ::ranges::to<std::string>(),
        result.aux() | ::ranges::views::transform(idx_to_string) |
            ::ranges::views::join(", "sv) | ::ranges::to<std::string>(),
        result.expression());
  }
};

#endif  // SEQUANT_EXTERNAL_FORMAT_SUPPORT_HPP
