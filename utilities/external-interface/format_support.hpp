#ifndef SEQUANT_EXTERNAL_FORMAT_SUPPORT_HPP
#define SEQUANT_EXTERNAL_FORMAT_SUPPORT_HPP

#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <range/v3/core.hpp>
#include <range/v3/view/join.hpp>
#include <range/v3/view/transform.hpp>

#include <format>
#include <ranges>
#include <string>
#include <string_view>

// Expr
template <>
struct std::formatter<sequant::Expr> : std::formatter<std::string_view> {
  template <typename FormatContext>
  auto format(const sequant::Expr &expr, FormatContext &ctx) const
      -> decltype(ctx.out()) {
    if (expr.is<sequant::Tensor>()) {
      return std::format_to(ctx.out(), "{}", expr.as<sequant::Tensor>());
    } else if (expr.is<sequant::Constant>()) {
      return std::format_to(ctx.out(), "{}", expr.as<sequant::Constant>());
    } else if (expr.is<sequant::Variable>()) {
      return std::format_to(ctx.out(), "{}", expr.as<sequant::Variable>());
    } else if (expr.is<sequant::Sum>()) {
      return std::format_to(ctx.out(), "{}", expr.as<sequant::Sum>());
    } else if (expr.is<sequant::Product>()) {
      return std::format_to(ctx.out(), "{}", expr.as<sequant::Product>());
    } else {
      return std::format_to(ctx.out(), "<Unknown expression type>)");
    }
  }
};

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

// Tensor
template <>
struct std::formatter<sequant::Tensor> : std::formatter<std::string_view> {
  template <typename FormatContext>
  auto format(const sequant::Tensor &tensor, FormatContext &ctx) const
      -> decltype(ctx.out()) {
    auto idx_to_string = [](const sequant::Index &idx) {
      return std::format("{}", idx);
    };

    using namespace std::literals;

    return std::format_to(
        ctx.out(), "{}[{};{};{}]", sequant::toUtf8(tensor.label()),
        tensor.bra() | ::ranges::views::transform(idx_to_string) |
            ::ranges::views::join(", "sv) | ::ranges::to<std::string>(),
        tensor.ket() | ::ranges::views::transform(idx_to_string) |
            ::ranges::views::join(", "sv) | ::ranges::to<std::string>(),
        tensor.aux() | ::ranges::views::transform(idx_to_string) |
            ::ranges::views::join(", "sv) | ::ranges::to<std::string>());
  }
};

// Variable
template <>
struct std::formatter<sequant::Variable> : std::formatter<std::string_view> {
  template <typename FormatContext>
  auto format(const sequant::Variable &variable, FormatContext &ctx) const
      -> decltype(ctx.out()) {
    return std::format_to(ctx.out(), "{}", sequant::toUtf8(variable.label()));
  }
};

// rational
template <>
struct std::formatter<sequant::rational> : std::formatter<std::string_view> {
  template <typename FormatContext>
  auto format(const sequant::rational &number, FormatContext &ctx) const
      -> decltype(ctx.out()) {
    std::stringstream sstream;
    sstream << number;
    return std::format_to(ctx.out(), "{}", sstream.str());
  }
};

// Complex
template <>
struct std::formatter<sequant::Complex<sequant::rational> >
    : std::formatter<std::string_view> {
  template <typename FormatContext>
  auto format(const sequant::Complex<sequant::rational> &number,
              FormatContext &ctx) const -> decltype(ctx.out()) {
    if (number.imag().is_zero()) {
      return std::format_to(ctx.out(), "{}", number.real());
    } else if (number.real().is_zero()) {
      return std::format_to(ctx.out(), "{}i", number.imag());
    } else if (number.imag() < 0) {
      // We need this intermediate step to prevent any template-magic types from
      // appearing due to the multiplication.
      decltype(number.imag()) imag = -number.imag();
      return std::format_to(ctx.out(), "({} - {}i)", number.real(), imag);
    } else {
      return std::format_to(ctx.out(), "({} + {}i)", number.real(),
                            number.imag());
    }
  }
};

// Constant
template <>
struct std::formatter<sequant::Constant> : std::formatter<std::string_view> {
  template <typename FormatContext>
  auto format(const sequant::Constant &constant, FormatContext &ctx) const
      -> decltype(ctx.out()) {
    return std::format_to(ctx.out(), "{}", constant.value());
  }
};

// Product
template <>
struct std::formatter<sequant::Product> : std::formatter<std::string_view> {
  template <typename FormatContext>
  auto format(const sequant::Product &product, FormatContext &ctx) const
      -> decltype(ctx.out()) {
    auto factors = product.factors() |
                   ::ranges::views::transform([](const sequant::ExprPtr &expr) {
                     return std::format("{}", *expr);
                   }) |
                   ::ranges::views::join(' ') | ::ranges::to<std::string>();

    if (product.scalar().is_identity()) {
      return std::format_to(ctx.out(), "{}", factors);
    } else {
      return std::format_to(ctx.out(), "{} {}", product.scalar(), factors);
    }
  }
};

// Sum
template <>
struct std::formatter<sequant::Sum> : std::formatter<std::string_view> {
  template <typename FormatContext>
  auto format(const sequant::Sum &sum, FormatContext &ctx) const
      -> decltype(ctx.out()) {
    using namespace std::literals;

    return std::format_to(
        ctx.out(), "{}",
        sum.summands() |
            ::ranges::views::transform([](const sequant::ExprPtr &expr) {
              return std::format("{}", *expr);
            }) |
            ::ranges::views::join("\n+ "sv) | ::ranges::to<std::string>());
  }
};

// ExprPtr
template <>
struct std::formatter<sequant::ExprPtr> : std::formatter<std::string_view> {
  template <typename FormatContext>
  auto format(const sequant::ExprPtr &expr, FormatContext &ctx) const
      -> decltype(ctx.out()) {
    return std::format_to(ctx.out(), "{}", *expr);
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
      return std::format_to(ctx.out(), "{} =\n{}", label, result.expression());
    }

    auto idx_to_string = [](const sequant::Index &idx) {
      return std::format("{}", idx);
    };

    return std::format_to(
        ctx.out(), "{}[{};{};{}] =\n {}", label,
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
