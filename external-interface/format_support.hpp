#ifndef SEQUANT_EXTERNAL_FORMAT_SUPPORT_HPP
#define SEQUANT_EXTERNAL_FORMAT_SUPPORT_HPP

#include <spdlog/fmt/bundled/format.h>
#include <spdlog/fmt/bundled/ranges.h>

#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/result_expr.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <string_view>

// Index
template<> struct fmt::formatter< sequant::Index > : fmt::formatter< std::string_view > {
	template< typename FormatContext >
	auto format(const sequant::Index &idx, FormatContext &ctx) const -> decltype(ctx.out()) {
		if (idx.has_proto_indices()) {
			return fmt::format_to(ctx.out(), "{}", sequant::toUtf8(idx.full_label()));
		}

		return fmt::format_to(ctx.out(), "{}", sequant::toUtf8(idx.label()));
	}
};

// Tensor
template<> struct fmt::formatter< sequant::Tensor > : fmt::formatter< std::string_view > {
	template< typename FormatContext >
	auto format(const sequant::Tensor &tensor, FormatContext &ctx) const -> decltype(ctx.out()) {
		return fmt::format_to(ctx.out(), "{}[{};{};{}]", sequant::toUtf8(tensor.label()), fmt::join(tensor.bra(), ", "),
						 fmt::join(tensor.ket(), ", "), fmt::join(tensor.aux(), ", "));
	}
};

// Variable
template<> struct fmt::formatter< sequant::Variable > : fmt::formatter< std::string_view > {
	template< typename FormatContext >
	auto format(const sequant::Variable &variable, FormatContext &ctx) const -> decltype(ctx.out()) {
		return fmt::format_to(ctx.out(), "{}", sequant::toUtf8(variable.label()));
	}
};

// rational
template<> struct fmt::formatter< sequant::rational > : fmt::formatter< std::string_view > {
	template< typename FormatContext >
	auto format(const sequant::rational &number, FormatContext &ctx) const -> decltype(ctx.out()) {
		std::stringstream sstream;
		sstream << number;
		return fmt::format_to(ctx.out(), "{}", sstream.str());
	}
};

// Complex
template<> struct fmt::formatter< sequant::Complex< sequant::rational > > : fmt::formatter< std::string_view > {
	template< typename FormatContext >
	auto format(const sequant::Complex< sequant::rational > &number, FormatContext &ctx) const -> decltype(ctx.out()) {
		if (number.imag().is_zero()) {
			return fmt::format_to(ctx.out(), "{}", number.real());
		} else if (number.real().is_zero()) {
			return fmt::format_to(ctx.out(), "{}i", number.imag());
		} else if (number.imag() < 0) {
			// We need this intermediate step to prevent any template-magic types from appearing
			// due to the multiplication.
			decltype(number.imag()) imag = -number.imag();
			return fmt::format_to(ctx.out(), "({} - {}i)", number.real(), imag);
		} else {
			return fmt::format_to(ctx.out(), "({} + {}i)", number.real(), number.imag());
		}
	}
};

// Constant
template<> struct fmt::formatter< sequant::Constant > : fmt::formatter< std::string_view > {
	template< typename FormatContext >
	auto format(const sequant::Constant &constant, FormatContext &ctx) const -> decltype(ctx.out()) {
		return fmt::format_to(ctx.out(), "{}", constant.value());
	}
};

// Product
template<> struct fmt::formatter< sequant::Product > : fmt::formatter< std::string_view > {
	template< typename FormatContext >
	auto format(const sequant::Product &product, FormatContext &ctx) const -> decltype(ctx.out()) {
		if (product.scalar().is_identity()) {
			return fmt::format_to(ctx.out(), "{}", fmt::join(product.factors(), " "));
		} else {
			return fmt::format_to(ctx.out(), "{} {}", product.scalar(), fmt::join(product.factors(), " "));
		}
	}
};

// Sum
template<> struct fmt::formatter< sequant::Sum > : fmt::formatter< std::string_view > {
	template< typename FormatContext >
	auto format(const sequant::Sum &sum, FormatContext &ctx) const -> decltype(ctx.out()) {
		return fmt::format_to(ctx.out(), "{}", fmt::join(sum.summands(), "\n+ "));
	}
};

// Expr
template<> struct fmt::formatter< sequant::Expr > : fmt::formatter< std::string_view > {
	template< typename FormatContext >
	auto format(const sequant::Expr &expr, FormatContext &ctx) const -> decltype(ctx.out()) {
		if (expr.is< sequant::Tensor >()) {
			return fmt::format_to(ctx.out(), "{}", expr.as< sequant::Tensor >());
		} else if (expr.is< sequant::Constant >()) {
			return fmt::format_to(ctx.out(), "{}", expr.as< sequant::Constant >());
		} else if (expr.is< sequant::Variable >()) {
			return fmt::format_to(ctx.out(), "{}", expr.as< sequant::Variable >());
		} else if (expr.is< sequant::Sum >()) {
			return fmt::format_to(ctx.out(), "{}", expr.as< sequant::Sum >());
		} else if (expr.is< sequant::Product >()) {
			return fmt::format_to(ctx.out(), "{}", expr.as< sequant::Product >());
		} else {
			return fmt::format_to(ctx.out(), "<Unknown expression type>)");
		}
	}
};

// ExprPtr
template<> struct fmt::formatter< sequant::ExprPtr > : fmt::formatter< std::string_view > {
	template< typename FormatContext >
	auto format(const sequant::ExprPtr &expr, FormatContext &ctx) const -> decltype(ctx.out()) {
		return fmt::format_to(ctx.out(), "{}", *expr);
	}
};

// ResultExpr
template<> struct fmt::formatter< sequant::ResultExpr > : fmt::formatter< std::string_view > {
	template< typename FormatContext >
	auto format(const sequant::ResultExpr &result, FormatContext &ctx) const -> decltype(ctx.out()) {
		std::string label = result.has_label() ? sequant::toUtf8(result.label()) : "?";

		if (result.bra().empty() && result.ket().empty() && result.aux().empty()) {
			return fmt::format_to(ctx.out(), "{} =\n{}", label, result.expression());
		}

		return fmt::format_to(ctx.out(), "{}[{};{};{}] =\n{}", label, fmt::join(result.bra(), ", "), fmt::join(result.ket(), ", "),
						 fmt::join(result.aux(), ", "), result.expression());
	}
};

#endif // SEQUANT_EXTERNAL_FORMAT_SUPPORT_HPP
