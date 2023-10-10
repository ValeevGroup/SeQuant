//
// Created by Robert Adam on 2023-09-20
//

#include "SeQuant/core/parse_expr.hpp"
#include "ast.hpp"
#include "ast_conversions.hpp"
#include "semantic_actions.hpp"

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/tensor.hpp>

#define BOOST_SPIRIT_X3_UNICODE
#include <boost/core/demangle.hpp>
#include <boost/spirit/home/x3.hpp>
#include <boost/spirit/home/x3/support/utility/error_reporting.hpp>
#include <boost/variant.hpp>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <functional>
#include <iostream>
#include <memory>
#include <type_traits>

namespace sequant {

ParseError::ParseError(std::size_t offset, std::size_t length,
                       std::string message)
    : std::runtime_error(std::move(message)), offset(offset), length(length) {}

namespace x3 = boost::spirit::x3;

namespace parse {

struct NumberRule;
struct VariableRule;
struct TensorRule;
struct ProductRule;
struct SumRule;
struct ExprRule;
struct IndexLabelRule;
struct IndexRule;
struct IndexGroupRule;

// Types
x3::rule<NumberRule, ast::Number> number{"Number"};
x3::rule<VariableRule, ast::Variable> variable{"Variable"};
x3::rule<TensorRule, ast::Tensor> tensor{"Tensor"};

// Expression structure
x3::rule<ProductRule, ast::Product> product{"Product"};
x3::rule<SumRule, ast::Sum> sum{"Sum"};
x3::rule<ExprRule, ast::Sum> expr{"Expression"};

// Auxiliaries
x3::rule<struct NameRule, std::wstring> name{"Name"};
x3::rule<IndexLabelRule, ast::IndexLabel> index_label{"IndexLabel"};
x3::rule<IndexRule, ast::Index> index{"Index"};
x3::rule<IndexGroupRule, ast::IndexGroups> index_groups{"IndexGroups"};

// clang-format off
auto word_components = x3::unicode::alnum | '_';
// A name begins with a letter, then can container letters, digits and
// underscores, but can not end with an underscore (to not confuse the parser
// with tensors á la t_{…}^{…}.
auto name_def         = x3::lexeme[
							x3::unicode::alpha >> -( *(word_components >> &word_components) >> x3::unicode::alnum )
						];

auto number_def       = x3::double_ >> -('/' >> x3::double_);

auto variable_def     = name;

auto index_label_def  = x3::lexeme[
	   						+x3::unicode::alpha >> -x3::lit('_') >> x3::uint_
						];

auto index_def        = x3::lexeme[
							index_label >> -('<' >> index_label % ',' >> ">")
						];

auto index_groups_def =   L"_{" > -(index % ',') > L"}^{" > -(index % ',') > L"}" >> x3::attr(false)
                        | L"^{" > -(index % ',') > L"}_{" > -(index % ',') > L"}" >> x3::attr(true)
                        |  '{'  > -(index % ',') > ';'    > -(index % ',') >  '}' >> x3::attr(false);

auto tensor_def       = x3::lexeme[
							name >> x3::skip[index_groups] >> -(':' >> x3::upper)
						];

auto nullary          = number | tensor | variable;

auto grouped          = '(' > sum > ')' | nullary;

auto product_def      = grouped % -x3::lit('*');

auto first_addend     = (('-' >> x3::attr(-1) | -x3::lit('+') >> x3::attr(1)) >> product)[actions::process_addend{}];

auto addend           = (('+' >> x3::attr(1) | '-' >> x3::attr(-1)) > product)[actions::process_addend{}];

auto sum_def          = first_addend >> *addend;

auto expr_def         = -sum > x3::eoi;
// clang-format on

BOOST_SPIRIT_DEFINE(name, number, variable, index_label, index, index_groups,
                    tensor, product, sum, expr);

struct position_cache_tag;
struct error_handler_tag;

namespace helpers {

struct annotate_position {
  template <typename T, typename Iterator, typename Context>
  void on_success(const Iterator &first, const Iterator &last, T &ast,
                  const Context &ctx) {
    auto &position_cache =
        boost::spirit::x3::get<position_cache_tag>(ctx).get();
    position_cache.annotate(ast, first, last);
  }
};

struct error_handler {
  template <typename Iterator, typename Exception, typename Context>
  x3::error_handler_result on_error(const Iterator &first, const Iterator &last,
                                    const Exception &e, const Context &ctx) {
    auto &error_handler = x3::get<error_handler_tag>(ctx).get();
    error_handler(e.where(), boost::core::demangle(e.which().data()));
    return x3::error_handler_result::fail;
  }
};

}  // namespace helpers

struct NumberRule : helpers::annotate_position, helpers::error_handler {};
struct VariableRule : helpers::annotate_position, helpers::error_handler {};
struct TensorRule : helpers::annotate_position, helpers::error_handler {};
struct ProductRule : helpers::annotate_position, helpers::error_handler {};
struct SumRule : helpers::annotate_position, helpers::error_handler {};
struct ExprRule : helpers::annotate_position, helpers::error_handler {};
struct IndexLabelRule : helpers::annotate_position, helpers::error_handler {};
struct IndexRule : helpers::annotate_position, helpers::error_handler {};
struct IndexGroupRule : helpers::annotate_position, helpers::error_handler {};

}  // namespace parse

template <typename Iterator>
struct ErrorHandler {
  Iterator begin;

  ErrorHandler(Iterator begin) : begin(std::move(begin)) {}

  void operator()(Iterator where, std::string expected) const {
    std::size_t offset = std::distance(begin, where);
    throw ParseError(offset, 1,
                     std::string("Parse failure at offset ") +
                         std::to_string(offset) + ": expected " + expected);
  }
};

ExprPtr parse_expr(std::wstring_view input, Symmetry default_symmetry) {
  using iterator_type = decltype(input)::iterator;
  x3::position_cache<std::vector<iterator_type>> positions(input.begin(),
                                                           input.end());

  ErrorHandler<iterator_type> error_handler(input.begin());

  parse::ast::Sum ast;

  const auto parser = x3::with<parse::error_handler_tag>(
      std::ref(error_handler))[x3::with<parse::position_cache_tag>(
      std::ref(positions))[parse::expr]];

  auto start = input.begin();
  try {
    bool success =
        x3::phrase_parse(start, input.end(), parser, x3::unicode::space, ast);

    if (!success) {
      // Normally, this shouldn't happen as any error should itself throw a
      // ParseError already
      throw ParseError(0, input.size(),
                       "Parsing was unsuccessful for an unknown reason");
    }
    if (start != input.end()) {
      // This should also not happen as the parser requires matching EOI
      throw ParseError(std::distance(input.begin(), start),
                       std::distance(start, input.end()),
                       "Couldn't parse the entire input");
    }
  } catch (const boost::spirit::x3::expectation_failure<iterator_type> &e) {
    std::wcout << "Caught expectation_failure\nwhere: " << e.where()
               << "\nwhat: " << e.what() << "\nwhich: " << e.which().data()
               << std::endl;
    throw;
  }

  return parse::transform::ast_to_expr(ast, positions, input.begin(),
                                       default_symmetry);
}

}  // namespace sequant
