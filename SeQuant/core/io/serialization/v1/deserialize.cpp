//
// Created by Robert Adam on 2023-09-20
//

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/io/serialization/v1/ast.hpp>
#include <SeQuant/core/io/serialization/v1/ast_conversions.hpp>
#include <SeQuant/core/io/serialization/v1/semantic_actions.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/string.hpp>

#define BOOST_SPIRIT_X3_UNICODE
#include <boost/core/demangle.hpp>
#include <boost/spirit/home/x3.hpp>
#include <boost/spirit/home/x3/support/utility/error_reporting.hpp>
#include <boost/variant.hpp>

#include <utf8/unchecked.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <iterator>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sequant::io::serialization::v1 {

namespace x3 = boost::spirit::x3;

/// The parser operates on code points obtained by decoding the (UTF-8 encoded)
/// input on the fly; the underlying input position is accessible via base()
using iterator_type =
    utf8::unchecked::iterator<std::string_view::const_iterator>;

namespace parse {

struct NumberRule;
struct VariableRule;
struct TensorRule;
struct PowerRule;
struct ProductRule;
struct SumRule;
struct ExprRule;
struct ResultExprRule;
struct IndexLabelRule;
struct ProtoLabelRule;
struct IndexRule;
struct IndexGroupRule;
struct SymmetrySpecRule;

// Types
x3::rule<NumberRule, ast::Number> number{"Number"};
x3::rule<VariableRule, ast::Variable> variable{"Variable"};
x3::rule<TensorRule, ast::Tensor> tensor{"Tensor"};
x3::rule<PowerRule, ast::Power> power{"Power"};

// Expression structure
x3::rule<ProductRule, ast::Product> product{"Product"};
x3::rule<SumRule, ast::Sum> sum{"Sum"};
x3::rule<ExprRule, ast::Sum> expr{"Expression"};
x3::rule<ResultExprRule, ast::ResultExpr> resultExpr{"ResultExpr"};

// Auxiliaries
x3::rule<struct NameRule, std::string> name{"Name"};
x3::rule<struct BaseKeyRule, std::string> base_key{"BaseKey"};
x3::rule<IndexLabelRule, ast::IndexLabel> index_label{"IndexLabel"};
x3::rule<ProtoLabelRule, ast::ProtoLabel> proto_label{"ProtoLabel"};
x3::rule<IndexRule, ast::Index> index{"Index"};
x3::rule<IndexGroupRule, ast::IndexGroups> index_groups{"IndexGroups"};
x3::rule<SymmetrySpecRule, ast::SymmetrySpec> symmetry_spec{"SymmetrySpec"};

auto to_char_type = [](auto c) {
  return static_cast<x3::unicode::char_type::char_type>(c);
};

/// semantic action that sets the rule's attribute to the (UTF-8 encoded) input
/// that has been matched by the subject parser (which has to be wrapped in
/// x3::raw)
auto set_matched_input = [](auto &ctx) {
  const auto &range = x3::_attr(ctx);
  x3::_val(ctx) = std::string(range.begin().base(), range.end().base());
};

// clang-format off
auto word_components = x3::unicode::alnum
                       | x3::char_('_') | x3::unicode::char_(U'⁔')
                       | x3::unicode::char_(to_char_type(0x0303)) // combining tilde
                       | x3::unicode::char_(to_char_type(0x0302)) // combining circumflex/hat
                       // Superscript and Subscript block
                       | (x3::unicode::char_(to_char_type(0x2070), to_char_type(0x209F)) - x3::unicode::unassigned)
                       // These are defined in the Latin-1 Supplement block and thus need to be listed explicitly
                       | x3::unicode::char_(U'¹') | x3::unicode::char_(U'²') | x3::unicode::char_(U'³')
                       // Arrow block
                       | (x3::unicode::char_(to_char_type(0x2190), to_char_type(0x21FF)) - x3::unicode::unassigned);
// A name begins with a letter, then can container letters, digits and
// underscores, but can not end with an underscore (to not confuse the parser
// with tensors á la t_{…}^{…}.
auto name_def         = x3::raw[x3::lexeme[
                            x3::unicode::alpha >> -( *(word_components >> &word_components) >> (word_components - '_') )
                        ]][set_matched_input];

auto number_def       = x3::double_ >> -('/' >> x3::double_);

auto variable_def     = x3::lexeme[name >> -(x3::lit('^') >> '*' >> x3::attr(true))];

// letters, ⁺/⁻, combining diacritics (e.g. μ̃, f̌), arrows (e.g. ↑/↓) and
// primes (e.g. α')
auto base_key_chars   = +(  x3::unicode::alpha | x3::unicode::char_(U'⁺') | x3::unicode::char_(U'⁻')
                          | x3::unicode::char_(to_char_type(0x0300), to_char_type(0x036F))
                          | (x3::unicode::char_(to_char_type(0x2190), to_char_type(0x21FF)) - x3::unicode::unassigned)
                          | x3::unicode::char_('\'')
                         );

auto base_key_def     = x3::raw[base_key_chars][set_matched_input];

auto index_label_def  = x3::lexeme[
                               base_key >> -x3::lit('_') >> x3::uint_
                        ];

// the basis instance suffix of a domain ("<protos;N>") or a domainless proto
// ("i_1<;3>"); overflow fails the parse rather than silently wrapping
auto inst              = x3::int_parser<IndexBasis::instance_type>();

auto proto_label_def   = index_label >> -('<' >> x3::lit(';') >> inst >> '>');

// a domain is <protos>, <protos;N> or <;N>; the empty domain <> is rejected
// in to_index
auto index_def        = x3::lexeme[
                            index_label >> -x3::skip['<' >> -(proto_label % ',') >> -(';' >> inst) >> '>']
                        ];

const std::vector<ast::Index> noIndices;
// n.b. multi-character string literals have to be (UTF-32) Unicode literals as
// they are matched against the decoded code points of the input
auto index_groups_def =   x3::unicode::lit(U"_{") > -(index % ',') > x3::unicode::lit(U"}^{") > -(index % ',')  > '}' >> x3::attr(noIndices) >> x3::attr(false)
                        | x3::unicode::lit(U"^{") > -(index % ',') > x3::unicode::lit(U"}_{") > -(index % ',')  > '}' >> x3::attr(noIndices) >> x3::attr(true)
                        |  '{'  > -(index % ',') > -( ';' > -(index % ',')) > -(';' > -(index % ','))     >  '}'  >> x3::attr(false);

auto symmetry_spec_def= x3::lexeme[
                         ':' >> x3::upper >> -('-' >> x3::upper) >> -('-' >> x3::upper)
                        ];

auto tensor_def       = x3::lexeme[
                            name >> x3::skip[index_groups] >> -(symmetry_spec)
                        ];

// TODO(power): per comments on PR #513, promote `^` to a binary operator (with higher precedence than *) and then reject unsupported cases while traversing the AST.
// Two forms are accepted: b^(e) and (b^(e))^*
auto power_core       = (('(' >> (number | variable) >> ')') | number | variable)
                        >> x3::lit('^') >> '(' >> number >> ')';
auto power_def        = (power_core >> x3::attr(false))
                        | ('(' >> power_core >> ')' >> x3::unicode::lit(U"^*") >> x3::attr(true));

auto nullary          = number | tensor | variable;

auto grouped          = power | '(' > sum > ')' | nullary;

auto product_def      = grouped % -x3::lit('*');

auto first_addend     = (('-' >> x3::attr(-1) | -x3::lit('+') >> x3::attr(1)) >> product)[actions::process_addend{}];

auto addend           = (('+' >> x3::attr(1) | '-' >> x3::attr(-1)) > product)[actions::process_addend{}];

auto sum_def          = first_addend >> *addend;

auto expr_def         = -sum > x3::eoi;

auto resultExpr_def       = (tensor | variable) > ('=' | x3::unicode::lit(U"->")) >> expr;
// clang-format on

BOOST_SPIRIT_DEFINE(name, base_key, number, variable, index_label, proto_label,
                    index, index_groups, tensor, power, product, sum, expr,
                    symmetry_spec, resultExpr);

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
  x3::error_handler_result on_error(const Iterator & /*first*/,
                                    const Iterator & /*last*/,
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
struct PowerRule : helpers::annotate_position, helpers::error_handler {};
struct ProductRule : helpers::annotate_position, helpers::error_handler {};
struct SumRule : helpers::annotate_position, helpers::error_handler {};
struct ExprRule : helpers::annotate_position, helpers::error_handler {};
struct ResultRule : helpers::annotate_position, helpers::error_handler {};
struct IndexLabelRule : helpers::annotate_position, helpers::error_handler {};
struct ProtoLabelRule : helpers::annotate_position, helpers::error_handler {};
struct IndexRule : helpers::annotate_position, helpers::error_handler {};
struct IndexGroupRule : helpers::annotate_position, helpers::error_handler {};
struct SymmetrySpecRule : helpers::annotate_position, helpers::error_handler {};

template <typename Iterator>
struct ErrorHandler {
  Iterator begin;

  ErrorHandler(Iterator begin) : begin(std::move(begin)) {}

  void operator()(Iterator where, std::string expected) const {
    std::size_t offset = std::distance(begin.base(), where.base());
    throw SerializationError(offset, 1,
                             std::string("Parse failure at offset ") +
                                 std::to_string(offset) + ": expected " +
                                 expected);
  }
};

template <typename AST, typename StartRule, typename PositionCache>
AST parse(const StartRule &start, std::string_view input,
          PositionCache &positions) {
  const iterator_type input_begin(input.begin());
  const iterator_type input_end(input.end());

  ErrorHandler<iterator_type> error_handler(input_begin);

  AST ast;

  const auto parser = x3::with<error_handler_tag>(std::ref(
      error_handler))[x3::with<position_cache_tag>(std::ref(positions))[start]];

  auto begin = input_begin;
  try {
    bool success =
        x3::phrase_parse(begin, input_end, parser, x3::unicode::space, ast);

    if (!success) {
      // Normally, this shouldn't happen as any error should itself throw a
      // SerializationError already
      throw SerializationError(
          0, input.size(), "Parsing was unsuccessful for an unknown reason");
    }
    if (begin != input_end) {
      // This should also not happen as the parser requires matching EOI
      throw SerializationError(std::distance(input.begin(), begin.base()),
                               std::distance(begin.base(), input.end()),
                               "Couldn't parse the entire input");
    }
  } catch ([[maybe_unused]] const boost::spirit::x3::expectation_failure<
           iterator_type> &e) {
    // std::cout << "Caught expectation_failure\nwhere: " << e.where()
    //           << "\nwhat: " << e.what() << "\nwhich: " << e.which().data()
    //           << std::endl;
    throw;
  }

  return ast;
}

}  // namespace parse

bool is_base_key(std::string_view label) {
  if (!is_valid_utf8(label)) return false;
  auto begin = iterator_type(label.begin());
  const auto end = iterator_type(label.end());
  return x3::parse(begin, end, parse::base_key_chars) && begin == end;
}

transform::DefaultSymmetries to_default_symms(
    const DeserializationOptions &options) {
  // unspecified symmetries default to the active Context's; the Hermiticity is
  // resolved against each tensor's base_field downstream, as in the Tensor ctor
  const Context &ctx = get_default_context();
  transform::DefaultSymmetries symms{ctx.deserialization_symmetry(),
                                     ctx.deserialization_hermiticity(),
                                     ctx.deserialization_column_symmetry()};

  if (options.def_perm_symm.has_value()) {
    std::get<0>(symms) = options.def_perm_symm.value();
  }
  if (options.def_braket_symm.has_value()) {
    std::get<1>(symms) = options.def_braket_symm.value();
  }
  if (options.def_col_symm.has_value()) {
    std::get<2>(symms) = options.def_col_symm.value();
  }
  if (std::get<0>(symms) != Symmetry::Nonsymm) {
    // antisymmetry/symmetric bra and ket imply particle symmetry
    std::get<2>(symms) = ColumnSymmetry::Symm;
  }

  return symms;
}

/// @throw SerializationError if @p input is not valid UTF-8
void validate_input(std::string_view input) {
  if (!is_valid_utf8(input)) {
    const auto invalid_pos = std::distance(
        input.begin(), utf8::find_invalid(input.begin(), input.end()));
    throw SerializationError(
        invalid_pos, 1,
        "Invalid UTF-8 at offset " + std::to_string(invalid_pos));
  }
}

#define SEQUANT_DESERIALIZATION_SPECIALIZATION(ExprType, StartRule, ASTType, \
                                               transformFunc)                \
  template <>                                                                \
  ExprType from_string<ExprType>(std::string_view input,                     \
                                 const DeserializationOptions &options) {    \
    validate_input(input);                                                   \
    x3::position_cache<std::vector<iterator_type>> positions(                \
        iterator_type(input.begin()), iterator_type(input.end()));           \
    auto ast = parse::parse<ASTType>(StartRule, input, positions);           \
                                                                             \
    return transformFunc(ast, positions, iterator_type(input.begin()),       \
                         to_default_symms(options));                         \
  }

SEQUANT_DESERIALIZATION_SPECIALIZATION(ResultExpr, parse::resultExpr,
                                       ast::ResultExpr,
                                       transform::ast_to_result)
SEQUANT_DESERIALIZATION_SPECIALIZATION(ExprPtr, parse::expr, ast::Sum,
                                       transform::ast_to_expr)

}  // namespace sequant::io::serialization::v1
