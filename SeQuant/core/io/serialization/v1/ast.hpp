//
// Created by Robert Adam on 2023-09-21
//

#ifndef SEQUANT_CORE_PARSE_V1_AST_HPP
#define SEQUANT_CORE_PARSE_V1_AST_HPP

#define BOOST_SPIRIT_X3_UNICODE
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/optional.hpp>
#include <boost/spirit/home/x3.hpp>
#include <boost/spirit/home/x3/support/ast/position_tagged.hpp>
#include <boost/variant.hpp>

#include <string>
#include <variant>
#include <vector>

namespace sequant::io::serialization::v1::ast {

struct IndexLabel : boost::spirit::x3::position_tagged {
  std::wstring label;
  unsigned int id;

  IndexLabel(std::wstring label = {}, unsigned int id = {})
      : label(std::move(label)), id(id) {}
};

struct Index : boost::spirit::x3::position_tagged {
  IndexLabel label;
  std::vector<IndexLabel> protoLabels;

  Index(IndexLabel label = {}, std::vector<IndexLabel> protoLabels = {})
      : label(std::move(label)), protoLabels(std::move(protoLabels)) {}
};

struct Number : boost::spirit::x3::position_tagged {
  double numerator;
  double denominator;
  /// whether an `i` abutted the digits, i.e. the literal sits on the
  /// imaginary axis (`2i`, `1/2i`)
  bool imaginary;

  Number(double numerator = {}, double denominator = 1, bool imaginary = false)
      : numerator(numerator), denominator(denominator), imaginary(imaginary) {}
};

struct Variable : boost::spirit::x3::position_tagged {
  /// the name as written, a trailing conjugation mark included
  std::wstring name;

  /// @note not explicit, and the struct is not a fusion sequence: the
  ///       variable rule's attribute is the bare name, which x3 assigns
  ///       through this conversion
  Variable(std::wstring name = {}) : name(std::move(name)) {}
};

struct IndexGroups : boost::spirit::x3::position_tagged {
  std::vector<Index> bra;
  std::vector<Index> ket;
  std::vector<Index> auxiliaries;
  bool reverse_bra_ket;

  IndexGroups(std::vector<Index> bra = {}, std::vector<Index> ket = {},
              std::vector<Index> auxiliaries = {}, bool reverse_bra_ket = {})
      : bra(std::move(bra)),
        ket(std::move(ket)),
        auxiliaries(std::move(auxiliaries)),
        reverse_bra_ket(reverse_bra_ket) {}
};

struct SymmetrySpec : boost::spirit::x3::position_tagged {
  static constexpr char unspecified = '\0';
  char perm_symm = unspecified;
  char braket_symm = unspecified;
  char column_symm = unspecified;
  char conjugation_parity = unspecified;
};

// represents AbstractTensor, i.e. Tensor or NormalOperator
struct Tensor : boost::spirit::x3::position_tagged {
  /// the name as written, trailing state marks (`⁺`, `꙳`) included; the
  /// conversion to Tensor splits them off
  std::wstring name;
  IndexGroups indices;
  boost::optional<SymmetrySpec> symmetry;

  Tensor(std::wstring name = {}, IndexGroups indices = {},
         boost::optional<SymmetrySpec> symmetry = {})
      : name(std::move(name)),
        indices(std::move(indices)),
        symmetry(std::move(symmetry)) {}
};

struct Product;
struct Sum;
struct RealImagPart;

struct Power : boost::spirit::x3::position_tagged {
  boost::variant<Number, Variable> base;
  Number exponent;
  bool conjugated;

  Power() noexcept = default;
  Power(boost::variant<Number, Variable> base, Number exponent,
        bool conjugated = false)
      : base(std::move(base)),
        exponent(std::move(exponent)),
        conjugated(conjugated) {}
};

using NullaryValue =
    boost::variant<Number, Tensor, Variable, Power, Product, Sum, RealImagPart>;

struct Product : boost::spirit::x3::position_tagged {
  std::vector<NullaryValue> factors;

  Product() noexcept = default;

  template <typename T>
  Product(T value);

  Product(std::vector<NullaryValue> factors);

  // Required to use as a container
  using value_type = decltype(factors)::value_type;
};

struct Sum : boost::spirit::x3::position_tagged {
  std::vector<Product> summands;

  Sum() noexcept = default;

  Sum(std::vector<Product> summands);

  // Required to use as a container
  using value_type = decltype(summands)::value_type;
};

template <typename T>
Product::Product(T value) : factors({std::move(value)}) {}

/// the real or the imaginary part of a wrapped expression, spelled `Re[...]`
/// or `Im[...]`
struct RealImagPart : boost::spirit::x3::position_tagged {
  /// false spells `Re`, true spells `Im`
  bool imaginary;
  Sum inner;

  RealImagPart() noexcept : imaginary(false) {}
  RealImagPart(bool imaginary, Sum inner)
      : imaginary(imaginary), inner(std::move(inner)) {}
};

struct ResultExpr : boost::spirit::x3::position_tagged {
  std::variant<Tensor, Variable> lhs;
  Sum rhs;

  ResultExpr(Variable variable = {}, Sum expr = {})
      : lhs(std::move(variable)), rhs(std::move(expr)) {}

  ResultExpr(Tensor tensor, Sum expr)
      : lhs(std::move(tensor)), rhs(std::move(expr)) {}
};

}  // namespace sequant::io::serialization::v1::ast

BOOST_FUSION_ADAPT_STRUCT(sequant::io::serialization::v1::ast::IndexLabel,
                          label, id);
BOOST_FUSION_ADAPT_STRUCT(sequant::io::serialization::v1::ast::Index, label,
                          protoLabels);
BOOST_FUSION_ADAPT_STRUCT(sequant::io::serialization::v1::ast::Number,
                          numerator, denominator, imaginary);
BOOST_FUSION_ADAPT_STRUCT(sequant::io::serialization::v1::ast::IndexGroups, bra,
                          ket, auxiliaries, reverse_bra_ket);
BOOST_FUSION_ADAPT_STRUCT(sequant::io::serialization::v1::ast::SymmetrySpec,
                          perm_symm, braket_symm, column_symm,
                          conjugation_parity);
BOOST_FUSION_ADAPT_STRUCT(sequant::io::serialization::v1::ast::Tensor, name,
                          indices, symmetry);
BOOST_FUSION_ADAPT_STRUCT(sequant::io::serialization::v1::ast::Power, base,
                          exponent, conjugated);

BOOST_FUSION_ADAPT_STRUCT(sequant::io::serialization::v1::ast::RealImagPart,
                          imaginary, inner);

BOOST_FUSION_ADAPT_STRUCT(sequant::io::serialization::v1::ast::Product,
                          factors);
BOOST_FUSION_ADAPT_STRUCT(sequant::io::serialization::v1::ast::Sum, summands);
BOOST_FUSION_ADAPT_STRUCT(sequant::io::serialization::v1::ast::ResultExpr, lhs,
                          rhs);

#endif  // SEQUANT_CORE_PARSE_AST_V1_HPP
