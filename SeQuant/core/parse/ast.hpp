//
// Created by Robert Adam on 2023-09-21
//

#ifndef SEQUANT_CORE_PARSE_AST_HPP
#define SEQUANT_CORE_PARSE_AST_HPP

#define BOOST_SPIRIT_X3_UNICODE
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/optional.hpp>
#include <boost/spirit/home/x3.hpp>
#include <boost/spirit/home/x3/support/ast/position_tagged.hpp>
#include <boost/variant.hpp>

#include <cstdint>
#include <string>
#include <vector>

namespace sequant::parse::ast {

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

  Number(double numerator = {}, double denominator = 1)
      : numerator(numerator), denominator(denominator) {}
};

struct Variable : boost::spirit::x3::position_tagged {
  std::wstring name;
  bool conjugated;

  Variable(std::wstring name = {}, bool conjugated = false)
      : name(std::move(name)), conjugated(conjugated) {}

  // Required to use as a container
  using value_type = decltype(name)::value_type;
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
  char particle_symm = unspecified;
};

struct Tensor : boost::spirit::x3::position_tagged {
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

using NullaryValue = boost::variant<Number, Tensor, Variable, Product, Sum>;

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

Product::Product(std::vector<NullaryValue> factors)
    : factors(std::move(factors)) {}

Sum::Sum(std::vector<Product> summands) : summands(std::move(summands)) {}

}  // namespace sequant::parse::ast

BOOST_FUSION_ADAPT_STRUCT(sequant::parse::ast::IndexLabel, label, id);
BOOST_FUSION_ADAPT_STRUCT(sequant::parse::ast::Index, label, protoLabels);
BOOST_FUSION_ADAPT_STRUCT(sequant::parse::ast::Number, numerator, denominator);
BOOST_FUSION_ADAPT_STRUCT(sequant::parse::ast::Variable, name, conjugated);
BOOST_FUSION_ADAPT_STRUCT(sequant::parse::ast::IndexGroups, bra, ket,
                          auxiliaries, reverse_bra_ket);
BOOST_FUSION_ADAPT_STRUCT(sequant::parse::ast::SymmetrySpec, perm_symm,
                          braket_symm, particle_symm);
BOOST_FUSION_ADAPT_STRUCT(sequant::parse::ast::Tensor, name, indices, symmetry);

BOOST_FUSION_ADAPT_STRUCT(sequant::parse::ast::Product, factors);
BOOST_FUSION_ADAPT_STRUCT(sequant::parse::ast::Sum, summands);

#endif  // SEQUANT_CORE_PARSE_AST_HPP
