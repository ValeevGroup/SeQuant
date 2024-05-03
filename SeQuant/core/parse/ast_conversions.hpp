//
// Created by Robert Adam on 2023-09-22
//

#ifndef SEQUANT_CORE_PARSE_AST_CONVERSIONS_HPP
#define SEQUANT_CORE_PARSE_AST_CONVERSIONS_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/parse/ast.hpp>
#include <SeQuant/core/parse_expr.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <boost/variant.hpp>

#include <algorithm>
#include <string>
#include <tuple>

namespace sequant::parse::transform {

template <typename AST, typename PositionCache, typename Iterator>
std::tuple<std::size_t, std::size_t> get_pos(const AST &ast,
                                             const PositionCache &cache,
                                             const Iterator &begin) {
  const auto range = cache.position_of(ast);

  return {std::distance(begin, range.begin()),
          std::distance(range.begin(), range.end())};
}

template <typename PositionCache, typename Iterator>
Index to_index(const parse::ast::Index &index,
               const PositionCache &position_cache, const Iterator &begin) {
  container::vector<Index> protoIndices;
  protoIndices.reserve(index.protoLabels.size());

  for (const parse::ast::IndexLabel &current : index.protoLabels) {
    try {
      std::wstring label = current.label + L"_" + std::to_wstring(current.id);
      IndexSpace space = IndexSpace::instance(label);
      protoIndices.push_back(Index(std::move(label), std::move(space)));
    } catch (const IndexSpace::bad_key &) {
      auto [offset, length] = get_pos(current, position_cache, begin);
      throw ParseError(offset, length,
                       "Unknown index space '" + toUtf8(current.label) +
                           "' in proto index specification");
    } catch (const std::invalid_argument &e) {
      auto [offset, length] = get_pos(current, position_cache, begin);
      throw ParseError(offset, length,
                       "Invalid index '" + toUtf8(current.label) + "_" +
                           std::to_string(current.id) + ": " + e.what());
    }
  }

  try {
    std::wstring label =
        index.label.label + L"_" + std::to_wstring(index.label.id);
    IndexSpace space = IndexSpace::instance(label);
    return Index(std::move(label), std::move(space), std::move(protoIndices));
  } catch (const IndexSpace::bad_key &e) {
    auto [offset, length] = get_pos(index.label, position_cache, begin);
    throw ParseError(offset, length,
                     "Unknown index space '" + toUtf8(index.label.label) +
                         "' in index specification");
  } catch (const std::invalid_argument &e) {
    auto [offset, length] = get_pos(index.label, position_cache, begin);
    throw ParseError(offset, length,
                     "Invalid index '" + toUtf8(index.label.label) + "_" +
                         std::to_string(index.label.id) + ": " + e.what());
  }
}

template <typename PositionCache, typename Iterator>
std::tuple<container::vector<Index>, container::vector<Index>> make_indices(
    const parse::ast::IndexGroups &groups, const PositionCache &position_cache,
    const Iterator &begin) {
  container::vector<Index> braIndices;
  container::vector<Index> ketIndices;

  static_assert(std::is_same_v<decltype(groups.bra), decltype(groups.ket)>,
                "Types for bra and ket indices must be equal for pointer "
                "aliasing to work");

  const auto *bra = &groups.bra;
  const auto *ket = &groups.ket;
  if (groups.reverse_bra_ket) {
    bra = &groups.ket;
    ket = &groups.bra;
  }

  braIndices.reserve(bra->size());
  ketIndices.reserve(ket->size());

  for (const parse::ast::Index &current : *bra) {
    braIndices.push_back(to_index(current, position_cache, begin));
  }
  for (const parse::ast::Index &current : *ket) {
    ketIndices.push_back(to_index(current, position_cache, begin));
  }

  return {std::move(braIndices), std::move(ketIndices)};
}

template <typename Iterator>
Symmetry to_symmetry(char c, std::size_t offset, const Iterator &begin,
                     Symmetry default_symmetry) {
  if (c == parse::ast::Tensor::unspecified_symmetry) {
    return default_symmetry;
  }

  switch (c) {
    case 'A':
    case 'a':
      return Symmetry::antisymm;
    case 'S':
    case 's':
      return Symmetry::symm;
    case 'N':
    case 'n':
      return Symmetry::nonsymm;
  }

  throw ParseError(offset, 1,
                   std::string("Invalid symmetry specifier '") + c + "'");
}

template <typename PositionCache, typename Iterator>
Constant to_constant(const parse::ast::Number &number,
                     const PositionCache &position_cache,
                     const Iterator &begin) {
  if (static_cast<std::int64_t>(number.numerator) == number.numerator &&
      static_cast<std::int64_t>(number.denominator) == number.denominator) {
    // Integer fraction
    return Constant(
        ::sequant::rational(static_cast<std::int64_t>(number.numerator),
                            static_cast<std::int64_t>(number.denominator)));
  } else {
    // Construct from floating point value
    return Constant(::sequant::rational(number.numerator / number.denominator));
  }
}

template <typename PositionCache, typename Iterator>
ExprPtr ast_to_expr(const parse::ast::Product &product,
                    const PositionCache &position_cache, const Iterator &begin,
                    Symmetry default_symmetry);
template <typename PositionCache, typename Iterator>
ExprPtr ast_to_expr(const parse::ast::Sum &sum,
                    const PositionCache &position_cache, const Iterator &begin,
                    Symmetry default_symmetry);

template <typename PositionCache, typename Iterator>
ExprPtr ast_to_expr(const parse::ast::NullaryValue &value,
                    const PositionCache &position_cache, const Iterator &begin,
                    Symmetry default_symmetry) {
  struct Transformer {
    std::reference_wrapper<const PositionCache> position_cache;
    std::reference_wrapper<const Iterator> begin;
    std::reference_wrapper<Symmetry> default_symmetry;

    ExprPtr operator()(const parse::ast::Product &product) const {
      return ast_to_expr<PositionCache>(product, position_cache.get(),
                                        begin.get(), default_symmetry);
    }

    ExprPtr operator()(const parse::ast::Sum &sum) const {
      return ast_to_expr<PositionCache>(sum, position_cache.get(), begin.get(),
                                        default_symmetry);
    }

    ExprPtr operator()(const parse::ast::Tensor &tensor) const {
      auto [braIndices, ketIndices] =
          make_indices(tensor.indices, position_cache.get(), begin.get());

      auto [offset, length] =
          get_pos(tensor, position_cache.get(), begin.get());

      return ex<Tensor>(tensor.name, std::move(braIndices),
                        std::move(ketIndices),
                        to_symmetry(tensor.symmetry, offset + length - 1,
                                    begin.get(), default_symmetry));
    }

    ExprPtr operator()(const parse::ast::Variable &variable) const {
      if (variable.conjugated) {
        return ex<Variable>(variable.name + L"^*");
      } else {
        return ex<Variable>(variable.name);
      }
    }

    ExprPtr operator()(const parse::ast::Number &number) const {
      return ex<Constant>(
          to_constant(number, position_cache.get(), begin.get()));
    }
  };

  return boost::apply_visitor(
      Transformer{std::ref(position_cache), std::ref(begin),
                  std::ref(default_symmetry)},
      value);
}

template <typename T, typename... Ts>
bool holds_alternative(const boost::variant<Ts...> &v) noexcept {
  return boost::get<T>(&v) != nullptr;
}

template <typename PositionCache, typename Iterator>
ExprPtr ast_to_expr(const parse::ast::Product &product,
                    const PositionCache &position_cache, const Iterator &begin,
                    Symmetry default_symmetry) {
  if (product.factors.empty()) {
    // This shouldn't happen
    assert(false);
    throw std::runtime_error(
        "ast_to_expr: Reached supposed-to-be unreachable code");
  }

  if (product.factors.size() == 1) {
    return ast_to_expr(product.factors.front(), position_cache, begin,
                       default_symmetry);
  }

  std::vector<ExprPtr> factors;
  factors.reserve(product.factors.size());
  Constant prefactor(1);

  // We perform constant folding
  for (const parse::ast::NullaryValue &value : product.factors) {
    if (holds_alternative<parse::ast::Number>(value)) {
      prefactor *= to_constant(boost::get<parse::ast::Number>(value),
                               position_cache, begin);
    } else {
      factors.push_back(
          ast_to_expr(value, position_cache, begin, default_symmetry));
    }
  }

  if (factors.empty()) {
    // Only constants
    return ex<Constant>(std::move(prefactor));
  }

  if (factors.size() == 1 && prefactor.value() == 1) {
    return factors.front();
  }

  return ex<Product>(prefactor.value(), std::move(factors),
                     Product::Flatten::No);
}

template <typename PositionCache, typename Iterator>
ExprPtr ast_to_expr(const parse::ast::Sum &sum,
                    const PositionCache &position_cache, const Iterator &begin,
                    Symmetry default_symmetry) {
  if (sum.summands.empty()) {
    return {};
  }
  if (sum.summands.size() == 1) {
    return ast_to_expr(sum.summands.front(), position_cache, begin,
                       default_symmetry);
  }

  std::vector<ExprPtr> summands;
  summands.reserve(sum.summands.size());
  std::transform(
      sum.summands.begin(), sum.summands.end(), std::back_inserter(summands),
      [&](const parse::ast::Product &product) {
        return ast_to_expr(product, position_cache, begin, default_symmetry);
      });

  return ex<Sum>(std::move(summands));
}

}  // namespace sequant::parse::transform

#endif  // SEQUANT_CORE_PARSE_AST_CONVERSIONS_HPP
