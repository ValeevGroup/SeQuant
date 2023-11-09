//
// Created by Robert Adam on 2023-09-21
//

#ifndef SEQUANT_CORE_PARSE_SEMANTIC_ACTIONS_HPP
#define SEQUANT_CORE_PARSE_SEMANTIC_ACTIONS_HPP

#include "ast.hpp"

#include <SeQuant/core/utility/string.hpp>

#define BOOST_SPIRIT_X3_UNICODE
#include <boost/spirit/home/x3.hpp>
#include <boost/variant/apply_visitor.hpp>

#include <type_traits>
#include <vector>

namespace sequant::parse::actions {

namespace x3 = boost::spirit::x3;

/// Auxiliary type traits
namespace {

// ctx_tags
template <typename Context>
struct ctx_tags {
  using attribute_type = std::remove_reference_t<
      std::remove_cv_t<decltype(x3::_attr(std::declval<Context>()))>>;
  using target_type = std::remove_reference_t<
      std::remove_cv_t<decltype(x3::_val(std::declval<Context>()))>>;
};

template <typename Context>
using ctx_attribute_t = typename ctx_tags<Context>::attribute_type;
template <typename Context>
using ctx_target_t = typename ctx_tags<Context>::target_type;

}  // namespace

// Helper function to get the current value out of a variant
template <typename Variant>
auto get_val(Variant &variant) {
  return boost::apply_visitor([](auto &element) { return element; }, variant);
}

struct process_addend {
  template <typename Context>
  void operator()(Context &ctx) const {
    static_assert(std::is_same_v<ctx_target_t<Context>, ast::Sum>);

    using boost::fusion::at_c;

    int factor;
    if constexpr (std::is_same_v<std::remove_cv_t<std::remove_reference_t<
                                     decltype(at_c<0>(x3::_attr(ctx)))>>,
                                 int>) {
      factor = at_c<0>(x3::_attr(ctx));
    } else {
      factor = get_val(at_c<0>(x3::_attr(ctx)));
    }

    ast::Product &prod = at_c<1>(x3::_attr(ctx));

    if (factor != 1) {
      prod.factors.insert(prod.factors.begin(), ast::Number(factor));
    }

    x3::_val(ctx).summands.emplace_back(std::move(prod));
  }
};

}  // namespace sequant::parse::actions

#endif  // SEQUANT_CORE_PARSE_SEMANTIC_ACTIONS_HPP
