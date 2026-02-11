#ifndef SEQUANT_DOMAIN_MBPT_DETAIL_CONCEPTS_HPP
#define SEQUANT_DOMAIN_MBPT_DETAIL_CONCEPTS_HPP

#include <SeQuant/core/index.hpp>

#include <ranges>

namespace sequant::mbpt::detail {

/// Concept for a range of a range of Index objects
template <typename T>
concept index_group_range =
    std::ranges::range<T> &&
    std::ranges::range<std::ranges::range_value_t<T>> &&
    std::convertible_to<
        std::ranges::range_value_t<std::ranges::range_value_t<T>>, Index>;

}  // namespace sequant::mbpt::detail

#endif  // SEQUANT_DOMAIN_MBPT_DETAIL_CONCEPTS_HPP
