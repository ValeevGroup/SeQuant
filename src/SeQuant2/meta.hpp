//
// Created by Eduard Valeyev on 3/29/18.
//

#ifndef SEQUANT2_META_HPP
#define SEQUANT2_META_HPP

#include <type_traits>
#include <memory>

namespace sequant2 {
namespace meta {

template<typename T>
struct type_printer;

template<class T>
struct is_shared_ptr
    : std::false_type {
};
template<class T>
struct is_shared_ptr<std::shared_ptr<T> >
    : std::true_type {
};
template<class T> static constexpr bool is_shared_ptr_v = is_shared_ptr<T>::value;

template<typename T, typename = std::void_t<>>
struct is_less_than_comparable : public std::false_type {};

template<typename T>
struct is_less_than_comparable<T, std::void_t<decltype(std::declval<const T &>() < std::declval<const T &>())>>
    : public std::true_type {
};

template<typename T> static constexpr bool is_less_than_comparable_v = is_less_than_comparable<T>::value;

}  // namespae meta
}  // namespace sequant2

#endif //SEQUANT2_META_HPP
