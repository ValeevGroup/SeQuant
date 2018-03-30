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

template <class T>
struct is_shared_ptr
    : std::false_type {};
template <class T>
struct is_shared_ptr<std::shared_ptr<T> >
: std::true_type {};
template <class T> using is_shared_ptr_v = is_shared_ptr<T>::value;

}
}

#endif //SEQUANT2_META_HPP
