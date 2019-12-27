//
// Created by Bimal Gaudel on 12/20/19.
//

#ifndef SEQUANT_FACTORIZER_HPP
#define SEQUANT_FACTORIZER_HPP

#include "path_tree.hpp"

#include <SeQuant/core/container.hpp>

#include <cstddef>
#include <string>

namespace sequant {
namespace factorize {
namespace detail {

  void path_scanner(sequant::container::svector<std::shared_ptr<PathTree>>&);

}  // namespace detail
}  // namespace factorize
}  // namespace sequant

#endif  // SEQUANT_FACTORIZER_HPP
