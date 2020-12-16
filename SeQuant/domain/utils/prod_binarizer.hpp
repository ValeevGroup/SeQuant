#ifndef SEQUANT_UTILS_PROD_BINARIZER_HPP
#define SEQUANT_UTILS_PROD_BINARIZER_HPP

#include "eval_expr.hpp"

namespace sequant::utils {

/**
 * Convert eval_sequence of positions of factors to a binary_expr node of
 * eval_expr.
 */
struct prod_binarizer {
 private:
  const Product& prod;

 public:
  prod_binarizer(const Product& p) : prod{p} {};

  eval_expr operator()(size_t pos) const {
    return eval_expr{prod.factors().at(pos)->as<Tensor>()};
  };

  eval_expr operator()(const eval_expr& expr1, const eval_expr& expr2) const {
    return eval_expr{expr1, expr2};
  }
};

}  // namespace sequant::utils

#endif  // SEQUANT_UTILS_PROD_BINARIZER_HPP
