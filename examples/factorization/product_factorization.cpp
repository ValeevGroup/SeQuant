//
// Created by Bimal Gaudel on 12/22/19.
//

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/evaluate/eval_tree.hpp>
#include <SeQuant/domain/factorize/factorizer.hpp>

#include "../sequant_setup.hpp"

#include <iostream>
#include <limits>
#include <memory>
#include <ostream>
#include <string>

#include <chrono>      // for seeding
#include <functional>  // for std::bind
#include <random>      // for std::mt19937

// drops A and P tensors
// sets scalars to 1
// expr is a Product of Tensors'
auto trimmed_prod = [](const auto& expr) {
  auto result = std::make_shared<Product>();
  for (const auto& xpr : *expr) {
    if (auto label = xpr->template as<Tensor>().label();
        label == L"A" || label == L"P")
      continue;
    else
      result->append(1, xpr);
  }
  return result;
};

// trim summands of expr
// expr is a Sum of Products' of Tensors'
auto trimmed_sum = [](const auto& expr) {
  auto result = std::make_shared<Sum>();
  for (const auto& xpr : *expr) result->append(trimmed_prod(xpr));
  return result;
};

// expr1 and expr2 are Products' Tensors'
void print_factors(const ExprPtr& expr1, const ExprPtr& expr2) {
  std::wcout << "$" << expr1->to_latex() << "$  vs  $" << expr2->to_latex()
             << "$  :  ";
  auto [fact1, fact2] = factorize::factorize_pair(expr1, expr2);
  std::wcout << "$" << fact1->to_latex() << "$  and $" << fact2->to_latex()
             << "$\n";
}

// expr is a Sum of Products' of Tensors'
void print_factors(const ExprPtr& expr) {
  for (auto ii = 0; ii < expr->size(); ++ii) {
    for (auto jj = ii + 1; jj < expr->size(); ++jj) {
      std::wcout << "term " << ii + 1 << " vs " << jj + 1 << "\n";
      print_factors(expr->at(ii), expr->at(jj));
      std::wcout << "\n";
    }
    std::wcout << std::endl;
  }
}

int main() {
  // global sequant setup...
  std::setlocale(LC_ALL, "en_US.UTF-8");
  std::wcout.precision(std::numeric_limits<double>::max_digits10);
  std::wcerr.precision(std::numeric_limits<double>::max_digits10);
  std::wostream::sync_with_stdio(true);
  std::wostream::sync_with_stdio(true);
  std::wcout.imbue(std::locale("en_US.UTF-8"));
  std::wcerr.imbue(std::locale("en_US.UTF-8"));
  std::wostream::sync_with_stdio(true);
  std::wostream::sync_with_stdio(true);
  sequant::detail::OpIdRegistrar op_id_registrar;

  sequant::mbpt::set_default_convention();

  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());
  Logger::get_instance().wick_stats = false;

  // get positions in a container where different 'kind' of tensors begin
  // container: {f_ov, f_ov, t_oo, t_oovv, g_oovv, g_oovv}
  // output:    {(0, 2), (2, 3), (3, 4), (4, 6)}
  const auto parts_indices = [](const auto& container) {
    using evaluate::EvalTree;  // for hashing tensors by their types
    using pos_type = std::size_t;

    container::svector<pos_type> indices;
    indices.push_back(0);
    for (auto ii = 0; ii < container.size();) {
      auto lead_hvalue = EvalTree(container.at(ii)).hash_value();
      auto jj = ii + 1;
      for (; jj < container.size(); ++jj) {
        auto trail_hvalue = EvalTree(container.at(jj)).hash_value();
        if (lead_hvalue != trail_hvalue) break;
      }
      ii = jj;
      indices.push_back(ii);
    }
    container::svector<std::tuple<pos_type, pos_type>> result;
    for (auto ii = 0; ii < indices.size() - 1; ++ii)
      result.push_back(std::make_tuple(indices.at(ii), indices.at(ii + 1)));
    return result;
  };

  // lambda to get the common tensors in container_t1 and container_t2
  auto common_tensors = [](const auto& container_t1, const auto& container_t2) {
    container::set<ExprPtr> common_t1, common_t2;

    for (const auto& t1 : container_t1)
      for (const auto& t2 : container_t2) {
        //
        // NOTE: As an example: t_{i j}^{a b} has the same
        // hash value as t_{a b}^{i j}. To hash such expressions
        // differently, use EvalTree(expr, false).
        //
        if (evaluate::EvalTree(t1).hash_value() ==
            evaluate::EvalTree(t2).hash_value()) {
          common_t1.insert(t1);
          common_t2.insert(t2);
        }
      }

    // convert set to vector
    auto exprptr_vec = [](const auto& container) {
      container::svector<ExprPtr> result;
      result.reserve(container.size());
      for (auto&& xpr : container) result.push_back(xpr);
      return result;
    };

    return std::make_tuple(exprptr_vec(common_t1), exprptr_vec(common_t2));
  };  // lambda common_tensors

  // CC equations
  auto cc_r = cceqvec{2, 2}(true, true, true, true);
  auto ampl = 2;
  auto expr = trimmed_sum(cc_r[ampl]);

  /*   for (auto ii = 0; ii < expr->size(); ++ii) */
  /*     for (auto jj = ii + 1; jj < expr->size(); ++jj) { */
  /*       auto pos1 = ii; */
  /*       auto pos2 = jj; */

  /*       auto term1 = trimmed_prod(expr->at(pos1)); */
  /*       auto term2 = trimmed_prod(expr->at(pos2)); */

  /*       std::wcout << "term #" << pos1 + 1                 // */
  /*                  << " = " << term1->to_latex() << "\n"   // */
  /*                  << "term #" << pos2 + 1                 // */
  /*                  << " = " << term2->to_latex() << "\n";  // */

  /*       auto [tensorsA, tensorsB] = common_tensors(*term1, *term2); */

  /*       std::wcout << "tensorsA.size() = " << tensorsA.size() */
  /*                  << "    tensorsB.size() = " << tensorsB.size() << "\n"; */

  /*       auto parts_term1 = parts_indices(tensorsA); */
  /*       auto parts_term2 = parts_indices(tensorsB); */
  /*       std::wcout << "parts_term1.size() = "        // */
  /*                  << parts_term1.size()             // */
  /*                  << "\nparts_term2.size() = "      // */
  /*                  << parts_term2.size() << "\n\n";  // */
  /*     } */

  std::wcout << R"(\section*{Equation})"
             << "\n"
             << R"(\begin{dmath*})"
             << "\n"
             << expr->to_latex() << "\n"
             << R"(\end{dmath*})" << std::endl;
  std::wcout << R"(\section*{Factors})"
             << "\n";
  print_factors(expr);

  return 0;
}
