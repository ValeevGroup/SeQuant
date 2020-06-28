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

// get a sequant Tensor made out of specs
// specs -> {label, b1, ..., b(n/2), k1, ..., k(n/2)}
// eg. {"g", "i_1", "i_2", "a_1", "a_2"}
auto make_tensor_expr =
    [](const sequant::container::svector<std::string>& specs) {
      // only equal bra-ket ranks are expected
      assert((specs.size() > 2) && (specs.size() % 2 != 0));
      std::wstring label = std::wstring(specs[0].begin(), specs[0].end());
      sequant::container::svector<sequant::Index, 4> bra_indices, ket_indices;
      for (auto i = 1; i < specs.size(); ++i) {
        if (i <= specs.size() / 2)
          bra_indices.push_back(
              sequant::Index(std::wstring(specs[i].begin(), specs[i].end())));
        else
          ket_indices.push_back(
              sequant::Index(std::wstring(specs[i].begin(), specs[i].end())));
      }
      return std::make_shared<sequant::Tensor>(label, bra_indices, ket_indices,
                                               Symmetry::antisymm,
                                               BraKetSymmetry::conjugate);
    };

// drops A and S tensors
// sets scalars to 1
// expr is a Product of Tensors'
auto trimmed_prod = [](const auto& expr) {
  auto result = std::make_shared<Product>();
  for (const auto& xpr : *expr) {
    if (auto label = xpr->template as<Tensor>().label();
        label == L"A" || label == L"S")
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

  ExprPtr expr{nullptr};

  //
  // CC equations
  auto cc_r = cceqvec{2, 2}(true, true, true, true);
  auto ampl = 2;
  expr = trimmed_sum(cc_r[ampl]);

  auto prod1 = std::make_shared<Product>();
  prod1->append(1, make_tensor_expr({"g", "i_3", "i_4", "i_1", "a_3"}));
  prod1->append(1, make_tensor_expr({"t", "a_3", "i_2"}));
  prod1->append(1, make_tensor_expr({"t", "a_1", "a_2", "i_3", "i_4"}));

  auto prod2 = std::make_shared<Product>();
  prod2->append(1, make_tensor_expr({"g", "i_3", "i_4", "i_1", "a_3"}));
  prod2->append(1, make_tensor_expr({"t", "a_1", "i_3"}));
  prod2->append(1, make_tensor_expr({"t", "a_2", "i_4"}));
  prod2->append(1, make_tensor_expr({"t", "a_3", "i_2"}));

  auto sum = std::make_shared<Sum>();
  // sum->append(prod1);
  // sum->append(prod2);
  sum->append(expr->at(12));
  sum->append(expr->at(25));

  std::wcout << R"(\section*{Equation})"
             << "\n"
             << R"(\begin{dmath*})"
             << "\n"
             << expr->to_latex()
             << "\n"
             // << sum->to_latex() << "\n"
             << R"(\end{dmath*})" << std::endl;

  std::wcout << R"(\section*{Factors})"
             << "\n";
  print_factors(expr);
  // print_factors(sum);

  return 0;
}
