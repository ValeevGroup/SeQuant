//
// created by Bimal Gaudel on Apr 18, 2020
//

#include "catch.hpp"

#include <SeQuant/domain/evaluate/eval_tree.hpp>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/tensor.hpp>

#include <tiledarray.h>

#include <clocale>
#include <iostream>
#include <memory>
#include <random>

using namespace sequant;
using namespace sequant::evaluate;

// initialize MADWorld
int argc = 1;
char* args[]{};
char** argv = {args};
auto& world = TA::initialize(argc, argv);

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
      return std::make_shared<sequant::Tensor>(label, bra_indices, ket_indices);
    };

TEST_CASE("CONSTRUCTOR TESTS", "[eval_tree]") {
  SECTION("Testing construction") {
    auto t = make_tensor_expr({"t", "i_1", "i_2", "a_1", "a_2"});
    auto g = make_tensor_expr({"g", "i_1", "i_2", "a_1", "a_2"});

    REQUIRE_NOTHROW(EvalTree(t));
    REQUIRE_NOTHROW(EvalTree(t, true));

    REQUIRE_NOTHROW(std::make_shared<Product>(Product({t, g})));
    REQUIRE_NOTHROW(std::make_shared<Sum>(Sum({t, g})));
  }
}

TEST_CASE("EVALUATIONS TESTS", "[eval_tree]") {
  using DTensorType = TA::TArrayD;
  using ContextMapType =
      sequant::container::map<HashType, std::shared_ptr<DTensorType>>;

  const size_t nocc = 10;
  const size_t nvirt = 20;

  TA::TiledRange tr_ov{{0, nocc}, {0, nvirt}};
  TA::TiledRange tr_oovv{{0, nocc}, {0, nocc}, {0, nvirt}, {0, nvirt}};

  auto T_ov = std::make_shared<DTensorType>(world, tr_ov);
  auto T_oovv = std::make_shared<DTensorType>(world, tr_oovv);
  auto G_oovv = std::make_shared<DTensorType>(world, tr_oovv);

  T_ov->fill_random();
  T_oovv->fill_random();
  G_oovv->fill_random();

  SECTION("Testing Sum type evaluations") {
    //
    auto t = make_tensor_expr({"t", "i_1", "i_2", "a_1", "a_2"});
    auto g = make_tensor_expr({"g", "i_1", "i_2", "a_1", "a_2"});

    DTensorType manual_result;
    manual_result("i, j, a, b") =
        (*T_oovv)("i, j, a, b") + (*G_oovv)("i, j, a, b");

    ContextMapType context;
    context.insert(
        ContextMapType::value_type(EvalTree(t).hash_value(), T_oovv));
    context.insert(
        ContextMapType::value_type(EvalTree(g).hash_value(), G_oovv));

    auto expr = std::make_shared<Sum>(Sum({g, t}));
    auto tree = EvalTree(expr);
    auto eval_result = tree.evaluate(context);

    auto manual_norm =
        std::sqrt(manual_result("0,1,2,3").dot(manual_result("0,1,2,3")));
    auto eval_norm =
        std::sqrt(eval_result("0,1,2,3").dot(eval_result("0,1,2,3")));

    REQUIRE(manual_norm == Approx(eval_norm));

    // sum by permutation test
    g = make_tensor_expr({"g", "i_1", "i_2", "a_2", "a_1"});
    context.insert(
        ContextMapType::value_type(EvalTree(g).hash_value(), G_oovv));

    manual_result("i, j, a, b") =
        (*T_oovv)("i, j, a, b") + (*G_oovv)("i, j, b, a");
    manual_norm =
        std::sqrt(manual_result("0,1,2,3").dot(manual_result("0,1,2,3")));

    expr = std::make_shared<Sum>(Sum({g, t}));
    tree = EvalTree(expr);
    eval_result = tree.evaluate(context);
    eval_norm = std::sqrt(eval_result("0,1,2,3").dot(eval_result("0,1,2,3")));

    REQUIRE(manual_norm == Approx(eval_norm));
  }

  SECTION("Testing product type evaluation") {
    auto t = make_tensor_expr({"t", "i_1", "a_1"});
    auto g = make_tensor_expr({"g", "i_1", "i_2", "a_1", "a_2"});

    DTensorType manual_result;
    manual_result("j,b") = (*T_ov)("i,a") * (*G_oovv)("i,j,a,b");

    ContextMapType context;
    context.insert(ContextMapType::value_type(EvalTree(t).hash_value(), T_ov));
    context.insert(
        ContextMapType::value_type(EvalTree(g).hash_value(), G_oovv));

    auto expr = std::make_shared<Product>(Product({t, g}));
    auto tree = EvalTree(expr);
    auto eval_result = tree.evaluate(context);

    auto manual_norm =
        std::sqrt(manual_result("0,1").dot(manual_result("0,1")));
    auto eval_norm = std::sqrt(eval_result("0,1").dot(eval_result("0,1")));

    REQUIRE(manual_norm == Approx(eval_norm));
  }

  SECTION("Testing antisymmetrization evaluation") {
    auto t = make_tensor_expr({"t", "i_1", "i_2", "a_1", "a_2"});
    auto A = make_tensor_expr({"A", "i_1", "i_2", "a_1", "a_2"});

    ContextMapType context;
    context.insert(
        ContextMapType::value_type(EvalTree(t).hash_value(), T_oovv));

    DTensorType manual_result;
    manual_result("i,j,a,b") = (*T_oovv)("i,j,a,b") - (*T_oovv)("i,j,b,a") +
                               (*T_oovv)("j,i,b,a") - (*T_oovv)("j,i,a,b");

    auto expr = std::make_shared<Product>(Product({A, t}));
    auto tree = EvalTree(expr);
    auto eval_result = tree.evaluate(context);

    auto manual_norm =
        std::sqrt(manual_result("0,1,2,3").dot(manual_result("0,1,2,3")));

    auto eval_norm =
        std::sqrt(eval_result("0,1,2,3").dot(eval_result("0,1,2,3")));

    REQUIRE(manual_norm == Approx(eval_norm));
  }
}
