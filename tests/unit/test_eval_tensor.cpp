//
// created by Bimal Gaudel on Mar 11, 2020
//

#include <clocale>
#include <iostream>
#include <memory>
#include <random>

#include <tiledarray.h>

#include "../../src/SeQuant/core/expr.hpp"
#include "../../src/SeQuant/core/op.hpp"
#include "../../src/SeQuant/core/runtime.hpp"
#include "../../src/SeQuant/core/space.hpp"
#include "../../src/SeQuant/core/tensor.hpp"
#include "../../src/SeQuant/core/utility.hpp"
#include "../../src/SeQuant/domain/evaluate/eval_context.hpp"
#include "../../src/SeQuant/domain/mbpt/convention.hpp"
#include "catch.hpp"

#include "../../examples/sequant_setup.hpp"
#include "../../src/SeQuant/domain/evaluate/eval_tensor.hpp"
#include "../../src/SeQuant/domain/evaluate/eval_tensor_builder.hpp"

using namespace sequant;
using namespace std;

using namespace sequant::evaluate;
using DTensorType = TA::TArrayD;
using ContextMapType =
    sequant::container::map<sequant::ExprPtr, std::shared_ptr<DTensorType>>;

// initialize MADWorld
int argc = 1;
char* args[]{};
char** argv = {args};
auto& world = TA::initialize(argc, argv);

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

// context map builder for evaluating the eval_tensor
auto builder = EvalTensorBuilder<DTensorType>();

const size_t nocc = 10;
const size_t nvirt = 20;

TA::TiledRange tr_oo{{0, nocc}, {0, nocc}};
TA::TiledRange tr_ov{{0, nocc}, {0, nvirt}};
TA::TiledRange tr_vv{{0, nvirt}, {0, nvirt}};
TA::TiledRange tr_oooo{{0, nocc}, {0, nocc}, {0, nocc}, {0, nocc}};
TA::TiledRange tr_ooov{{0, nocc}, {0, nocc}, {0, nocc}, {0, nvirt}};
TA::TiledRange tr_oovv{{0, nocc}, {0, nocc}, {0, nvirt}, {0, nvirt}};
TA::TiledRange tr_ovov{{0, nocc}, {0, nvirt}, {0, nocc}, {0, nvirt}};
TA::TiledRange tr_ovvv{{0, nocc}, {0, nvirt}, {0, nvirt}, {0, nvirt}};
TA::TiledRange tr_vvvv{{0, nvirt}, {0, nvirt}, {0, nvirt}, {0, nvirt}};

TEST_CASE("EVAL_TENSOR_CONSTRUCTOR_TESTS", "[eval_tensor]") {
  SECTION("Intermediate construction") {
    EvalTensorIntermediate<DTensorType> evt_imed;
    REQUIRE(!evt_imed.is_leaf());
    REQUIRE(evt_imed.get_operation() == Operation::INVALID);
    REQUIRE(evt_imed.get_left_tensor() == nullptr);
    REQUIRE(evt_imed.get_right_tensor() == nullptr);
  }

  SECTION("Leaf construction") {
    EvalTensorLeaf<DTensorType> evt_leaf(nullptr);
    REQUIRE(evt_leaf.is_leaf());
  }

  // SECTION("Scalars in the evaluation tree") {
  //   auto visitor = [](const EvalTensor<DTensorType>& node) {
  //       std::wcout << "scalar found: " << node.get_scalar() << std::endl;
  //   };
  //   // global sequant setup...
  //   sequant::detail::OpIdRegistrar op_id_registrar;
  //   sequant::mbpt::set_default_convention();
  //   TensorCanonicalizer::register_instance(
  //       std::make_shared<DefaultTensorCanonicalizer>());
  //   auto cc_r = cceqvec(2, 2)(true, true, true, true);

  //   auto builder = EvalTensorBuilder<DTensorType>();
  //   auto r1_tree = builder.build_tree(cc_r[1]);
  //   r1_tree->visit(visitor);
  // }
}

TEST_CASE("EVAL_TENSOR_EVALUATE_TESTS", "[eval_tensor_builder]") {
  // global sequant setup...
  sequant::detail::OpIdRegistrar op_id_registrar;
  sequant::mbpt::set_default_convention();

  using sequant::Product;
  using sequant::Sum;

  // creating some random tensors
  auto T_ov = std::make_shared<DTensorType>(world, tr_ov);
  auto T_oovv = std::make_shared<DTensorType>(world, tr_oovv);
  auto G_oovv = std::make_shared<DTensorType>(world, tr_oovv);
  //

  T_ov->fill_random();
  T_oovv->fill_random();
  G_oovv->fill_random();

  SECTION("Testing Sum type evaluation") {
    auto t = make_tensor_expr({"t", "i_1", "i_2", "a_1", "a_2"});
    auto g = make_tensor_expr({"g", "i_1", "i_2", "a_1", "a_2"});

    ContextMapType context;
    context.insert(ContextMapType::value_type(t, T_oovv));
    context.insert(ContextMapType::value_type(g, G_oovv));
    auto ev_context = EvalContext(context);

    DTensorType manual_sum;
    manual_sum("i, j, a, b") =
        (*G_oovv)("i, j, a, b") + (*T_oovv)("i, j, a, b");

    auto expr = std::make_shared<Sum>(Sum({g, t}));
    auto tree = builder.build_tree(expr);
    auto eval_sum = tree->evaluate(ev_context.get_map());

    auto manual_norm =
        std::sqrt(manual_sum("i,j,a,b").dot(manual_sum("i,j,a,b")));

    auto eval_norm = std::sqrt(eval_sum("i,j,a,b").dot(eval_sum("i,j,a,b")));

    REQUIRE(manual_norm == Approx(eval_norm));

    // sum by permutation test
    t = make_tensor_expr({"t", "i_1", "i_2", "a_1", "a_2"});
    g = make_tensor_expr({"g", "i_1", "i_2", "a_2", "a_1"});
    context.clear();
    context.insert(ContextMapType::value_type(t, T_oovv));
    context.insert(ContextMapType::value_type(g, G_oovv));
    ev_context = EvalContext(context);
    expr = std::make_shared<Sum>(Sum({g, t}));
    tree = builder.build_tree(expr);

    eval_sum = tree->evaluate(ev_context.get_map());

    manual_sum("i, j, a, b") =
        (*G_oovv)("i, j, a, b") + (*T_oovv)("i, j, b, a");
    manual_norm = std::sqrt(manual_sum("i,j,a,b").dot(manual_sum("i,j,a,b")));

    eval_norm = std::sqrt(eval_sum("i,j,a,b").dot(eval_sum("i,j,a,b")));

    REQUIRE(manual_norm == Approx(eval_norm));
  }

  SECTION("Testing Product type evaluation") {
    auto t = make_tensor_expr({"t", "i_1", "a_1"});
    auto g = make_tensor_expr({"g", "i_1", "i_2", "a_1", "a_2"});

    ContextMapType context;
    context.insert(ContextMapType::value_type(t, T_ov));
    context.insert(ContextMapType::value_type(g, G_oovv));
    auto ev_context = EvalContext(context);

    DTensorType manual_prod;
    manual_prod("j,b") = (*T_ov)("i,a") * (*G_oovv)("i,j,a,b");

    auto expr = std::make_shared<Product>(Product({t, g}));
    auto tree = builder.build_tree(expr);
    auto eval_prod = tree->evaluate(ev_context.get_map());

    auto manual_norm = std::sqrt(manual_prod("j,b").dot(manual_prod("j,b")));

    auto eval_norm = std::sqrt(eval_prod("j,b").dot(eval_prod("j,b")));

    REQUIRE(manual_norm == Approx(eval_norm));
  }

  SECTION("Testing antisymmetrization evaluation") {
    auto t = make_tensor_expr({"t", "i_1", "i_2", "a_1", "a_2"});
    auto A = make_tensor_expr({"A", "i_1", "i_2", "a_1", "a_2"});

    ContextMapType context;
    context.insert(ContextMapType::value_type(t, T_oovv));
    auto ev_context = EvalContext(context);

    DTensorType manual_result;
    manual_result("i,j,a,b") = (*T_oovv)("i,j,a,b") - (*T_oovv)("i,j,b,a") +
                               (*T_oovv)("j,i,b,a") - (*T_oovv)("j,i,a,b");

    auto manual_norm =
        std::sqrt(manual_result("i,j,a,b").dot(manual_result("i,j,a,b")));

    auto expr = std::make_shared<Product>(Product({A, t}));
    auto tree = builder.build_tree(expr);
    auto eval_result = tree->evaluate(ev_context.get_map());

    auto eval_norm =
        std::sqrt(eval_result("i,j,a,b").dot(eval_result("i,j,a,b")));

    REQUIRE(manual_norm == Approx(eval_norm));
  }

  SECTION("Testing missing data tensor") {
    auto seq_tensor_bad = make_tensor_expr({"t", "a_1", "i_1", "a_2", "a_3"});
    auto seq_tensor_good = make_tensor_expr({"t", "i_1", "a_1", "a_2", "a_3"});
    ContextMapType context;
    context.insert(ContextMapType::value_type(seq_tensor_good, T_ov));
    auto ev_context = EvalContext(context);

    auto expr = seq_tensor_bad;
    auto tree = builder.build_tree(expr);

    REQUIRE_THROWS_AS(tree->evaluate(ev_context.get_map()), std::logic_error);
  }

  //SECTION("Bra and ket indices") {
  //  auto visitor = [](const EvalTensor<DTensorType>& evtensor) {
  //    std::wcout << "Indices: ";
  //    for (const auto& idx : evtensor.get_indices())
  //      std::wcout << idx.to_latex();
  //    std::wcout << std::endl;
  //  };
  //  auto t = make_tensor_expr({"t", "i_1", "i_2", "a_1", "a_2"});
  //  auto g = make_tensor_expr({"g", "i_1", "i_2", "a_1", "a_2"});
  //  auto expr = std::make_shared<Sum>(Sum({g, t}));
  //  auto tree = builder.build_tree(expr);
  //  tree->visit(visitor);
  //  std::wcout << "digraph G {\n";
  //  std::wcout << tree->to_digraph();
  //  std::wcout << "}\n";
  //}
  // TA::finalize();
}
