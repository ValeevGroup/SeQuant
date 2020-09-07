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

TEST_CASE("HASH_VALUES_TEST", "[eval_tree]") {
  auto f_oo = make_tensor_expr({"f", "i_1", "i_2"});
  auto f_ov = make_tensor_expr({"f", "i_1", "a_1"});
  auto f_vv = make_tensor_expr({"f", "a_1", "a_2"});

  auto g_oooo = make_tensor_expr({"g", "i_1", "i_2", "i_3", "i_4"});
  auto g_ooov = make_tensor_expr({"g", "i_1", "i_2", "i_3", "a_1"});
  auto g_oovv = make_tensor_expr({"g", "i_1", "i_2", "a_1", "a_2"});
  auto g_ovvv = make_tensor_expr({"g", "i_1", "a_1", "a_2", "a_3"});
  auto g_vvvv = make_tensor_expr({"g", "a_1", "a_2", "a_3", "a_4"});

  auto t_ov = make_tensor_expr({"t", "i_1", "a_1"});
  auto t_oovv = make_tensor_expr({"t", "i_1", "i_2", "a_1", "a_2"});

  auto A_ov = make_tensor_expr({"A", "i_1", "a_1"});
  auto A_oovv = make_tensor_expr({"A", "i_1", "i_2", "a_1", "a_2"});

  container::svector<decltype(f_oo)> tensor_exprs{
      f_oo,   f_ov,   f_vv, g_oooo, g_ooov, g_oovv,
      g_ovvv, g_vvvv, t_ov, t_oovv, A_ov,   A_oovv};

  SECTION("Allowed swapped vs not allowed swapped") {
    auto g_oovv_perm = make_tensor_expr({"g", "a_1", "a_2", "i_1", "i_2"});

    REQUIRE(EvalTree(g_oovv, true).hash_value() ==
            EvalTree(g_oovv_perm, true).hash_value());
    REQUIRE(EvalTree(g_oovv, false).hash_value() !=
            EvalTree(g_oovv_perm, false).hash_value());
  }

  SECTION("Unique hash values") {
    for (auto ii = 0; ii < tensor_exprs.size() - 1; ++ii) {
      REQUIRE(EvalTree(tensor_exprs.at(ii), true).hash_value() ==
              EvalTree(tensor_exprs.at(ii), false).hash_value());
      // param for swapping or unswapping braket indices
      // doesn't matter for above tensors as they are already
      // written in a canonicalized form

      for (auto jj = ii + 1; jj < tensor_exprs.size(); ++jj) {
        // confirm each tensor is uniquely identified by
        // its hash value
        REQUIRE(EvalTree(tensor_exprs.at(ii)).hash_value() !=
                EvalTree(tensor_exprs.at(jj)).hash_value());
      }
    }
  }

  SECTION("Order agnostic cases") {
    // iter by omitting the last two: A_ov and A_oovv tensors
    for (auto ii = 0; ii < tensor_exprs.size() - 3; ++ii) {
      for (auto jj = ii + 1; jj < tensor_exprs.size() - 2; ++jj) {
        auto prodA = std::make_shared<Product>(
            Product{tensor_exprs.at(ii), tensor_exprs.at(jj)});
        auto prodB = std::make_shared<Product>(
            Product{tensor_exprs.at(jj), tensor_exprs.at(ii)});
        // T1.T2 == T2.T1
        REQUIRE(EvalTree(prodA).hash_value() == EvalTree(prodB).hash_value());
      }
    }

    //
    // T1.(T2.T3) = (T2.T3).T1
    //
    auto prodA = std::make_shared<Product>(Product{g_oooo});
    auto prodB = std::make_shared<Product>(Product{t_oovv, g_vvvv});

    // append without flattening
    auto prodX = std::make_shared<Product>(*prodA);
    prodX->append(prodB);  // prodX: f_ov *  (g_oovv, g_vvvv)

    auto prodY = std::make_shared<Product>(*prodB);
    prodY->append(prodA);  // prodY: (g_oovv, g_vvvv) * f_ov

    REQUIRE(EvalTree(prodX).hash_value() == EvalTree(prodY).hash_value());
  }
}

TEST_CASE("VISITOR_TEST", "[eval_tree]") {
  auto g_oooo = make_tensor_expr({"g", "i_1", "i_2", "i_3", "i_4"});
  auto g_vvvv = make_tensor_expr({"g", "a_1", "a_2", "a_3", "a_4"});
  auto t_oovv = make_tensor_expr({"t", "i_1", "i_2", "a_1", "a_2"});
  auto prodA = std::make_shared<Product>(Product{g_oooo, t_oovv, g_vvvv});

  container::svector<decltype(g_oooo->bra().begin()->label())>
      reaped_idx_labels;

  auto visitor = [&reaped_idx_labels](const auto& node) {
    for (const auto& idx : node->indices())
      reaped_idx_labels.push_back(idx.label());
  };

  auto tree = EvalTree(prodA);
  // tree:
  //                         a_3 a_4
  //                        Q
  //                      /  i_3 i_4
  //                     /          \
  //                    /            \
  //                   /              \
  //                  /                \
  //            a_1 a_2                  a_3 a_4
  //           Q                        G
  //         /  i_3 i_4                  a_1 a_2
  //        /           \
  //       /             \
  //   i_3 i_4            \  a_1 a_2
  //  G                     T
  //   i_1 i_2               i_1 i_2
  //
  //
  tree.visit(visitor);

  decltype(reaped_idx_labels) expected_idx_labels{
      L"i_3", L"i_4", L"a_3", L"a_4", L"i_3", L"i_4", L"a_1",
      L"a_2", L"i_1", L"i_2", L"i_3", L"i_4", L"i_1", L"i_2",
      L"a_1", L"a_2", L"a_1", L"a_2", L"a_3", L"a_4"};

  REQUIRE(size(reaped_idx_labels) == size(expected_idx_labels));
  for (size_t c = 0; c != size(reaped_idx_labels); ++c) {
    REQUIRE(reaped_idx_labels[c] == *(begin(expected_idx_labels) + c));
  }
}

TEST_CASE("UNSWAP_BK_LABELS_TEST", "[eval_tree]") {
  auto G_vvoo = make_tensor_expr({"g", "a_1", "a_2", "i_1", "i_2"});
  container::svector<decltype(G_vvoo->bra().begin()->label())>
      reaped_idx_labels;

  // reap index labels from leaf nodes
  // and store them into reaped_idx_labels vector
  auto reaper = [&reaped_idx_labels](const auto& node) {
    if (node->is_leaf()) {
      for (const auto& idx : node->indices())
        reaped_idx_labels.push_back(idx.label());
    }
  };

  auto tree = EvalTree(G_vvoo, true);

  // G_vvoo was used with swap allowed
  // so, the bra-ket labels will be swapped as they should be for G_vvoo

  reaped_idx_labels.clear();
  tree.visit(reaper);

  {
    auto expected_idx_labels = {L"i_1", L"i_2", L"a_1", L"a_2"};
    REQUIRE(size(reaped_idx_labels) == size(expected_idx_labels));
    for (size_t c = 0; c != size(reaped_idx_labels); ++c) {
      REQUIRE(reaped_idx_labels[c] == *(begin(expected_idx_labels) + c));
    }
  }

  // let's unswap
  tree.swap_labels(G_vvoo);

  reaped_idx_labels.clear();
  tree.visit(reaper);

  {
    auto expected_idx_labels = {L"a_1", L"a_2", L"i_1", L"i_2"};
    REQUIRE(size(reaped_idx_labels) == size(expected_idx_labels));
    for (size_t c = 0; c != size(reaped_idx_labels); ++c) {
      REQUIRE(reaped_idx_labels[c] == *(begin(expected_idx_labels) + c));
    }
  }
}

TEST_CASE("OPS_COUNT TESTS", "[eval_tree]") {
  const size_t nocc = 10;
  const size_t nvirt = 20;

  container::map<IndexSpace::TypeAttr, size_t> space_size;
  space_size.insert(
      decltype(space_size)::value_type(IndexSpace::active_occupied, nocc));
  space_size.insert(
      decltype(space_size)::value_type(IndexSpace::active_unoccupied, nvirt));

  SECTION("Testing operation counts for product") {
    auto t = make_tensor_expr({"t", "i_1", "a_1"});
    auto g = make_tensor_expr({"g", "i_1", "i_2", "a_1", "a_2"});
    auto expr = std::make_shared<Product>(Product({t, g}));
    auto tree = EvalTree(expr);

    // by looking at 't' and 'g' tensors
    // we see the ops count should be
    // (num of active_occupied squared) * (num of active_unoccupied squared)

    REQUIRE(tree.ops_count(space_size) == nocc * nocc * nvirt * nvirt);
  }

  SECTION("Testing operation counts for sum") {
    auto t = make_tensor_expr({"t", "i_1", "a_1"});
    auto f = make_tensor_expr({"f", "i_2", "a_2"});
    auto g = make_tensor_expr({"g", "i_1", "i_2", "a_1", "a_2"});

    auto left_sumand = std::make_shared<Product>(Product{t, f});
    auto& right_sumand = g;

    auto expr = std::make_shared<Sum>(Sum{left_sumand, right_sumand});

    auto tree = EvalTree(expr);

    REQUIRE(tree.ops_count(space_size) == nocc * nocc * nvirt * nvirt);
  }
}

TEST_CASE("Digraph generation", "[eval_tree]") {
  auto prod1 = std::make_shared<Product>(
      Product(1 / 16.0, {make_tensor_expr({"A", "i_1", "i_2", "a_1", "a_2"}),
                         make_tensor_expr({"g", "i_3", "i_4", "a_3", "a_4"}),
                         make_tensor_expr({"t", "a_1", "a_2", "i_3", "i_4"}),
                         make_tensor_expr({"t", "a_3", "a_4", "i_1", "i_2"})}));

  auto prod2 = std::make_shared<Product>(
      Product(-1., {make_tensor_expr({"A", "i_1", "i_2", "a_1", "a_2"}),
                    make_tensor_expr({"g", "i_3", "a_1", "a_3", "a_4"}),
                    make_tensor_expr({"t", "a_3", "i_1"}),
                    make_tensor_expr({"t", "a_2", "a_4", "i_2", "i_3"})}));

  auto sum = std::make_shared<Sum>(Sum{prod1, prod2});

  EvalTree(sum).digraph(std::wcout);
}

TEST_CASE("EVALUATIONS TESTS", "[eval_tree]") {
  using DTensorType = TA::TArrayD;
  using ContextMapType =
      sequant::container::map<HashType, std::shared_ptr<DTensorType>>;

  const size_t nocc = 10;
  const size_t nvirt = 20;

  TA::TiledRange tr_ov{{0, nocc}, {0, nvirt}};
  TA::TiledRange tr_oovv{{0, nocc}, {0, nocc}, {0, nvirt}, {0, nvirt}};

  auto& world = TA::get_default_world();
  auto tnsr_T_ov = std::make_shared<DTensorType>(world, tr_ov);
  auto tnsr_T_oovv = std::make_shared<DTensorType>(world, tr_oovv);
  auto tnsr_G_oovv = std::make_shared<DTensorType>(world, tr_oovv);

  tnsr_T_ov->fill_random();
  tnsr_T_oovv->fill_random();
  tnsr_G_oovv->fill_random();

  SECTION("Testing Sum type evaluations") {
    //
    auto t = make_tensor_expr({"t", "i_1", "i_2", "a_1", "a_2"});
    auto g = make_tensor_expr({"g", "i_1", "i_2", "a_1", "a_2"});

    DTensorType manual_result;
    manual_result("i, j, a, b") =
        (*tnsr_T_oovv)("i, j, a, b") + (*tnsr_G_oovv)("i, j, a, b");

    ContextMapType context;
    context.insert(
        ContextMapType::value_type(EvalTree(t).hash_value(), tnsr_T_oovv));
    context.insert(
        ContextMapType::value_type(EvalTree(g).hash_value(), tnsr_G_oovv));

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
        ContextMapType::value_type(EvalTree(g).hash_value(), tnsr_G_oovv));

    manual_result("i, j, a, b") =
        (*tnsr_T_oovv)("i, j, a, b") + (*tnsr_G_oovv)("i, j, b, a");
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
    manual_result("j,b") = (*tnsr_T_ov)("i,a") * (*tnsr_G_oovv)("i,j,a,b");

    ContextMapType context;
    context.insert(
        ContextMapType::value_type(EvalTree(t).hash_value(), tnsr_T_ov));
    context.insert(
        ContextMapType::value_type(EvalTree(g).hash_value(), tnsr_G_oovv));

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
        ContextMapType::value_type(EvalTree(t).hash_value(), tnsr_T_oovv));

    DTensorType manual_result;
    manual_result("i,j,a,b") =
        (*tnsr_T_oovv)("i,j,a,b") - (*tnsr_T_oovv)("i,j,b,a") +
        (*tnsr_T_oovv)("j,i,b,a") - (*tnsr_T_oovv)("j,i,a,b");

    auto expr = std::make_shared<Product>(Product({A, t}));
    auto tree = EvalTree(expr);
    auto eval_result = tree.evaluate(context);

    auto manual_norm =
        std::sqrt(manual_result("0,1,2,3").dot(manual_result("0,1,2,3")));

    auto eval_norm =
        std::sqrt(eval_result("0,1,2,3").dot(eval_result("0,1,2,3")));

    REQUIRE(manual_norm == Approx(eval_norm));
  }

  SECTION("Testing symmetrization evaluation") {
    TA::TiledRange tr_ooovvv{{0, nocc},  {0, nocc},  {0, nocc},
                             {0, nvirt}, {0, nvirt}, {0, nvirt}};
    auto tnsr_T_ooovvv = std::make_shared<DTensorType>(world, tr_ooovvv);
    tnsr_T_ooovvv->fill_random();

    auto t_ov = make_tensor_expr({"t", "i_1", "a_1"});
    auto t_oovv = make_tensor_expr({"t", "i_1", "i_2", "a_1", "a_2"});
    auto t_ooovvv =
        make_tensor_expr({"t", "i_1", "i_2", "i_3", "a_1", "a_2", "a_3"});

    auto S_ov = make_tensor_expr({"S", "i_1", "i_2"});
    auto S_oovv = make_tensor_expr({"S", "i_1", "i_2", "a_1", "a_2"});
    auto S_ooovvv =
        make_tensor_expr({"S", "i_1", "i_2", "i_3", "a_1", "a_2", "a_3"});

    ContextMapType context;
    context.insert(
        ContextMapType::value_type(EvalTree(t_ov).hash_value(), tnsr_T_ov));
    context.insert(
        ContextMapType::value_type(EvalTree(t_oovv).hash_value(), tnsr_T_oovv));
    context.insert(  //
        ContextMapType::value_type(EvalTree(t_ooovvv).hash_value(),
                                   tnsr_T_ooovvv));

    DTensorType manual_result, eval_result;
    auto expr = std::make_shared<Product>(Product({S_ov, t_ov}));
    eval_result = EvalTree(expr).evaluate(context);
    manual_result("0,1") = (*tnsr_T_ov)("0,1");
    auto eval_norm = std::sqrt(eval_result("0,1").dot(eval_result("0,1")));
    auto manual_norm =
        std::sqrt(manual_result("0,1").dot(manual_result("0,1")));
    REQUIRE(manual_norm == Approx(eval_norm));

    expr = std::make_shared<Product>(Product({S_oovv, t_oovv}));
    eval_result = EvalTree(expr).evaluate(context);
    manual_result("0,1,2,3") =
        (*tnsr_T_oovv)("0,1,2,3") + (*tnsr_T_oovv)("1,0,3,2");
    eval_norm = std::sqrt(eval_result("0,1,2,3").dot(eval_result("0,1,2,3")));
    manual_norm =
        std::sqrt(manual_result("0,1,2,3").dot(manual_result("0,1,2,3")));
    REQUIRE(manual_norm == Approx(eval_norm));

    expr = std::make_shared<Product>(Product({S_ooovvv, t_ooovvv}));
    expr->scale(1.0 / 12.0);
    eval_result = EvalTree(expr).evaluate(context);
    manual_result("0,1,2,3,4,5") =
        (*tnsr_T_ooovvv)("0,2,1,3,5,4") + (*tnsr_T_ooovvv)("0,1,2,3,4,5") +
        (*tnsr_T_ooovvv)("1,0,2,4,3,5") + (*tnsr_T_ooovvv)("1,2,0,4,5,3") +
        (*tnsr_T_ooovvv)("2,0,1,5,3,4") + (*tnsr_T_ooovvv)("2,1,0,5,4,3");
    manual_result("0,1,2,3,4,5") = 1.0 / 12.0 * manual_result("0,1,2,3,4,5");
    eval_norm =
        std::sqrt(eval_result("0,1,2,3,4,5").dot(eval_result("0,1,2,3,4,5")));
    manual_norm = std::sqrt(
        manual_result("0,1,2,3,4,5").dot(manual_result("0,1,2,3,4,5")));
    REQUIRE(manual_norm == Approx(eval_norm));
  }

  SECTION("Testing missing data tensor") {
    auto seq_tensor_bad = make_tensor_expr({"t", "a_1", "i_1", "a_2", "a_3"});
    auto seq_tensor_good = make_tensor_expr({"t", "i_1", "a_1", "a_2", "a_3"});
    ContextMapType context;
    context.insert(ContextMapType::value_type(
        EvalTree(seq_tensor_good).hash_value(), tnsr_T_ov));

    auto expr = seq_tensor_bad;
    auto tree = EvalTree(expr);

    REQUIRE_THROWS_AS(tree.evaluate(context), std::logic_error);
  }

  SECTION("Testing mixed evaluations") {
    auto prod1 = std::make_shared<Product>(Product(
        1 / 16.0, {make_tensor_expr({"A", "i_1", "i_2", "a_1", "a_2"}),
                   make_tensor_expr({"g", "i_3", "i_4", "a_3", "a_4"}),
                   make_tensor_expr({"t", "a_1", "a_2", "i_3", "i_4"}),
                   make_tensor_expr({"t", "a_3", "a_4", "i_1", "i_2"})}));

    auto prod2 = std::make_shared<Product>(
        Product(-1., {make_tensor_expr({"A", "i_1", "i_2", "a_1", "a_2"}),
                      make_tensor_expr({"g", "i_3", "a_1", "a_3", "a_4"}),
                      make_tensor_expr({"t", "a_3", "i_1"}),
                      make_tensor_expr({"t", "a_2", "a_4", "i_2", "i_3"})}));

    auto sum = std::make_shared<Sum>(Sum{prod1, prod2});

    // g_oovv, t_oovv, t_ov are already declared
    TA::TiledRange tr_ovvv{{0, nocc}, {0, nvirt}, {0, nvirt}, {0, nvirt}};
    auto G_ovvv = std::make_shared<DTensorType>(world, tr_ovvv);
    G_ovvv->fill_random();

    // building context for evaluation
    ContextMapType context;

    auto& g_oovv_seq = prod1->as<Product>().factor(1);
    auto& t_oovv_seq = prod1->as<Product>().factor(2);

    auto& g_ovvv_seq = prod2->as<Product>().factor(1);
    auto& t_ov_seq = prod2->as<Product>().factor(2);

    context.insert(ContextMapType::value_type(EvalTree(g_oovv_seq).hash_value(),
                                              tnsr_G_oovv));

    context.insert(
        ContextMapType::value_type(EvalTree(g_ovvv_seq).hash_value(), G_ovvv));

    context.insert(ContextMapType::value_type(EvalTree(t_oovv_seq).hash_value(),
                                              tnsr_T_oovv));

    context.insert(
        ContextMapType::value_type(EvalTree(t_ov_seq).hash_value(), tnsr_T_ov));

    auto visitor = [](const auto& node) {
      if (node->is_leaf())
        std::wcout << "leaf node..\n";
      else
        std::wcout << "internal node..\n";

      for (const auto& idx : node->indices()) std::wcout << idx.label() << " ";
      std::wcout << std::endl;
    };

    auto tree = EvalTree(sum);

    // tree.visit(visitor);
    auto eval_result = tree.evaluate(context);
    auto eval_norm =
        std::sqrt(eval_result("0,1,2,3").dot(eval_result("0,1,2,3")));

    // manual evaluation
    DTensorType before_antisym;
    before_antisym("i_1, i_2, a_1, a_2") =
        1 / 16. * (*tnsr_G_oovv)("i_3, i_4, a_3, a_4") *
            (*tnsr_T_oovv)("i_3, i_4, a_1, a_2") *
            (*tnsr_T_oovv)("i_1, i_2, a_3, a_4") -
        1 * (*G_ovvv)("i_3, a_1, a_3, a_4") * (*tnsr_T_ov)("i_1, a_3") *
            (*tnsr_T_oovv)("i_2, i_3, a_2, a_4");
    DTensorType manual_result;
    manual_result("i,j,a,b") =
        1 * before_antisym("i,j,a,b") - 1 * before_antisym("i,j,b,a") +
        1 * before_antisym("j,i,b,a") - 1 * before_antisym("j,i,a,b");

    auto manual_norm =
        std::sqrt(manual_result("0,1,2,3").dot(manual_result("0,1,2,3")));

    REQUIRE(eval_norm == Approx(manual_norm));
  }
}
