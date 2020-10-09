#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/factorize/eval_sequence.hpp>
#include <SeQuant/domain/factorize/ops_count.hpp>
#include <cmath>

#include "catch.hpp"

TEST_CASE("TEST_OPS_COUNT", "[ops_count]") {
  using namespace sequant;

  auto g_oovv = ex<Tensor>(
      Tensor{L"g", {L"i_3", L"i_4"}, {L"a_3", L"a_4"}, Symmetry::antisymm});
  auto t_oovv_1 = ex<Tensor>(
      Tensor{L"t", {L"a_1", L"a_2"}, {L"i_3", L"i_4"}, Symmetry::antisymm});
  auto t_oovv_2 = ex<Tensor>(
      Tensor{L"t", {L"a_3", L"a_4"}, {L"i_1", L"i_2"}, Symmetry::antisymm});

  // ops count will be computed on the following product.
  auto prod = ex<Product>(Product{g_oovv, t_oovv_1, t_oovv_2});

  using namespace sequant::factorize;
  auto tree1 = rooted_tree{0};
  tree1.children.push_back(rooted_tree{1});
  tree1.children.push_back(rooted_tree{2});
  // tree1 == ((0 1) 2)
  // manually derived formula
  auto ops_manual_tree1 = [](size_t nocc,
                             size_t nvirt) -> OpsCalcResult::ops_type {
    return 4 * std::pow(nocc, 2) * std::pow(nvirt, 4);
  };

  auto tree2 = rooted_tree{0};
  tree2.children.push_back(rooted_tree{2});
  tree2.children.push_back(rooted_tree{1});
  // tree2 == ((0 2) 1)
  // manually derived formula
  auto ops_manual_tree2 = [](size_t nocc,
                             size_t nvirt) -> OpsCalcResult::ops_type {
    return 4 * std::pow(nocc, 4) * std::pow(nvirt, 2);
  };

  auto tree3 = rooted_tree{1};
  tree3.children.push_back(rooted_tree{2});
  tree3.children.push_back(rooted_tree{0});
  // tree3 == ((1 2) 0)
  // manually derived formula
  auto ops_manual_tree3 = [](size_t nocc,
                             size_t nvirt) -> OpsCalcResult::ops_type {
    return 4 * std::pow(nocc, 4) * std::pow(nvirt, 4);
  };

  REQUIRE(ops_manual_tree1(2, 3) == ops_count(prod, tree1, 2, 3).flops);
  REQUIRE(ops_manual_tree1(2, 2) == ops_count(prod, tree1, 2, 2).flops);
  REQUIRE(ops_manual_tree1(3, 2) == ops_count(prod, tree1, 3, 2).flops);

  REQUIRE(ops_manual_tree2(2, 3) == ops_count(prod, tree2, 2, 3).flops);
  REQUIRE(ops_manual_tree2(2, 2) == ops_count(prod, tree2, 2, 2).flops);
  REQUIRE(ops_manual_tree2(3, 2) == ops_count(prod, tree2, 3, 2).flops);

  REQUIRE(ops_manual_tree3(2, 3) == ops_count(prod, tree3, 2, 3).flops);
  REQUIRE(ops_manual_tree3(2, 2) == ops_count(prod, tree3, 2, 2).flops);
  REQUIRE(ops_manual_tree3(3, 2) == ops_count(prod, tree3, 3, 2).flops);
}
