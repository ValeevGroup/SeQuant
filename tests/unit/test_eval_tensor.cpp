//
// created by Bimal Gaudel on Mar 11, 2020
//

#include "SeQuant/core/expr.hpp"
#include "SeQuant/core/tensor.hpp"
#include "catch.hpp"

#include <iostream>
#include <memory>

#include "../../src/SeQuant/domain/evaluate/eval_tensor.hpp"
#include "../../src/SeQuant/domain/evaluate/eval_tensor_builder.hpp"

TEST_CASE("EVAL_TENSOR_TESTS", "[eval_tensor]") {
  using namespace sequant::evaluate;
  SECTION("intermediate construction") {
    EvalTensorIntermediate evt_imed;
    REQUIRE(!evt_imed.is_leaf());
    REQUIRE(evt_imed.get_operation() == Operation::INVALID);
    REQUIRE(evt_imed.get_left_tensor() == nullptr);
    REQUIRE(evt_imed.get_right_tensor() == nullptr);
  }
  SECTION("leaf construction") {
    EvalTensorLeaf evt_leaf;
    REQUIRE(evt_leaf.is_leaf());
  }
}

TEST_CASE("EVAL_TENSOR_BUILDER_TESTS", "[eval_tensor_builder]") {
  using namespace sequant::evaluate;
  using sequant::Product;
  using sequant::Tensor;
  auto builder = EvalTensorBuilder();
  SECTION("construction") { REQUIRE(builder.get_eval_tree() == nullptr); }
  SECTION("Testing eval_tensor_builder methods") {
    auto t1 = std::make_shared<Tensor>(Tensor{L"t", {L"i_1"}, {L"a_1"}});
    auto t2 = std::make_shared<Tensor>(Tensor{L"t", {L"i_2"}, {L"a_2"}});
    auto g1 = std::make_shared<Tensor>(
        Tensor{L"g", {L"i_2", L"i_3"}, {L"a_2", L"a_4"}});
    auto prod = std::make_shared<Product>(Product({t1, t2, g1}));
    std::wcout << "Formed product \n";
    std::wcout << prod->to_latex() << std::endl;
    builder.build_from_product(prod);
    // std::wcout << "printing the graph..\n";
    // std::wcout << builder.get_eval_tree()->to_digraph() << std::endl;
  }
}
