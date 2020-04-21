#include "catch.hpp"

#include <SeQuant/domain/evaluate/eval_tree.hpp>
#include <SeQuant/domain/evaluate/factorizer.hpp>

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

TEST_CASE("TN_FACTORIZE_BASED_ON_EVAL_TREE", "[tensor_network]") {
  SECTION("Testing largest common sub-expressions: Product") {
    // forming two tensor products
    auto prodA = std::make_shared<Product>(
        Product({make_tensor_expr({"t", "i_1", "i_2", "a_1", "a_2"}),
                 make_tensor_expr({"g", "i_1", "i_3", "a_1", "a_3"}),
                 make_tensor_expr({"f", "i_2", "a_2"})}));

    auto prodB = std::make_shared<Product>(
        Product({make_tensor_expr({"t", "i_4", "i_6", "a_4", "a_6"}),
                 make_tensor_expr({"g", "i_4", "i_8", "a_4", "a_8"}),
                 make_tensor_expr({"f", "i_8", "a_8"})}));
    // Note how each tensor in prodA is equivalent to a tensor in prodB
    // (equivalent in the sense that both represent the same data tensor)
    // however, prodA and prodB are not equivalent products
    // that is beacuse even though, the nodes of the two tensor networks (prodA
    // and prodB) are the same, their edges differ.
    //
    // The 't' and 'g' tensors are labelled in such a way that the network of
    // t<-->g in both products have the same edges, and thus t<-->g is the
    // subtensor network common to prodA and prodB.

    auto [subfactorA, subfactorB] = largest_common_subnet(prodA, prodB);
    REQUIRE(subfactorA == container::svector<size_t>{0, 1});
    REQUIRE(subfactorB == container::svector<size_t>{0, 1});

    // prodC is the same as prodB except the position of 'f' and 'g' tensors are
    // swapped this shows that absolute order of the constiuents that make up
    // the common subnet in the pair of tensor networks
    // doesn't need to be the same
    auto prodC = std::make_shared<Product>(
        Product({make_tensor_expr({"t", "i_4", "i_6", "a_4", "a_6"}),
                 make_tensor_expr({"f", "i_8", "a_8"}),
                 make_tensor_expr({"g", "i_4", "i_8", "a_4", "a_8"})}));
    auto [subfactorX, subfactorY] = largest_common_subnet(prodA, prodC);
    REQUIRE(subfactorX == container::svector<size_t>{0, 1});
    REQUIRE(subfactorY == container::svector<size_t>{0, 2});

    // lets completely reverse the order of the tensors in prodB
    // shows that the relative order of the constiuents of the common subnet
    // doesn't matter
    prodC = std::make_shared<Product>(
        Product({make_tensor_expr({"f", "i_8", "a_8"}),
                 make_tensor_expr({"g", "i_4", "i_8", "a_4", "a_8"}),
                 make_tensor_expr({"t", "i_4", "i_6", "a_4", "a_6"})}));

    auto [subfactorXX, subfactorYY] = largest_common_subnet(prodA, prodC);
    REQUIRE(subfactorXX == container::svector<size_t>{0, 1});
    REQUIRE(subfactorYY == container::svector<size_t>{2, 1});

    // no common subnet
    prodA = std::make_shared<Product>(
        Product({make_tensor_expr({"f", "i_8", "a_8"})}));
    prodB = std::make_shared<Product>(
        Product({make_tensor_expr({"g", "i_4", "i_8", "a_4", "a_8"})}));

    auto [subfactorA_null, subfactorB_null] =
        largest_common_subnet(prodA, prodB);

    REQUIRE(subfactorA_null.empty());
    REQUIRE(subfactorB_null.empty());
  }

  SECTION("Testing largest common sub-expressions: Sum") {
    // forming two sums
    auto sumA = std::make_shared<Sum>(
        Sum({make_tensor_expr({"t", "i_1", "i_2", "a_1", "a_2"}),
             make_tensor_expr({"g", "i_1", "i_2", "a_1", "a_4"})}));

    auto sumB = std::make_shared<Sum>(
        Sum({make_tensor_expr({"t", "i_4", "i_6", "a_4", "a_6"}),
             make_tensor_expr({"g", "i_4", "i_6", "a_4", "a_6"})}));

    auto [summandA, summandB] = largest_common_subnet(sumA, sumB);

    REQUIRE(summandA == container::svector<size_t>{0, 1});
    REQUIRE(summandB == container::svector<size_t>{0, 1});

    // forming sums whose summands are products
    auto prod1 = std::make_shared<Product>(
        Product({make_tensor_expr({"t", "i_1", "i_2", "a_1", "a_2"}),
                 make_tensor_expr({"g", "i_1", "i_3", "a_1", "a_3"}),
                 make_tensor_expr({"f", "i_2", "a_2"})}));

    auto prod2 = std::make_shared<Product>(
        Product({make_tensor_expr({"f", "i_3", "a_3"})}));

    auto prod3 = std::make_shared<Product>(
        Product({make_tensor_expr({"t", "i_2", "i_3", "a_2", "a_3"}),
                 make_tensor_expr({"f", "i_2", "a_2"})}));

    auto prod4 = std::make_shared<Product>(
        Product({make_tensor_expr({"g", "i_2", "i_3", "a_2", "a_3"}),
                 make_tensor_expr({"t", "i_2", "a_2"})}));

    sumA = std::make_shared<Sum>(Sum({prod1, prod2, prod3}));
    sumB = std::make_shared<Sum>(Sum({prod1, prod2, prod4}));

    auto expectedSubNet = std::make_shared<Sum>(Sum({prod1, prod2}));

    auto [summandAA, summandBB] = largest_common_subnet(sumA, sumB);

    REQUIRE(summandAA == container::svector<size_t>{0, 1});
    REQUIRE(summandBB == container::svector<size_t>{0, 1});

    // no common subnet
    sumA = std::make_shared<Sum>(Sum({prod1, prod2}));
    sumB = std::make_shared<Sum>(Sum({prod3, prod4}));

    auto [summandAnull, summandBnull] = largest_common_subnet(sumA, sumB);

    REQUIRE(summandAnull.empty());
    REQUIRE(summandBnull.empty());
  }

  SECTION("Testing largest common subfactor: Fails using EvalTree hash") {
    // forming two tensor products
    auto prodA = std::make_shared<Product>(
        Product({make_tensor_expr({"t", "i_1", "i_2", "a_1", "a_2"}),
                 make_tensor_expr({"g", "i_3", "i_4", "a_2", "a_4"}),
                 make_tensor_expr({"t", "i_3", "i_4", "a_3", "a_4"})}));

    auto prodB = std::make_shared<Product>(
        Product({make_tensor_expr({"t", "i_3", "i_4", "a_3", "a_4"}),
                 make_tensor_expr({"g", "i_3", "i_4", "a_2", "a_4"}),
                 make_tensor_expr({"t", "i_1", "i_2", "a_1", "a_2"})}));

    auto [subfactorA, subfactorB] = largest_common_subnet(prodA, prodB);
    // REQUIRE(subfactorA == container::svector<size_t>{0, 1, 2});
    // REQUIRE(subfactorB == container::svector<size_t>{2, 1, 0});
  }
}
