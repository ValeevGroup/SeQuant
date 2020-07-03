#include "catch.hpp"

#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/domain/evaluate/eval_tree.hpp>
#include <SeQuant/domain/factorize/factorizer.hpp>
#include <iomanip>
#include <ios>

using namespace sequant;
using namespace sequant::factorize;

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

TEST_CASE("ADJACENCY_MATRIX", "[factorizer]") {
  auto prod = std::make_shared<Product>(
      Product{make_tensor_expr({"g", "i_3", "i_4", "i_1", "a_3"}),
              make_tensor_expr({"t", "a_1", "i_3"}),
              make_tensor_expr({"t", "a_2", "i_4"}),
              make_tensor_expr({"t", "a_3", "i_2"})});
  auto adjMat = AdjacencyMatrix(
      prod, {Index{L"i_1"}, Index{L"i_2"}, Index{L"a_1"}, Index{L"a_2"}});

  size_t pos0 = 0, pos1 = 1, pos2 = 2, pos3 = 3;
  REQUIRE(adjMat.are_connected(pos0, pos1));
  REQUIRE(adjMat.are_connected(pos0, pos2));
  REQUIRE(adjMat.are_connected(pos0, pos3));

  REQUIRE(!adjMat.are_connected(pos1, pos2));
  REQUIRE(!adjMat.are_connected(pos1, pos3));
  REQUIRE(!adjMat.are_connected(pos2, pos3));
}

TEST_CASE("TENSOR_EXISTS", "[factorizer]") {
  auto G_oovv = make_tensor_expr({"g", "i_1", "i_2", "a_1", "a_2"});
  auto T_ov1 = make_tensor_expr({"t", "i_3", "a_2"});
  auto prod1 = std::make_shared<Product>(Product{G_oovv, T_ov1});

  auto T_ov2 = make_tensor_expr({"t", "i_4", "a_3"});
  auto F_oo = make_tensor_expr({"f", "i_1", "i_2"});
  auto prod2 = std::make_shared<Product>(Product{T_ov2, F_oo});

  auto prod = std::make_shared<Product>(Product{prod1, prod2});

  REQUIRE(tensor_exists(prod1, G_oovv));
  REQUIRE(tensor_exists(prod2, F_oo));
  REQUIRE(!tensor_exists(prod1, F_oo));
  REQUIRE(!tensor_exists(prod2, G_oovv));
  REQUIRE(tensor_exists(prod, T_ov1));
  REQUIRE(tensor_exists(prod, T_ov2));
}

TEST_CASE("FACTORIZER", "[factorizer]") {
  using std::endl;
  using std::wcout;

  auto print_adj_mat = [](const auto& mat, bool color = false) {
    for (auto ii = 0; ii < mat.num_verts(); ++ii) {
      for (auto jj = 0; jj < mat.num_verts(); ++jj) {
        if (color)
          std::wcout << std::setw(27) << mat.color(ii, jj);
        else
          std::wcout << mat.are_connected(ii, jj);
        std::wcout << "  ";
      }
      std::wcout << "\n";
    }
  };

  SECTION("Single common subnet") {
    auto prod1 = std::make_shared<Product>();
    prod1->append(1, make_tensor_expr({"g", "i_3", "i_4", "i_1", "a_3"}));
    prod1->append(1, make_tensor_expr({"t", "a_3", "i_2"}));
    prod1->append(1, make_tensor_expr({"t", "a_1", "a_2", "i_3", "i_4"}));

    auto prod2 = std::make_shared<Product>();
    prod2->append(1, make_tensor_expr({"g", "i_3", "i_4", "i_1", "a_3"}));
    prod2->append(1, make_tensor_expr({"t", "a_1", "i_3"}));
    prod2->append(1, make_tensor_expr({"t", "a_2", "i_4"}));
    prod2->append(1, make_tensor_expr({"t", "a_3", "i_2"}));

    std::wcout << "Factorization of single common subnet\n";

    auto [fact1, fact2] = factorize_pair(prod1, prod2);

    REQUIRE(evaluate::EvalTree(fact1->at(0)).hash_value() ==
            evaluate::EvalTree(fact2->at(0)).hash_value());

    prod1 = std::make_shared<Product>();
    prod1->append(1, make_tensor_expr({"g", "i_3", "i_4", "a_3", "a_4"}));
    prod1->append(1, make_tensor_expr({"t", "a_1", "a_2", "i_1", "i_3"}));
    prod1->append(1, make_tensor_expr({"t", "a_3", "a_4", "i_2", "i_4"}));

    prod2 = std::make_shared<Product>();
    prod2->append(1, make_tensor_expr({"g", "i_3", "i_4", "a_3", "a_4"}));
    prod2->append(1, make_tensor_expr({"t", "a_3", "i_3"}));
    prod2->append(1, make_tensor_expr({"t", "a_4", "i_1"}));
    prod2->append(1, make_tensor_expr({"t", "a_1", "a_2", "i_2", "i_4"}));

    auto [fact3, fact4] = factorize_pair(prod1, prod2);

    auto term1 = fact3->at(0)->clone();
    auto term2 = fact4->at(0)->clone();

    auto ext_indices1 = factorize::target_indices(prod1);
    auto ext_indices2 = factorize::target_indices(prod2);

    auto tn1 = TensorNetwork(*term1);
    auto tn2 = TensorNetwork(*term2);

    tn1.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), false,
                     &ext_indices1);
    tn2.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), false,
                     &ext_indices2);

    auto canon_term1 = std::make_shared<Product>();
    auto canon_term2 = std::make_shared<Product>();

    for (const auto& t : tn1.tensors())
      canon_term1->append(1, std::dynamic_pointer_cast<Expr>(t));
    for (const auto& t : tn2.tensors())
      canon_term2->append(1, std::dynamic_pointer_cast<Expr>(t));

    REQUIRE(evaluate::EvalTree(canon_term1).hash_value() ==
            evaluate::EvalTree(canon_term2).hash_value());
  }
}
