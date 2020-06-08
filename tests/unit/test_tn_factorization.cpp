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
      Product{make_tensor_expr({"g", "i_1", "i_2", "a_1", "a_2"}),
              make_tensor_expr({"t", "a_3", "a_4", "i_2", "i_3"}),
              make_tensor_expr({"f", "a_5", "i_1"})});

  auto adjMat = AdjacencyMatrix{prod};

  size_t pos0 = 0, pos1 = 1, pos2 = 2;
  REQUIRE(adjMat.are_connected(pos0, pos1));
  REQUIRE(adjMat.are_connected(pos0, pos2));
  REQUIRE(!adjMat.are_connected(pos1, pos2));
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
    auto F_oo = make_tensor_expr({"f", "i_1", "i_2"});
    auto G_oovv = make_tensor_expr({"g", "i_2", "i_3", "a_1", "a_2"});
    auto T_vo = make_tensor_expr({"t", "a_1", "i_1"});

    auto prod1 = std::make_shared<Product>();
    prod1->append(1, F_oo);
    prod1->append(1, G_oovv);
    prod1->append(1, T_vo);
    prod1->scale(1 / 12.0);

    auto prod2 = std::make_shared<Product>();
    prod2->append(1, F_oo);
    prod2->append(1, G_oovv);
    prod2->append(1, T_vo);
    prod2->append(1, make_tensor_expr({"t", "i_3", "a_2", "a_3", "a_4"}));
    prod2->scale(1 / 8.0);

    auto [fact1, fact2] = factorize_pair(prod1, prod2);

    auto common_subnet = std::make_shared<Product>(Product{F_oo, G_oovv, T_vo});

    auto fact1Expected = std::make_shared<Product>();
    fact1Expected->append(common_subnet);
    fact1Expected->scale(prod1->scalar());
    REQUIRE(*fact1Expected == *fact1);

    auto fact2Expected = std::make_shared<Product>();
    fact2Expected->append(common_subnet);
    fact2Expected->scale(prod2->scalar());

    fact2Expected->append(1, prod2->at(prod2->size() - 1));

    REQUIRE(*fact2Expected == *fact2);

    std::wcout << "Factorization of single common subnets\n";
    std::wcout << "prod1 = " << prod1->to_latex() << "\n";
    std::wcout << "prod2 = " << prod2->to_latex() << "\n";
    std::wcout << "after factorization..\n";
    std::wcout << "prod1 = " << fact1->to_latex() << "\n";
    std::wcout << "prod2 = " << fact2->to_latex() << "\n";
  }

  SECTION("Connected and unconnected common subnets") {
    auto prod1 = std::make_shared<Product>(
        Product{make_tensor_expr({"f", "i_1", "a_1"}),
                make_tensor_expr({"g", "a_1", "a_2", "i_1", "i_2"}),
                make_tensor_expr({"t", "i_2", "a_4"}),
                make_tensor_expr({"t", "a_4", "i_3"})});

    auto prod2 = std::make_shared<Product>(  //
        Product{prod1->at(0), prod1->at(1),
                make_tensor_expr({"t", "i_3", "a_4"}),
                make_tensor_expr({"t", "a_4", "i_4"})});

    auto [fact1, fact2] = factorize_pair(prod1, prod2);

    auto fact1Expected = std::make_shared<Product>();
    fact1Expected->append(
        std::make_shared<Product>(Product{prod1->at(0), prod1->at(1)}));
    fact1Expected->append(
        std::make_shared<Product>(Product{prod1->at(2), prod1->at(3)}));
    REQUIRE(*fact1 == *fact1Expected);

    auto fact2Expected = std::make_shared<Product>();
    fact2Expected->append(
        std::make_shared<Product>(Product{prod2->at(0), prod2->at(1)}));
    fact2Expected->append(
        std::make_shared<Product>(Product{prod2->at(2), prod2->at(3)}));
    REQUIRE(*fact2 == *fact2Expected);

    std::wcout
        << "Factorization of connected and disconnected common subnets\n";
    std::wcout << "prod1 = " << prod1->to_latex() << "\n";
    std::wcout << "prod2 = " << prod2->to_latex() << "\n";
    std::wcout << "after factorization..\n";
    std::wcout << "prod1 = " << fact1->to_latex() << "\n";
    std::wcout << "prod2 = " << fact2->to_latex() << "\n";
  }

  SECTION("Connected subnets with different colors") {
    auto prod1 = std::make_shared<Product>(
        Product{make_tensor_expr({"f", "i_5", "a_1"}),
                make_tensor_expr({"g", "a_1", "a_2", "i_1", "i_2"}),
                make_tensor_expr({"t", "a_4", "i_4"}),
                make_tensor_expr({"t", "i_1", "i_3", "a_3", "a_4"})});

    auto prod2 = std::make_shared<Product>(  //
        Product{prod1->at(0), prod1->at(1),
                make_tensor_expr({"t", "a_4", "i_4"}),
                make_tensor_expr({"t", "i_2", "i_3", "a_3", "a_4"})});

    auto [fact1, fact2] = factorize_pair(prod1, prod2);

    auto fact1Expected = std::make_shared<Product>();
    fact1Expected->append(
        std::make_shared<Product>(Product{prod1->at(0), prod1->at(1)}));
    fact1Expected->append(
        std::make_shared<Product>(Product{prod1->at(2), prod1->at(3)}));
    REQUIRE(*fact1 == *fact1Expected);

    auto fact2Expected = std::make_shared<Product>();
    fact2Expected->append(
        std::make_shared<Product>(Product{prod2->at(0), prod2->at(1)}));
    fact2Expected->append(
        std::make_shared<Product>(Product{prod2->at(2), prod2->at(3)}));
    REQUIRE(*fact2 == *fact2Expected);

    std::wcout << "Factorization of connected common subnets with different "
                  "connectivity\n";
    std::wcout << "prod1 = " << prod1->to_latex() << "\n";
    std::wcout << "prod2 = " << prod2->to_latex() << "\n";
    std::wcout << "after factorization..\n";
    std::wcout << "prod1 = " << fact1->to_latex() << "\n";
    std::wcout << "prod2 = " << fact2->to_latex() << "\n";
  }

  SECTION("Different factorizations based on the second param") {
    auto prod = std::make_shared<Product>(
        Product{make_tensor_expr({"f", "i_1", "a_1"}),
                make_tensor_expr({"g", "a_1", "a_2", "i_2", "i_3"}),
                make_tensor_expr({"t", "i_4", "a_2"}),
                make_tensor_expr({"t", "a_3", "a_4", "i_4", "i_5"})});
    auto prod1 = std::make_shared<Product>(  //
        Product{prod->at(0), prod->at(1),    //
                make_tensor_expr({"t", "a_3", "i_4"}),
                make_tensor_expr({"t", "i_4", "i_5", "a_4", "a_5"})});

    auto prod2 = std::make_shared<Product>(             //
        Product{prod->at(0), prod->at(1), prod->at(2),  //
                make_tensor_expr({"t", "a_3", "a_4", "a_5", "a_6"})});
    //
    // factorization of prod against prod1
    auto [x, y] = factorize_pair(prod, prod1);
    auto xExpected = std::make_shared<Product>();
    xExpected->append(
        std::make_shared<Product>(Product{prod->at(0), prod->at(1)}));
    xExpected->append(
        std::make_shared<Product>(Product{prod->at(2), prod->at(3)}));
    REQUIRE(*x == *xExpected);

    // TODO: Debug
    // factorization of prod against prod2
    /* auto [xx, yy] = factorize_pair(prod, prod2); */
    /* std::wcout << "prod  = " << prod->to_latex() << "\n" */
    /*            << "prod1 = " << prod1->to_latex() << "\n" */
    /*            << "prod2 = " << prod2->to_latex() << "\n"; */
    /* std::wcout << "xx = " << xx->to_latex() << "\n"; */
    /* std::wcout << "yy = " << yy->to_latex() << "\n"; */
    //
    /* auto tn = TensorNetwork(*std::make_shared<Product>( */
    /*     Product({prod->at(0), prod->at(1), prod->at(2)}))); */
    /* tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), false); */
    /* for (const auto& t : tn.tensors()) std::wcout << to_latex(*t) << " "; */
    /* std::wcout << "\n"; */
  }
}
