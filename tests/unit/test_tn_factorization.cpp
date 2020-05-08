#include "catch.hpp"

#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/domain/factorize/factorizer.hpp>

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

TEST_CASE("FACTORIZER", "[tensor_network]") {
  using std::endl;
  using std::wcout;

  auto prod1 = std::make_shared<Product>(
      Product{make_tensor_expr({"g", "i_1", "i_2", "a_1", "a_2"}),
              make_tensor_expr({"t", "i_1", "a_3"}),
              make_tensor_expr({"t", "i_2", "a_2"})});
  // make_tensor_expr({"f", "i_3", "i_4"})});

  auto prod2 = std::make_shared<Product>(
      Product{make_tensor_expr({"g", "i_1", "i_2", "a_1", "a_2"}),
              make_tensor_expr({"t", "i_3", "a_2"}),
              make_tensor_expr({"t", "i_4", "a_3"}),
              make_tensor_expr({"f", "i_1", "i_2"})});

std::wcout << "\nprod1 is:\n" << prod1->to_latex();
std::wcout << "\nprod2 is:\n" << prod2->to_latex() << std::endl;

  auto [x, y] = largest_common_subnet(prod1, prod2);

  REQUIRE(x.empty());
  REQUIRE(y.empty());
}
