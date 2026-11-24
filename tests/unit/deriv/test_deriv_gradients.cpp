#include <SeQuant/core/context.hpp>
#include <SeQuant/core/deriv/deriv.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/result_expr.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/spbasis_rotation_gradient.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "catch2_sequant.hpp"
#include "contexts.hpp"
#include "judge.hpp"

#include <algorithm>
#include <numeric>
#include <string>
#include <string_view>

using namespace sequant;

namespace {

/// @p e with every Â and Ŝ factor written out as the normalized projector it
/// stands for: the mean, over the permutations of the bra and of the ket
/// slots (Â, with their signs) or of the columns (Ŝ), of the rest of the term
/// with the indices permuted
ExprPtr projected(const ExprPtr& e) {
  auto expanded = e->clone();
  expand(expanded);
  ExprPtr result = ex<Constant>(0);
  const auto terms =
      expanded->is<Sum>()
          ? container::svector<ExprPtr>(expanded->as<Sum>().summands().begin(),
                                        expanded->as<Sum>().summands().end())
          : container::svector<ExprPtr>{expanded};
  for (const auto& term : terms) {
    if (!term->is<Product>()) {
      result = result + term;
      continue;
    }
    const auto& product = term->as<Product>();
    ExprPtr rest = ex<Constant>(product.scalar());
    const Tensor* projector = nullptr;
    for (const auto& f : product.factors()) {
      const bool is_projector =
          f->is<Tensor>() &&
          (f->as<Tensor>().label() == L"Â" || f->as<Tensor>().label() == L"Ŝ");
      if (is_projector && !projector)
        projector = &f->as<Tensor>();
      else
        rest = rest * f;
    }
    if (!projector) {
      result = result + term;
      continue;
    }
    const bool columns = projector->label() == L"Ŝ";
    const container::svector<Index> B(projector->bra().begin(),
                                      projector->bra().end()),
        K(projector->ket().begin(), projector->ket().end());
    container::svector<std::size_t> order_B(B.size()), order_K(K.size());
    std::iota(order_B.begin(), order_B.end(), std::size_t{0});
    std::iota(order_K.begin(), order_K.end(), std::size_t{0});
    const auto parity = [](const auto& order) {
      int inversions = 0;
      for (std::size_t k = 0; k != order.size(); ++k)
        for (std::size_t l = k + 1; l != order.size(); ++l)
          inversions += order[k] > order[l];
      return inversions % 2;
    };
    ExprPtr sum = ex<Constant>(0);
    int count = 0;
    do {
      do {
        const auto& order_ket = columns ? order_B : order_K;
        container::map<Index, Index> to;
        for (std::size_t k = 0; k != B.size(); ++k)
          to.emplace(B[k], B[order_B[k]]);
        for (std::size_t k = 0; k != K.size(); ++k)
          to.emplace(K[k], K[order_ket[k]]);
        const int sign =
            !columns && parity(order_B) != parity(order_K) ? -1 : 1;
        sum = sum + ex<Constant>(sign) * transform_expr(rest, to);
        ++count;
      } while (!columns &&
               std::next_permutation(order_K.begin(), order_K.end()));
    } while (std::next_permutation(order_B.begin(), order_B.end()));
    // a term may hold another projector
    result = result + projected(ex<Constant>(rational(1, count)) * sum);
  }
  expand(result);
  return result;
}

/// the rotation gradient of @p L by way of its block gradients: for each
/// block x of the leaves labelled as one of @p labels , and each bra or ket
/// slot of x over the space of @p p , the gradient at p in that slot times x
/// with @p q there, the other indices of x summed
ExprPtr via_block_gradients(const ExprPtr& L,
                            const container::svector<std::wstring>& labels,
                            const Index& p, const Index& q) {
  ExprPtr result = ex<Constant>(0);
  for (const auto& label : labels)
    for (const auto& [x, G] : deriv::block_gradients(L, label)) {
      // the evaluator applies a symmetrizer over free labels only
      const auto g = projected(G);
      const auto& t = x->as<Tensor>();
      const auto nb = t.bra_rank();
      for (std::size_t s = 0; s != nb + t.ket_rank(); ++s) {
        const auto u = s < nb ? t.bra()[s] : t.ket()[s - nb];
        if (u.basis() != p.basis()) continue;
        const container::map<Index, Index> u_to_p{{u, p}}, u_to_q{{u, q}};
        result = result + transform_expr(g, u_to_p) * transform_expr(x, u_to_q);
      }
    }
  expand(result);
  return result;
}

}  // namespace

TEST_CASE("deriv/gradients/blocks", "[deriv][gradients]") {
  using deriv::block_gradients;
  auto T = [](std::wstring_view s) { return deserialize(s); };

  SECTION("one entry per class of blocks") {
    // f{i;a} and f{a;i} are other blocks of a tensor without bra-ket symmetry
    const auto E =
        T(L"f{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N + "
          L"f{a_1;i_1}:N-N-N u{i_1;a_1}:N-N-N + "
          L"f{i_2;a_2}:N-N-N w{a_2;i_2}:N-N-N");
    const auto blocks = block_gradients(E, L"f");
    REQUIRE(blocks.size() == 2);
    const auto& x = blocks[0].variable->as<Tensor>();
    REQUIRE(x.label() == L"f");
    const container::map<Index, Index> to_x{{Index(L"a_1"), x.ket()[0]},
                                            {Index(L"i_1"), x.bra()[0]}};
    REQUIRE_THAT(blocks[0].gradient,
                 EquivalentTo(transform_expr(
                     T(L"t{a_1;i_1}:N-N-N + w{a_1;i_1}:N-N-N"), to_x)));
    REQUIRE(block_gradients(E, L"t").size() == 1);
    REQUIRE(block_gradients(E, L"z").empty());
  }

  SECTION("an antisymmetric block: the gradient of every element") {
    // E = 1/4 g t: dE = 1/4 Σ_J t(J) dg(J) over all J, so the gradient is
    // 1/4 t, antisymmetrized (t is antisymmetric already)
    const auto E = T(L"1/4 g{i_1,i_2;a_1,a_2}:A-N-S t{a_1,a_2;i_1,i_2}:A-N-S");
    const auto blocks = block_gradients(E, L"g");
    REQUIRE(blocks.size() == 1);
    const auto& x = blocks[0].variable->as<Tensor>();
    const container::map<Index, Index> to_x{{Index(L"i_1"), x.bra()[0]},
                                            {Index(L"i_2"), x.bra()[1]},
                                            {Index(L"a_1"), x.ket()[0]},
                                            {Index(L"a_2"), x.ket()[1]}};
    deriv::test::require_same(
        blocks[0].gradient,
        transform_expr(T(L"1/4 t{a_1,a_2;i_1,i_2}:A-N-S"), to_x),
        {x.bra()[0], x.bra()[1], x.ket()[0], x.ket()[1]}, {.projectors = true});
  }

  // Euler's identity for E homogeneous of @p degree in the leaves labelled
  // @p label : the gradients contracted with their blocks are degree times E
  const auto euler = [](const ExprPtr& E, std::wstring_view label,
                        std::size_t entries, int degree = 1,
                        const deriv::test::Data& data = {}) {
    const auto blocks = block_gradients(E, label);
    REQUIRE(blocks.size() == entries);
    ExprPtr sum = ex<Constant>(0);
    for (const auto& [x, G] : blocks) sum = sum + projected(G) * x;
    deriv::test::require_same(sum, ex<Constant>(degree) * E, {}, data);
  };

  SECTION("blocks that hold the same elements are one entry") {
    // a bra-ket symmetric tensor of unequal ranks: the exchange is no
    // operation, and the exchanged block still holds the same elements
    euler(T(L"W{a_1,a_2;i_1}:N-N-N g{i_1;a_1,a_2}:N-S-N + "
            L"V{i_1;a_1,a_2}:N-N-N g{a_1,a_2;i_1}:N-S-N"),
          L"g", 1);
    // equal ranks, antisymmetric: mixed blocks and a quadratic function
    euler(T(L"W{a_1,i_2;i_1,a_2}:N-N-N g{i_1,a_2;a_1,i_2}:A-N-S"), L"g", 1);
    euler(T(L"g{i_1,i_2;i_3,a_1}:A-N-S g{i_3,a_1;i_1,i_2}:A-N-S"), L"g", 2, 2);
    {
      const auto scope = deriv::test::real_spaces();
      // the evaluator's data of these symbols is real, as their field is
      const deriv::test::Data real{
          .real = {L"g", L"W", L"V", L"t", L"w", L"f"}};
      euler(T(L"W{a_1,a_2;i_1}:N-N-N g{i_1;a_1,a_2}:N-C-N + "
              L"V{i_1;a_1,a_2}:N-N-N g{a_1,a_2;i_1}:N-C-N"),
            L"g", 1, 1, real);
      // over a real field t⁺ holds the elements of t
      euler(T(L"w{a_1;i_1}:N-N-N t⁺{i_1;a_1}:N-N-N"), L"t", 1, 1, real);
      euler(T(L"f{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N + "
              L"w{a_1;i_1}:N-N-N t⁺{i_1;a_1}:N-N-N"),
            L"t", 1, 1, real);
    }
    // over a complex field t⁺ is another variable
    REQUIRE(block_gradients(T(L"w{a_1;i_1}:N-N-N t⁺{i_1;a_1}:N-N-N"), L"t")
                .empty());
  }

  SECTION("leaves held through a definition") {
    const container::svector<ResultExpr> defs{
        ResultExpr(Tensor(L"D", bra{L"i_7"}, ket{L"a_7"}, Symmetry::Nonsymm),
                   T(L"2 f{i_7;a_7}:N-N-N"))};
    const auto blocks =
        block_gradients(T(L"D{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N + "
                          L"f{i_1;i_2}:N-N-N U{i_2;i_1}:N-N-N"),
                        L"f", defs);
    REQUIRE(blocks.size() == 2);
    // the block met first is that of E itself, then that of the definition
    const auto& x = blocks[1].variable->as<Tensor>();
    const container::map<Index, Index> to_x{{Index(L"i_1"), x.bra()[0]},
                                            {Index(L"a_1"), x.ket()[0]}};
    REQUIRE_THAT(blocks[1].gradient,
                 EquivalentTo(transform_expr(T(L"2 t{a_1;i_1}:N-N-N"), to_x)));
  }

  SECTION("rejections") {
    REQUIRE_THROWS_WITH(
        block_gradients(T(L"f{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N + "
                          L"f{p_1;a_1}:N-N-N u{a_1;p_1}:N-N-N"),
                        L"f"),
        "deriv::block_gradients: two blocks overlap without being one block");
    REQUIRE_THROWS_WITH(
        block_gradients(T(L"f{i_1;a_1<i_1>}:N-N-N t{a_1<i_1>;i_1}:N-N-N"),
                        L"f"),
        "deriv::block_gradients: a leaf with the label has protoindices");
    REQUIRE_THROWS_WITH(block_gradients(nullptr, L"f"),
                        "deriv::block_gradients: null argument");
  }
}

TEST_CASE("deriv/gradients/rotation", "[deriv][gradients]") {
  using deriv::test::require_same;
  using mbpt::spbasis_rotation_gradient;
  auto T = [](std::wstring_view s) { return deserialize(s); };
  const container::svector<std::wstring> fg{L"f", L"g"};

  SECTION("the Hartree-Fock energy: the Fock operator") {
    const auto E = T(L"h{i_1;i_1}:N-N-N + 1/2 g{i_1,i_2;i_1,i_2}:A-N-S");
    const Index i(L"i_5"), a(L"a_5");
    require_same(spbasis_rotation_gradient(E, {L"h", L"g"}, i, a),
                 T(L"h{a_5;i_5}:N-N-N + g{a_5,i_1;i_5,i_1}:A-N-S + "
                   L"h{i_5;a_5}:N-N-N + g{i_5,i_1;a_5,i_1}:A-N-S"),
                 {i, a}, {});
    // a virtual p: no slot of the energy runs over it
    REQUIRE(spbasis_rotation_gradient(E, {L"h", L"g"}, a, i)->is_zero());
    // only the leaves named rotate
    require_same(spbasis_rotation_gradient(E, {L"h"}, i, a),
                 T(L"h{a_5;i_5}:N-N-N + h{i_5;a_5}:N-N-N"), {i, a}, {});
  }

  SECTION("a leaf spelled as the adjoint rotates too") {
    const Index i(L"i_5"), a(L"a_5");
    require_same(
        spbasis_rotation_gradient(T(L"f⁺{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N + "
                                    L"f{i_1;a_1}:N-N-N u{a_1;i_1}:N-N-N"),
                                  {L"f"}, i, a),
        T(L"f⁺{a_5;a_1}:N-N-N t{a_1;i_5}:N-N-N + "
          L"f{a_5;a_1}:N-N-N u{a_1;i_5}:N-N-N"),
        {i, a}, {});
  }

  SECTION("a Hadamard index and an auxiliary slot") {
    const Index i(L"i_5"), a(L"a_5");
    require_same(
        spbasis_rotation_gradient(
            T(L"B{i_1;a_1;x_1}:N-N-N C{a_1;i_1;x_1}:N-N-N D{;;x_1}:N-N-N"),
            {L"B"}, i, a),
        T(L"B{a_5;a_1;x_1}:N-N-N C{a_1;i_5;x_1}:N-N-N D{;;x_1}:N-N-N"), {i, a},
        {});
  }

  SECTION("the CCSD Lagrangian: the block gradients contracted") {
    // two routes to one array: the slots of the leaves substituted, and the
    // gradients by the blocks of f and g contracted with the rotated blocks
    for (const bool traced : {false, true}) {
      mbpt::CC cc{2};
      auto L = simplify(mbpt::op::ref_av((1 + mbpt::op::Λ(2)) * cc.hbar(),
                                         {.connect = cc.hbar_connections()}));
      if (traced) L = simplify(mbpt::closed_shell_spintrace(L, {}));
      const Index i(L"i_55"), j(L"i_56"), a(L"a_55"), b(L"a_56");
      for (const auto& [p, q] : {std::pair{i, a}, std::pair{a, i},
                                 std::pair{i, j}, std::pair{a, b}}) {
        const auto X = spbasis_rotation_gradient(L, fg, p, q);
        REQUIRE_FALSE(X->is_zero());
        require_same(X, via_block_gradients(L, fg, p, q), {p, q}, {});
      }
    }
  }

  SECTION("rejections") {
    const auto E = T(L"f{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N");
    const Index i(L"i_5"), a(L"a_5");
    const auto X = [&](const ExprPtr& e, const Index& p, const Index& q,
                       std::string text) {
      REQUIRE_THROWS_WITH(spbasis_rotation_gradient(e, {L"f"}, p, q),
                          "mbpt::spbasis_rotation_gradient: " + text);
    };
    X(nullptr, i, a, "null argument");
    X(E, i, i, "p and q are one index");
    X(E, Index(L"i_1"), a, "the function uses the label of p or q");
    X(E, Index(L"p_5"), a,
      "a slot of a leaf to rotate overlaps the space of p without including "
      "it");
    X(T(L"f{i_1;a_1}:N-N-N"), i, a, "the function is not a scalar");
    X(T(L"f{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N w{;a_2}:N-N-N"), i, a,
      "the function is not a scalar");
    X(T(L"Â{i_1;a_1}:A-N-S f{i_1;a_1}:N-N-N"), i, a,
      "the function is not a scalar");
    X(T(L"f{i_1;i_2}:N-N-N Y{a_1<i_1>;a_2<i_2>}:N-N-N "
        L"Z{a_2<i_2>;a_1<i_1>}:N-N-N"),
      i, Index(L"i_6"),
      "a slot of a leaf to rotate has protoindices or is a protoindex");
    X(pow(T(L"f{i_1;a_1}:N-N-N"), 2) * T(L"t{a_1;i_1}:N-N-N"), i, a,
      "a leaf to rotate stands under a TensorPower");
    X(mbpt::op::H(1), i, a, "the function holds operators");
  }
}
