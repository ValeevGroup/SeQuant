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

/// the rotation gradient of @p L by way of its derivatives by the tensors: for
/// each block x of the leaves labelled as one of @p labels , and each bra or
/// ket slot of x over the space of @p p , the gradient at p in that slot times
/// x with @p q there, the other indices of x summed
ExprPtr via_tensor_derivatives(const ExprPtr& L,
                               const container::svector<std::wstring>& labels,
                               const Index& p, const Index& q) {
  ExprPtr result = ex<Constant>(0);
  for (const auto& label : labels)
    for (const auto& [x, G] : deriv::tensor_derivatives(L, label)) {
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

TEST_CASE("deriv/gradients/tensors", "[deriv][gradients]") {
  using deriv::tensor_derivatives;
  auto T = [](std::wstring_view s) { return deserialize(s); };

  SECTION("one entry per class of blocks") {
    // f{i;a} and f{a;i} are other blocks of a tensor without bra-ket symmetry
    const auto E =
        T(L"f{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N + "
          L"f{a_1;i_1}:N-N-N u{i_1;a_1}:N-N-N + "
          L"f{i_2;a_2}:N-N-N w{a_2;i_2}:N-N-N");
    const auto blocks = tensor_derivatives(E, L"f");
    REQUIRE(blocks.size() == 2);
    const auto& x = blocks[0].variable->as<Tensor>();
    REQUIRE(x.label() == L"f");
    const container::map<Index, Index> to_x{{Index(L"a_1"), x.ket()[0]},
                                            {Index(L"i_1"), x.bra()[0]}};
    REQUIRE_THAT(blocks[0].derivative,
                 EquivalentTo(transform_expr(
                     T(L"t{a_1;i_1}:N-N-N + w{a_1;i_1}:N-N-N"), to_x)));
    REQUIRE(tensor_derivatives(E, L"t").size() == 1);
    REQUIRE(tensor_derivatives(E, L"z").empty());
  }

  SECTION("an antisymmetric block: the gradient of every element") {
    // E = 1/4 g t: dE = 1/4 Σ_J t(J) dg(J) over all J, so the gradient is
    // 1/4 t, antisymmetrized (t is antisymmetric already)
    const auto E = T(L"1/4 g{i_1,i_2;a_1,a_2}:A-N-S t{a_1,a_2;i_1,i_2}:A-N-S");
    const auto blocks = tensor_derivatives(E, L"g");
    REQUIRE(blocks.size() == 1);
    const auto& x = blocks[0].variable->as<Tensor>();
    const container::map<Index, Index> to_x{{Index(L"i_1"), x.bra()[0]},
                                            {Index(L"i_2"), x.bra()[1]},
                                            {Index(L"a_1"), x.ket()[0]},
                                            {Index(L"a_2"), x.ket()[1]}};
    deriv::test::require_same(
        blocks[0].derivative,
        transform_expr(T(L"1/4 t{a_1,a_2;i_1,i_2}:A-N-S"), to_x),
        {x.bra()[0], x.bra()[1], x.ket()[0], x.ket()[1]}, {.projectors = true});
  }

  // Euler's identity for E homogeneous of @p degree in the leaves labelled
  // @p label : the gradients contracted with their blocks are degree times E
  const auto euler = [](const ExprPtr& E, std::wstring_view label,
                        std::size_t entries, int degree = 1,
                        const deriv::test::Data& data = {}) {
    const auto blocks = tensor_derivatives(E, label);
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
    REQUIRE(tensor_derivatives(T(L"w{a_1;i_1}:N-N-N t⁺{i_1;a_1}:N-N-N"), L"t")
                .empty());
  }

  SECTION("leaves held through a definition") {
    const container::svector<ResultExpr> defs{
        ResultExpr(Tensor(L"D", bra{L"i_7"}, ket{L"a_7"}, Symmetry::Nonsymm),
                   T(L"2 f{i_7;a_7}:N-N-N"))};
    const auto blocks =
        tensor_derivatives(T(L"D{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N + "
                             L"f{i_1;i_2}:N-N-N U{i_2;i_1}:N-N-N"),
                           L"f", defs);
    REQUIRE(blocks.size() == 2);
    // the block met first is that of E itself, then that of the definition
    const auto& x = blocks[1].variable->as<Tensor>();
    const container::map<Index, Index> to_x{{Index(L"i_1"), x.bra()[0]},
                                            {Index(L"a_1"), x.ket()[0]}};
    REQUIRE_THAT(blocks[1].derivative,
                 EquivalentTo(transform_expr(T(L"2 t{a_1;i_1}:N-N-N"), to_x)));
  }

  SECTION("rejections") {
    REQUIRE_THROWS_WITH(
        tensor_derivatives(T(L"f{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N + "
                             L"f{p_1;a_1}:N-N-N u{a_1;p_1}:N-N-N"),
                           L"f"),
        "deriv::tensor_derivatives: two blocks overlap without being one "
        "block");
    REQUIRE_THROWS_WITH(
        tensor_derivatives(T(L"f{i_1;a_1<i_1>}:N-N-N t{a_1<i_1>;i_1}:N-N-N"),
                           L"f"),
        "deriv::tensor_derivatives: a leaf with the label has protoindices");
    REQUIRE_THROWS_WITH(tensor_derivatives(nullptr, L"f"),
                        "deriv::tensor_derivatives: null argument");
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

  SECTION("the CCSD Lagrangian: the derivatives by the tensors contracted") {
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
        require_same(X, via_tensor_derivatives(L, fg, p, q), {p, q}, {});
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
      "a slot of a leaf that rotates overlaps the space of p without "
      "including it");
    X(T(L"Â{i_1;a_1}:A-N-S f{i_1;a_1}:N-N-N"), i, a,
      "the function holds an Â or Ŝ");
    X(T(L"f{i_1;i_2}:N-N-N Y{a_1<i_1>;a_2<i_2>}:N-N-N "
        L"Z{a_2<i_2>;a_1<i_1>}:N-N-N"),
      i, Index(L"i_6"),
      "a slot of a leaf that rotates has protoindices or is a protoindex");
    X(pow(T(L"f{i_1;a_1}:N-N-N"), 2) * T(L"t{a_1;i_1}:N-N-N"), i, a,
      "a leaf that rotates stands under a TensorPower");
    X(mbpt::op::H(1), i, a, "the function holds operators");
    X(T(L"f{i_1;i_1}:N-N-N n{;;i_1}:N-N-N"), i, a,
      "a slot of a leaf that rotates holds an index that stands in more than "
      "two slots");

    const auto with = [&](std::wstring_view definition, std::string text) {
      const container::svector<ResultExpr> defs{
          ResultExpr(Tensor(L"f", bra{L"p_7"}, ket{L"p_8"}, Symmetry::Nonsymm),
                     T(definition))};
      REQUIRE_THROWS_WITH(spbasis_rotation_gradient(E, {L"g"}, i, a, defs),
                          "mbpt::spbasis_rotation_gradient: " + text);
    };
    with(L"h{p_7;p_8}:N-N-N + g{p_7,i_7;p_8,i_7}:A-N-S",
         "a tensor of the expression of a definition holds an index of the "
         "head and does not rotate");
    with(
        L"g{p_7,i_7;p_8,i_7}:A-N-S + f{p_7;i_7}:N-N-N g{i_7,i_8;p_8,i_8}:A-N-S",
        "a definition depends on itself");
    with(L"g{p_7,i_7;p_8,i_7}:A-N-S + w{p_7;p_8}:N-N-N^(2)",
         "a tensor of the expression of a definition holds an index of the "
         "head and does not rotate");
  }

  // the Fock operator as a function of the basis through the occupied
  // orbitals
  const container::svector<ResultExpr> fock{
      ResultExpr(Tensor(L"f", bra{L"p_7"}, ket{L"p_8"}, Symmetry::Nonsymm),
                 T(L"h{p_7;p_8}:N-N-N + g{p_7,i_7;p_8,i_7}:A-N-S"))};
  // @p e with every f written out as its definition
  const auto written_out = [&](const ExprPtr& e) {
    auto expanded = e->clone();
    expand(expanded);
    ExprPtr result = ex<Constant>(0);
    for (const auto& term : expanded->is<Sum>()
                                ? container::svector<ExprPtr>(
                                      expanded->as<Sum>().summands().begin(),
                                      expanded->as<Sum>().summands().end())
                                : container::svector<ExprPtr>{expanded}) {
      container::set<std::wstring> taken;
      for (const auto& i : get_used_indices(term)) taken.emplace(i.label());
      ExprPtr product =
          ex<Constant>(term->is<Product>() ? term->as<Product>().scalar() : 1);
      for (const auto& f : term->is<Product>()
                               ? container::svector<ExprPtr>(
                                     term->as<Product>().factors().begin(),
                                     term->as<Product>().factors().end())
                               : container::svector<ExprPtr>{term}) {
        const bool is_f = f->is<Tensor>() && f->as<Tensor>().label() == L"f";
        auto g =
            is_f ? deriv::detail::instantiate(fock[0], f->as<Tensor>(), taken)
                 : f;
        for (const auto& i : get_used_indices(g)) taken.emplace(i.label());
        product = product * g;
      }
      result = result + product;
    }
    expand(result);
    return result;
  };

  SECTION("a leaf with a definition: the condition for the Fock operator") {
    // Σ z f: the Fock operator rotates in its own slots, and through the
    // occupied orbitals of its definition
    const auto L = T(L"z{a_1;i_1}:N-N-N f{i_1;a_1}:N-N-N");
    const Index i(L"i_5"), a(L"a_5");
    const container::svector<std::wstring> hg{L"h", L"g"};
    require_same(spbasis_rotation_gradient(L, hg, i, a, fock),
                 T(L"z{a_1;i_5}:N-N-N f{a_5;a_1}:N-N-N + "
                   L"z{a_1;i_1}:N-N-N g{i_1,a_5;a_1,i_5}:A-N-S + "
                   L"z{a_1;i_1}:N-N-N g{i_1,i_5;a_1,a_5}:A-N-S"),
                 {i, a}, {});
    REQUIRE_THROWS_WITH(
        spbasis_rotation_gradient(T(L"f⁺{a_1;i_1}:N-N-N t{i_1;a_1}:N-N-N"), hg,
                                  i, a, fock),
        "mbpt::spbasis_rotation_gradient: a leaf is spelled as the adjoint of "
        "one with a definition");
    // the same array as that of the function with f written out
    for (const auto& [p, q] :
         {std::pair{i, a}, std::pair{a, i}, std::pair{i, Index(L"i_6")},
          std::pair{a, Index(L"a_6")}})
      require_same(written_out(spbasis_rotation_gradient(L, hg, p, q, fock)),
                   spbasis_rotation_gradient(written_out(L), hg, p, q), {p, q},
                   {});
  }

  SECTION("the CCSD Lagrangian with the Fock operator defined") {
    mbpt::CC cc{2};
    const auto L = simplify(mbpt::op::ref_av(
        (1 + mbpt::op::Λ(2)) * cc.hbar(), {.connect = cc.hbar_connections()}));
    const container::svector<std::wstring> hg{L"h", L"g"};
    const Index i(L"i_55"), a(L"a_55");
    for (const auto& [p, q] : {std::pair{i, a}, std::pair{a, i}}) {
      const auto X = spbasis_rotation_gradient(L, hg, p, q, fock);
      require_same(written_out(X),
                   spbasis_rotation_gradient(written_out(L), hg, p, q), {p, q},
                   {});
      // the occupied orbitals of the definition rotate with an occupied p
      // only: then there is more than at a fixed Fock operator
      REQUIRE(simplify(X - spbasis_rotation_gradient(L, fg, p, q))->is_zero() ==
              (p == a));
    }
  }

  SECTION("free indices, and the second derivative") {
    // a slot that holds a free index: its Kronecker delta with p
    const Index i(L"i_5"), a(L"a_5"), j(L"i_6"), b(L"a_6");
    REQUIRE_THAT(
        spbasis_rotation_gradient(T(L"f{i_1;a_1}:N-N-N"), {L"f"}, i, a),
        EquivalentTo(T(L"f{a_5;a_1}:N-N-N δ{i_1;i_5}:N-S-S")));
    // L = Σ h{i;i}: X{q;p} = h{q;p} + h{p;q}, and with the slots of q held
    // the slots of p rotate: ∂²L/∂κ{q;p}∂κ{s;r} = δ{p;r} (h{q;s} + h{s;q})
    const auto L = T(L"h{i_1;i_1}:N-N-N");
    const auto X = spbasis_rotation_gradient(L, {L"h"}, i, a);
    REQUIRE_THAT(X, EquivalentTo(T(L"h{a_5;i_5}:N-N-N + h{i_5;a_5}:N-N-N")));
    REQUIRE_THAT(spbasis_rotation_gradient(X, {L"h"}, j, b, {}, {a}),
                 EquivalentTo(T(L"h{a_5;a_6}:N-N-N δ{i_6;i_5}:N-S-S + "
                                L"h{a_6;a_5}:N-N-N δ{i_5;i_6}:N-S-S")));
  }
}
