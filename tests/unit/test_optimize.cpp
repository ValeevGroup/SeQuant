#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/optimize/common_subexpression_elimination.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/optimize/single_term.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <cstddef>
#include <initializer_list>
#include <memory>

sequant::ExprPtr extract(sequant::ExprPtr expr,
                         std::initializer_list<size_t> const& idxs) {
  using namespace sequant;
  ExprPtr result = expr;
  for (auto s : idxs) result = result->at(s);
  return result;
}

// number of Tensor leaves in a (binarized) expression tree
size_t count_tensor_leaves(sequant::ExprPtr const& expr) {
  using namespace sequant;
  size_t n = 0;
  expr->visit([&n](auto const& x) { n += x->template is<Tensor>() ? 1 : 0; },
              /*atoms_only=*/true);
  return n;
}

TEST_CASE("optimize", "[optimize]") {
  using namespace sequant;

  // for optimization tests, need to specify index space sizes, so make a clone
  // of the context
  {
    auto ctx_resetter =
        set_scoped_default_context(get_default_context().clone());
    auto reg = get_default_context().mutable_index_space_registry();

    {
      auto occ = reg->retrieve_ptr(L"i");
      auto uocc = reg->retrieve_ptr(L"a");
      auto aux = reg->retrieve_ptr(L"x");
      REQUIRE(occ);
      REQUIRE(uocc);
      REQUIRE(aux);
      occ->approximate_size(10);
      uocc->approximate_size(100);
      aux->approximate_size(4);
      REQUIRE(uocc->approximate_size() == 100);
    }

    auto single_term_opt = [](Product const& prod) {
      return opt::single_term_opt(prod, [](Index const& ix) {
        // null space contributes x1 to the size
        auto sz = ix.nonnull() ? ix.space().approximate_size() : 1;
        return sz;
      });
    };

    auto parse_expr_antisymm = [](auto const& xpr) {
      return deserialize(xpr, {.def_perm_symm = Symmetry::Antisymm});
    };

    SECTION("Single term optimization") {
      const auto prod1 = parse_expr_antisymm(
                             L"g_{i3,i4}^{a3,a4}"     // T1
                             " * t_{a1,a2}^{i3,i4}"   // T2
                             " * t_{a3,a4}^{i1,i2}")  // T3
                             ->as<Product>();
      //
      // Cost of evaluation prod1:
      //
      // ((T1 * T2) * T3)  : 2 * O^2 * V^4  best if nocc > nvirt
      //
      // this is the one we want to find
      // ((T1 * T3) * T2)  : 2 * O^4 * V^2  best if nvirt > nocc
      //
      // (T1 * (T2 * T3))  : 2 * O^4 * V^4  worst sequence of evaluation
      //

      const auto res1 = single_term_opt(prod1);

      REQUIRE(extract(res1, {0, 0}) == prod1.at(0));
      REQUIRE(extract(res1, {0, 1}) == prod1.at(2));
      REQUIRE(extract(res1, {1}) == prod1.at(1));

      const auto prod2 = parse_expr_antisymm(
                             L"   g_{i3,i4}^{a3,a4}"
                             L" * t_{a3,a4}^{i1,i2}"
                             L" * t_{a1}^{i3}"
                             L" * t_{a2}^{i4}")
                             ->as<Product>();

      const auto res2 = single_term_opt(prod2);

      REQUIRE(extract(res2, {0, 0, 0}) == prod2.at(0));
      REQUIRE(extract(res2, {0, 0, 1}) == prod2.at(1));
      REQUIRE(extract(res2, {0, 1}) == prod2.at(2));
      REQUIRE(extract(res2, {1}) == prod2.at(3));

      const auto prod3 = parse_expr_antisymm(
                             L""                   //
                             " g_{i3,i4}^{a3,a4}"  //
                             " t_{a1}^{i3}"        //
                             " t_{a2}^{i4}"        //
                             " t_{a3,a4}^{i1,i2}"  //
                             )
                             ->as<Product>();
      auto res3 = single_term_opt(prod3);

      REQUIRE(extract(res3, {0, 0, 0}) == prod3.at(0));
      REQUIRE(extract(res3, {0, 0, 1}) == prod3.at(3));
      REQUIRE(extract(res3, {0, 1}) == prod3.at(1));
      REQUIRE(extract(res3, {1}) == prod3.at(2));

      //
      // single-term optimization when a dot product occurs in the tensor
      // network
      // ========================

      auto prod4 =
          parse_expr_antisymm(L"1/4 λ{i1;a1} g{i2,i3;a2,a3} t{a2,a3;i2,i3}")
              ->as<Product>();
      auto res4 = single_term_opt(prod4);

      REQUIRE(extract(res4, {0}) == prod4.at(0));
      REQUIRE(extract(res4, {1, 0}) == prod4.at(1));
      REQUIRE(extract(res4, {1, 1}) == prod4.at(2));

      auto prod5 =
          parse_expr_antisymm(L"x{i1,i2;a3,a4} y{a1,a2;i1,i2} z{a3,a4;a1,a2}")
              ->as<Product>();
      auto res5 = single_term_opt(prod5);
      REQUIRE(extract(res5, {0, 0}) == prod5.at(0));
      REQUIRE(extract(res5, {0, 1}) == prod5.at(2));
      REQUIRE(extract(res5, {1}) == prod5.at(1));

      //
      // single-term optimization when sequant::Variables appear in a product
      //
      auto prod6 = deserialize(
                       L"α * β * γ * "
                       "g_{i3,i4}^{a3,a4}"      // T1
                       " * t_{a1,a2}^{i3,i4}"   // T2
                       " * t_{a3,a4}^{i1,i2}")  // T3
                       ->as<Product>();
      auto res6 = single_term_opt(prod6);

      // this is the one we want to find
      // α * β * γ * ((T1 * T3) * T2)  : 2 * O^4 * V^2  best if nvirt > nocc
      REQUIRE(extract(res6, {0}) == prod6.at(0));
      REQUIRE(extract(res6, {1}) == prod6.at(1));
      REQUIRE(extract(res6, {2}) == prod6.at(2));
      REQUIRE(extract(res6, {3, 0}) == prod6.at(3));
      REQUIRE(extract(res6, {3, 1}) == prod6.at(5));
      REQUIRE(extract(res6, {4}) == prod6.at(4));

      //
      // single-term optimization including tensors with auxiliary indices
      //
      auto prod7 = deserialize(
                       L"DF{a_1;a_3;x_1} "  // T1
                       "DF{a_2;i_1;x_1} "   // T2
                       "t{a_3;i_2}"         // T3
                       )
                       ->as<Product>();
      auto res7 = single_term_opt(prod7);

      // this is the one we want to find
      // (T1 T3) T2: V^2 O^1 A^1 + V^2 O^2 A^1 best if nvirt > nocc and nvirt >
      // nact
      REQUIRE(extract(res7, {0, 0}) == prod7.at(0));
      REQUIRE(extract(res7, {0, 1}) == prod7.at(2));
      REQUIRE(extract(res7, {1}) == prod7.at(1));

      auto prod8 =
          deserialize(
              L"T1{i_1;i_2;x_1,x_2,x_3,x_4} T2{i_2;i_1;x_5,x_6,x_7,x_8} "
              L"T3{i_3;;x_1,x_2,x_3,x_4} T4{i_4;;x_5,x_6,x_7,x_8}")
              ->as<Product>();
      auto res8 = single_term_opt(prod8);

      // this is the one we want to find
      // (T1 T3)(T2 T4)
      REQUIRE(extract(res8, {0, 0}) == prod8.at(0));
      REQUIRE(extract(res8, {0, 1}) == prod8.at(2));
      REQUIRE(extract(res8, {1, 0}) == prod8.at(1));
      REQUIRE(extract(res8, {1, 1}) == prod8.at(3));
    }

    SECTION("Single term optimization: n_replay volatility weighting") {
      using namespace sequant;

      // PPL-shaped motif, fully contracted to a scalar:
      //   A = g_{i1,a1}^{x1}   (persistent integral)
      //   B = g_{i2,a2}^{x1}   (persistent integral)
      //   t = t_{a1,a2}^{i1,i2} (VOLATILE amplitude)
      // sizes: i=10 (O), a=100 (V), x=4 (X).
      //
      // (A*B)*t : build I=A*B over x  -> {i1,a1,i2,a2}  cost O^2 V^2 X
      // (persistent)
      //           then I*t            -> scalar         cost O^2 V^2 (volatile)
      // (A*t)*B : build J=A*t over i1,a1 -> {x,i2,a2}   cost O^2 V^2 X
      // (VOLATILE)
      //           then J*B               -> scalar      cost X O V (volatile)
      //
      // n_replay=1  : (A*t)*B wins (O^2 V^2 X + X O V  <  O^2 V^2 X + O^2 V^2)
      //               => t buried in an inner volatile intermediate.
      // n_replay=10 : (A*B)*t wins (persistent build counted once; the only
      //               x10 term is the cheap O^2 V^2 final step)
      //               => t contracted LAST, persistent integral formed first.
      auto idxsz = [](Index const& ix) {
        return ix.nonnull() ? ix.space().approximate_size() : std::size_t{1};
      };

      auto prod = parse_expr_antisymm(
                      L"g_{i1,a1}^{x1}"
                      L" * g_{i2,a2}^{x1}"
                      L" * t_{a1,a2}^{i1,i2}")
                      ->as<Product>();

      auto is_t = [](Tensor const& t) { return t.label() == L"t"; };

      OptimizeOptions base;
      base.idx_to_extent = idxsz;

      // baseline: predicate set but n_replay==1 => weight is 1 everywhere =>
      // reverts to current behavior. (The empty-predicate no-op is checked
      // separately via opts_off below.)
      auto opts1 = base;
      opts1.is_volatile_leaf = is_t;
      opts1.n_replay = 1;
      auto res1 = optimize(ex<Product>(prod), opts1);

      auto opts10 = base;
      opts10.is_volatile_leaf = is_t;
      opts10.n_replay = 10;
      auto res10 = optimize(ex<Product>(prod), opts10);

      // a bare top-level t leaf means t was contracted last (persistent-first)
      auto top_has_bare_t = [](ExprPtr const& e) {
        if (!e->is<Product>()) return false;
        for (auto const& c : *e)
          if (c->is<Tensor>() && c->as<Tensor>().label() == L"t") return true;
        return false;
      };

      // weighting flips the chosen factorization
      REQUIRE(res1 != res10);
      // n_replay=1 reproduces today's behavior: t buried in an inner
      // intermediate
      REQUIRE_FALSE(top_has_bare_t(res1));
      // n_replay=10: persistent g*g built first, t contracted last
      REQUIRE(top_has_bare_t(res10));

      // empty predicate => weighting off => identical to n_replay=1 regardless
      auto opts_off = base;
      opts_off.n_replay = 10;  // ignored: predicate empty
      auto res_off = optimize(ex<Product>(prod), opts_off);
      REQUIRE(res_off == res1);
    }

    SECTION("Ensure single-value sums/products are not discarded") {
      auto sum = ex<Sum>();
      sum->as<Sum>().append(
          ex<Product>(ExprPtrList{deserialize(L"f{a_1;i_1}")}));
      REQUIRE(sum->as<Sum>().summand(0).as<Product>().factors().size() == 1);
      auto optimized = optimize(sum);
      REQUIRE(optimized->is<Sum>());
      REQUIRE(optimized->as<Sum>().summands().size() == 1);
      REQUIRE(sum->as<Sum>().summand(0).as<Product>().factors().size() == 1);
    }

    SECTION("Non-covariant indices") {
      auto uocc = reg->retrieve_ptr(L"a");
      auto aux = reg->retrieve_ptr(L"x");
      auto const aux_sz = aux->approximate_size();
      aux->approximate_size(3 * uocc->approximate_size());

      auto const G_abcd_thc =
          deserialize(L"X{a1;;x1} X{;a2;x1} Y{;;x1,x2} X{a3;;x2} X{;a4;x2}")
              ->as<Product>();
      auto const G_abcd_thc_opt =
          deserialize(
              L"((X{a1;;x1} X{;a2;x1}) Y{;;x1,x2})(X{a3;;x2} X{;a4;x2})")
              ->as<Product>();
      REQUIRE(single_term_opt(G_abcd_thc)->as<Product>() == G_abcd_thc_opt);

      auto const GT_abij_thc = deserialize(
                                   L"X{a1;;x1} X{;a2;x1} Y{;;x1,x2} X{a3;;x2} "
                                   L"X{;a4;x2} T{a2,a4;i1,i2}")
                                   ->as<Product>();
      auto const GT_abij_thc_opt = deserialize(
                                       L"(((X{a1;;x1} X{;a2;x1}) Y{;;x1,x2}) ( "
                                       L"X{;a4;x2} T{a2,a4;i1,i2} )) X{a3;;x2}")
                                       ->as<Product>();
      REQUIRE(single_term_opt(GT_abij_thc)->as<Product>() == GT_abij_thc_opt);

      aux->approximate_size(aux_sz);
    }

    SECTION("OptimizeOptions: cost metric and reorder knobs") {
      auto const prod = parse_expr_antisymm(
          L"g_{i3,i4}^{a3,a4} t_{a1,a2}^{i3,i4} t_{a3,a4}^{i1,i2}");

      // both metrics must binarize the 3-tensor product into a binary tree:
      // a 2-factor top product whose leaves are the 3 original tensors
      for (auto opt_for : {OptFor::Flops, OptFor::Memsize}) {
        CAPTURE(static_cast<int>(opt_for));
        auto res = optimize(prod, OptimizeOptions{.opt_for = opt_for});
        REQUIRE(res->is<Product>());
        REQUIRE(res->as<Product>().factors().size() == 2);
        REQUIRE(count_tensor_leaves(res) == 3);
      }

      // reorder knob: a two-summand sum is optimized either way, and the
      // optimize() default (reorder) matches an explicit Reorder request
      auto const sum = parse_expr_antisymm(
          L"g_{i3,i4}^{a3,a4} t_{a1,a2}^{i3,i4} t_{a3,a4}^{i1,i2}"
          L" + g_{i3,i4}^{a3,a4} t_{a3,a4}^{i1,i2} t_{a1}^{i3} t_{a2}^{i4}");
      REQUIRE(sum->is<Sum>());

      auto no_reorder =
          optimize(sum, OptimizeOptions{.reorder = ReorderSum::NoReorder});
      auto reorder =
          optimize(sum, OptimizeOptions{.reorder = ReorderSum::Reorder});
      REQUIRE(no_reorder->is<Sum>());
      REQUIRE(reorder->is<Sum>());
      REQUIRE(no_reorder->as<Sum>().size() == sum->as<Sum>().size());
      REQUIRE(reorder->as<Sum>().size() == sum->as<Sum>().size());
      // default options == explicit Reorder
      REQUIRE(*optimize(sum) == *reorder);
    }

    SECTION("Parallel optimization of summands matches sequential") {
      // exercise optimize_impl(..., parallel_outer=true): a multi-summand sum
      // optimized concurrently must yield the same result as single-threaded.
      auto const sum = parse_expr_antisymm(
          L"g_{i3,i4}^{a3,a4} t_{a1,a2}^{i3,i4} t_{a3,a4}^{i1,i2}"
          L" + g_{i3,i4}^{a3,a4} t_{a3,a4}^{i1,i2} t_{a1}^{i3} t_{a2}^{i4}"
          L" + g_{i3,i4}^{a3,a4} t_{a1}^{i3} t_{a2}^{i4} t_{a3,a4}^{i1,i2}");
      REQUIRE(sum->is<Sum>());
      REQUIRE(sum->as<Sum>().size() > 1);

      auto const nthreads_save = num_threads();
      struct ThreadGuard {
        int n;
        ~ThreadGuard() { set_num_threads(n); }
      } guard{nthreads_save};

      set_num_threads(1);
      auto const seq = optimize(sum);
      set_num_threads(4);
      auto const par = optimize(sum);

      REQUIRE(*seq == *par);
    }
  }

  SECTION("CSE") {
    auto ctx_resetter =
        set_scoped_default_context(get_default_context().clone());
    IndexSpaceRegistry registry;
    registry.add("a", 0b001);
    registry.add("i", 0b010);
    registry.add("u", 0b100);
    *get_default_context().mutable_index_space_registry() = registry;

    auto binarizer = [](auto&& expr) {
      // CSE drives binarize() on subexpressions for hash-equivalence
      // detection; positional head is irrelevant here.
      SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
      return binarize(expr);
      SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    };

    auto collect_as_expr = [](auto&& expressions) {
      std::vector<ResultExpr> actual;
      for (const auto& current : expressions) {
        if (current->is_tensor()) {
          actual.emplace_back(current->expr()->template as<Tensor>(),
                              to_expr(current));
        } else {
          REQUIRE(current->is_scalar());
          actual.emplace_back(current->expr()->template as<Variable>(),
                              to_expr(current));
        }
      }

      return actual;
    };

    auto parse_inputs = [](auto&& inputs) {
      std::vector<EvalNode<EvalExpr>> expressions;
      for (const std::wstring& current : inputs) {
        expressions.push_back(binarize(deserialize<ResultExpr>(
            current, {.def_perm_symm = Symmetry::Nonsymm,
                      .def_braket_symm = BraKetSymmetry::Nonsymm,
                      .def_col_symm = ColumnSymmetry::Nonsymm})));
      }
      return expressions;
    };

    auto parse_expected = [](auto&& outputs) {
      std::vector<ResultExpr> expected;
      for (const std::wstring& current : outputs) {
        expected.push_back(deserialize<ResultExpr>(
            current, {.def_perm_symm = Symmetry::Nonsymm,
                      .def_braket_symm = BraKetSymmetry::Nonsymm,
                      .def_col_symm = ColumnSymmetry::Nonsymm}));
      }
      return expected;
    };

    SECTION("standard") {
      for (bool force_hash_collisions : {false, true}) {
        CAPTURE(force_hash_collisions);

        for (const auto& [inputs, outputs] :
             std::vector<std::pair<std::vector<std::wstring>,
                                   std::vector<std::wstring>>>{
                 // Basic example with only scalars
                 {{L"R1 = (A B) C", L"R2 = D (A B)"},
                  {L"CSE1 = A B", L"R1 = CSE1 C", L"R2 = D CSE1"}},
                 // Test case in which the same intermediate is reused but
                 // requires different indexing
                 {{L"R{a1,a3;i2,i3} = 2 GAM0{a1,a3;a4,a5} T2g{a4,a5;i2,i3} - "
                   L"GAM0{a1,a3;a4,a5} T2g{a4,a5;i3,i2}"},
                  {L"CSE1{;;a3,a1,i2,i3} = GAM0{a1,a3;a4,a5} T2g{a4,a5;i2,i3}",
                   L"R{a1,a3;i2,i3} = 2 CSE1{;;a3,a1,i2,i3} - "
                   L"CSE1{;;a3,a1,i3,i2}"}},
                 // Scalar CSE with proto-index tensors
                 {{L"R1 = (f{i1;a1<i1>} t{a1<i1>;i1}) A",
                   L"R2 = B (f{i1;a1<i1>} t{a1<i1>;i1})"},
                  {L"CSE1 = f{i1;a1<i1>} t{a1<i1>;i1}", L"R1 = CSE1 A",
                   L"R2 = B CSE1"}},
                 // Tensor CSE with proto-index tensors: reused
                 // contraction with different external indexing
                 {{L"R{i1;i2} = 2 g{i1;a1<i1>} t{a1<i1>;i2} - "
                   L"g{i2;a1<i2>} t{a1<i2>;i1}"},
                  {L"CSE1{;;i1,i2} = g{i1;a1<i1>} t{a1<i1>;i2}",
                   L"R{i1;i2} = 2 CSE1{;;i1,i2} - CSE1{;;i2,i1}"}},
                 // ToT CSE: the intermediate itself has proto-indexed
                 // indices (tensor-of-tensor)
                 {{L"R1{i1;i2} = (g{i1;a1} C{a1;a1<i1>}) h{a1<i1>;i2}",
                   L"R2{i1;i2} = (g{i1;a1} C{a1;a1<i1>}) k{a1<i1>;i2}"},
                  {L"CSE1{;;a1<i1>,i1} = g{i1;a1} C{a1;a1<i1>}",
                   L"R1{i1;i2} = CSE1{;;a1<i1>,i1} h{a1<i1>;i2}",
                   L"R2{i1;i2} = CSE1{;;a1<i1>,i1} k{a1<i1>;i2}"}},
                 // In this case it is important that the computation of the
                 // subexpression isn't simply thrown at the beginning of the
                 // expression list as it depends on B, which has to be computed
                 // first.
                 {{L"B = K J", L"R = (A B) C + (A B) D"},
                  {L"B = K J", L"CSE1 = A B", L"R = CSE1 C + CSE1 D"}},
                 // CSE in the presence of bra-ket symmetry
                 {{L"R2{u2,a1;u1,i1} = -2 f{u3;u4}:N-S Y{u2,u3;u1,u5} "
                   L"t{a1,u5;i1,u4} + f{u3;u4}:N-S Y{u2,u4;u5,u1} "
                   L"t{a1,u5;i1,u3} "
                   L"+ f{u3;u4}:N-S Y{u2,u4;u1,u5} t{a1,u5;u3,i1}"},
                  {L"CSE1{;;u4,u2,u1,u5} = f{u3;u4}:N-S Y{u2,u3;u1,u5}",
                   L"R2{u2,a1;u1,i1} = -2 CSE1{;;u4,u2,u1,u5} t{a1,u5;i1,u4}"
                   L" + CSE1{;;u3,u2,u5,u1} t{a1,u5;i1,u3}"
                   L" + CSE1{;;u3,u2,u1,u5} t{a1,u5;u3,i1}"}},
             }) {
          CAPTURE(inputs);

          std::vector<EvalNode<EvalExpr>> expressions = parse_inputs(inputs);
          const std::vector<ResultExpr> expected = parse_expected(outputs);

          if (force_hash_collisions) {
            // This code path makes all hashes be computed to be zero and hence
            // every pair of objects will yield a hash collision which need to
            // be dealt with by using proper comparison operators.
            static constexpr bool force_collisions = true;
            opt::eliminate_common_subexpressions<
                decltype(expressions), decltype(binarizer), force_collisions>(
                expressions, binarizer);
          } else {
            opt::eliminate_common_subexpressions(expressions, binarizer);
          }

          REQUIRE(collect_as_expr(expressions) == expected);
        }
      }
    }
    SECTION("batch indices") {
      const opt::CSEOptions<EvalNode<EvalExpr>> opts = {.batch_indices = {
                                                            "i5",
                                                            "i6",
                                                            "a5",
                                                            "a6",
                                                        }};

      for (const auto& [inputs, outputs] : std::vector<
               std::pair<std::vector<std::wstring>, std::vector<std::wstring>>>{
               // Can't eliminate any CSE due to differences in batching indices
               {{L"R = (A{i1;i5} B{i5;i2}) C{i2;i1} + "
                 L"(A{i1;i6} B{i6;i2}) D{i2;i1}"},
                {L"R = (A{i1;i5} B{i5;i2}) C{i2;i1} + "
                 L"(A{i1;i6} B{i6;i2}) D{i2;i1}"}},
               {{L"R = (A{i1;i5} B{i5;i2}) C{i2;i1} + "
                 L"(A{i1;i3} B{i3;i2}) D{i2;i1}"},
                {L"R = (A{i1;i5} B{i5;i2}) C{i2;i1} + "
                 L"(A{i1;i3} B{i3;i2}) D{i2;i1}"}},
               // Can eliminate if batched index is same
               {{L"R = (A{i1;i5} B{i5;i2}) C{i2;i1} + "
                 L"(A{i3;i5} B{i5;i4}) D{i4;i3}"},
                {L"CSE1{;;i1,i2} = A{i1;i5} B{i5;i2}",
                 L"R = CSE1{;;i1,i2} C{i2;i1} + "
                 L"CSE1{;;i3,i4} D{i4;i3}"}},
           }) {
        CAPTURE(inputs);

        std::vector<EvalNode<EvalExpr>> expressions = parse_inputs(inputs);
        const std::vector<ResultExpr> expected = parse_expected(outputs);

        opt::eliminate_common_subexpressions(expressions, binarizer, opts);

        REQUIRE(collect_as_expr(expressions) == expected);
      }
    }
  }

  SECTION("Single term optimization with CSE") {
    auto ctx_resetter =
        set_scoped_default_context(get_default_context().clone());
    auto reg = get_default_context().mutable_index_space_registry();
    mbpt::add_df_spaces(reg);
    mbpt::add_pao_spaces(reg);
    mbpt::add_ao_spaces(reg);
    // i 10
    // a 40
    // μ̃ 50
    // Κ 90
    for (auto&& [k, v] :
         std::initializer_list<std::pair<std::wstring_view, size_t>>{
             {L"i", 10}, {L"a", 40}, {L"μ̃", 50}, {L"Κ", 90}}) {
      reg->retrieve_ptr(k)->approximate_size(v);
    }

    auto single_term_opt = [](Product const& prod, bool cse = true) {
      return opt::single_term_opt(
          prod,
          [](Index const& ix) {
            // null space contributes x1 to the size
            auto sz = ix.nonnull() ? ix.space().approximate_size() : 1;
            return sz;
          },
          /*subnet_cse=*/cse);
    };

    auto prod9 =
        deserialize("X{i1;a1} X{i2;a2} Y{a2;i3} Y{a1;i4}")->as<Product>();
    auto res9 = single_term_opt(prod9);
    auto res9_no_cse = single_term_opt(prod9, false);
    // this is the one we want to find
    // (X Y) (X Y)
    REQUIRE(extract(res9, {0, 0}) == prod9.at(0));
    REQUIRE(extract(res9, {0, 1}) == prod9.at(3));
    REQUIRE(extract(res9, {1, 0}) == prod9.at(1));
    REQUIRE(extract(res9, {1, 1}) == prod9.at(2));

    // take a look at res9_no_cse for a result with subnet_cse disabled
    // should give the same result in this case as it's already optimal
    REQUIRE(extract(res9_no_cse, {0, 0}) == prod9.at(0));
    REQUIRE(extract(res9_no_cse, {0, 1}) == prod9.at(3));
    REQUIRE(extract(res9_no_cse, {1, 0}) == prod9.at(1));
    REQUIRE(extract(res9_no_cse, {1, 1}) == prod9.at(2));

    SECTION("CSE effect on optimization result") {
      auto ctx_resetter =
          set_scoped_default_context(get_default_context().clone());
      auto reg = get_default_context().mutable_index_space_registry();
      // Use sizes that make the unbalanced tree better without CSE,
      // but the balanced tree better with CSE.
      // Balanced: ( (X1 Y1) (X2 Y2) )
      // Cost(X1*Y1) = size(i)*size(a)*size(j) = 12*10*12 = 1440.
      // Cost(Inter) = 12^3 = 1728.
      // Total no-CSE: 2*1440 + 1728 = 4608.
      // Total CSE: 1440 + 1728 = 3168.
      // Unbalanced: ( ( (X1 Y1) X2 ) Y2 )
      // Cost(X1*Y1) = 12*10*12 = 1440.
      // Cost((X1*Y1)*X2) = size(i)*size(i)*size(a) = 12*12*10 = 1440.
      // Cost(...) * Y2 = 12*10*12 = 1440.
      // Total Unbalanced: 1440 + 1440 + 1440 = 4320.
      // 3168 < 4320 < 4608.
      reg->retrieve_ptr(L"i")->approximate_size(12);
      reg->retrieve_ptr(L"a")->approximate_size(10);

      auto single_term_opt = [](Product const& prod, bool cse) {
        return opt::single_term_opt(
            prod,
            [](Index const& ix) {
              return ix.nonnull() ? ix.space().approximate_size() : 1;
            },
            cse);
      };

      // X{i1;a1} Y{a1;i2} X{i2;a2} Y{a2;i3}
      auto prod =
          deserialize(L"X{i1;a1} Y{a1;i2} X{i2;a2} Y{a2;i3}")->as<Product>();

      auto res_cse = single_term_opt(prod, true);
      auto res_no_cse = single_term_opt(prod, false);

      // With CSE: Balanced tree
      REQUIRE(res_cse->as<Product>().factors().size() == 2);
      REQUIRE(res_cse->at(0)->is<Product>());
      REQUIRE(res_cse->at(1)->is<Product>());

      // Without CSE: Unbalanced tree
      bool is_unbalanced =
          (res_no_cse->at(0)->is<Tensor>() || res_no_cse->at(1)->is<Tensor>());
      REQUIRE(is_unbalanced);
    }

    SECTION("subnet_cse flows through OptimizeOptions") {
      auto ctx_resetter =
          set_scoped_default_context(get_default_context().clone());
      auto reg = get_default_context().mutable_index_space_registry();
      // Same sizing trick as the section above: CSE prefers balanced,
      // no-CSE prefers unbalanced.
      reg->retrieve_ptr(L"i")->approximate_size(12);
      reg->retrieve_ptr(L"a")->approximate_size(10);

      auto idx_to_extent = [](Index const& ix) -> std::size_t {
        return ix.nonnull() ? ix.space().approximate_size() : 1;
      };

      auto prod =
          deserialize(L"X{i1;a1} Y{a1;i2} X{i2;a2} Y{a2;i3}")->as<Product>();
      auto expr = ex<Product>(prod);

      auto res_cse =
          optimize(expr, OptimizeOptions{.subnet_cse = SubnetCSE::Enable,
                                         .idx_to_extent = idx_to_extent});
      auto res_no_cse =
          optimize(expr, OptimizeOptions{.subnet_cse = SubnetCSE::Disable,
                                         .idx_to_extent = idx_to_extent});

      // With CSE: balanced tree -- both children are Products.
      REQUIRE(res_cse->is<Product>());
      REQUIRE(res_cse->as<Product>().factors().size() == 2);
      REQUIRE(res_cse->at(0)->is<Product>());
      REQUIRE(res_cse->at(1)->is<Product>());

      // Without CSE: unbalanced tree -- at least one child is a bare Tensor.
      REQUIRE(res_no_cse->is<Product>());
      REQUIRE(res_no_cse->as<Product>().factors().size() == 2);
      bool is_unbalanced =
          res_no_cse->at(0)->is<Tensor>() || res_no_cse->at(1)->is<Tensor>();
      REQUIRE(is_unbalanced);

      // Default OptimizeOptions => subnet_cse Disable => same as no-CSE shape.
      auto res_default =
          optimize(expr, OptimizeOptions{.idx_to_extent = idx_to_extent});
      REQUIRE(res_default->is<Product>());
      REQUIRE(res_default->as<Product>().factors().size() == 2);
      bool default_is_unbalanced =
          res_default->at(0)->is<Tensor>() || res_default->at(1)->is<Tensor>();
      REQUIRE(default_is_unbalanced);
    }
  }

  /// verify that space changes did not leak
  auto reg_check = get_default_context().index_space_registry();
  auto uocc_check = reg_check->retrieve_ptr(L"a");
  REQUIRE(uocc_check);
  REQUIRE(uocc_check->approximate_size() == 10);
}
