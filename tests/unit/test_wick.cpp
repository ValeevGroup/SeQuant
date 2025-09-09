//
// Created by Eduard Valeyev on 3/23/18.
//

#include "./gwt.hpp"

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/latex.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/utility/debug.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/nodiscard.hpp>
#include <SeQuant/core/utility/timer.hpp>
#include <SeQuant/core/wick.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <catch2/catch_test_macros.hpp>
#include "catch2_sequant.hpp"

#include <range/v3/all.hpp>

#include <algorithm>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace sequant {
struct WickAccessor {};

template <>
template <>
struct WickTheorem<Statistics::FermiDirac>::access_by<WickAccessor> {
  auto compute_nontensor_wick(WickTheorem<Statistics::FermiDirac>& wick) {
    return wick.compute_nontensor_wick(false);
  }
};

auto compute_nontensor_wick(WickTheorem<Statistics::FermiDirac>& wick) {
  return WickTheorem<Statistics::FermiDirac>::access_by<WickAccessor>{}
      .compute_nontensor_wick(wick);
}

}  // namespace sequant

#if 1
TEST_CASE("wick", "[algorithms][wick]") {
  using namespace sequant;

  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());
  Index::reset_tmp_index();

  SECTION("Op contractions") {
    REQUIRE(FWickTheorem::can_contract(fann(L"i_1"), fcre(L"i_2"),
                                       Vacuum::Physical));
    REQUIRE(!FWickTheorem::can_contract(fcre(L"i_1"), fcre(L"i_2"),
                                        Vacuum::Physical));
    REQUIRE(!FWickTheorem::can_contract(fcre(L"i_1"), fann(L"i_2"),
                                        Vacuum::Physical));
    REQUIRE(!FWickTheorem::can_contract(fann(L"i_1"), fann(L"i_2"),
                                        Vacuum::Physical));

    REQUIRE(!FWickTheorem::can_contract(fann(L"i_1"), fcre(L"i_2"),
                                        Vacuum::SingleProduct));
    REQUIRE(!FWickTheorem::can_contract(fcre(L"i_1"), fcre(L"i_2"),
                                        Vacuum::SingleProduct));
    REQUIRE(FWickTheorem::can_contract(fcre(L"i_1"), fann(L"i_2"),
                                       Vacuum::SingleProduct));
    REQUIRE(!FWickTheorem::can_contract(fann(L"i_1"), fann(L"i_2"),
                                        Vacuum::SingleProduct));

    REQUIRE(!FWickTheorem::can_contract(fann(L"i_1"), fcre(L"a_2"),
                                        Vacuum::SingleProduct));
    REQUIRE(!FWickTheorem::can_contract(fcre(L"i_1"), fcre(L"a_2"),
                                        Vacuum::SingleProduct));
    REQUIRE(!FWickTheorem::can_contract(fcre(L"i_1"), fann(L"a_2"),
                                        Vacuum::SingleProduct));
    REQUIRE(!FWickTheorem::can_contract(fann(L"i_1"), fann(L"a_2"),
                                        Vacuum::SingleProduct));

    REQUIRE(!FWickTheorem::can_contract(fann(L"a_1"), fcre(L"i_2"),
                                        Vacuum::SingleProduct));
    REQUIRE(!FWickTheorem::can_contract(fcre(L"a_1"), fcre(L"i_2"),
                                        Vacuum::SingleProduct));
    REQUIRE(!FWickTheorem::can_contract(fcre(L"a_1"), fann(L"i_2"),
                                        Vacuum::SingleProduct));
    REQUIRE(!FWickTheorem::can_contract(fann(L"a_1"), fann(L"i_2"),
                                        Vacuum::SingleProduct));

    REQUIRE(FWickTheorem::can_contract(fann(L"a_1"), fcre(L"a_2"),
                                       Vacuum::SingleProduct));
    REQUIRE(!FWickTheorem::can_contract(fcre(L"a_1"), fcre(L"a_2"),
                                        Vacuum::SingleProduct));
    REQUIRE(!FWickTheorem::can_contract(fcre(L"a_1"), fann(L"a_2"),
                                        Vacuum::SingleProduct));
    REQUIRE(!FWickTheorem::can_contract(fann(L"a_1"), fann(L"a_2"),
                                        Vacuum::SingleProduct));

    REQUIRE(!FWickTheorem::can_contract(fann(L"p_1"), fcre(L"i_2"),
                                        Vacuum::SingleProduct));
    REQUIRE(!FWickTheorem::can_contract(fcre(L"p_1"), fcre(L"i_2"),
                                        Vacuum::SingleProduct));
    REQUIRE(FWickTheorem::can_contract(fcre(L"p_1"), fann(L"i_2"),
                                       Vacuum::SingleProduct));
    REQUIRE(!FWickTheorem::can_contract(fann(L"p_1"), fann(L"i_2"),
                                        Vacuum::SingleProduct));

    REQUIRE(BWickTheorem::can_contract(bann(L"i_1"), bcre(L"i_2"),
                                       Vacuum::Physical));
    REQUIRE(!BWickTheorem::can_contract(bcre(L"i_1"), bcre(L"i_2"),
                                        Vacuum::Physical));
    REQUIRE(!BWickTheorem::can_contract(bcre(L"i_1"), bann(L"i_2"),
                                        Vacuum::Physical));
    REQUIRE(!BWickTheorem::can_contract(bann(L"i_1"), bann(L"i_2"),
                                        Vacuum::Physical));
  }

  SECTION("constructors") {
    REQUIRE_NOTHROW(FWickTheorem{FNOperatorSeq{}});
    REQUIRE_NOTHROW(BWickTheorem{BNOperatorSeq{}});

    {
      auto opseq1 = ex<FNOperatorSeq>(FNOperator(cre({L"i_1"}), ann({L"i_2"})),
                                      FNOperator(cre({L"i_3"}), ann({L"i_4"})),
                                      FNOperator(cre({L"i_5"}), ann({L"i_6"})));
      REQUIRE_NOTHROW(FWickTheorem{opseq1});
      auto wick1 = FWickTheorem{opseq1};

      SEQUANT_PRAGMA_CLANG(diagnostic push)
      SEQUANT_PRAGMA_CLANG(diagnostic ignored "-Wdeprecated-declarations")
      SEQUANT_PRAGMA_GCC(diagnostic push)
      SEQUANT_PRAGMA_GCC(diagnostic ignored "-Wdeprecated-declarations")

      if (get_default_context().spbasis() == SPBasis::Spinor) {
        REQUIRE_NOTHROW(wick1.spinfree(false));
        REQUIRE_THROWS_AS(wick1.spinfree(true), std::invalid_argument);
      }
      if (get_default_context().spbasis() == SPBasis::Spinfree) {
        REQUIRE_NOTHROW(wick1.spinfree(true));
        REQUIRE_THROWS_AS(wick1.spinfree(false), std::invalid_argument);
      }

      SEQUANT_PRAGMA_GCC(diagnostic pop)
      SEQUANT_PRAGMA_CLANG(diagnostic pop)
    }

  }  // SECTION("constructors")

  SECTION("physical vacuum") {
    constexpr Vacuum V = Vacuum::Physical;
    auto ctx = get_default_context();
    ctx.set(sequant::mbpt::make_sr_spaces());
    ctx.set(V);
    auto raii_tmp = set_scoped_default_context(ctx);

    auto switch_to_spinfree_context = detail::NoDiscard([&]() {
      auto context_sf = get_default_context();
      context_sf.set(SPBasis::Spinfree);
      return set_scoped_default_context(context_sf);
    });

    // number operator
    {
      {
        auto opseq1 = ex<FNOperatorSeq>(FNOperator(cre({L"i_1"}), ann()),
                                        FNOperator(cre(), ann({L"i_2"})));
        auto wick1 = FWickTheorem{opseq1};
        REQUIRE_NOTHROW(wick1.compute());
        // full contractions = null (N is already in normal form)
        auto full_contractions = FWickTheorem{opseq1}.compute();
        REQUIRE(full_contractions->is<Constant>());
        REQUIRE(full_contractions->as<Constant>().value<int>() == 0);
        // partial contractions = N
        auto partial_contractions =
            FWickTheorem{opseq1}.full_contractions(false).compute();
        // std::wcout << "partial_contractions=" <<
        // to_latex(partial_contractions)
        // << std::endl;
        REQUIRE(partial_contractions->is<Product>());
        REQUIRE(partial_contractions->as<Product>().size() == 1);
      }
      {
        auto opseq1 = ex<BNOperatorSeq>(BNOperator(cre({L"i_1"}), ann()),
                                        BNOperator(cre(), ann({L"i_2"})));
        auto wick1 = BWickTheorem{opseq1};
        REQUIRE_NOTHROW(wick1.compute());
        // full contractions = null
        auto full_contractions = BWickTheorem{opseq1}.compute();
        REQUIRE(full_contractions->is<Constant>());
        REQUIRE(full_contractions->as<Constant>().value<int>() == 0);
        // partial contractions = N
        auto partial_contractions =
            BWickTheorem{opseq1}.full_contractions(false).compute();
        // std::wcout << "partial_contractions=" <<
        // to_latex(partial_contractions)
        // << std::endl;
        REQUIRE(partial_contractions->is<Product>());
        REQUIRE(partial_contractions->as<Product>().size() == 1);
      }
    }

    // hole number operator
    {
      {
        auto opseq1 = ex<FNOperatorSeq>(FNOperator(cre({}), ann({L"i_1"})),
                                        FNOperator(cre({L"i_2"}), ann({})));
        auto wick1 = FWickTheorem{opseq1};
        REQUIRE_NOTHROW(wick1.compute());
        // full contractions = delta
        auto full_contractions = FWickTheorem{opseq1}.compute();
        REQUIRE(full_contractions->is<Product>());
        REQUIRE(full_contractions->as<Product>().size() == 1);
        // partial contractions = delta - N
        auto partial_contractions =
            FWickTheorem{opseq1}.full_contractions(false).compute();
        // std::wcout << "partial_contractions=" <<
        // to_latex(partial_contractions) << std::endl;
        REQUIRE(partial_contractions->is<Sum>());
        REQUIRE(partial_contractions->as<Sum>().size() == 2);
        REQUIRE(
            to_latex(partial_contractions) ==
            L"{ \\bigl({{s^{{i_2}}_{{i_1}}}} - {{a^{{i_2}}_{{i_1}}}}\\bigr) }");
      }
      {
        auto opseq1 = ex<BNOperatorSeq>(BNOperator(cre({}), ann({L"i_1"})),
                                        BNOperator(cre({L"i_2"}), ann({})));
        auto wick1 = BWickTheorem{opseq1};
        REQUIRE_NOTHROW(wick1.compute());
        // full contractions = delta
        auto full_contractions = BWickTheorem{opseq1}.compute();
        REQUIRE(full_contractions->is<Product>());
        REQUIRE(full_contractions->as<Product>().size() == 1);
        // partial contractions = delta + N
        auto partial_contractions =
            BWickTheorem{opseq1}.full_contractions(false).compute();
        // std::wcout << "partial_contractions=" <<
        // to_latex(partial_contractions) << std::endl;
        REQUIRE(partial_contractions->is<Sum>());
        REQUIRE(partial_contractions->as<Sum>().size() == 2);
        REQUIRE(
            to_latex(partial_contractions) ==
            L"{ \\bigl({{s^{{i_2}}_{{i_1}}}} + {{b^{{i_2}}_{{i_1}}}}\\bigr) }");
      }
    }

    // general string of creators/annihilators, from Nick Mayhall
    {
      // sequence of individual ops
      const auto ops = {fann("i_1"), fcre("i_2"), fcre("i_3"),
                        fann("i_4"), fann("i_5"), fcre("i_1")};
      REQUIRE_NOTHROW(FWickTheorem(ops));
      FWickTheorem w(ops);
      REQUIRE_NOTHROW(w.full_contractions(false).compute());
      auto partial_contractions =
          FWickTheorem{ops}.full_contractions(false).compute();
      REQUIRE(to_latex(partial_contractions) ==
              L"{ \\bigl( - {{s^{{i_2}}_{{i_1}}}{a^{{i_3}{i_1}}_{{i_4}{i_5}}}} "
              L"- {{s^{{i_2}}_{{i_1}}}{s^{{i_1}}_{{i_4}}}{a^{{i_3}}_{{i_5}}}} "
              L"+ {{s^{{i_2}}_{{i_1}}}{s^{{i_1}}_{{i_5}}}{a^{{i_3}}_{{i_4}}}} "
              L"+ {{s^{{i_3}}_{{i_1}}}{a^{{i_2}{i_1}}_{{i_4}{i_5}}}} + "
              L"{{s^{{i_3}}_{{i_1}}}{s^{{i_1}}_{{i_4}}}{a^{{i_2}}_{{i_5}}}} - "
              L"{{s^{{i_3}}_{{i_1}}}{s^{{i_1}}_{{i_5}}}{a^{{i_2}}_{{i_4}}}} - "
              L"{{s^{{i_1}}_{{i_1}}}{a^{{i_2}{i_3}}_{{i_4}{i_5}}}} + "
              L"{{s^{{i_1}}_{{i_4}}}{a^{{i_2}{i_3}}_{{i_1}{i_5}}}} - "
              L"{{s^{{i_1}}_{{i_5}}}{a^{{i_2}{i_3}}_{{i_1}{i_4}}}} + "
              L"{{a^{{i_2}{i_3}{i_1}}_{{i_1}{i_4}{i_5}}}}\\bigr) }");
    }

    // three 1-body operators
    {
      auto opseq1 = ex<FNOperatorSeq>(FNOperator(cre({L"i_1"}), ann({L"i_2"})),
                                      FNOperator(cre({L"i_3"}), ann({L"i_4"})),
                                      FNOperator(cre({L"i_5"}), ann({L"i_6"})));
      auto wick1 = FWickTheorem{opseq1};
      REQUIRE_NOTHROW(wick1.compute());
      auto result = FWickTheorem{opseq1}.compute();
      REQUIRE(result->is<Constant>());
      REQUIRE(result->as<Constant>().value<int>() == 0);
    }

    // two 2-body operators
    {
      auto opseq =
          ex<FNOperatorSeq>(FNOperator(cre({}), ann({L"i_1", L"i_2"})),
                            FNOperator(cre({L"i_3", L"i_4"}), ann({})));
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.compute());
      auto result = wick.compute();
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 2);
    }

    // two 3-body operators
    {
      auto opseq =
          ex<FNOperatorSeq>(FNOperator(cre({}), ann({L"i_1", L"i_2", L"i_3"})),
                            FNOperator(cre({L"i_4", L"i_5", L"i_6"}), ann({})));
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.compute());
      auto result = wick.compute();
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 6);
    }

    // two 4-body operators
    {
      auto opseq = ex<FNOperatorSeq>(
          FNOperator(cre({}), ann({L"i_1", L"i_2", L"i_3", L"i_4"})),
          FNOperator(cre({L"i_5", L"i_6", L"i_7", L"i_8"}), ann({})));
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.compute());
      auto result = wick.compute();
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 24);
    }

    // 1/2 * 1 * 1/2 body ops, full contraction
    {
      auto opseq = ex<FNOperatorSeq>(FNOperator(cre({}), ann({L"i_1"})),
                                     FNOperator(cre({L"i_2"}), ann({L"i_3"})),
                                     FNOperator(cre({L"i_4"}), ann({})));
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.compute());
      auto result = wick.compute();
      REQUIRE(result->is<Product>());
      REQUIRE(result->size() == 2);  // product of 2 terms
    }

    // 1/2 * 1 * 1/2 body ops, partial contraction
    {
      auto opseq = ex<FNOperatorSeq>(FNOperator(cre({}), ann({L"i_1"})),
                                     FNOperator(cre({L"i_2"}), ann({L"i_3"})),
                                     FNOperator(cre({L"i_4"}), ann({})));
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.full_contractions(false).compute());
      auto result = wick.full_contractions(false).compute();
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 5);  // sum of 4 terms
      REQUIRE(to_latex(result) ==
              L"{ \\bigl( - {{s^{{i_2}}_{{i_1}}}{a^{{i_4}}_{{i_3}}}} + "
              L"{{s^{{i_2}}_{{i_1}}}{s^{{i_4}}_{{i_3}}}} + "
              L"{{s^{{i_4}}_{{i_1}}}{a^{{i_2}}_{{i_3}}}} - "
              L"{{s^{{i_4}}_{{i_3}}}{a^{{i_2}}_{{i_1}}}} - "
              L"{{a^{{i_2}{i_4}}_{{i_3}{i_1}}}}\\bigr) }");
    }

    // three 1-body operators, partial contraction
    {
      auto opseq1 = ex<FNOperatorSeq>(FNOperator(cre({L"i_1"}), ann({L"i_2"})),
                                      FNOperator(cre({L"i_3"}), ann({L"i_4"})),
                                      FNOperator(cre({L"i_5"}), ann({L"i_6"})));
      auto wick1 = FWickTheorem{opseq1};
      REQUIRE_NOTHROW(wick1.full_contractions(false).compute());
      auto result = FWickTheorem{opseq1}.full_contractions(false).compute();
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 5);
      REQUIRE(to_latex(result) ==
              L"{ \\bigl({{s^{{i_3}}_{{i_2}}}{a^{{i_1}{i_5}}_{{i_4}{i_6}}}} + "
              L"{{s^{{i_3}}_{{i_2}}}{s^{{i_5}}_{{i_4}}}{a^{{i_1}}_{{i_6}}}} + "
              L"{{s^{{i_5}}_{{i_2}}}{a^{{i_1}{i_3}}_{{i_6}{i_4}}}} + "
              L"{{s^{{i_5}}_{{i_4}}}{a^{{i_1}{i_3}}_{{i_2}{i_6}}}} + "
              L"{{a^{{i_1}{i_3}{i_5}}_{{i_2}{i_4}{i_6}}}}\\bigr) }");
    }

    // two 2-body operators, partial contraction: Eq. 9b of DOI 10.1063/1.474405
    {
      auto opseq = ex<FNOperatorSeq>(
          FNOperator(cre({L"i_1", L"i_2"}), ann({L"i_3", L"i_4"})),
          FNOperator(cre({L"i_5", L"i_6"}), ann({L"i_7", L"i_8"})));
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.full_contractions(false).compute());
      auto result = wick.full_contractions(false).compute();
      auto result_latex = to_latex(result);
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 7);
      REQUIRE(
          result_latex ==
          L"{ "
          L"\\bigl({{s^{{i_5}}_{{i_4}}}{a^{{i_1}{i_2}{i_6}}_{{i_3}{i_7}{i_8}}}}"
          L" + "
          L"{{s^{{i_5}}_{{i_4}}}{s^{{i_6}}_{{i_3}}}{a^{{i_1}{i_2}}_{{i_8}{i_7}}"
          L"}} + {{s^{{i_6}}_{{i_4}}}{a^{{i_1}{i_2}{i_5}}_{{i_3}{i_8}{i_7}}}} "
          L"+ "
          L"{{s^{{i_6}}_{{i_4}}}{s^{{i_5}}_{{i_3}}}{a^{{i_1}{i_2}}_{{i_7}{i_8}}"
          L"}} + {{s^{{i_5}}_{{i_3}}}{a^{{i_1}{i_2}{i_6}}_{{i_7}{i_4}{i_8}}}} "
          L"+ {{s^{{i_6}}_{{i_3}}}{a^{{i_1}{i_2}{i_5}}_{{i_8}{i_4}{i_7}}}} + "
          L"{{a^{{i_1}{i_2}{i_5}{i_6}}_{{i_3}{i_4}{i_7}{i_8}}}}\\bigr) }");

      // we use this example in the paper, both as bare ops and contracted with
      // scalars to make all indices dummy and the topological optimizations
      // kicking in ... see if the paper's exlanation of how
      // topology works is correct
      {
        auto expr =
            ex<Tensor>(L"g", bra{L"i_1", L"i_2"}, ket{L"i_3", L"i_4"},
                       Symmetry::Antisymm) *
            ex<FNOperator>(cre({L"i_1", L"i_2"}), ann({L"i_3", L"i_4"})) *
            ex<Tensor>(L"g", bra{L"i_5", L"i_6"}, ket{L"i_7", L"i_8"},
                       Symmetry::Antisymm) *
            ex<FNOperator>(cre({L"i_5", L"i_6"}), ann({L"i_7", L"i_8"}));
        auto wick = FWickTheorem{expr};

        // first with topology utilization
        REQUIRE_NOTHROW(
            wick.full_contractions(false).use_topology(true).compute());
        wick.stats().reset();
        auto result =
            wick.full_contractions(false).use_topology(true).compute();
        result->rapid_canonicalize();
        REQUIRE(result->is<Sum>());
        REQUIRE(result->size() == 3);
        REQUIRE(wick.stats().num_attempted_contractions == 2);
        REQUIRE(wick.stats().num_useful_contractions == 3);
        auto result_latex_w_topology = to_latex(result);

        // now without topology
        wick = FWickTheorem{expr};
        REQUIRE_NOTHROW(
            wick.full_contractions(false).use_topology(false).compute());
        wick.stats().reset();
        result = wick.full_contractions(false).use_topology(false).compute();
        result->rapid_canonicalize();
        REQUIRE(result->is<Sum>());
        REQUIRE(result->size() == 3);
        REQUIRE(wick.stats().num_attempted_contractions == 6);
        REQUIRE(wick.stats().num_useful_contractions == 8);
        auto result_latex_wo_topology = to_latex(result);
        REQUIRE(result_latex_wo_topology == result_latex_w_topology);
      }

      // if Wick's theorem's result is in "canonical" (columns-matching-inputs
      // ... this is what Kutzelnigg calls generalized Wick's theorem) it works
      // same for spinor and spinfree basis for physical vacuum
      auto raii_tmp = switch_to_spinfree_context();
      REQUIRE_NOTHROW(wick.full_contractions(false).compute());
      auto result_sf = wick.full_contractions(false).compute();
      auto result_sf_latex = to_latex(result_sf);
      // std::wcout << "result_sf: " << to_latex(result_sf) << std::endl;
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 7);
      auto result_sf_latex_with_E_replaced_by_a =
          result_sf_latex | ranges::views::replace(L'E', L'a') |
          ranges::to<std::wstring>();
      REQUIRE(result_latex == result_sf_latex_with_E_replaced_by_a);
    }

  }  // SECTION("physical vacuum")

  SECTION("fermi vacuum") {
    // default vacuum is already spin-orbital Fermi vacuum

    auto switch_to_spinfree_context = detail::NoDiscard([&]() {
      auto context_sf = get_default_context();
      context_sf.set(SPBasis::Spinfree);
      return set_scoped_default_context(context_sf);
    });

    // two (pure qp) 1-body operators
    {
      auto opseq = ex<FNOperatorSeq>(FNOperator(cre({L"i_1"}), ann({L"a_1"})),
                                     FNOperator(cre({L"a_2"}), ann({L"i_2"})));
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.compute());
      auto result = wick.compute();
      REQUIRE(result->is<Product>());
      REQUIRE(result->size() == 2);  // product of 2 terms

      // spin-free result is simply twice the spin-orbital result
      auto raii_tmp = switch_to_spinfree_context();
      auto result_sf = wick.compute();
      REQUIRE(simplify(result_sf - result * ex<Constant>(2)) ==
              ex<Constant>(0));
    }

    // two (pure qp) N-nonconserving 2-body operators
    {
      auto opseq =
          ex<FNOperatorSeq>(FNOperator(cre({L"i_1", L"i_2"}), ann({L"a_1"})),
                            FNOperator(cre({L"a_2"}), ann({L"i_3", L"i_4"})));
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.compute());
      auto result = wick.compute();
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 2);
    }

    // two general 1-body operators
    {
      auto opseq = ex<FNOperatorSeq>(FNOperator(cre({L"p_1"}), ann({L"p_2"})),
                                     FNOperator(cre({L"p_3"}), ann({L"p_4"})));
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.compute());
      auto result = wick.compute();
      REQUIRE(result->is<Product>());
      REQUIRE(to_latex(result) ==
              L"{{s^{{m_{104}}}_{{m_{105}}}}{\\delta^{{m_{105}}}_{{p_4}}}{"
              L"\\delta^{{p_1}}_{{m_{104}}}}{s^{{e_{107}}}_{{e_{106}}}}{"
              L"\\delta^{{e_{106}}}_{{p_2}}}{\\delta^{{p_3}}_{{e_{107}}}}}");
    }
    // two general 1-body operators, partial contractions: Eq. 21a of
    // DOI 10.1063/1.474405
    {
      auto opseq = ex<FNOperatorSeq>(FNOperator(cre({L"p_1"}), ann({L"p_2"})),
                                     FNOperator(cre({L"p_3"}), ann({L"p_4"})));
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.full_contractions(false).compute());
      auto result = wick.full_contractions(false).compute();
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 4);
      REQUIRE(
          4 ==
          GWT({1, 1}, /* full_contractions_only = */ false).result().size());
      REQUIRE(to_latex(result) ==
              L"{ \\bigl( - "
              L"{{s^{{m_{114}}}_{{m_{115}}}}{\\delta^{{m_{115}}}_{{p_4}}}{"
              L"\\delta^{{p_1}}_{{m_{114}}}}{\\tilde{a}^{{p_3}}_{{p_2}}}} + "
              L"{{s^{{m_{114}}}_{{m_{115}}}}{\\delta^{{m_{115}}}_{{p_4}}}{"
              L"\\delta^{{p_1}}_{{m_{114}}}}{s^{{e_{117}}}_{{e_{116}}}}{"
              L"\\delta^{{e_{116}}}_{{p_2}}}{\\delta^{{p_3}}_{{e_{117}}}}} + "
              L"{{s^{{e_{119}}}_{{e_{118}}}}{\\delta^{{e_{118}}}_{{p_2}}}{"
              L"\\delta^{{p_3}}_{{e_{119}}}}{\\tilde{a}^{{p_1}}_{{p_4}}}} + "
              L"{{\\tilde{a}^{{p_1}{p_3}}_{{p_2}{p_4}}}}\\bigr) }");
    }

    // two (pure qp) 2-body operators
    {
      auto opseq = ex<FNOperatorSeq>(
          FNOperator(cre({L"i_1", L"i_2"}), ann({L"a_1", L"a_2"})),
          FNOperator(cre({L"a_3", L"a_4"}), ann({L"i_3", L"i_4"})));
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.compute());
      auto result = wick.compute();
      auto result_latex = to_latex(result);
      // std::wcout << "<" << to_latex(opseq) << "> = " << result_latex <<
      // std::endl;
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 4);
      REQUIRE(result_latex ==
              L"{ "
              L"\\bigl({{s^{{i_1}}_{{i_4}}}{s^{{i_2}}_{{i_3}}}{s^{{a_3}}_{{a_2}"
              L"}}{s^{{a_4}}_{{a_1}}}} - "
              L"{{s^{{i_1}}_{{i_4}}}{s^{{i_2}}_{{i_3}}}{s^{{a_4}}_{{a_2}}}{s^{{"
              L"a_3}}_{{a_1}}}} - "
              L"{{s^{{i_1}}_{{i_3}}}{s^{{i_2}}_{{i_4}}}{s^{{a_3}}_{{a_2}}}{s^{{"
              L"a_4}}_{{a_1}}}} + "
              L"{{s^{{i_1}}_{{i_3}}}{s^{{i_2}}_{{i_4}}}{s^{{a_4}}_{{a_2}}}{s^{{"
              L"a_3}}_{{a_1}}}}\\bigr) }");

      // in spin-free result first and fourth terms are multiplied by 4, second
      // and third terms multiplied by 2
      auto raii_tmp = switch_to_spinfree_context();
      auto result_sf = wick.compute();
      auto result_sf_latex = to_latex(result_sf);
      //     std::wcout << "<" << to_latex(opseq) << "> = " <<
      //     to_latex(result_sf)
      //     << std::endl;
      REQUIRE(result_sf->is<Sum>());
      REQUIRE(result_sf->size() == 4);
      REQUIRE(result_sf_latex ==
              L"{ "
              L"\\bigl({{{4}}{s^{{i_1}}_{{i_4}}}{s^{{i_2}}_{{i_3}}}{s^{{a_3}}_{"
              L"{a_2}}}{s^{{a_4}}_{{a_1}}}} - "
              L"{{{2}}{s^{{i_1}}_{{i_4}}}{s^{{i_2}}_{{i_3}}}{s^{{a_4}}_{{a_2}}}"
              L"{s^{{a_3}}_{{a_1}}}} - "
              L"{{{2}}{s^{{i_1}}_{{i_3}}}{s^{{i_2}}_{{i_4}}}{s^{{a_3}}_{{a_2}}}"
              L"{s^{{a_4}}_{{a_1}}}} + "
              L"{{{4}}{s^{{i_1}}_{{i_3}}}{s^{{i_2}}_{{i_4}}}{s^{{a_4}}_{{a_2}}}"
              L"{s^{{a_3}}_{{a_1}}}}\\bigr) }");
    }
    // two (pure qp) 3-body operators
    {
      auto opseq = ex<FNOperatorSeq>(FNOperator(cre({L"i_1", L"i_2", L"i_3"}),
                                                ann({L"a_1", L"a_2", L"a_3"})),
                                     FNOperator(cre({L"a_4", L"a_5", L"a_6"}),
                                                ann({L"i_4", L"i_5", L"i_6"})));
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.compute());
      auto result = wick.compute();
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 36);
    }

    // one general 1-body operator + one general 2-body operator, partial
    // contraction: Eq. 9 of DOI 10.1063/1.474405
    {
      auto opseq = ex<FNOperatorSeq>(
          FNOperator(cre({L"p_1"}), ann({L"p_2"})),
          FNOperator(cre({L"p_3", L"p_4"}), ann({L"p_5", L"p_6"})));
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.full_contractions(false).compute());
      auto result = wick.full_contractions(false).compute();
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 9);
      REQUIRE(
          9 ==
          GWT({1, 2}, /* full_contractions_only = */ false).result().size());
    }

    // two general 2-body operators
    {
      auto opseq = ex<FNOperatorSeq>(
          FNOperator(cre({L"p_1", L"p_2"}), ann({L"p_3", L"p_4"})),
          FNOperator(cre({L"p_5", L"p_6"}), ann({L"p_7", L"p_8"})));
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.compute());
      auto result = wick.compute();
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 4);
      REQUIRE(4 == GWT({2, 2}).result().size());
    }
    // two general 2-body operators, partial contractions: Eqs. 22 of
    // DOI 10.1063/1.474405
    {
      auto opseq = ex<FNOperatorSeq>(
          FNOperator(cre({L"p_1", L"p_2"}), ann({L"p_3", L"p_4"})),
          FNOperator(cre({L"p_5", L"p_6"}), ann({L"p_7", L"p_8"})));
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.full_contractions(false).compute());
      auto result = wick.full_contractions(false).compute();
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 49);  // the MK paper only gives 47 terms,
                                      // misses the 2 double-hole contractions
      REQUIRE(
          49 ==
          GWT({2, 2}, /* full_contractions_only = */ false).result().size());
    }
    // one general 2-body operator and one 2-body excitation operator
    {
      auto opseq = ex<FNOperatorSeq>(
          FNOperator(cre({L"p_1", L"p_2"}), ann({L"p_3", L"p_4"})),
          FNOperator(cre({L"a_3", L"a_4"}), ann({L"i_3", L"i_4"})));
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.compute());
      auto result = wick.compute();
      std::wcout << "<" << to_latex(opseq) << "> = " << to_latex(result)
                 << std::endl;
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 4);
      REQUIRE(
          to_latex(result) ==
          L"{ "
          L"\\bigl({{s^{{p_1}}_{{i_4}}}{s^{{p_2}}_{{i_3}}}{s^{{a_3}}_{{p_4}}}"
          L"{s^{{a_4}}_{{p_3}}}} - "
          L"{{s^{{p_1}}_{{i_4}}}{s^{{p_2}}_{{i_3}}}{s^{{a_4}}_{{p_4}}}{s^{{a_"
          L"3}}_{{p_3}}}} - "
          L"{{s^{{p_1}}_{{i_3}}}{s^{{p_2}}_{{i_4}}}{s^{{a_3}}_{{p_4}}}{s^{{a_"
          L"4}}_{{p_3}}}} + "
          L"{{s^{{p_1}}_{{i_3}}}{s^{{p_2}}_{{i_4}}}{s^{{a_4}}_{{p_4}}}{s^{{a_"
          L"3}}_{{p_3}}}}\\bigr) }");

      // in spin-free result first and fourth terms are multiplied by 4, second
      // and third terms multiplied by 2
      auto raii_tmp = switch_to_spinfree_context();
      auto result_sf = wick.compute();
      auto result_sf_latex = to_latex(result_sf);
      //     std::wcout << "<" << to_latex(opseq) << "> = " <<
      //     to_latex(result_sf)
      //     << std::endl;
      REQUIRE(result_sf->is<Sum>());
      REQUIRE(result_sf->size() == 4);
      REQUIRE(
          to_latex(result_sf) ==
          L"{ "
          L"\\bigl({{{4}}{s^{{p_1}}_{{i_4}}}{s^{{p_2}}_{{i_3}}}{s^{{a_3}}_{{"
          L"p_4}}}{s^{{a_4}}_{{p_3}}}} - "
          L"{{{2}}{s^{{p_1}}_{{i_4}}}{s^{{p_2}}_{{i_3}}}{s^{{a_4}}_{{p_4}}}{"
          L"s^{{a_3}}_{{p_3}}}} - "
          L"{{{2}}{s^{{p_1}}_{{i_3}}}{s^{{p_2}}_{{i_4}}}{s^{{a_3}}_{{p_4}}}{"
          L"s^{{a_4}}_{{p_3}}}} + "
          L"{{{4}}{s^{{p_1}}_{{i_3}}}{s^{{p_2}}_{{i_4}}}{s^{{a_4}}_{{p_4}}}{"
          L"s^{{a_3}}_{{p_3}}}}\\bigr) }");
    }

    // two general 3-body operators
    {
      auto opseq =
          ex<FNOperatorSeq>(FNOperator(cre({L"p_1", L"p_2", L"p_3"}),
                                       ann({L"p_4", L"p_5", L"p_6"})),
                            FNOperator(cre({L"p_7", L"p_8", L"p_9"}),
                                       ann({L"p_10", L"p_11", L"p_12"})));
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.compute());
      auto result = wick.compute();
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 36);
      REQUIRE(36 == GWT({3, 3}).result().size());
    }

    // two N-nonconserving operators
    {
      auto opseq = ex<FNOperatorSeq>(
          FNOperator(cre({L"p_1", L"p_2", L"p_3"}), ann({L"p_4", L"p_5"})),
          FNOperator(cre({L"p_7", L"p_8"}), ann({L"p_10", L"p_11", L"p_12"})));
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.compute());
      auto result = wick.compute();
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 12);
    }

    // more N-nonconserving operators
    {
      auto input =
          ex<FNOperator>(cre({L"i_1"}), ann(L"a_3", L"a_4")) *
          (ex<Constant>(rational{1, 4}) *
           ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"},
                      Symmetry::Antisymm) *
           ex<FNOperator>(cre({L"p_1", L"p_2"}), ann({L"p_3", L"p_4"}))) *
          ex<FNOperator>(cre({L"a_2"}), ann({}));
      auto wick = FWickTheorem{input};
      wick.set_external_indices(IndexList{L"i_1", L"a_3", L"a_4", L"a_2"})
          .use_topology(true);
      ExprPtr result;
      REQUIRE_NOTHROW(result = wick.compute());
      // std::wcout << "result = " << to_latex(result) << std::endl;
      REQUIRE(to_latex(result) == L"{{-}{\\bar{g}^{{a_2}{i_1}}_{{a_4}{a_3}}}}");
      canonicalize(result, {.method = CanonicalizationMethod::Rapid});
      REQUIRE(to_latex(result) == L"{{-}{\\bar{g}^{{i_1}{a_2}}_{{a_3}{a_4}}}}");
    }

    // general string of creators/annihilators, from Nick Mayhall
    {
      // sequence of individual ops
      const auto ops = {fann("p_1"), fcre("p_2"), fcre("p_3"),
                        fann("p_4"), fann("p_5"), fcre("p_1")};
      REQUIRE_NOTHROW(FWickTheorem(ops));
      FWickTheorem w(ops);
      REQUIRE_NOTHROW(w.full_contractions(false).compute());
      auto partial_contractions =
          FWickTheorem{ops}.full_contractions(false).compute();
      REQUIRE(partial_contractions.is<Sum>());
      REQUIRE(partial_contractions.as<Sum>().size() == 34);
      // REQUIRE(
      //       to_latex(partial_contractions) ==
      //       L"{ \\bigl( -
      //       {{s^{{e_{133}}}_{{p_1}}}{s^{{p_2}}_{{e_{133}}}}{\\tilde{a}^{{p_3}{p_1}}_{{p_4}{p_5}}}}
      //       -
      //       {{s^{{e_{133}}}_{{p_1}}}{s^{{p_2}}_{{e_{133}}}}{s^{{p_3}}_{{m_{134}}}}{s^{{m_{134}}}_{{p_4}}}{\\tilde{a}^{{p_1}}_{{p_5}}}}
      //       +
      //       {{s^{{e_{133}}}_{{p_1}}}{s^{{p_2}}_{{e_{133}}}}{s^{{p_3}}_{{m_{134}}}}{s^{{m_{134}}}_{{p_4}}}{s^{{e_{135}}}_{{p_5}}}{s^{{p_1}}_{{e_{135}}}}}
      //       +
      //       {{s^{{e_{133}}}_{{p_1}}}{s^{{p_2}}_{{e_{133}}}}{s^{{p_3}}_{{m_{136}}}}{s^{{m_{136}}}_{{p_5}}}{\\tilde{a}^{{p_1}}_{{p_4}}}}
      //       -
      //       {{s^{{e_{133}}}_{{p_1}}}{s^{{p_2}}_{{e_{133}}}}{s^{{p_3}}_{{m_{136}}}}{s^{{m_{136}}}_{{p_5}}}{s^{{e_{137}}}_{{p_4}}}{s^{{p_1}}_{{e_{137}}}}}
      //       -
      //       {{s^{{e_{133}}}_{{p_1}}}{s^{{p_2}}_{{e_{133}}}}{s^{{e_{138}}}_{{p_4}}}{s^{{p_1}}_{{e_{138}}}}{\\tilde{a}^{{p_3}}_{{p_5}}}}
      //       +
      //       {{s^{{e_{133}}}_{{p_1}}}{s^{{p_2}}_{{e_{133}}}}{s^{{e_{139}}}_{{p_5}}}{s^{{p_1}}_{{e_{139}}}}{\\tilde{a}^{{p_3}}_{{p_4}}}}
      //       +
      //       {{s^{{e_{140}}}_{{p_1}}}{s^{{p_3}}_{{e_{140}}}}{\\tilde{a}^{{p_2}{p_1}}_{{p_4}{p_5}}}}
      //       +
      //       {{s^{{e_{140}}}_{{p_1}}}{s^{{p_3}}_{{e_{140}}}}{s^{{p_2}}_{{m_{141}}}}{s^{{m_{141}}}_{{p_4}}}{\\tilde{a}^{{p_1}}_{{p_5}}}}
      //       -
      //       {{s^{{e_{140}}}_{{p_1}}}{s^{{p_3}}_{{e_{140}}}}{s^{{p_2}}_{{m_{141}}}}{s^{{m_{141}}}_{{p_4}}}{s^{{e_{142}}}_{{p_5}}}{s^{{p_1}}_{{e_{142}}}}}
      //       -
      //       {{s^{{e_{140}}}_{{p_1}}}{s^{{p_3}}_{{e_{140}}}}{s^{{p_2}}_{{m_{143}}}}{s^{{m_{143}}}_{{p_5}}}{\\tilde{a}^{{p_1}}_{{p_4}}}}
      //       +
      //       {{s^{{e_{140}}}_{{p_1}}}{s^{{p_3}}_{{e_{140}}}}{s^{{p_2}}_{{m_{143}}}}{s^{{m_{143}}}_{{p_5}}}{s^{{e_{144}}}_{{p_4}}}{s^{{p_1}}_{{e_{144}}}}}
      //       +
      //       {{s^{{e_{140}}}_{{p_1}}}{s^{{p_3}}_{{e_{140}}}}{s^{{e_{145}}}_{{p_4}}}{s^{{p_1}}_{{e_{145}}}}{\\tilde{a}^{{p_2}}_{{p_5}}}}
      //       -
      //       {{s^{{e_{140}}}_{{p_1}}}{s^{{p_3}}_{{e_{140}}}}{s^{{e_{146}}}_{{p_5}}}{s^{{p_1}}_{{e_{146}}}}{\\tilde{a}^{{p_2}}_{{p_4}}}}
      //       -
      //       {{s^{{e_{147}}}_{{p_1}}}{s^{{p_1}}_{{e_{147}}}}{\\tilde{a}^{{p_2}{p_3}}_{{p_4}{p_5}}}}
      //       -
      //       {{s^{{e_{147}}}_{{p_1}}}{s^{{p_1}}_{{e_{147}}}}{s^{{p_2}}_{{m_{148}}}}{s^{{m_{148}}}_{{p_4}}}{\\tilde{a}^{{p_3}}_{{p_5}}}}
      //       -
      //       {{s^{{e_{147}}}_{{p_1}}}{s^{{p_1}}_{{e_{147}}}}{s^{{p_2}}_{{m_{148}}}}{s^{{m_{148}}}_{{p_4}}}{s^{{p_3}}_{{m_{149}}}}{s^{{m_{149}}}_{{p_5}}}}
      //       +
      //       {{s^{{e_{147}}}_{{p_1}}}{s^{{p_1}}_{{e_{147}}}}{s^{{p_2}}_{{m_{150}}}}{s^{{m_{150}}}_{{p_5}}}{\\tilde{a}^{{p_3}}_{{p_4}}}}
      //       +
      //       {{s^{{e_{147}}}_{{p_1}}}{s^{{p_1}}_{{e_{147}}}}{s^{{p_2}}_{{m_{150}}}}{s^{{m_{150}}}_{{p_5}}}{s^{{p_3}}_{{m_{151}}}}{s^{{m_{151}}}_{{p_4}}}}
      //       +
      //       {{s^{{e_{147}}}_{{p_1}}}{s^{{p_1}}_{{e_{147}}}}{s^{{p_3}}_{{m_{152}}}}{s^{{m_{152}}}_{{p_4}}}{\\tilde{a}^{{p_2}}_{{p_5}}}}
      //       -
      //       {{s^{{e_{147}}}_{{p_1}}}{s^{{p_1}}_{{e_{147}}}}{s^{{p_3}}_{{m_{153}}}}{s^{{m_{153}}}_{{p_5}}}{\\tilde{a}^{{p_2}}_{{p_4}}}}
      //       -
      //       {{s^{{p_2}}_{{m_{154}}}}{s^{{m_{154}}}_{{p_4}}}{\\tilde{a}^{{p_3}{p_1}}_{{p_1}{p_5}}}}
      //       +
      //       {{s^{{p_2}}_{{m_{154}}}}{s^{{m_{154}}}_{{p_4}}}{s^{{p_3}}_{{m_{155}}}}{s^{{m_{155}}}_{{p_5}}}{\\tilde{a}^{{p_1}}_{{p_1}}}}
      //       +
      //       {{s^{{p_2}}_{{m_{154}}}}{s^{{m_{154}}}_{{p_4}}}{s^{{e_{156}}}_{{p_5}}}{s^{{p_1}}_{{e_{156}}}}{\\tilde{a}^{{p_3}}_{{p_1}}}}
      //       +
      //       {{s^{{p_2}}_{{m_{157}}}}{s^{{m_{157}}}_{{p_5}}}{\\tilde{a}^{{p_3}{p_1}}_{{p_1}{p_4}}}}
      //       -
      //       {{s^{{p_2}}_{{m_{157}}}}{s^{{m_{157}}}_{{p_5}}}{s^{{p_3}}_{{m_{158}}}}{s^{{m_{158}}}_{{p_4}}}{\\tilde{a}^{{p_1}}_{{p_1}}}}
      //       -
      //       {{s^{{p_2}}_{{m_{157}}}}{s^{{m_{157}}}_{{p_5}}}{s^{{e_{159}}}_{{p_4}}}{s^{{p_1}}_{{e_{159}}}}{\\tilde{a}^{{p_3}}_{{p_1}}}}
      //       +
      //       {{s^{{p_3}}_{{m_{160}}}}{s^{{m_{160}}}_{{p_4}}}{\\tilde{a}^{{p_2}{p_1}}_{{p_1}{p_5}}}}
      //       -
      //       {{s^{{p_3}}_{{m_{160}}}}{s^{{m_{160}}}_{{p_4}}}{s^{{e_{161}}}_{{p_5}}}{s^{{p_1}}_{{e_{161}}}}{\\tilde{a}^{{p_2}}_{{p_1}}}}
      //       -
      //       {{s^{{p_3}}_{{m_{162}}}}{s^{{m_{162}}}_{{p_5}}}{\\tilde{a}^{{p_2}{p_1}}_{{p_1}{p_4}}}}
      //       +
      //       {{s^{{p_3}}_{{m_{162}}}}{s^{{m_{162}}}_{{p_5}}}{s^{{e_{163}}}_{{p_4}}}{s^{{p_1}}_{{e_{163}}}}{\\tilde{a}^{{p_2}}_{{p_1}}}}
      //       +
      //       {{s^{{e_{164}}}_{{p_4}}}{s^{{p_1}}_{{e_{164}}}}{\\tilde{a}^{{p_2}{p_3}}_{{p_1}{p_5}}}}
      //       -
      //       {{s^{{e_{165}}}_{{p_5}}}{s^{{p_1}}_{{e_{165}}}}{\\tilde{a}^{{p_2}{p_3}}_{{p_1}{p_4}}}}
      //       + {{\\tilde{a}^{{p_2}{p_3}{p_1}}_{{p_1}{p_4}{p_5}}}}\\bigr) }");
    }

    // odd number of ops -> full contraction is 0
    {
      auto opseq = ex<FNOperatorSeq>(
          FNOperator(cre({L"p_1", L"p_2"}), ann({L"p_4", L"p_5"})),
          FNOperator(cre({L"p_7", L"p_8"}), ann({L"p_10", L"p_11", L"p_12"})));
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.compute());
      auto result = wick.compute();
      REQUIRE(result->is<Constant>());
      REQUIRE(result->as<Constant>().value<int>() == 0);
    }

    // 4-body ^ 4-body
    SECTION("wick(4^4)") {
      auto opseq = ex<FNOperatorSeq>(
          FNOperator(cre({L"p_1", L"p_2", L"p_3", L"p_4"}),
                     ann({L"p_5", L"p_6", L"p_7", L"p_8"})),
          FNOperator(cre({L"p_21", L"p_22", L"p_23", L"p_24"}),
                     ann({L"p_25", L"p_26", L"p_27", L"p_28"})));
      auto wick = FWickTheorem{opseq};
      auto result = wick.compute(true);
      REQUIRE(result->is<Constant>());
      REQUIRE(result->as<Constant>().value<int>() == 576);
      REQUIRE(576 == GWT({4, 4}).result().size());
    }

    // three general 1-body operators
    {
      auto opseq = ex<FNOperatorSeq>(FNOperator(cre({L"p_1"}), ann({L"p_2"})),
                                     FNOperator(cre({L"p_3"}), ann({L"p_4"})),
                                     FNOperator(cre({L"p_5"}), ann({L"p_6"})));
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.compute());
      auto result = wick.compute();
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 2);
      REQUIRE(2 == GWT({1, 1, 1}).result().size());
    }

    // 4 general 1-body operators
    {
      auto opseq = ex<FNOperatorSeq>(FNOperator(cre({L"p_1"}), ann({L"p_2"})),
                                     FNOperator(cre({L"p_3"}), ann({L"p_4"})),
                                     FNOperator(cre({L"p_5"}), ann({L"p_6"})),
                                     FNOperator(cre({L"p_7"}), ann({L"p_8"})));
      auto ext_indices = make_indices<std::vector<Index>>(WstrList{
          L"p_1", L"p_2", L"p_3", L"p_4", L"p_5", L"p_6", L"p_7", L"p_8"});
      auto wick1 = FWickTheorem{opseq};
      auto result1 = wick1.set_external_indices(ext_indices).compute();
      REQUIRE(result1->is<Sum>());
      REQUIRE(result1->size() == 9);
      REQUIRE(9 == GWT({1, 1, 1, 1}).result().size());
      auto wick2 = FWickTheorem{opseq};
      auto result2 = wick2.set_external_indices(ext_indices)
                         .set_nop_connections({{1, 2}, {1, 3}})
                         .compute();
      REQUIRE(result2->is<Sum>());
      REQUIRE(result2->size() == 2);
    }

    // 4-body ^ 2-body ^ 2-body
    {
      auto opseq = ex<FNOperatorSeq>(
          FNOperator(cre({L"p_1", L"p_2", L"p_3", L"p_4"}),
                     ann({L"p_5", L"p_6", L"p_7", L"p_8"})),
          FNOperator(cre({L"p_9", L"p_10"}), ann({L"p_11", L"p_12"})),
          FNOperator(cre({L"p_13", L"p_14"}), ann({L"p_15", L"p_16"})));
      auto wick = FWickTheorem{opseq};
      auto result = wick.compute();
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 576);
      REQUIRE(576 == GWT({4, 2, 2}).result().size());
    }

    // 2-body ^ 2-body ^ 2-body
    {
      auto opseq = ex<FNOperatorSeq>(
          FNOperator(cre({L"p_1", L"p_2"}), ann({L"p_5", L"p_6"})),
          FNOperator(cre({L"p_9", L"p_10"}), ann({L"p_11", L"p_12"})),
          FNOperator(cre({L"p_17", L"p_18"}), ann({L"p_19", L"p_20"})));
      auto wick = FWickTheorem{opseq};
      auto result = wick.compute();
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 80);
      REQUIRE(80 == GWT({2, 2, 2}).result().size());
    }

    // 2-body ^ 2-body ^ 2-body ^ 2-body
    SECTION("wick(2^2^2^2)") {
      auto opseq = ex<FNOperatorSeq>(
          FNOperator(cre({L"p_1", L"p_2"}), ann({L"p_5", L"p_6"})),
          FNOperator(cre({L"p_9", L"p_10"}), ann({L"p_11", L"p_12"})),
          FNOperator(cre({L"p_13", L"p_14"}), ann({L"p_15", L"p_16"})),
          FNOperator(cre({L"p_17", L"p_18"}), ann({L"p_19", L"p_20"})));
      auto wick = FWickTheorem{opseq};
      auto result = wick.compute(true);
      REQUIRE(result->is<Constant>());
      REQUIRE(result->as<Constant>().value<int>() == 4752);
    }
#ifndef SEQUANT_SKIP_LONG_TESTS
    REQUIRE(4752 == GWT({2, 2, 2, 2}).result().size());
#endif

#ifndef SEQUANT_SKIP_LONG_TESTS

// set to 1 to use GWT to count terms for large contractions
#define SEQUANT_EXPENSIVELY_VALIDATE_LONG_TESTS 0

    // 4-body ^ 2-body ^ 4-body
    SECTION("wick(4^2^4)") {
      auto opseq = ex<FNOperatorSeq>(
          FNOperator(cre({L"p_1", L"p_2", L"p_3", L"p_4"}),
                     ann({L"p_5", L"p_6", L"p_7", L"p_8"})),
          FNOperator(cre({L"p_9", L"p_10"}), ann({L"p_11", L"p_12"})),
          FNOperator(cre({L"p_21", L"p_22", L"p_23", L"p_24"}),
                     ann({L"p_25", L"p_26", L"p_27", L"p_28"})));
      auto wick = FWickTheorem{opseq};
      auto result = wick.compute(true);
      REQUIRE(result->is<Constant>());
      REQUIRE(result->as<Constant>().value<int>() == 50688);
      if (SEQUANT_EXPENSIVELY_VALIDATE_LONG_TESTS)
        REQUIRE(50688 == GWT({4, 2, 4}).result().size());
    }

// comment out to run superlong tests
#define SEQUANT_SKIP_SUPERLONG_TESTS

#ifndef SEQUANT_SKIP_SUPERLONG_TESTS
    // 4-body ^ 2-body ^ 2-body ^ 2-body
    SECTION("wick(4^2^2^2)") {
      auto opseq = ex<FNOperatorSeq>(
          FNOperator(cre({L"p_1", L"p_2", L"p_3", L"p_4"}),
                     ann({L"p_5", L"p_6", L"p_7", L"p_8"})),
          FNOperator(cre({L"p_9", L"p_10"}), ann({L"p_11", L"p_12"})),
          FNOperator(cre({L"p_13", L"p_14"}), ann({L"p_15", L"p_16"})),
          FNOperator(cre({L"p_17", L"p_18"}), ann({L"p_19", L"p_20"})));
      auto wick = FWickTheorem{opseq};
      auto result = wick.compute(true);
      REQUIRE(result->is<Constant>());
      REQUIRE(result->as<Constant>().value<int>() == 117504);
      if (SEQUANT_EXPENSIVELY_VALIDATE_LONG_TESTS)
        REQUIRE(117504 == GWT({4, 2, 2, 2}).result().size());
    }

    // 3-body ^ 2-body ^ 2-body ^ 3-body
    SECTION("wick(3^2^2^3)") {
      auto opseq = ex<FNOperatorSeq>(
          FNOperator(cre({L"p_1", L"p_2", L"p_3"}),
                     ann({L"p_5", L"p_6", L"p_7"})),
          FNOperator(cre({L"p_9", L"p_10"}), ann({L"p_11", L"p_12"})),
          FNOperator(cre({L"p_13", L"p_14"}), ann({L"p_15", L"p_16"})),
          FNOperator(cre({L"p_17", L"p_18", L"p_19"}),
                     ann({L"p_20", L"p_21", L"p_22"})));
      auto wick = FWickTheorem{opseq};
      auto result = wick.compute(true);
      REQUIRE(result->is<Constant>());
      REQUIRE(result->as<Constant>().value<int>() == 202320);
      if (SEQUANT_EXPENSIVELY_VALIDATE_LONG_TESTS)
        REQUIRE(202320 == GWT({3, 2, 2, 3}).result().size());
    }

    // 4-body ^ 4-body ^ 4-body
    SECTION("wick(4^4^4)") {
      auto opseq = ex<FNOperatorSeq>(
          FNOperator(cre({L"p_1", L"p_2", L"p_3", L"p_4"}),
                     ann({L"p_5", L"p_6", L"p_7", L"p_8"})),
          FNOperator(cre({L"p_11", L"p_12", L"p_13", L"p_14"}),
                     ann({L"p_15", L"p_16", L"p_17", L"p_18"})),
          FNOperator(cre({L"p_21", L"p_22", L"p_23", L"p_24"}),
                     ann({L"p_25", L"p_26", L"p_27", L"p_28"})));
      auto wick = FWickTheorem{opseq};
      auto result = wick.compute(true);
      REQUIRE(result->is<Constant>());
      REQUIRE(result->as<Constant>().value<int>() == 4783104);
      if (SEQUANT_EXPENSIVELY_VALIDATE_LONG_TESTS)
        REQUIRE(4783104 == GWT({4, 4, 4}).result().size());
    }
#endif

#endif

#if 0
    // impossible: 4-body ^ 4-body ^ 4-body ^ 4-body ^ 4-body ^ 4-body
    {
      auto opseq =
          ex<FNOperatorSeq>(FNOperator(cre({L"p_1", L"p_2", L"p_3", L"p_4"}), ann({L"p_5", L"p_6", L"p_7", L"p_8"})),
                         FNOperator(cre({L"p_11", L"p_12", L"p_13", L"p_14"}), ann({L"p_15", L"p_16", L"p_17", L"p_18"})),
                         FNOperator(cre({L"p_21", L"p_22", L"p_23", L"p_24"}), ann({L"p_25", L"p_26", L"p_27", L"p_28"})),
                         FNOperator(cre({L"p_31", L"p_32", L"p_33", L"p_34"}), ann({L"p_35", L"p_36", L"p_37", L"p_38"})),
                         FNOperator(cre({L"p_41", L"p_42", L"p_43", L"p_44"}), ann({L"p_45", L"p_46", L"p_47", L"p_48"})),
                         FNOperator(cre({L"p_51", L"p_52", L"p_53", L"p_54"}), ann({L"p_55", L"p_56", L"p_57", L"p_58"}))
                        );
      auto wick = FWickTheorem{opseq};
      auto result = wick.compute(true);
    }
#endif
  }  // SECTION("fermi vacuum")

  SECTION("Expression Reduction") {
    constexpr Vacuum V = Vacuum::SingleProduct;
    // default vacuum is already spin-orbital Fermi vacuum

    auto switch_to_spinfree_context = detail::NoDiscard([&]() {
      auto context_sf = get_default_context();
      context_sf.set(SPBasis::Spinfree);
      return set_scoped_default_context(context_sf);
    });

    // 2-body ^ 2-body
    SECTION("wick(H2*T2)") {
      auto opseq = ex<FNOperatorSeq>(
          FNOperator(cre({L"p_1", L"p_2"}), ann({L"p_3", L"p_4"})),
          FNOperator(cre({L"a_4", L"a_5"}), ann({L"i_4", L"i_5"})));
      auto wick = FWickTheorem{opseq};
      auto wick_result = wick.compute();
      REQUIRE(wick_result->is<Sum>());
      REQUIRE(wick_result->size() == 4);

      // multiply tensor factors and expand
      auto wick_result_2 = ex<Tensor>(L"g", bra{L"p_1", L"p_2"},
                                      ket{L"p_3", L"p_4"}, Symmetry::Antisymm) *
                           ex<Tensor>(L"t", bra{L"a_4", L"a_5"},
                                      ket{L"i_4", L"i_5"}, Symmetry::Antisymm) *
                           wick_result;
      expand(wick_result_2);
      REQUIRE(to_latex(wick_result_2) ==
              L"{ "
              L"\\bigl({{\\bar{g}^{{p_3}{p_4}}_{{p_1}{p_2}}}{\\bar{t}^{{i_4}{i_"
              L"5}}_{{a_4}{a_"
              L"5}}}{s^{{p_1}}_{{i_5}}}{s^{{p_2}}_{{i_4}}}{s^{{a_4}}_{{p_4}}}{"
              L"s^{{a_5}}_{{p_3}}}} - {"
              L"{\\bar{g}^{{p_3}{p_4}}_{{p_1}{p_2}}}{\\bar{t}^{{i_4}{i_5}}_{{a_"
              L"4}{a_5}}}{s^{{"
              L"p_1}}_{{i_5}}}{s^{{p_2}}_{{i_4}}}{s^{{a_5}}_{{p_4}}}{s^{{a_4}}_"
              L"{{p_3}}}} - {"
              L"{\\bar{g}^{{p_3}{p_4}}_{{p_1}{p_2}}}{\\bar{t}^{{i_4}{i_5}}_{{a_"
              L"4}{a_5}}}{s^{{"
              L"p_1}}_{{i_4}}}{s^{{p_2}}_{{i_5}}}{s^{{a_4}}_{{p_4}}}{s^{{a_5}}_"
              L"{{p_3}}}} + "
              L"{{\\bar{g}^{{p_3}{p_4}}_{{p_1}{p_2}}}{\\bar{t}^{{i_4}{i_5}}_{{"
              L"a_4}{a_5}}}{s^{"
              L"{p_1}}_{{i_4}}}{s^{{p_2}}_{{i_5}}}{s^{{a_5}}_{{p_4}}}{s^{{a_4}}"
              L"_{{p_3}}}}\\bigr) }");
      wick.reduce(wick_result_2);
      rapid_simplify(wick_result_2);
      TensorCanonicalizer::register_instance(
          std::make_shared<DefaultTensorCanonicalizer>());
      canonicalize(wick_result_2, {.method = CanonicalizationMethod::Complete});
      rapid_simplify(wick_result_2);

      std::wcout << L"H2*T2 = " << to_latex(wick_result_2) << std::endl;
      // std::wcout << L"H2*T2 = " << to_wolfram(wick_result_2) << std::endl;
      REQUIRE(to_latex(wick_result_2) ==
              L"{{{4}}"
              L"{\\bar{g}^{{a_1}{a_2}}_{{i_1}{i_2}}}{\\bar{t}^{{i_1}{i_2}}_{{a_"
              L"1}{a_2}}}}");

      // spin-free case will produce 2 terms
      {
        auto raii_tmp = switch_to_spinfree_context();

        auto wick = FWickTheorem{opseq};
        auto wick_result = wick.compute();
        REQUIRE(wick_result->is<Sum>());
        REQUIRE(wick_result->size() == 4);

        // multiply tensor factors and expand
        auto wick_result_2 =
            ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"},
                       Symmetry::Nonsymm) *
            ex<Tensor>(L"t", bra{L"a_4", L"a_5"}, ket{L"i_4", L"i_5"},
                       Symmetry::Nonsymm) *
            wick_result;
        expand(wick_result_2);
        REQUIRE(wick_result_2->size() == 4);  // still 4 terms

        wick.reduce(wick_result_2);
        rapid_simplify(wick_result_2);
        TensorCanonicalizer::register_instance(
            std::make_shared<DefaultTensorCanonicalizer>());
        canonicalize(wick_result_2);
        rapid_simplify(wick_result_2);
        REQUIRE(wick_result_2->size() == 2);  // now 2 terms

        std::wcout << L"spinfree H2*T2 = " << to_latex(wick_result_2)
                   << std::endl;
        REQUIRE_THAT(wick_result_2,
                     EquivalentTo("-4 * g{i1,i2;a1,a2} t{a1,a2;i2,i1} + 8 "
                                  "g{i1,i2;a1,a2} t{a1,a2;i1,i2}"));
      }
    }

    // 2-body ^ 1-body ^ 1-body, with/without using topology
    SECTION("wick(H2*T1*T1)") {
      for (auto&& use_nop_partitions : {false}) {
        for (auto&& use_op_partitions : {true, false}) {
          std::wostringstream oss;
          oss << "use_op_partitions=" << use_op_partitions << "}: H2*T1*T1 = ";

          auto opseq = ex<FNOperatorSeq>(
              FNOperator(cre({L"p_1", L"p_2"}), ann({L"p_3", L"p_4"})),
              FNOperator(cre({L"a_4"}), ann({L"i_4"})),
              FNOperator(cre({L"a_5"}), ann({L"i_5"})));
          auto wick = FWickTheorem{opseq};
          wick.use_topology(use_nop_partitions || use_op_partitions);
          // if (use_nop_partitions) wick.set_nop_partitions({{1, 2}});
          if (use_op_partitions) wick.set_op_partitions({{0, 1}, {2, 3}});
          auto wick_result = wick.compute();
          // print(oss.str() + L" (nopseq only) ", wick_result);
          if (use_op_partitions) {
            REQUIRE(wick_result->is<Product>());
            REQUIRE(wick_result->size() == 4 /* factors */);
          } else {
            REQUIRE(wick_result->is<Sum>());
            REQUIRE(wick_result->size() == 4 /* summands */);
          }

          // multiply tensor factors and expand
          auto wick_result_2 =
              ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"},
                         Symmetry::Antisymm) *
              ex<Tensor>(L"t", bra{L"a_4"}, ket{L"i_4"}, Symmetry::Antisymm) *
              ex<Tensor>(L"t", bra{L"a_5"}, ket{L"i_5"}, Symmetry::Antisymm) *
              wick_result;
          expand(wick_result_2);
          wick.reduce(wick_result_2);
          rapid_simplify(wick_result_2);
          TensorCanonicalizer::register_instance(
              std::make_shared<DefaultTensorCanonicalizer>(
                  std::vector<Index>{}));
          canonicalize(wick_result_2,
                       {.method = CanonicalizationMethod::Complete});
          rapid_simplify(wick_result_2);

          // print(oss.str(), wick_result_2);
          REQUIRE(wick_result_2->size() == 3 /* factors */);
          REQUIRE(to_latex(wick_result_2) ==
                  L"{{{4}}"
                  L"{\\bar{g}^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}}_{{a_1}}}{"
                  L"t^{{i_"
                  L"2}}_{{a_"
                  L"2}}}}");
        }  // use_op_partitions
      }    // use_nop_partitions
    }

    // 2-body ^ 1-body ^ 2-body with dependent (PNO) indices
    SECTION("wick(P2*H1*T2)") {
      auto opseq =
          ex<FNOperatorSeq>(FNOperator(cre({L"i_1", L"i_2"}),
                                       ann({Index(L"a_1", {L"i_1", L"i_2"}),
                                            Index(L"a_2", {L"i_1", L"i_2"})}),
                                       V),
                            FNOperator(cre({L"p_1"}), ann({L"p_2"})),
                            FNOperator(cre({Index(L"a_3", {L"i_3", L"i_4"}),
                                            Index(L"a_4", {L"i_3", L"i_4"})}),
                                       ann({L"i_3", L"i_4"})));
      auto wick = FWickTheorem{opseq};
      auto wick_result = wick.compute();
      REQUIRE(wick_result->is<Sum>());
      REQUIRE(wick_result->size() == 16);

      // multiply tensor factors and expand
      auto wick_result_2 =
          ex<Tensor>(L"A", bra{L"i_1", L"i_2"},
                     ket{Index{L"a_1", {L"i_1", L"i_2"}},
                         Index{L"a_2", {L"i_1", L"i_2"}}},
                     Symmetry::Antisymm) *
          ex<Tensor>(L"f", bra{L"p_1"}, ket{L"p_2"}, Symmetry::Antisymm) *
          ex<Tensor>(L"t",
                     bra{Index{L"a_3", {L"i_3", L"i_4"}},
                         Index{L"a_4", {L"i_3", L"i_4"}}},
                     ket{L"i_3", L"i_4"}, Symmetry::Antisymm) *
          wick_result;
      expand(wick_result_2);
      wick.reduce(wick_result_2);
      rapid_simplify(wick_result_2);
      TensorCanonicalizer::register_instance(
          std::make_shared<DefaultTensorCanonicalizer>());
      canonicalize(wick_result_2);
      rapid_simplify(wick_result_2);
    }

    // 2=body ^ 2-body ^ 2-body ^ 2-body with dependent (PNO) indices
    SECTION("wick(P2*H2*T2*T2)") {
      for (auto&& use_nop_partitions : {false}) {
        for (auto&& use_op_partitions : {true, false}) {
          std::wostringstream oss;
          oss << "use_{nop,op}_partitions={" << use_nop_partitions << ","
              << use_op_partitions << "}: P2*H2*T2*T2(PNO) = ";

          auto opseq = ex<FNOperatorSeq>(
              FNOperator(cre({L"i_1", L"i_2"}),
                         ann({Index(L"a_1", {L"i_1", L"i_2"}),
                              Index(L"a_2", {L"i_1", L"i_2"})}),
                         V),
              FNOperator(cre({L"p_1", L"p_2"}), ann({L"p_3", L"p_4"})),
              FNOperator(cre({Index(L"a_3", {L"i_3", L"i_4"}),
                              Index(L"a_4", {L"i_3", L"i_4"})}),
                         ann({L"i_3", L"i_4"})),
              FNOperator(cre({Index(L"a_5", {L"i_5", L"i_6"}),
                              Index(L"a_6", {L"i_5", L"i_6"})}),
                         ann({L"i_5", L"i_6"})));
          auto wick = FWickTheorem{opseq};
          wick.set_nop_connections({{1, 2}, {1, 3}}).use_topology(true);

          if (use_nop_partitions) wick.set_nop_partitions({{2, 3}});
          if (use_op_partitions)
            wick.set_op_partitions({{0, 1},
                                    {2, 3},
                                    {4, 5},
                                    {6, 7},
                                    {8, 9},
                                    {10, 11},
                                    {12, 13},
                                    {14, 15}});
          auto wick_result = wick.compute();
          REQUIRE(wick_result->is<Sum>());
          if (use_op_partitions) {
            REQUIRE(wick_result->size() == 7);
          } else {
            REQUIRE(wick_result->size() == 544);
          }

          // multiply tensor factors and expand
          auto wick_result_2 =
              ex<Constant>(rational{1, 256}) *
              ex<Tensor>(L"A", bra{L"i_1", L"i_2"},
                         ket{Index{L"a_1", {L"i_1", L"i_2"}},
                             Index{L"a_2", {L"i_1", L"i_2"}}},
                         Symmetry::Antisymm) *
              ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"},
                         Symmetry::Antisymm) *
              ex<Tensor>(L"t",
                         bra{Index{L"a_3", {L"i_3", L"i_4"}},
                             Index{L"a_4", {L"i_3", L"i_4"}}},
                         ket{L"i_3", L"i_4"}, Symmetry::Antisymm) *
              ex<Tensor>(L"t",
                         bra{Index{L"a_5", {L"i_5", L"i_6"}},
                             Index{L"a_6", {L"i_5", L"i_6"}}},
                         ket{L"i_5", L"i_6"}, Symmetry::Antisymm) *
              wick_result;
          expand(wick_result_2);
          wick.reduce(wick_result_2);
          rapid_simplify(wick_result_2);
          TensorCanonicalizer::register_instance(
              std::make_shared<DefaultTensorCanonicalizer>());
          canonicalize(wick_result_2);
          canonicalize(wick_result_2);
          canonicalize(wick_result_2);
          rapid_simplify(wick_result_2);

          //        std::wcout << oss.str() << to_latex_align(wick_result_2, 20)
          //                   << std::endl;
          REQUIRE(wick_result_2->is<Sum>());
          REQUIRE(wick_result_2->size() == 4);
        }  // use_op_partitions
      }    // use_nop_partitions
    }

#if 1
    // 3-body ^ 2-body ^ 2-body ^ 3-body
    SECTION("wick(P3*H2*T2*T3)") {
      constexpr bool connected_only = true;
      constexpr bool topology = true;
      auto P3 = ex<Constant>(rational{1, 36}) *
                ex<Tensor>(L"A", bra{L"i_1", L"i_2", L"i_3"},
                           ket{L"a_1", L"a_2", L"a_3"}, Symmetry::Antisymm) *
                ex<FNOperator>(cre{L"i_1", L"i_2", L"i_3"},
                               ann{L"a_1", L"a_2", L"a_3"});
      auto H2 = ex<Constant>(rational{1, 4}) *
                ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"},
                           Symmetry::Antisymm) *
                ex<FNOperator>(cre{L"p_1", L"p_2"}, ann{L"p_3", L"p_4"});
      auto T2 = ex<Constant>(rational{1, 4}) *
                ex<Tensor>(L"t", bra{L"a_4", L"a_5"}, ket{L"i_4", L"i_5"},
                           Symmetry::Antisymm) *
                ex<FNOperator>(cre{L"a_4", L"a_5"}, ann{L"i_4", L"i_5"});
      auto T3 = ex<Constant>(rational{1, 36}) *
                ex<Tensor>(L"t", bra{L"a_6", L"a_7", L"a_8"},
                           ket{L"i_6", L"i_7", L"i_8"}, Symmetry::Antisymm) *
                ex<FNOperator>(cre{L"a_6", L"a_7", L"a_8"},
                               ann{L"i_6", L"i_7", L"i_8"});
      FWickTheorem wick{P3 * H2 * T2 * T3};
      wick.use_topology(topology);
      if (connected_only) wick.set_nop_connections({{1, 2}, {1, 3}});
      auto wick_result = wick.compute();

      std::wcout << "P3*H2*T2*T3 = " << to_latex_align(wick_result, 20)
                 << std::endl;
      REQUIRE(wick_result->is<Sum>());
      REQUIRE(
          wick_result->size() ==
          (connected_only ? 7 : 9));  // 9 = 2 disconnected + 7 connected terms
    }
#endif

    // example with "diagonal" operator from Nick Mayhall
    {
      auto _ = set_scoped_default_context(
          {.index_space_registry_shared_ptr = mbpt::make_min_sr_spaces(),
           .vacuum = Vacuum::SingleProduct,
           .braket_symmetry = BraKetSymmetry::Symm,
           .spbasis = SPBasis::Spinor});

      // sequence of individual ops
      const auto ops = {fann("p_1"), fcre("p_2"), fcre("p_3"),
                        fann("p_5"), fann("p_4"), fcre("p_1")};
      REQUIRE_NOTHROW(FWickTheorem(ops));
      FWickTheorem w(ops);
      REQUIRE_NOTHROW(w.full_contractions(false).compute());
      auto wresult = FWickTheorem{ops}.full_contractions(false).compute();

      auto result0 =
          wresult * ex<Constant>(ratio(1, 4)) *
          ex<Tensor>(L"v", bra{L"p_2", L"p_3"}, ket{L"p_4", L"p_5"},
                     Symmetry::Antisymm) *
          ex<Tensor>(L"w", bra{}, ket{}, aux{L"p_1"}, Symmetry::Nonsymm);
      // std::wcout << "before expand: op = " << to_latex(op) << std::endl;
      expand(result0);
      // std::wcout << "after expand: op = " << to_latex(op) << std::endl;
      w.reduce(result0);
      // std::wcout << "after reduce: op = " << to_latex(op) << std::endl;
      simplify(result0, {{.named_indices = IndexList{}}});
      // sequant::wprintf(L"after simplify: op = ", to_latex_align(result0, 3),
      // L"\n");

      auto Ld_H2_L =
          ex<Constant>(ratio(1, 4)) *
          ex<Tensor>(L"v", bra{L"p_2", L"p_3"}, ket{L"p_4", L"p_5"},
                     Symmetry::Antisymm) *
          ex<Tensor>(L"w", bra{}, ket{}, aux{L"p_1"}, Symmetry::Nonsymm) *
          fannx("p_1") * fcrex("p_2") * fcrex("p_3") * fannx("p_5") *
          fannx("p_4") * fcrex("p_1");
      auto result1 = FWickTheorem{Ld_H2_L}.full_contractions(false).compute();
      simplify(result1, {{.named_indices = IndexList{}}});
      // sequant::wprintf(to_latex_align(Ld_H2_L), L" = \n",
      // to_latex_align(result1, 0, 2), L"\n");
      REQUIRE(result0 == result1);

      auto Ld_H2N_L =
          ex<Constant>(ratio(1, 4)) *
          ex<Tensor>(L"v", bra{L"p_2", L"p_3"}, ket{L"p_4", L"p_5"},
                     Symmetry::Antisymm) *
          ex<Tensor>(L"w", bra{}, ket{}, aux{L"p_1"}, Symmetry::Nonsymm) *
          fannx("p_1") *
          ex<FNOperator>(cre({L"p_2", L"p_3"}), ann({L"p_4", L"p_5"})) *
          fcrex("p_1");
      auto result2 = FWickTheorem{Ld_H2N_L}.full_contractions(false).compute();
      simplify(result2, {{.method = CanonicalizationMethod::Complete,
                          .named_indices = IndexList{}}});
      sequant::wprintf(to_latex_align(Ld_H2N_L), L" = \n",
                       to_latex_align(result2, 0, 2), L"\n");
      REQUIRE(result2.as<Sum>().size() == 5);
      REQUIRE(result2.to_latex() ==
              L"{ \\bigl( - "
              L"{{{\\frac{1}{4}}}{\\tensor*{\\tilde{a}}{*^{p_3}_{p_1}*^{p_4}_{"
              L"p_2}*^{p_5}_{p_5}}}{\\tensor*{\\bar{v}}{*^{p_3}_{p_1}*^{p_4}_{"
              L"p_2}}}{\\tensor*{w}{}[{p_5}]}} + "
              L"{{\\tensor*{\\tilde{a}}{*^{p_2}_{p_1}}}{\\tensor*{\\bar{v}}{*^{"
              L"a_1}_{a_1}*^{p_2}_{p_1}}}{\\tensor*{w}{}[{a_1}]}} - "
              L"{{{\\frac{1}{2}}}{\\tensor*{\\tilde{a}}{*^{p_2}_{a_1}*^{p_3}_{"
              L"p_1}}}{\\tensor*{\\bar{v}}{*^{a_1}_{p_2}*^{p_1}_{p_3}}}{"
              L"\\tensor*{w}{}[{a_1}]}} + "
              L"{{{\\frac{1}{4}}}{\\tensor*{\\tilde{a}}{*^{p_3}_{p_1}*^{p_4}_{"
              L"p_2}}}{\\tensor*{\\bar{v}}{*^{p_1}_{p_3}*^{p_2}_{p_4}}}{"
              L"\\tensor*{w}{}[{a_1}]}} - "
              L"{{{\\frac{1}{2}}}{\\tensor*{\\tilde{a}}{*^{a_1}_{p_1}*^{p_3}_{"
              L"p_2}}}{\\tensor*{\\bar{v}}{*^{a_1}_{p_1}*^{p_3}_{p_2}}}{"
              L"\\tensor*{w}{}[{a_1}]}}\\bigr) }");
    }
  }
}
#endif
