//
// Created by Eduard Valeyev on 3/20/18.
//

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>

#include <SeQuant/core/op.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

TEST_CASE("Op", "[elements]") {
  using namespace sequant;

  SECTION("constructors") {
    REQUIRE_NOTHROW(FOp{});

    REQUIRE_NOTHROW(FOp(Index(L"i_1"), Action::create));
    FOp o1(Index(L"i_1"), Action::create);
    REQUIRE(o1.statistics == Statistics::FermiDirac);
    REQUIRE(o1.index() == Index(L"i_1"));
    REQUIRE(o1.action() == Action::create);

    {
      REQUIRE_NOTHROW(FOp(L"i_1", Action::create));
      FOp o1(L"i_1", Action::create);
      REQUIRE(o1.statistics == Statistics::FermiDirac);
      REQUIRE(o1.index() == Index(L"i_1"));
      REQUIRE(o1.action() == Action::create);
    }

    REQUIRE_NOTHROW(bann(Index(L"i_2")));
    auto o2 = bann(Index(L"i_2"));
    REQUIRE(o2.statistics == Statistics::BoseEinstein);
    REQUIRE(o2.index() == Index(L"i_2"));
    REQUIRE(o2.action() == Action::annihilate);

    REQUIRE_NOTHROW(bcre(L"i_3"));
    auto o3 = bcre(L"i_3");
    REQUIRE(o3.statistics == Statistics::BoseEinstein);
    REQUIRE(o3.index() == Index(L"i_3"));
    REQUIRE(o3.action() == Action::create);

    REQUIRE_NOTHROW(fann(L"i_4", {L"a_1", L"a_2"}));
    auto o4 = fann(L"i_4", {L"a_1", L"a_2"});
    REQUIRE(o4.statistics == Statistics::FermiDirac);
    REQUIRE(o4.index() == Index(L"i_4", {L"a_1", L"a_2"}));
    REQUIRE(o4.action() == Action::annihilate);

    REQUIRE_NOTHROW(FOperator{fcre(L"i_1"), fann(L"i_1")});
    auto oper1 = FOperator{fcre(L"i_1"), fann(L"i_1")};
    REQUIRE(oper1.statistics == Statistics::FermiDirac);
    REQUIRE(oper1.size() == 2);
    REQUIRE(oper1[0] == fcre(L"i_1"));
    REQUIRE(oper1[1] == fann(L"i_1"));

    /////////////////// N-conserving normal op
    // 1 body
    REQUIRE_NOTHROW(FNOperator(cre({L"i_1"}), ann({L"a_1"})));
    auto nop1 = FNOperator(cre({L"i_1"}), ann({L"a_1"}));
    REQUIRE(nop1.creators().size() == 1);
    REQUIRE(nop1.annihilators().size() == 1);
    REQUIRE(nop1.creators()[0] == fcre(L"i_1"));
    REQUIRE(nop1.annihilators()[0] == fann(L"a_1"));

    // 2 body
    REQUIRE_NOTHROW(FNOperator(cre({Index{L"i_1"}, Index{L"i_2"}}),
                               ann({Index{L"a_1", {L"i_1", L"i_2"}},
                                    Index{L"a_2", {L"i_1", L"i_2"}}})));
    auto nop2 = FNOperator(cre({Index{L"i_1"}, Index{L"i_2"}}),
                           ann({Index{L"a_1", {L"i_1", L"i_2"}},
                                Index{L"a_2", {L"i_1", L"i_2"}}}));
    REQUIRE(nop2.creators().size() == 2);
    REQUIRE(nop2.annihilators().size() == 2);
    REQUIRE(nop2.creators()[0] == fcre(L"i_1"));
    REQUIRE(nop2.creators()[1] == fcre(L"i_2"));
    REQUIRE(nop2.annihilators()[0] == fann(Index{L"a_1", {L"i_1", L"i_2"}}));
    REQUIRE(nop2.annihilators()[1] == fann(Index{L"a_2", {L"i_1", L"i_2"}}));

    /////////////////// N-nonconserving normal op
    REQUIRE_NOTHROW(FNOperator(cre({L"i_1"}), ann({})));
    auto nop3 = FNOperator(cre({L"i_1"}), ann({}));
    REQUIRE(nop3.creators().size() == 1);
    REQUIRE(nop3.annihilators().size() == 0);
    REQUIRE(nop3.creators()[0] == fcre(L"i_1"));
    REQUIRE_NOTHROW(FNOperator(cre({}), ann({Index{L"a_1", {L"i_1"}}})));
    auto nop4 = FNOperator(cre({}), ann({Index{L"a_1", {L"i_1"}}}));
    REQUIRE(nop4.creators().size() == 0);
    REQUIRE(nop4.annihilators().size() == 1);
    REQUIRE(nop4.annihilators()[0] == fann(Index{L"a_1", {L"i_1"}}));

    /////////////////// normal op sequence
    // can include N-conserving ...
    REQUIRE_NOTHROW(FNOperatorSeq({FNOperator(cre({L"i_1"}), ann({L"i_2"})),
                                   FNOperator(cre({L"i_3"}), ann({L"i_4"})),
                                   FNOperator(cre({L"i_5"}), ann({L"i_6"}))}));
    auto nopseq1 = FNOperatorSeq({FNOperator(cre({L"i_1"}), ann({L"i_2"})),
                                  FNOperator(cre({L"i_3"}), ann({L"i_4"})),
                                  FNOperator(cre({L"i_5"}), ann({L"i_6"}))});
    REQUIRE(nopseq1.size() == 3);
    REQUIRE(nopseq1[0] == FNOperator(cre({L"i_1"}), ann({L"i_2"})));
    REQUIRE(nopseq1[1] == FNOperator(cre({L"i_3"}), ann({L"i_4"})));
    REQUIRE(nopseq1[2] == FNOperator(cre({L"i_5"}), ann({L"i_6"})));
    REQUIRE(nopseq1.at(0) == FNOperator(cre({L"i_1"}), ann({L"i_2"})));
    REQUIRE(nopseq1.at(1) == FNOperator(cre({L"i_3"}), ann({L"i_4"})));
    REQUIRE(nopseq1.at(2) == FNOperator(cre({L"i_5"}), ann({L"i_6"})));

    // ... N-nonconserving normal ops ...
    auto nopseq2 = FNOperatorSeq({FNOperator(cre({}), ann({L"i_1"}))});
    REQUIRE(nopseq2.size() == 1);
    REQUIRE(nopseq2[0] == FNOperator(cre({}), ann({L"i_1"})));

    // but all vacua must match
    REQUIRE_THROWS(FNOperatorSeq(
        {FNOperator(cre({L"i_1"}), ann({L"i_2"}), Vacuum::Physical),
         FNOperator(cre({L"i_3"}), ann({L"i_4"}), Vacuum::SingleProduct),
         FNOperator(cre({L"i_5"}), ann({L"i_6"}))}));
  }

  SECTION("adjoint") {
    auto o1 = adjoint(FOp(Index(L"i_1"), Action::create));
    REQUIRE(o1.statistics == Statistics::FermiDirac);
    REQUIRE(o1.index() == Index(L"i_1"));
    REQUIRE(o1.action() == Action::annihilate);
    o1.adjoint();
    REQUIRE(o1.action() == Action::create);
    o1.adjoint();
    o1.adjoint();
    REQUIRE(o1.action() == Action::create);

    auto oper1 = adjoint(FOperator{fcre(L"i_1"), fann(L"i_2")});
    REQUIRE(oper1.statistics == Statistics::FermiDirac);
    REQUIRE(oper1.size() == 2);
    REQUIRE(oper1[0] == fcre(L"i_2"));
    REQUIRE(oper1[1] == fann(L"i_1"));

    auto nop2 = adjoint(
        FNOperator(cre({Index{L"i_1"}}),
                   ann({Index{L"a_1", {L"i_1"}}, Index{L"a_2", {L"i_1"}}})));
    REQUIRE(nop2.creators().size() == 2);
    REQUIRE(nop2.annihilators().size() == 1);
    REQUIRE(nop2.annihilators()[0] == fann(L"i_1"));
    REQUIRE(nop2.creators()[0] == fcre(Index{L"a_1", {L"i_1"}}));
    REQUIRE(nop2.creators()[1] == fcre(Index{L"a_2", {L"i_1"}}));

    auto nopseq1 =
        adjoint(FNOperatorSeq({FNOperator(cre({L"i_1"}), ann({L"i_2"})),
                               FNOperator(cre({L"i_3"}), ann({L"i_4"})),
                               FNOperator(cre({L"i_5"}), ann({L"i_6"}))}));
    REQUIRE(nopseq1.size() == 3);
    REQUIRE(nopseq1[2] == adjoint(FNOperator(cre({L"i_1"}), ann({L"i_2"}))));
    REQUIRE(nopseq1[1] == adjoint(FNOperator(cre({L"i_3"}), ann({L"i_4"}))));
    REQUIRE(nopseq1[0] == adjoint(FNOperator(cre({L"i_5"}), ann({L"i_6"}))));
  }

  SECTION("conversion") {
    auto nop1 = FNOperator(cre({Index{L"i_1"}, Index{L"i_2"}}),
                           ann({Index{L"a_1", {L"i_1", L"i_2"}},
                                Index{L"a_2", {L"i_1", L"i_2"}}}));
    REQUIRE_NOTHROW(static_cast<FOperator>(nop1));
    FOperator op1(static_cast<FOperator>(nop1));
    REQUIRE(op1.size() == 4);
    REQUIRE(op1[0] == fcre(L"i_1"));
    REQUIRE(op1[1] == fcre(L"i_2"));
    REQUIRE(op1[2] == fann(Index{L"a_2", {L"i_1", L"i_2"}}));
    REQUIRE(op1[3] == fann(Index{L"a_1", {L"i_1", L"i_2"}}));
  }

  SECTION("quasiparticle character") {
    {
      constexpr const Vacuum V = Vacuum::SingleProduct;

      REQUIRE(!is_qpcreator(fcre(L"i_1"), V));
      REQUIRE(is_qpcreator(fcre(L"a_1"), V));
      REQUIRE(is_qpcreator(fcre(L"p_1"), V));
      REQUIRE(!is_pure_qpcreator(fcre(L"i_1"), V));
      REQUIRE(is_pure_qpcreator(fcre(L"a_1"), V));
      REQUIRE(!is_pure_qpcreator(fcre(L"p_1"), V));

      REQUIRE(is_qpannihilator(fcre(L"i_1"), V));
      REQUIRE(!is_qpannihilator(fcre(L"a_1"), V));
      REQUIRE(is_qpannihilator(fcre(L"p_1"), V));
      REQUIRE(is_pure_qpannihilator(fcre(L"i_1"), V));
      REQUIRE(!is_pure_qpannihilator(fcre(L"a_1"), V));
      REQUIRE(!is_pure_qpannihilator(fcre(L"p_1"), V));

      REQUIRE(is_qpcreator(fann(L"i_1"), V));
      REQUIRE(!is_qpcreator(fann(L"a_1"), V));
      REQUIRE(is_qpcreator(fann(L"p_1"), V));
      REQUIRE(is_pure_qpcreator(fann(L"i_1"), V));
      REQUIRE(!is_pure_qpcreator(fann(L"a_1"), V));
      REQUIRE(!is_pure_qpcreator(fann(L"p_1"), V));

      REQUIRE(!is_qpannihilator(fann(L"i_1"), V));
      REQUIRE(is_qpannihilator(fann(L"a_1"), V));
      REQUIRE(is_qpannihilator(fann(L"p_1"), V));
      REQUIRE(!is_pure_qpannihilator(fann(L"i_1"), V));
      REQUIRE(is_pure_qpannihilator(fann(L"a_1"), V));
      REQUIRE(!is_pure_qpannihilator(fann(L"p_1"), V));
      auto isr = get_default_context().index_space_registry();
      REQUIRE(!qpannihilator_space(fann(L"i_1"), V));
      REQUIRE(qpannihilator_space(fcre(L"i_1"), V) == isr->retrieve(L"i_1"));
      REQUIRE(!qpannihilator_space(fcre(L"a_1"), V));
      REQUIRE(qpannihilator_space(fann(L"a_1"), V) == isr->retrieve(L"a_1"));
      REQUIRE(qpannihilator_space(fcre(L"p_1"), V) == isr->retrieve(L"m_1"));
      REQUIRE(qpannihilator_space(fann(L"p_1"), V) == isr->retrieve(L"e_1"));
    }
    {
      constexpr const Vacuum V = Vacuum::Physical;

      REQUIRE(is_qpcreator(fcre(L"i_1"), V));
      REQUIRE(is_qpcreator(fcre(L"a_1"), V));
      REQUIRE(is_qpcreator(fcre(L"p_1"), V));
      REQUIRE(is_pure_qpcreator(fcre(L"i_1"), V));
      REQUIRE(is_pure_qpcreator(fcre(L"a_1"), V));
      REQUIRE(is_pure_qpcreator(fcre(L"p_1"), V));

      REQUIRE(!is_qpcreator(fann(L"i_1"), V));
      REQUIRE(!is_qpcreator(fann(L"a_1"), V));
      REQUIRE(!is_qpcreator(fann(L"p_1"), V));
      REQUIRE(!is_pure_qpcreator(fann(L"i_1"), V));
      REQUIRE(!is_pure_qpcreator(fann(L"a_1"), V));
      REQUIRE(!is_pure_qpcreator(fann(L"p_1"), V));

      REQUIRE(!is_qpannihilator(fcre(L"i_1"), V));
      REQUIRE(!is_qpannihilator(fcre(L"a_1"), V));
      REQUIRE(!is_qpannihilator(fcre(L"p_1"), V));
      REQUIRE(!is_pure_qpannihilator(fcre(L"i_1"), V));
      REQUIRE(!is_pure_qpannihilator(fcre(L"a_1"), V));
      REQUIRE(!is_pure_qpannihilator(fcre(L"p_1"), V));

      REQUIRE(is_qpannihilator(fann(L"i_1"), V));
      REQUIRE(is_qpannihilator(fann(L"a_1"), V));
      REQUIRE(is_qpannihilator(fann(L"p_1"), V));
      REQUIRE(is_pure_qpannihilator(fann(L"i_1"), V));
      REQUIRE(is_pure_qpannihilator(fann(L"a_1"), V));
      REQUIRE(is_pure_qpannihilator(fann(L"p_1"), V));
    }
  }

  SECTION("hashing") {
    REQUIRE_NOTHROW(hash_value(FOp{}));
    REQUIRE_NOTHROW(hash_value(BOp{}));
    auto fc1 = fcre(L"i_1");
    auto fc1_copy = fcre(L"i_1");
    auto fa1 = fann(L"i_1");
    auto bc1 = bcre(L"i_1");
    auto ba1 = bann(L"i_1");
    auto fc2 = fcre(L"i_2");
    auto fa2 = fann(L"i_2");
    auto fc2_34 = fcre(L"i_2", {L"i_3", L"i_4"});
    REQUIRE_NOTHROW(hash_value(fc1));
    REQUIRE_NOTHROW(hash_value(fa1));
    REQUIRE_NOTHROW(hash_value(bc1));
    REQUIRE_NOTHROW(hash_value(ba1));
    REQUIRE_NOTHROW(hash_value(fc2_34));
    REQUIRE(hash_value(fc1) != hash_value(FOp{}));
    REQUIRE(hash_value(fc1) != hash_value(BOp{}));
    REQUIRE(hash_value(fc1) == hash_value(fc1_copy));
    REQUIRE(hash_value(fc1) != hash_value(fa1));
    REQUIRE(hash_value(fc1) != hash_value(fc2));
    REQUIRE(hash_value(fc1) !=
            hash_value(bc1));  // hash is depends on statistics
    REQUIRE(hash_value(fa1) !=
            hash_value(ba1));  // hash is depends on statistics
    REQUIRE(hash_value(fc1) == hash_value(adjoint(fa1)));
  }

  SECTION("hug") {
    auto nop1 = FNOperator(
        cre({Index{L"i_1"}, Index{L"i_2"}, Index{L"i_3"}, Index{L"a_1"}}),
        ann({Index{L"a_1", {L"i_1", L"i_2"}}, Index{L"a_2", {L"i_1", L"i_2"}},
             Index{L"a_1", {L"i_1", L"i_3"}}, Index{L"a_4"}}));
    REQUIRE_NOTHROW(nop1.hug());
    auto& hug1 = nop1.hug();
    REQUIRE(hug1->num_edges() == 8);
    REQUIRE(hug1->num_groups() == 5);
    REQUIRE(hug1->num_nonempty_groups() == 5);
    using group_idxs = decltype(hug1->group(0).second);
    REQUIRE(hug1->group(0).second == group_idxs{0, 1, 2});
    REQUIRE(hug1->group(1).second == group_idxs{0, 1, 2});
    REQUIRE(hug1->group(2).second == group_idxs{0, 1, 2});
    REQUIRE(hug1->group(3).second == group_idxs{3});
    REQUIRE(hug1->group(4).second == group_idxs{4});
    REQUIRE(hug1->group(5).second == group_idxs{5});
    REQUIRE(hug1->group(6).second == group_idxs{6, 7});
    REQUIRE(hug1->group(7).second == group_idxs{6, 7});
    REQUIRE(hug1->group_at(0).second == group_idxs{0, 1, 2});
    REQUIRE(hug1->group_at(1).second == group_idxs{3});
    REQUIRE(hug1->group_at(2).second == group_idxs{4});
    REQUIRE(hug1->group_at(3).second == group_idxs{5});
    REQUIRE(hug1->group_at(4).second == group_idxs{6, 7});
    REQUIRE(hug1->group_size(0) == 3);
    REQUIRE(hug1->group_size(1) == 3);
    REQUIRE(hug1->group_size(2) == 3);
    REQUIRE(hug1->group_size(3) == 1);
    REQUIRE(hug1->group_size(4) == 1);
    REQUIRE(hug1->group_size(5) == 1);
    REQUIRE(hug1->group_size(6) == 2);
    REQUIRE(hug1->group_size(7) == 2);

    // erase
    REQUIRE_NOTHROW(hug1->erase(0, fcre({L"i_1"})));
    REQUIRE(hug1->num_edges() == 7);
    REQUIRE(hug1->num_groups() == 5);
    REQUIRE(hug1->num_nonempty_groups() == 5);
    REQUIRE(hug1->group(0).second == group_idxs{0, 1});
    REQUIRE(hug1->group(1).second == group_idxs{0, 1});
    REQUIRE(hug1->group(2).second == group_idxs{2});
    REQUIRE(hug1->group(3).second == group_idxs{3});
    REQUIRE(hug1->group(4).second == group_idxs{4});
    REQUIRE(hug1->group(5).second == group_idxs{5, 6});
    REQUIRE(hug1->group(6).second == group_idxs{5, 6});
    REQUIRE(hug1->group_at(0).second == group_idxs{0, 1});
    REQUIRE(hug1->group_at(1).second == group_idxs{2});
    REQUIRE(hug1->group_at(2).second == group_idxs{3});
    REQUIRE(hug1->group_at(3).second == group_idxs{4});
    REQUIRE(hug1->group_at(4).second == group_idxs{5, 6});

    // another erase
    REQUIRE_NOTHROW(hug1->erase(2, fcre({L"a_17"})));
    REQUIRE(hug1->num_edges() == 6);
    REQUIRE(hug1->num_groups() == 5);
    REQUIRE(hug1->num_nonempty_groups() == 4);
    REQUIRE(hug1->group(0).second == group_idxs{0, 1});
    REQUIRE(hug1->group(1).second == group_idxs{0, 1});
    REQUIRE(hug1->group(2).second == group_idxs{2});
    REQUIRE(hug1->group(3).second == group_idxs{3});
    REQUIRE(hug1->group(4).second == group_idxs{4, 5});
    REQUIRE(hug1->group(5).second == group_idxs{4, 5});
    REQUIRE(hug1->group_at(0).second == group_idxs{0, 1});
    REQUIRE(hug1->group_at(1).second == group_idxs{});
    REQUIRE(hug1->group_at(2).second == group_idxs{2});
    REQUIRE(hug1->group_at(3).second == group_idxs{3});
    REQUIRE(hug1->group_at(4).second == group_idxs{4, 5});

    // now insert an index in new group
    REQUIRE_NOTHROW(hug1->insert(1, fcre({L"p_17"})));
    REQUIRE(hug1->num_edges() == 7);
    REQUIRE(hug1->num_groups() == 6);
    REQUIRE(hug1->num_nonempty_groups() == 5);
    REQUIRE(hug1->group(0).second == group_idxs{0, 2});
    REQUIRE(hug1->group(1).second == group_idxs{1});
    REQUIRE(hug1->group(2).second == group_idxs{0, 2});
    REQUIRE(hug1->group(3).second == group_idxs{3});
    REQUIRE(hug1->group(4).second == group_idxs{4});
    REQUIRE(hug1->group(5).second == group_idxs{5, 6});
    REQUIRE(hug1->group(6).second == group_idxs{5, 6});
    REQUIRE(hug1->group_at(0).second == group_idxs{0, 2});
    REQUIRE(hug1->group_at(1).second == group_idxs{});
    REQUIRE(hug1->group_at(2).second == group_idxs{3});
    REQUIRE(hug1->group_at(3).second == group_idxs{4});
    REQUIRE(hug1->group_at(4).second == group_idxs{5, 6});
    REQUIRE(hug1->group_at(5).second == group_idxs{1});

    // now insert an index in old group
    REQUIRE_NOTHROW(hug1->insert(6, fcre({L"i_17"})));
    REQUIRE(hug1->num_edges() == 8);
    REQUIRE(hug1->num_groups() == 6);
    REQUIRE(hug1->num_nonempty_groups() == 5);
    REQUIRE(hug1->group(0).second == group_idxs{0, 2, 6});
    REQUIRE(hug1->group(1).second == group_idxs{1});
    REQUIRE(hug1->group(2).second == group_idxs{0, 2, 6});
    REQUIRE(hug1->group(3).second == group_idxs{3});
    REQUIRE(hug1->group(4).second == group_idxs{4});
    REQUIRE(hug1->group(5).second == group_idxs{5, 7});
    REQUIRE(hug1->group(6).second == group_idxs{0, 2, 6});
    REQUIRE(hug1->group(7).second == group_idxs{5, 7});
    REQUIRE(hug1->group_at(0).second == group_idxs{0, 2, 6});
    REQUIRE(hug1->group_at(1).second == group_idxs{});
    REQUIRE(hug1->group_at(2).second == group_idxs{3});
    REQUIRE(hug1->group_at(3).second == group_idxs{4});
    REQUIRE(hug1->group_at(4).second == group_idxs{5, 7});
    REQUIRE(hug1->group_at(5).second == group_idxs{1});

    // inserting past the end is no-no
    REQUIRE_THROWS(hug1->insert(9, fcre({L"i_17"})));
  }

  SECTION("latex") {
    FOp o1(Index(L"i_1"), Action::create);
    REQUIRE(to_latex(o1) == L"{a^{\\dagger}_{i_1}}");

    BOp o2(Index(L"a_1"), Action::create);
    REQUIRE(to_latex(o2) == L"{b^{\\dagger}_{a_1}}");

    auto oper0 = FOperator{};
    REQUIRE(to_latex(oper0) == L"{}");

    auto oper1 = FOperator{fcre(L"i_1"), fann(L"i_1")};
    REQUIRE(to_latex(oper1) == L"{{a^{\\dagger}_{i_1}}{a_{i_1}}}");

    auto nop1 = FNOperator(cre({L"i_1", L"i_2"}), ann({L"a_1", L"a_2"}),
                           Vacuum::SingleProduct);
    REQUIRE(to_latex(nop1) == L"{\\tilde{a}^{{i_1}{i_2}}_{{a_1}{a_2}}}");

    auto nop2 = FNOperator(
        cre({Index{L"i_1"}, Index{L"i_2"}}),
        ann({Index{L"a_1", {L"i_1", L"i_2"}}, Index{L"a_2", {L"i_1", L"i_2"}}}),
        Vacuum::SingleProduct);
    REQUIRE(
        to_latex(nop2) ==
        L"{\\tilde{a}^{{i_1}{i_2}}_{{a_1^{{i_1}{i_2}}}{a_2^{{i_1}{i_2}}}}}");

    auto nop3 =
        FNOperator(cre({L"i_1", L"i_2"}), ann({L"a_2"}), Vacuum::SingleProduct);
    REQUIRE(to_latex(nop3) ==
            L"{\\tilde{a}^{{i_1}{i_2}}_{\\textvisiblespace\\,{a_2}}}");

    auto nop4 =
        FNOperator(cre({L"i_2"}), ann({L"a_1", L"a_2"}), Vacuum::SingleProduct);
    REQUIRE(to_latex(nop4) ==
            L"{\\tilde{a}^{\\textvisiblespace\\,{i_2}}_{{a_1}{a_2}}}");

    auto nop5 = FNOperator(cre({L"i_1"}), ann({}), Vacuum::SingleProduct);
    REQUIRE(to_latex(nop5) == L"{\\tilde{a}^{{i_1}}}");

    auto nop6 = FNOperator(cre({}), ann({L"a_1"}), Vacuum::SingleProduct);
    REQUIRE(to_latex(nop6) == L"{\\tilde{a}_{{a_1}}}");

    auto nopseq1 = FNOperatorSeq({nop1, nop2});
    REQUIRE(to_latex(nopseq1) ==
            L"{{\\tilde{a}^{{i_1}{i_2}}_{{a_1}{a_2}}}{\\tilde{a}^{{i_1}{i_2}}_{"
            L"{a_1^{{i_1}{i_2}}}{a_2^{{i_1}{i_2}}}}}}");
  }

  SECTION("commutativity") {
    auto nop1 = FNOperator(cre({L"i_1"}), ann({L"a_1"}));
    auto nop2 = FNOperator(cre({Index{L"i_1"}, Index{L"i_2"}}),
                           ann({Index{L"a_1", {L"i_1", L"i_2"}},
                                Index{L"a_2", {L"i_1", L"i_2"}}}));

    REQUIRE(nop1.commutes_with(nop1));
    REQUIRE(nop1.commutes_with(nop2));
    REQUIRE(nop2.commutes_with(nop1));
    REQUIRE(nop2.commutes_with(nop2));
    REQUIRE(!nop1.commutes_with(adjoint(FNOperator(nop1))));
    REQUIRE(!nop1.commutes_with(adjoint(FNOperator(nop2))));
    REQUIRE(!nop2.commutes_with(adjoint(FNOperator(nop1))));
    REQUIRE(!nop2.commutes_with(adjoint(FNOperator(nop2))));
    REQUIRE(adjoint(FNOperator(nop1)).commutes_with(adjoint(FNOperator(nop2))));
    REQUIRE(adjoint(FNOperator(nop2)).commutes_with(adjoint(FNOperator(nop1))));
  }

}  // TEST_CASE("Op")
