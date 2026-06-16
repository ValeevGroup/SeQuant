//
// Created by Eduard Valeyev on 3/23/18.
//

#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/tag.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <cstddef>
#include <initializer_list>
#include <map>
#include <string>
#include <string_view>
#include <type_traits>

TEST_CASE("tensor", "[elements]") {
  using namespace sequant;

  SECTION("constructors") {
    REQUIRE_NOTHROW(Tensor{});
    auto t1 = Tensor{};
    REQUIRE(!t1);
    REQUIRE(t1.bra_rank() == 0);
    REQUIRE(t1.ket_rank() == 0);
    REQUIRE(t1.aux_rank() == 0);
    REQUIRE(t1.bra_net_rank() == 0);
    REQUIRE(t1.ket_net_rank() == 0);
    REQUIRE(t1.rank() == 0);
    REQUIRE(t1.symmetry() == Symmetry::Nonsymm);
    REQUIRE(t1.braket_symmetry() == BraKetSymmetry::Nonsymm);
    REQUIRE(t1.column_symmetry() == ColumnSymmetry::Nonsymm);
    REQUIRE(t1.label() == L"");

    REQUIRE_NOTHROW(Tensor(L"F", bra{L"i_1"}, ket{L"i_1"}));
    auto t2 = Tensor(L"F", bra{L"i_1"}, ket{L"i_1"});
    REQUIRE(t2);
    REQUIRE(t2.bra_rank() == 1);
    REQUIRE(t2.ket_rank() == 1);
    REQUIRE(t2.aux_rank() == 0);
    REQUIRE(t2.rank() == 1);
    REQUIRE(t2.const_braketaux().size() == 2);
    REQUIRE(t2.const_slots().size() == 2);
    REQUIRE(t2.num_slots() == 2);
    REQUIRE(t2.num_indices() == 2);
    REQUIRE(ranges::distance(t2.const_braketaux_indices().begin(),
                             t2.const_braketaux_indices().end()) == 2);
    REQUIRE(ranges::distance(t2.const_indices().begin(),
                             t2.const_indices().end()) == 2);
    REQUIRE(t2.symmetry() == Symmetry::Nonsymm);
    // bare ctor: braket symmetry is derived from the default Context's
    // Hermiticity (NonHermitian here) -> Nonsymm; column symmetry follows the
    // Context default (Symm in the test context)
    REQUIRE(t2.braket_symmetry() == BraKetSymmetry::Nonsymm);
    REQUIRE(t2.column_symmetry() == ColumnSymmetry::Symm);
    REQUIRE(t2.label() == L"F");

    // bra/kets of different rank
    REQUIRE_NOTHROW(Tensor(L"N", bra{L"i_1"}, ket{}, aux{L"a_1"}));
    auto t3 = Tensor(L"N", bra{L"i_1"}, ket{}, aux{L"a_1"});
    REQUIRE(t3);
    REQUIRE(t3.bra_rank() == 1);
    REQUIRE(t3.ket_rank() == 0);
    REQUIRE(t3.aux_rank() == 1);
    REQUIRE(t3.bra_net_rank() == 1);
    REQUIRE(t3.ket_net_rank() == 0);
    REQUIRE_THROWS(t3.rank());
    REQUIRE(t3.const_braketaux().size() == 2);
    REQUIRE(t3.const_slots().size() == 2);
    REQUIRE(t3.num_slots() == 2);
    REQUIRE(t3.num_indices() == 2);
    REQUIRE(ranges::distance(t3.const_braketaux_indices().begin(),
                             t3.const_braketaux_indices().end()) == 2);
    REQUIRE(ranges::distance(t3.const_indices().begin(),
                             t3.const_indices().end()) == 2);
    REQUIRE(t3.symmetry() == Symmetry::Nonsymm);
    REQUIRE(t3.braket_symmetry() == BraKetSymmetry::Nonsymm);
    REQUIRE(t3.column_symmetry() == ColumnSymmetry::Symm);
    REQUIRE(t3.label() == L"N");

    REQUIRE_NOTHROW(Tensor(L"g", bra{Index{L"i_1"}, Index{L"i_2"}},
                           ket{Index{L"i_3"}, Index{L"i_4"}},
                           aux{Index{L"i_5"}}, Symmetry::Nonsymm,
                           BraKetSymmetry::Symm, ColumnSymmetry::Nonsymm));
    auto t4 = Tensor(L"g", bra{Index{L"i_1"}, Index{L"i_2"}},
                     ket{Index{L"i_3"}, Index{L"i_4"}}, aux{Index{L"i_5"}},
                     Symmetry::Nonsymm, BraKetSymmetry::Symm,
                     ColumnSymmetry::Nonsymm);
    REQUIRE(t4);
    REQUIRE(t4.bra_rank() == 2);
    REQUIRE(t4.ket_rank() == 2);
    REQUIRE(t4.aux_rank() == 1);
    REQUIRE(t4.rank() == 2);
    REQUIRE(ranges::size(t4.const_braketaux()) == 5);
    REQUIRE(ranges::size(t4.const_slots()) == 5);
    REQUIRE(t4.num_slots() == 5);
    REQUIRE(t4.symmetry() == Symmetry::Nonsymm);
    REQUIRE(t4.braket_symmetry() == BraKetSymmetry::Symm);
    REQUIRE(t4.column_symmetry() == ColumnSymmetry::Nonsymm);
    REQUIRE(t4.label() == L"g");

    SECTION("null indices") {
      // null indices ok in asymmetric bra or ket
      REQUIRE_NOTHROW(Tensor(L"N", bra{L"i_2", L"", L"i_3"},
                             ket{L"", L"i_1", L""}, aux{}, Symmetry::Nonsymm));
      REQUIRE_NOTHROW(Tensor(L"N", bra{L"", L"i_1", L""},
                             ket{L"i_2", L"", L"i_3"}, aux{},
                             Symmetry::Nonsymm));

      REQUIRE_NOTHROW(
          Tensor(L"N", bra{L""}, ket{L"i_1"}, aux{}, Symmetry::Nonsymm));
      Tensor t(L"N", bra{L""}, ket{L"i_1"}, aux{}, Symmetry::Nonsymm);
      REQUIRE(t.bra_rank() == 0);
      REQUIRE(t.ket_rank() == 1);
      REQUIRE(t.bra_net_rank() == 0);
      REQUIRE(t.ket_net_rank() == 1);
      REQUIRE(t.num_slots() == 1);
      REQUIRE(t.num_indices() == 1);

      // in fact slots of asymmetric particle-symmetric tensors are kept in
      // canonical order
      REQUIRE(Tensor(L"N", bra{L""}, ket{L"i_1"}, aux{}, Symmetry::Nonsymm) ==
              Tensor(L"N", bra{}, ket{L"i_1"}, aux{}, Symmetry::Nonsymm));
      REQUIRE(Tensor(L"N", bra{L"i_2", L"", L"i_3"}, ket{L"", L"i_1", L""},
                     aux{}, Symmetry::Nonsymm) ==
              Tensor(L"N", bra{L"i_2", L"i_3"}, ket{L"", L"", L"i_1"}, aux{},
                     Symmetry::Nonsymm));
      REQUIRE(Tensor(L"N", bra{L"", L"i_1", L""}, ket{L"i_2", L"", L"i_3"},
                     aux{}, Symmetry::Nonsymm) ==
              Tensor(L"N", bra{L"i_1", L"", L""}, ket{L"", L"i_2", L"i_3"},
                     aux{}, Symmetry::Nonsymm));
      REQUIRE(
          Tensor(L"N", bra{L"", L"i_1", L"", L"i_4"},
                 ket{L"i_2", L"", L"i_3", L"i_5"}, aux{}, Symmetry::Nonsymm) ==
          Tensor(L"N", bra{L"i_4", L"i_1", L"", L""},
                 ket{L"i_5", L"", L"i_2", L"i_3"}, aux{}, Symmetry::Nonsymm));
      // in fact unnecessary null indices are dropped in canonicalization
      REQUIRE(
          Tensor(L"N", bra{L"", L"i_1", L"", L"i_4"},
                 ket{L"i_2", L"", L"i_3", L"i_5"}, aux{}, Symmetry::Nonsymm) ==
          Tensor(L"N", bra{L"i_4", L"i_1"}, ket{L"i_5", L"", L"i_2", L"i_3"},
                 aux{}, Symmetry::Nonsymm));
      Tensor t5(L"N", bra{L"", L"i_1", L"", L"i_4"},
                ket{L"i_2", L"", L"i_3", L"i_5"}, aux{}, Symmetry::Nonsymm);
      REQUIRE(t5.bra_rank() == 2);
      REQUIRE(t5.bra_net_rank() == 2);
      REQUIRE(t5.bra()[0] == L"i_4");
      REQUIRE(t5.bra()[1] == L"i_1");
      REQUIRE(t5.ket_rank() == 4);
      REQUIRE(t5.ket_net_rank() == 3);
      REQUIRE(t5.ket()[0] == L"i_5");
      REQUIRE(t5.ket()[1].nonnull() == false);
      REQUIRE(t5.ket()[2] == L"i_2");
      REQUIRE(t5.ket()[3] == L"i_3");

      // but asymmetric particle-NONsymmetric tensors are kept in given order
      REQUIRE(Tensor(L"N", bra{L""}, ket{L"i_1"}, aux{}, Symmetry::Nonsymm,
                     BraKetSymmetry::Conjugate, ColumnSymmetry::Nonsymm) !=
              Tensor(L"N", bra{}, ket{L"i_1"}, aux{}, Symmetry::Nonsymm,
                     BraKetSymmetry::Conjugate, ColumnSymmetry::Nonsymm));
      REQUIRE(Tensor(L"N", bra{L"", L"i_1", L"", L"i_4"},
                     ket{L"i_2", L"", L"i_3", L"i_5"}, aux{}, Symmetry::Nonsymm,
                     BraKetSymmetry::Conjugate, ColumnSymmetry::Nonsymm) !=
              Tensor(L"N", bra{L"i_4", L"i_1", L"", L""},
                     ket{L"i_5", L"", L"i_2", L"i_3"}, aux{}, Symmetry::Nonsymm,
                     BraKetSymmetry::Conjugate, ColumnSymmetry::Nonsymm));
      Tensor t6(L"N", bra{L"", L"i_1", L"", L"i_4"},
                ket{L"i_2", L"", L"i_3", L"i_5"}, aux{}, Symmetry::Nonsymm,
                BraKetSymmetry::Conjugate, ColumnSymmetry::Nonsymm);
      REQUIRE(t6.bra_rank() == 4);
      REQUIRE(t6.bra_net_rank() == 2);
      REQUIRE(t6.bra()[0].nonnull() == false);
      REQUIRE(t6.bra()[1] == L"i_1");
      REQUIRE(t6.bra()[2].nonnull() == false);
      REQUIRE(t6.bra()[3] == L"i_4");
      REQUIRE(t6.ket_rank() == 4);
      REQUIRE(t6.ket_net_rank() == 3);
      REQUIRE(t6.ket()[0] == L"i_2");
      REQUIRE(t6.ket()[1].nonnull() == false);
      REQUIRE(t6.ket()[2] == L"i_3");
      REQUIRE(t6.ket()[3] == L"i_5");
      REQUIRE(t6.num_slots() == 8);
      REQUIRE(t6.num_indices() == 5);

      // check errors
      if (sequant::assert_behavior() == sequant::AssertBehavior::Throw) {
        // no null indices in antisymmetric bra or ket
        REQUIRE_THROWS_AS(
            Tensor(L"N", bra{L""}, ket{L"i_1"}, aux{}, Symmetry::Antisymm),
            sequant::Exception);
        REQUIRE_THROWS_AS(
            Tensor(L"N", bra{L"i_1"}, ket{L""}, aux{}, Symmetry::Antisymm),
            sequant::Exception);

        // no null indices in symmetric bra or ket
        REQUIRE_THROWS_AS(
            Tensor(L"N", bra{L""}, ket{L"i_1"}, aux{}, Symmetry::Symm),
            sequant::Exception);
        REQUIRE_THROWS_AS(
            Tensor(L"N", bra{L"i_1"}, ket{L""}, aux{}, Symmetry::Symm),
            sequant::Exception);

        // no paired null indices in asymmetric bra/ket
        REQUIRE_THROWS_AS(Tensor(L"N", bra{L"i_2", L""}, ket{L"i_1", L""},
                                 aux{}, Symmetry::Nonsymm),
                          sequant::Exception);

        // no unpaired null indices in asymmetric bra/ket
        REQUIRE_THROWS_AS(Tensor(L"N", bra{L"i_2"}, ket{L"i_1", L""}, aux{},
                                 Symmetry::Nonsymm),
                          sequant::Exception);
        REQUIRE_THROWS_AS(Tensor(L"N", bra{L"i_2", L""}, ket{L"i_1"}, aux{},
                                 Symmetry::Nonsymm),
                          sequant::Exception);

        // no null aux indices
        REQUIRE_THROWS_AS(Tensor(L"N", bra{L"i_1"}, ket{}, aux{L""}),
                          sequant::Exception);
      }
    }

    SECTION("duplicate indices") {
      // no duplicate indices allowed in bra, in ket, or in aux
      // but duplicates OK across the bundles
      REQUIRE_NOTHROW(Tensor(L"N", bra{L"i_1"}, ket{L"i_1"}, aux{L"i_2"}));
      // null indices are ignored in duplicate checks
      REQUIRE_NOTHROW(Tensor(L"N", bra{L"", L"", L"i_1"},
                             ket{L"i_1", L"i_2", L""}, aux{L"i_3"}));
      if (sequant::assert_behavior() == sequant::AssertBehavior::Throw) {
        REQUIRE_THROWS_AS(Tensor(L"N", bra{L"i_1", L"i_1"}, ket{}, aux{}),
                          sequant::Exception);
        REQUIRE_THROWS_AS(Tensor(L"N", bra{}, ket{L"i_1", L"i_1"}, aux{}),
                          sequant::Exception);
        REQUIRE_THROWS_AS(Tensor(L"N", bra{}, ket{}, aux{L"i_1", L"i_1"}),
                          sequant::Exception);
      }
    }
  }  // SECTION("constructors")

  SECTION("index transformation") {
    auto t = Tensor(L"g", bra{Index{L"i_1"}, Index{L"i_2"}},
                    ket{Index{L"i_3"}, Index{L"i_4"}}, aux{Index{L"i_5"}},
                    Symmetry::Antisymm);
    // index renaming preserves slot pointers
    const auto* slot_b0 = &(t.bra()[0]);
    const auto* slot_b1 = &(t.bra()[1]);
    const auto* slot_k0 = &(t.ket()[0]);
    const auto* slot_k1 = &(t.ket()[1]);
    const auto b0 = t.bra()[0];
    const auto b1 = t.bra()[1];
    std::map<Index, Index> idxmap = {{Index{L"i_1"}, Index{L"i_2"}},
                                     {Index{L"i_2"}, Index{L"i_1"}}};
    REQUIRE(t.transform_indices(idxmap));
    REQUIRE(slot_b0 == &(t.bra()[0]));
    REQUIRE(slot_b1 == &(t.bra()[1]));
    REQUIRE(slot_k0 == &(t.ket()[0]));
    REQUIRE(slot_k1 == &(t.ket()[1]));
    REQUIRE(b1 == t.bra()[0]);
    REQUIRE(b0 == t.bra()[1]);

    REQUIRE(t.bra()[0].tag().has_value());
    const auto t_bra0_tag_value = t.bra()[0].tag().value<int>();
    REQUIRE(t_bra0_tag_value == 0);
    REQUIRE(t.bra()[1].tag().has_value());
    const auto t_bra1_tag_value = t.bra()[1].tag().value<int>();
    REQUIRE(t_bra1_tag_value == 0);
    REQUIRE(!t.ket()[0].tag().has_value());
    REQUIRE(!t.ket()[1].tag().has_value());
    REQUIRE(t == Tensor(L"g", bra{Index{L"i_2"}, Index{L"i_1"}},
                        ket{Index{L"i_3"}, Index{L"i_4"}}, aux{Index{L"i_5"}},
                        Symmetry::Antisymm));
    // tagged indices are protected, so no replacements the second goaround
    REQUIRE(!t.transform_indices(idxmap));
    t.reset_tags();
    REQUIRE(t.transform_indices(idxmap));
    REQUIRE(t == Tensor(L"g", bra{Index{L"i_1"}, Index{L"i_2"}},
                        ket{Index{L"i_3"}, Index{L"i_4"}}, aux{Index{L"i_5"}},
                        Symmetry::Antisymm));
    t.reset_tags();
    REQUIRE(!t.bra()[0].tag().has_value());
    REQUIRE(!t.bra()[1].tag().has_value());
    REQUIRE(!t.ket()[0].tag().has_value());
    REQUIRE(!t.ket()[1].tag().has_value());

    SECTION("proto indices") {
      Tensor tensor = deserialize(L"g{i2,a1<i1>;a2<i2>,i1}")->as<Tensor>();

      const auto* slot_b0 = &(tensor.bra()[0]);
      const auto* slot_b1 = &(tensor.bra()[1]);
      const auto* slot_k0 = &(tensor.ket()[0]);
      const auto* slot_k1 = &(tensor.ket()[1]);

      // Swap columns of g
      std::map<Index, Index> idxmap = {
          {Index{L"i_2"}, Index{L"a_1", {L"i_1"}}},
          {Index{L"a_1", {L"i_1"}}, Index{L"i_2"}},
          {Index{L"a_2", {L"i_2"}}, Index{L"i_1"}},
          {Index{L"i_1"}, Index{L"a_2", {L"i_2"}}},
      };

      const Tensor expected =
          deserialize(L"g{a1<i1>,i2;i1,a2<i2>}")->as<Tensor>();
      tensor.transform_indices(idxmap);

      REQUIRE(tensor == expected);

      // verify slot stability
      REQUIRE(slot_b0 == &(tensor.bra()[0]));
      REQUIRE(slot_b1 == &(tensor.bra()[1]));
      REQUIRE(slot_k0 == &(tensor.ket()[0]));
      REQUIRE(slot_k1 == &(tensor.ket()[1]));
    }
  }  // SECTION("index transformation")

  SECTION("hash") {
    auto t1 = Tensor(L"F", bra{L"i_1"}, ket{L"i_2"});
    REQUIRE_NOTHROW(hash_value(t1));
    size_t t1_hash = hash_value(t1);

    auto t2 = Tensor(L"F", bra{L"i_2"}, ket{L"i_1"});
    REQUIRE_NOTHROW(hash_value(t2));
    size_t t2_hash = hash_value(t2);

    REQUIRE(t1_hash != t2_hash);

    auto t3 = Tensor(L"F", bra{L"i_2"}, ket{L"i_1"}, aux{L"i_3"});
    REQUIRE_NOTHROW(hash_value(t3));
    size_t t3_hash = hash_value(t3);
    REQUIRE(t2_hash != t3_hash);
    REQUIRE(t1_hash != t3_hash);

  }  // SECTION("hash")

  SECTION("latex") {
    auto t1 = Tensor(L"F", bra{L"i_1"}, ket{L"i_2"});
    auto t2 = Tensor(L"F", bra{L"i_1"}, ket{L"i_2"}, aux{L"i_3"});
    auto t3 = Tensor(L"F", bra{L"i_1"}, ket{L"i_2"}, aux{L"i_3", L"i_4"});
    auto t4 = Tensor(L"F", bra{Index(L"i_1", {L"i_5", L"i_6"}), Index{}},
                     ket{L"", L"i_2"}, aux{L"i_3", L"i_4"}, Symmetry::Nonsymm);
    auto h1 = ex<Tensor>(L"F", bra{L"i_1"}, ket{L"i_2"}) *
              ex<FNOperator>(cre({L"i_1"}), ann({L"i_2"}));

    SECTION("default (brasub, naive) typesetting") {
      REQUIRE(to_latex(t1) == L"{F^{{i_2}}_{{i_1}}}");
      REQUIRE(to_latex(t2) == L"{F^{{i_2}}_{{i_1}}[{i_3}]}");
      REQUIRE(to_latex(t3) == L"{F^{{i_2}}_{{i_1}}[{i_3},{i_4}]}");
      REQUIRE(to_latex(t4) ==
              L"{F^{\\textvisiblespace{i_2}}_{{i_1^{{i_5}{i_6}}}"
              L"}[{i_3},{i_4}]}");
      REQUIRE(to_latex(h1) ==
              L"{{F^{{i_2}}_{{i_1}}}{\\tilde{a}^{{i_1}}_{{i_2}}}}");
    }

    SECTION("ketsub naive typesetting") {
      Context ctx = get_default_context();
      REQUIRE(get_default_context().braket_typesetting() ==
              BraKetTypesetting::KetSuper);
      ctx.set(BraKetTypesetting::KetSub);
      auto resetter = set_scoped_default_context(ctx);
      REQUIRE(get_default_context().braket_typesetting() ==
              BraKetTypesetting::KetSub);

      REQUIRE(to_latex(t1) == L"{F_{{i_2}}^{{i_1}}}");
      REQUIRE(to_latex(t2) == L"{F_{{i_2}}^{{i_1}}[{i_3}]}");
      REQUIRE(to_latex(t3) == L"{F_{{i_2}}^{{i_1}}[{i_3},{i_4}]}");
      REQUIRE(to_latex(t4) ==
              L"{F_{\\textvisiblespace{i_2}}^{{i_1^{{i_5}{i_6}}}"
              L"}[{i_3},{i_4}]}");
      REQUIRE(to_latex(h1) ==
              L"{{F_{{i_2}}^{{i_1}}}{\\tilde{a}_{{i_1}}^{{i_2}}}}");
    }

    SECTION("brasub tensor typesetting") {
      Context ctx = get_default_context();
      REQUIRE(get_default_context().braket_typesetting() ==
              BraKetTypesetting::BraSub);
      REQUIRE(get_default_context().braket_slot_typesetting() ==
              BraKetSlotTypesetting::Naive);
      ctx.set(BraKetSlotTypesetting::TensorPackage);
      auto resetter = set_scoped_default_context(ctx);
      REQUIRE(get_default_context().braket_typesetting() ==
              BraKetTypesetting::BraSub);
      REQUIRE(get_default_context().braket_slot_typesetting() ==
              BraKetSlotTypesetting::TensorPackage);

      REQUIRE(to_latex(t1) == L"{\\tensor*{F}{*^{i_2}_{i_1}}}");
      REQUIRE(to_latex(t2) == L"{\\tensor*{F}{*^{i_2}_{i_1}}[{i_3}]}");
      REQUIRE(to_latex(t3) == L"{\\tensor*{F}{*^{i_2}_{i_1}}[{i_3},{i_4}]}");
      REQUIRE(
          to_latex(t4) ==
          L"{\\tensor*{F}{*^{}_{i_1^{{i_5}{i_6}}}*^{i_2}_{}}[{i_3},{i_4}]}");
      REQUIRE(to_latex(h1) ==
              L"{{\\tensor*{F}{*^{i_2}_{i_1}}}{\\tensor*{\\tilde{a}}{*^{i_1}_{"
              L"i_2}}}}");
    }

    SECTION("ketsub tensor typesetting") {
      Context ctx = get_default_context();
      REQUIRE(get_default_context().braket_typesetting() ==
              BraKetTypesetting::BraSub);
      REQUIRE(get_default_context().braket_slot_typesetting() ==
              BraKetSlotTypesetting::Naive);
      ctx.set(BraKetTypesetting::KetSub);
      ctx.set(BraKetSlotTypesetting::TensorPackage);
      auto resetter = set_scoped_default_context(ctx);
      REQUIRE(get_default_context().braket_typesetting() ==
              BraKetTypesetting::KetSub);
      REQUIRE(get_default_context().braket_slot_typesetting() ==
              BraKetSlotTypesetting::TensorPackage);

      REQUIRE(to_latex(t1) == L"{\\tensor*{F}{*^{i_1}_{i_2}}}");
      REQUIRE(to_latex(t2) == L"{\\tensor*{F}{*^{i_1}_{i_2}}[{i_3}]}");
      REQUIRE(to_latex(t3) == L"{\\tensor*{F}{*^{i_1}_{i_2}}[{i_3},{i_4}]}");
      REQUIRE(
          to_latex(t4) ==
          L"{\\tensor*{F}{*^{i_1^{{i_5}{i_6}}}_{}*^{}_{i_2}}[{i_3},{i_4}]}");
      REQUIRE(to_latex(h1) ==
              L"{{\\tensor*{F}{*^{i_1}_{i_2}}}{\\tensor*{\\tilde{a}}{*^{i_2}_{"
              L"i_1}}}}");
    }

  }  // SECTION("latex")

  SECTION("adjoint") {
    auto f1 = Tensor(L"F", bra{L"i_1", L"i_2"}, ket{L"i_3", L"i_4"});
    REQUIRE_NOTHROW(f1.adjoint());
    // F is now non-Hermitian by default (braket Nonsymm), so its adjoint is
    // marked with the conjugation superscript
    REQUIRE(to_latex(f1) == L"{F⁺^{{i_1}{i_2}}_{{i_3}{i_4}}}");

    auto t1 = Tensor(L"t", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm,
                     BraKetSymmetry::Nonsymm);
    REQUIRE_NOTHROW(t1.adjoint());
    REQUIRE(to_latex(t1) == L"{t⁺^{{a_1}}_{{i_1}}}");
    t1.adjoint();
    REQUIRE(to_latex(t1) == L"{t^{{i_1}}_{{a_1}}}");

    auto h1 = ex<Tensor>(L"F", bra{L"i_1"}, ket{L"i_2"}) *
              ex<FNOperator>(cre{L"i_1"}, ann{L"i_2"});
    h1 = adjoint(h1);
    REQUIRE(to_latex(h1) ==
            L"{{\\tilde{a}^{{i_2}}_{{i_1}}}{F⁺^{{i_1}}_{{i_2}}}}");
    h1 = adjoint(h1);
    REQUIRE(to_latex(h1) ==
            L"{{F^{{i_2}}_{{i_1}}}{\\tilde{a}^{{i_1}}_{{i_2}}}}");

  }  // SECTION("adjoint")
}

TEST_CASE("tensor_hermiticity", "[elements]") {
  using namespace sequant;

  SECTION("field<->hermiticity helpers") {
    REQUIRE(to_braket_symmetry(Hermiticity::Hermitian, Field::Real) ==
            BraKetSymmetry::Symm);
    REQUIRE(to_braket_symmetry(Hermiticity::Hermitian, Field::Complex) ==
            BraKetSymmetry::Conjugate);
    REQUIRE(to_braket_symmetry(Hermiticity::NonHermitian, Field::Real) ==
            BraKetSymmetry::Nonsymm);
    REQUIRE(to_braket_symmetry(Hermiticity::NonHermitian, Field::Complex) ==
            BraKetSymmetry::Nonsymm);
    // AntiHermitian cannot yet be represented in BraKetSymmetry -> Nonsymm
    REQUIRE(to_braket_symmetry(Hermiticity::AntiHermitian, Field::Real) ==
            BraKetSymmetry::Nonsymm);

    REQUIRE(to_hermiticity(BraKetSymmetry::Symm) == Hermiticity::Hermitian);
    REQUIRE(to_hermiticity(BraKetSymmetry::Conjugate) ==
            Hermiticity::Hermitian);
    REQUIRE(to_hermiticity(BraKetSymmetry::Nonsymm) ==
            Hermiticity::NonHermitian);
  }

  SECTION("braket_symmetry is derived from Hermiticity + base_field") {
    // an index whose space is forced to a given Field (field is non-identity
    // IndexSpace metadata, so the index still denotes the same space)
    auto idx = [](std::wstring_view label, Field f) {
      Index i{label};
      IndexSpace space = i.space();
      space.field(f);
      return Index(space, i.ordinal());
    };
    // an explicitly non-Hermitian tensor is bra-ket nonsymmetric regardless of
    // the field
    {
      auto t = Tensor(L"t", bra{idx(L"i_1", Field::Real)},
                      ket{idx(L"a_1", Field::Real)}, Symmetry::Nonsymm,
                      Hermiticity::NonHermitian);
      REQUIRE(t.hermiticity() == Hermiticity::NonHermitian);
      REQUIRE(t.braket_symmetry() == BraKetSymmetry::Nonsymm);
    }
    // a Hermitian tensor over real spaces is bra-ket Symm ...
    {
      auto g = Tensor(L"g", bra{idx(L"i_1", Field::Real)},
                      ket{idx(L"i_2", Field::Real)}, Symmetry::Nonsymm,
                      Hermiticity::Hermitian);
      REQUIRE(g.hermiticity() == Hermiticity::Hermitian);
      REQUIRE(g.base_field() == Field::Real);
      REQUIRE(g.braket_symmetry() == BraKetSymmetry::Symm);
    }
    // ... and bra-ket Conjugate over complex spaces (the default)
    {
      auto g = Tensor(L"g", bra{idx(L"i_1", Field::Complex)},
                      ket{idx(L"i_2", Field::Complex)}, Symmetry::Nonsymm,
                      Hermiticity::Hermitian);
      REQUIRE(g.base_field() == Field::Complex);
      REQUIRE(g.braket_symmetry() == BraKetSymmetry::Conjugate);
    }
    // base_field ORs the spaces: any complex space makes the tensor complex
    {
      auto g = Tensor(L"g", bra{idx(L"i_1", Field::Real)},
                      ket{idx(L"i_2", Field::Complex)}, Symmetry::Nonsymm,
                      Hermiticity::Hermitian);
      REQUIRE(g.base_field() == Field::Complex);
      REQUIRE(g.braket_symmetry() == BraKetSymmetry::Conjugate);
    }
    // AntiHermitian is recorded but currently derives to Nonsymm
    {
      auto w = Tensor(L"w", bra{idx(L"i_1", Field::Real)},
                      ket{idx(L"i_2", Field::Real)}, Symmetry::Nonsymm,
                      Hermiticity::AntiHermitian);
      REQUIRE(w.hermiticity() == Hermiticity::AntiHermitian);
      REQUIRE(w.braket_symmetry() == BraKetSymmetry::Nonsymm);
    }
  }

  SECTION("legacy BraKetSymmetry ctor back-fills Hermiticity") {
    auto g = Tensor(L"g", bra{L"i_1"}, ket{L"i_2"}, Symmetry::Nonsymm,
                    BraKetSymmetry::Symm);
    REQUIRE(g.braket_symmetry() == BraKetSymmetry::Symm);
    REQUIRE(g.hermiticity() == Hermiticity::Hermitian);
    auto t = Tensor(L"t", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Nonsymm,
                    BraKetSymmetry::Nonsymm);
    REQUIRE(t.hermiticity() == Hermiticity::NonHermitian);
  }

  // The payoff: over real spaces a Hermitian integral is bra<->ket symmetric,
  // so g_{i1}^{i2} and the bra<->ket-swapped g_{i2}^{i1} canonicalize to the
  // same tensor and cancel; over complex spaces (Conjugate) they do not.
  SECTION("bra-ket-symmetric integral collapses swapped contraction") {
    auto idx = [](std::wstring_view label, Field f) {
      Index i{label};
      IndexSpace space = i.space();
      space.field(f);
      return Index(space, i.ordinal());
    };
    auto make_diff = [&idx](Field f) {
      auto a = ex<Tensor>(L"g", bra{idx(L"i_1", f)}, ket{idx(L"i_2", f)},
                          Symmetry::Nonsymm, Hermiticity::Hermitian);
      auto b = ex<Tensor>(L"g", bra{idx(L"i_2", f)}, ket{idx(L"i_1", f)},
                          Symmetry::Nonsymm, Hermiticity::Hermitian);
      ExprPtr diff = a - b;
      simplify(diff);
      return diff;
    };

    // real spaces -> Symm -> the swapped copies are equal and cancel
    REQUIRE(make_diff(Field::Real) == ex<Constant>(0));
    // complex spaces -> Conjugate -> the two orientations remain distinct
    REQUIRE_FALSE(make_diff(Field::Complex) == ex<Constant>(0));
  }

  SECTION("default symmetries are taken from the Context") {
    // unspecified symmetries are resolved against the active default Context;
    // braket symmetry is *derived* from the Context's default Hermiticity and
    // the tensor's base field (there is no braket-symmetry Context knob)
    auto ctx = get_default_context();
    ctx.set(Symmetry::Nonsymm)
        .set(Hermiticity::Hermitian)
        .set(ColumnSymmetry::Nonsymm);
    auto resetter = set_scoped_default_context(ctx);

    auto g = Tensor(L"g", bra{L"i_1"}, ket{L"i_2"});
    REQUIRE(g.symmetry() == Symmetry::Nonsymm);
    REQUIRE(g.hermiticity() == Hermiticity::Hermitian);
    REQUIRE(g.column_symmetry() == ColumnSymmetry::Nonsymm);
    REQUIRE(g.braket_symmetry() ==
            to_braket_symmetry(Hermiticity::Hermitian, g.base_field()));

    // explicitly-specified symmetries override the Context defaults
    auto t = Tensor(L"t", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm,
                    BraKetSymmetry::Nonsymm, ColumnSymmetry::Symm);
    REQUIRE(t.braket_symmetry() == BraKetSymmetry::Nonsymm);
    REQUIRE(t.hermiticity() == Hermiticity::NonHermitian);
    REQUIRE(t.column_symmetry() == ColumnSymmetry::Symm);
  }
}
