//
// Created by Eduard Valeyev on 3/23/18.
//

#include <catch2/catch_test_macros.hpp>

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/latex.hpp>
#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/tag.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <cstddef>
#include <initializer_list>
#include <map>
#include <string>
#include <string_view>
#include <type_traits>

TEST_CASE("Tensor", "[elements]") {
  using namespace sequant;

  SECTION("constructors") {
    REQUIRE_NOTHROW(Tensor{});
    auto t1 = Tensor{};
    REQUIRE(!t1);
    REQUIRE(t1.bra_rank() == 0);
    REQUIRE(t1.ket_rank() == 0);
    REQUIRE(t1.aux_rank() == 0);
    REQUIRE(t1.rank() == 0);
    REQUIRE(t1.symmetry() == Symmetry::invalid);
    REQUIRE(t1.braket_symmetry() == BraKetSymmetry::invalid);
    REQUIRE(t1.particle_symmetry() == ParticleSymmetry::invalid);
    REQUIRE(t1.label() == L"");

    REQUIRE_NOTHROW(Tensor(L"F", bra{L"i_1"}, ket{L"i_1"}));
    auto t2 = Tensor(L"F", bra{L"i_1"}, ket{L"i_1"});
    REQUIRE(t2);
    REQUIRE(t2.bra_rank() == 1);
    REQUIRE(t2.ket_rank() == 1);
    REQUIRE(t2.aux_rank() == 0);
    REQUIRE(t2.rank() == 1);
    REQUIRE(t2.const_indices().size() == 2);
    REQUIRE(t2.symmetry() == Symmetry::nonsymm);
    REQUIRE(t2.braket_symmetry() == BraKetSymmetry::conjugate);
    REQUIRE(t2.particle_symmetry() == ParticleSymmetry::symm);
    REQUIRE(t2.label() == L"F");

    REQUIRE_NOTHROW(Tensor(L"N", bra{L"i_1"}, ket{}, aux{L"a_1"}));
    auto t3 = Tensor(L"N", bra{L"i_1"}, ket{}, aux{L"a_1"});
    REQUIRE(t3);
    REQUIRE(t3.bra_rank() == 1);
    REQUIRE(t3.ket_rank() == 0);
    REQUIRE(t3.aux_rank() == 1);
    REQUIRE_THROWS(t3.rank());
    REQUIRE(t3.const_indices().size() == 2);
    REQUIRE(t3.symmetry() == Symmetry::nonsymm);
    REQUIRE(t3.braket_symmetry() == BraKetSymmetry::conjugate);
    REQUIRE(t3.particle_symmetry() == ParticleSymmetry::symm);
    REQUIRE(t3.label() == L"N");

    REQUIRE_NOTHROW(Tensor(L"g", bra{Index{L"i_1"}, Index{L"i_2"}},
                           ket{Index{L"i_3"}, Index{L"i_4"}},
                           aux{Index{L"i_5"}}, Symmetry::nonsymm,
                           BraKetSymmetry::symm, ParticleSymmetry::nonsymm));
    auto t4 = Tensor(L"g", bra{Index{L"i_1"}, Index{L"i_2"}},
                     ket{Index{L"i_3"}, Index{L"i_4"}}, aux{Index{L"i_5"}},
                     Symmetry::nonsymm, BraKetSymmetry::symm,
                     ParticleSymmetry::nonsymm);
    REQUIRE(t4);
    REQUIRE(t4.bra_rank() == 2);
    REQUIRE(t4.ket_rank() == 2);
    REQUIRE(t4.aux_rank() == 1);
    REQUIRE(t4.rank() == 2);
    REQUIRE(t4.const_indices().size() == 5);
    REQUIRE(t4.symmetry() == Symmetry::nonsymm);
    REQUIRE(t4.braket_symmetry() == BraKetSymmetry::symm);
    REQUIRE(t4.particle_symmetry() == ParticleSymmetry::nonsymm);
    REQUIRE(t4.label() == L"g");
  }  // SECTION("constructors")

  SECTION("index transformation") {
    auto t = Tensor(L"g", bra{Index{L"i_1"}, Index{L"i_2"}},
                    ket{Index{L"i_3"}, Index{L"i_4"}}, aux{Index{L"i_5"}},
                    Symmetry::antisymm);
    std::map<Index, Index> idxmap = {{Index{L"i_1"}, Index{L"i_2"}},
                                     {Index{L"i_2"}, Index{L"i_1"}}};
    REQUIRE(t.transform_indices(idxmap));
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
                        Symmetry::antisymm));
    // tagged indices are protected, so no replacements the second goaround
    REQUIRE(!t.transform_indices(idxmap));
    t.reset_tags();
    REQUIRE(t.transform_indices(idxmap));
    REQUIRE(t == Tensor(L"g", bra{Index{L"i_1"}, Index{L"i_2"}},
                        ket{Index{L"i_3"}, Index{L"i_4"}}, aux{Index{L"i_5"}},
                        Symmetry::antisymm));
    t.reset_tags();
    REQUIRE(!t.bra()[0].tag().has_value());
    REQUIRE(!t.bra()[1].tag().has_value());
    REQUIRE(!t.ket()[0].tag().has_value());
    REQUIRE(!t.ket()[1].tag().has_value());
  }  // SECTION("index transformation")

  SECTION("hash") {
    auto t1 = Tensor(L"F", bra{L"i_1"}, ket{L"i_2"});
    size_t t1_hash;
    REQUIRE_NOTHROW(t1_hash = hash_value(t1));

    auto t2 = Tensor(L"F", bra{L"i_2"}, ket{L"i_1"});
    size_t t2_hash;
    REQUIRE_NOTHROW(t2_hash = hash_value(t2));
    REQUIRE(t1_hash != t2_hash);

    auto t3 = Tensor(L"F", bra{L"i_2"}, ket{L"i_1"}, aux{L"i_3"});
    size_t t3_hash;
    REQUIRE_NOTHROW(t3_hash = hash_value(t3));
    REQUIRE(t2_hash != t3_hash);
    REQUIRE(t1_hash != t3_hash);

  }  // SECTION("hash")

  SECTION("latex") {
    auto t1 = Tensor(L"F", bra{L"i_1"}, ket{L"i_2"});
    REQUIRE(to_latex(t1) == L"{F^{{i_2}}_{{i_1}}}");

    auto t2 = Tensor(L"F", bra{L"i_1"}, ket{L"i_2"}, aux{L"i_3"});
    REQUIRE(to_latex(t2) == L"{F^{{i_2}}_{{i_1}}[{i_3}]}");

    auto t3 = Tensor(L"F", bra{L"i_1"}, ket{L"i_2"}, aux{L"i_3", L"i_4"});
    REQUIRE(to_latex(t3) == L"{F^{{i_2}}_{{i_1}}[{i_3},{i_4}]}");

    auto h1 = ex<Tensor>(L"F", bra{L"i_1"}, ket{L"i_2"}) *
              ex<FNOperator>(cre({L"i_1"}), ann({L"i_2"}));
    REQUIRE(to_latex(h1) ==
            L"{{F^{{i_2}}_{{i_1}}}{\\tilde{a}^{{i_1}}_{{i_2}}}}");

  }  // SECTION("latex")

  SECTION("adjoint") {
    auto f1 = Tensor(L"F", bra{L"i_1", L"i_2"}, ket{L"i_3", L"i_4"});
    REQUIRE_NOTHROW(f1.adjoint());
    REQUIRE(to_latex(f1) == L"{F^{{i_1}{i_2}}_{{i_3}{i_4}}}");

    auto t1 = Tensor(L"t", bra{L"a_1"}, ket{L"i_1"}, Symmetry::nonsymm,
                     BraKetSymmetry::nonsymm);
    REQUIRE_NOTHROW(t1.adjoint());
    REQUIRE(to_latex(t1) == L"{t⁺^{{a_1}}_{{i_1}}}");
    t1.adjoint();
    REQUIRE(to_latex(t1) == L"{t^{{i_1}}_{{a_1}}}");

    auto h1 = ex<Tensor>(L"F", bra{L"i_1"}, ket{L"i_2"}) *
              ex<FNOperator>(cre{L"i_1"}, ann{L"i_2"});
    h1 = adjoint(h1);
    REQUIRE(to_latex(h1) ==
            L"{{\\tilde{a}^{{i_2}}_{{i_1}}}{F^{{i_1}}_{{i_2}}}}");
    h1 = adjoint(h1);
    REQUIRE(to_latex(h1) ==
            L"{{F^{{i_2}}_{{i_1}}}{\\tilde{a}^{{i_1}}_{{i_2}}}}");

  }  // SECTION("adjoint")

}  // TEST_CASE("Tensor")
