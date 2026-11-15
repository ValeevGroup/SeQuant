//
// Created by Eduard Valeyev on 8/1/25.
//

#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include "SeQuant/core/expr.hpp"
#include "SeQuant/core/op.hpp"

#include <range/v3/algorithm/equal.hpp>
#include <range/v3/view/transform.hpp>

#include <array>
#include <span>

TEST_CASE("abstract_tensor", "[elements]") {
  using namespace sequant;

  std::map<SlotType, std::array<const wchar_t*, 3>> orig_idxs;
  std::array braidxs{L"p_1", L"p_2", L"p_3"};
  std::array ketidxs{L"p_4", L"p_5", L"p_6"};
  std::array auxidxs{L"p_7", L"p_8", L"p_9"};
  orig_idxs[SlotType::Bra] = braidxs;
  orig_idxs[SlotType::Ket] = ketidxs;
  orig_idxs[SlotType::Aux] = auxidxs;
  std::shared_ptr<AbstractTensor> tensor = std::make_shared<Tensor>(
      L"t", bra<decltype(braidxs)>{braidxs}, ket<decltype(ketidxs)>{ketidxs},
      aux<decltype(auxidxs)>{auxidxs});
  std::shared_ptr<AbstractTensor> nop = std::make_shared<FNOperator>(
      cre<decltype(ketidxs)>{ketidxs}, ann<decltype(braidxs)>{braidxs});

  SECTION("to_latex") {
    REQUIRE_NOTHROW(tensor->_to_latex());

    // can handle empty bra ket and aux bundles
    Tensor t("x", bra{}, ket{}, aux{});
    REQUIRE_NOTHROW(t.to_latex());
    auto _ = set_scoped_default_context(
        {.braket_slot_typesetting = BraKetSlotTypesetting::TensorPackage});
    REQUIRE_NOTHROW(t.to_latex());
  }

  SECTION("permute") {
    SECTION("1-index") {
      for (auto& tptr : {tensor, nop}) {
        auto& t = *tptr;

        // permute 1-index slots
        for (auto slottype : {SlotType::Bra, SlotType::Ket, SlotType::Aux}) {
          // skip aux for nop
          if (t._label() != L"t" && slottype == SlotType::Aux) continue;

          auto idxs = [&]() {
            switch (slottype) {
              case SlotType::Bra:
                return t._bra();
              case SlotType::Ket:
                return t._ket();
              case SlotType::Aux:
                return t._aux();
              default:
                abort();
            }
          };
          auto permute = [&](auto& p) {
            switch (slottype) {
              case SlotType::Bra:
                return t._permute_bra(p);
              case SlotType::Ket:
                return t._permute_ket(p);
              case SlotType::Aux:
                return t._permute_aux(p);
              default:
                abort();
            }
          };

          const std::array<std::size_t, 3> pdata{1, 0, 2};
          std::span p{pdata.data(), 3};
          // permutation moves indices around, i.e. no reallocation, pointers
          // are stable
          auto original_idxptrs =
              idxs() |
              ranges::views::transform([](const auto& idx) { return &idx; }) |
              ranges::to_vector;
          permute(p);
          for (auto i = 0; i != 3; ++i) {
            REQUIRE(idxs()[i].label() == orig_idxs[slottype][pdata[i]]);
          }
          auto idxptrs =
              idxs() |
              ranges::views::transform([](const auto& idx) { return &idx; }) |
              ranges::to_vector;
          REQUIRE(original_idxptrs == idxptrs);
        }
      }
    }  // SECTION("1-index")

    // permute 2-index (braket) slots
    SECTION("2-index") {
      for (auto& tptr : {tensor, nop}) {
        auto& t = *tptr;
        std::array<std::size_t, 3> pdata{1, 0, 2};
        std::span p{pdata.data(), 3};
        // permutation moves indices around, i.e. no reallocation, pointers are
        // stable
        auto original_idxptrs =
            t._braket() |
            ranges::views::transform([](const auto& idx) { return &idx; }) |
            ranges::to_vector;
        t._permute_columns(p);
        for (auto i = 0; i != 3; ++i) {
          REQUIRE(t._bra()[i].label() == braidxs[pdata[i]]);
          REQUIRE(t._ket()[i].label() == ketidxs[pdata[i]]);
        }
        auto idxptrs =
            t._braket() |
            ranges::views::transform([](const auto& idx) { return &idx; }) |
            ranges::to_vector;
        REQUIRE(original_idxptrs == idxptrs);
      }
    }  // SECTION("2-index")
  }

  SECTION("tensor_power") {
    auto D = ex<Tensor>(L"D", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                        TensorSymmetries{.perm = Symmetry::Symm,
                                         .braket = BraKetSymmetry::Symm,
                                         .column = ColumnSymmetry::Symm});
    auto p = ex<TensorPower>(D, rational{-1});
    AbstractTensor& at = p->as<AbstractTensor>();
    CHECK(ranges::equal(at._bra(), D->as<Tensor>().bra()));
    CHECK(ranges::equal(at._ket(), D->as<Tensor>().ket()));
    CHECK(at._bra_rank() == 2);
    CHECK(at._ket_rank() == 2);
    CHECK(at._aux_rank() == 0);
    CHECK(at._label() == L"D");
    CHECK(at._color() == D->as<AbstractTensor>()._color());
    CHECK(at._hash_value() == p->hash_value());
    const auto h0 = p->hash_value();
    // index replacement reaches the base and invalidates the power's hash
    container::map<Index, Index> repl{{Index{L"i_1"}, Index{L"i_3"}}};
    CHECK(transform_indices(at, repl));
    CHECK(p->hash_value() != h0);
    auto D3 = ex<Tensor>(L"D", bra{L"i_3", L"i_2"}, ket{L"a_1", L"a_2"},
                         TensorSymmetries{.perm = Symmetry::Symm,
                                          .braket = BraKetSymmetry::Symm,
                                          .column = ColumnSymmetry::Symm});
    CHECK(*p == TensorPower(D3, rational{-1}));
    // slot permutation through the mutable view, hash tracks it
    const auto h1 = p->hash_value();
    at._permute_bra(std::array<std::size_t, 2>{1, 0});
    CHECK(p->hash_value() != h1);
    CHECK(ranges::equal(at._bra(), std::vector{Index{L"i_2"}, Index{L"i_3"}}));
    // clone_shared is deep and round-trips to Expr (what
    // Product::canonicalize_impl does)
    auto cs = at._clone_shared();
    CHECK(std::dynamic_pointer_cast<Expr>(cs));
    CHECK(*std::dynamic_pointer_cast<Expr>(cs) == *p);
    CHECK(std::dynamic_pointer_cast<TensorPower>(cs)->base().get() !=
          p->as<TensorPower>().base().get());
  }
}
