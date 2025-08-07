//
// Created by Eduard Valeyev on 8/1/25.
//

#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include "SeQuant/core/op.hpp"
#include "SeQuant/core/tensor.hpp"

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
        const std::array<std::size_t, 3> pdata{1, 0, 2};
        std::span p{pdata.data(), 3};
        // permutation moves indices around, i.e. no reallocation, pointers are
        // stable
        auto original_idxptrs =
            t._braket() |
            ranges::views::transform([](const auto& idx) { return &idx; }) |
            ranges::to_vector;
        t._permute_braket(p);
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
}
