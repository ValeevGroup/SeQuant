// the PAO basis as a named instance of the particle space through the CSV-CCSD
// derivation: round trips, metadata and the committed fixture
#include <catch2/catch_test_macros.hpp>
#include "catch2_sequant.hpp"
#include "csv_test_utils.hpp"

#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/domain/mbpt/rules/csv.hpp>
#include <SeQuant/domain/mbpt/rules/df.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <algorithm>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

namespace sequant::tests::pao {
constexpr std::size_t n_pao = 300;

/// the CSV-CCSD registry with the PAO extent set; sized before a Context
/// adopts it, since a registry inside a Context is immutable
std::shared_ptr<IndexBasisRegistry> sized_registry() {
  auto isr = csv::csv_cc_registry();
  isr->approximate_size(L"μ̃", n_pao);
  return isr;
}

/// csv_cc_context over the sized registry, with tmp ordinals allowed in parsed
/// text
Context gate_context() {
  auto ctx = csv::csv_cc_context(sized_registry());
  ctx.set_first_dummy_index_ordinal(1000000);
  return ctx;
}

/// [E, R1, R2] in mpqc's order: CC{2}.t(2,0) with t granted 0, V1 spintrace,
/// DF, csv_transform(μ̃), flatten
std::vector<ExprPtr> derive() {
  auto core = set_scoped_default_context(csv::csv_cc_context(sized_registry()));
  auto mbpt_ctx = mbpt::set_scoped_default_mbpt_context(
      mbpt::Context({.csv = mbpt::CSV::Yes,
                     .op_registry_ptr = csv::granted_registry({{L"t", 0}})}));
  const auto& isr = get_default_context().index_basis_registry();
  Index::reset_tmp_index();
  std::vector<ExprPtr> out;
  for (auto const& eq : mbpt::CC{2}.t(2, 0)) {
    ExprPtr e = mbpt::closed_shell_CC_spintrace(
        eq, {.method = mbpt::BiorthogonalizationMethod::V1});
    e = mbpt::density_fit(e, isr->retrieve(L"Κ"), L"g", L"g");
    e = mbpt::csv_transform(e, isr->retrieve_basis(L"μ̃"));
    flatten(e);
    out.push_back(e);
  }
  return out;
}
}  // namespace sequant::tests::pao

TEST_CASE("pao-named-basis", "[mbpt][csv][basis][valgrind_skip]") {
  using namespace sequant;
  using namespace sequant::tests::pao;
  const auto eqs = derive();
  REQUIRE(eqs.size() == 3);

  const Context ctx = gate_context();
  const IndexBasis pao = ctx.index_basis_registry()->retrieve_basis(L"μ̃");
  auto scoped = set_scoped_default_context(ctx);

  SECTION("round trip and metadata") {
    for (auto const& e : eqs) CHECK(deserialize<ExprPtr>(serialize(e)) == e);
    CHECK(Index(L"μ̃_3") == Index(pao, 3));
    CHECK(Index(pao, 3).full_label() == L"μ̃_3");
    CHECK(Index(pao, 3).basis().base_key() == L"μ̃");
    CHECK(Index(pao, 3).basis().extent() == n_pao);
    CHECK(Index(pao, 3).basis().metric() == IndexSpaceMetric::General);
    // the bare (space, instance) pair that operator grants and the integral
    // projection build resolves to the basis as registered
    const IndexBasis bare{ctx.index_basis_registry()->retrieve(L"a"),
                          pao.basis_instance()};
    REQUIRE(bare.extent() != n_pao);
    CHECK(default_registry_resolved(bare) == pao);
    CHECK(default_registry_resolved(bare).extent() == n_pao);
    CHECK(ctx.index_basis_registry()->retrieve(L"a").approximate_size() !=
          n_pao);
    CHECK_THROWS_AS(ctx.index_basis_registry()->retrieve(L"μ̃"),
                    IndexBasisRegistry::not_a_space);
  }
  SECTION("the committed fixture: every μ̃ index is the named basis") {
    std::ifstream in(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                     "/data/csv_ccsd_doubles_residual_df.txt");
    std::string line;
    REQUIRE(std::getline(in, line));
    const ExprPtr fixture = deserialize<ExprPtr>(toUtf16(line));
    std::size_t n = 0;
    for (auto const& idx : get_used_indices(fixture))
      if (idx.basis().base_key() == L"μ̃") {
        ++n;
        CHECK(idx.basis() == pao);
        CHECK_FALSE(idx.has_proto_indices());
      }
    CHECK(n > 0);
  }
  SECTION("slot layout: the PAO basis sorts after every basis of a") {
    const IndexSpace a = ctx.index_basis_registry()->retrieve(L"a");
    CHECK((IndexBasis{a} < pao && IndexBasis{a, 0} < pao &&
           IndexBasis{a, 2} < pao && IndexBasis{a, 10} < pao));
    for (std::wstring text :
         {L"C{μ̃_1;a_1<i_1,i_2;0>}:N-C-S", L"C{μ̃_1;a_1<i_1>}:N-C-S"}) {
      const ExprPtr expr = deserialize<ExprPtr>(text);
      const EvalExpr node(expr->as<Tensor>());
      const auto slots = node.canon_indices();
      // the PAO leg is one slot, keyed as the named basis; the CSV leg and
      // its protos are the others
      CHECK(std::ranges::count_if(slots, [&pao](const Index& ix) {
              return ix.basis() == pao;
            }) == 1);
      CHECK(std::ranges::count_if(slots, [](const Index& ix) {
              return ix.has_proto_indices();
            }) == 1);
    }
  }
}
