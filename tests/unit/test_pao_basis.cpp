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
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace sequant::tests::pao {
constexpr std::size_t n_pao = 300;

/// the CSV-CCSD registry with the AO (and so the PAO) extent set; sized before
/// a Context adopts it, since a registry inside a Context is immutable
std::shared_ptr<IndexBasisRegistry> sized_registry(
    std::optional<Field> field = {}) {
  auto ibr = csv::csv_cc_registry();
  ibr->extent(L"μ", n_pao);
  if (field) {  // as mpqc::scoped_sequant_field
    std::vector<std::wstring> keys;
    for (auto const& s : *ibr) keys.push_back(s.base_key());
    for (auto const& k : keys)
      if (!k.empty()) ibr->field(k, *field);
  }
  return ibr;
}

/// csv_cc_context over the sized registry, with tmp ordinals allowed in parsed
/// text
Context gate_context(std::optional<Field> field = {}) {
  auto ctx = csv::csv_cc_context(sized_registry(field));
  ctx.set_first_dummy_index_ordinal(1000000);
  return ctx;
}

/// [E, R1, R2] in mpqc's order: CC{2}.t(2,0) with t granted 0, V1 spintrace,
/// DF, csv_transform(μ̃), flatten
std::vector<ExprPtr> derive(std::optional<Field> field = {}) {
  auto core =
      set_scoped_default_context(csv::csv_cc_context(sized_registry(field)));
  auto mbpt_ctx = mbpt::set_scoped_default_mbpt_context(
      mbpt::Context({.csv = mbpt::CSV::Yes,
                     .op_registry_ptr = csv::granted_registry({{L"t", 0}})}));
  const auto& ibr = get_default_context().index_basis_registry();
  Index::reset_tmp_index();
  std::vector<ExprPtr> out;
  for (auto const& eq : mbpt::CC{2}.t(2, 0)) {
    ExprPtr e = mbpt::closed_shell_CC_spintrace(
        eq, {.method = mbpt::BiorthogonalizationMethod::V1});
    e = mbpt::density_fit(e, ibr->retrieve_space(L"Κ"), L"g", L"g");
    e = mbpt::csv_transform(e, ibr->retrieve(L"μ̃"));
    flatten(e);
    out.push_back(e);
  }
  return out;
}
}  // namespace sequant::tests::pao

TEST_CASE("pao-named-basis", "[mbpt][csv][basis][valgrind_skip]") {
  using namespace sequant;
  using namespace sequant::tests::pao;
  const Context ctx = gate_context();
  const IndexBasis pao = ctx.index_basis_registry()->retrieve(L"μ̃");
  auto scoped = set_scoped_default_context(ctx);

  SECTION("round trip and metadata") {
    const auto eqs = derive();
    REQUIRE(eqs.size() == 3);
    for (auto const& e : eqs) CHECK(deserialize<ExprPtr>(serialize(e)) == e);
    CHECK(Index(L"μ̃_3") == Index(pao, 3));
    CHECK(Index(pao, 3).full_label() == L"μ̃_3");
    CHECK(Index(pao, 3).basis().base_key() == L"μ̃");
    CHECK(Index(pao, 3).basis().extent() == n_pao);
    CHECK(Index(pao, 3).basis().metric() == IndexBasisMetric::General);
    // the bare (space, instance) pair that operator grants and the integral
    // projection build resolves to the basis as registered
    const IndexBasis bare{ctx.index_basis_registry()->retrieve_space(L"a"),
                          pao.basis_instance()};
    REQUIRE(bare.extent() != n_pao);
    CHECK(default_registry_resolved(bare) == pao);
    CHECK(default_registry_resolved(bare).extent() == n_pao);
    CHECK(ctx.index_basis_registry()->retrieve_space(L"a").dimension() !=
          n_pao);
    CHECK_THROWS_AS(ctx.index_basis_registry()->retrieve_space(L"μ̃"),
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
  SECTION("bra-ket symmetry of every tensor in the derived equations") {
    // integrals, the overlap and the CSV coefficients C are Hermitian:
    // bra<->ket symmetric over a real field, conjugate over a complex one;
    // amplitudes and the spintrace symmetrizer Ŝ are not
    for (const Field field : {Field::Real, Field::Complex}) {
      CAPTURE(field == Field::Real ? "real" : "complex");
      const auto hermitian = field == Field::Real ? BraKetSymmetry::Symm
                                                  : BraKetSymmetry::Conjugate;
      std::size_t n_C = 0;
      for (auto const& eq : derive(field))
        eq->visit(
            [&](ExprPtr const& x) {
              if (!x->is<Tensor>()) return;
              auto const& tnsr = x->as<Tensor>();
              const std::wstring label(tnsr.label());
              CAPTURE(toUtf8(label));
              if (label == L"C") {
                ++n_C;
                // like the overlap
                CHECK(tnsr.column_symmetry() == ColumnSymmetry::Symm);
              }
              if (label == L"C" || label == L"g" || label == L"f" ||
                  label == reserved::overlap_label())
                CHECK(tnsr.braket_symmetry() == hermitian);
              else if (label == L"t" || label == reserved::symm_label())
                CHECK(tnsr.braket_symmetry() == BraKetSymmetry::Nonsymm);
              else
                FAIL_CHECK("unclassified tensor label");
            },
            /* atoms_only = */ true);
      CHECK(n_C > 0);
    }
  }
  SECTION("half-transforming either leg of a DF factor is one intermediate") {
    // Σ C{a;μ̃} g{μ̃;i;Κ} and Σ C{μ̃;a} g{i;μ̃;Κ} are one object when g and C
    // are bra-ket symmetric (real field); canonicalization must fold them
    auto scoped_real = set_scoped_default_context(gate_context(Field::Real));
    const IndexBasis real_pao =
        get_default_context().index_basis_registry()->retrieve(L"μ̃");
    auto half_transform = [&real_pao](std::wstring_view g) {
      ExprPtr e = mbpt::csv_transform(deserialize<ExprPtr>(g), real_pao);
      REQUIRE(e->is<Product>());
      return canonicalize(std::move(e));
    };
    const ExprPtr bra_leg = half_transform(L"g{a_1<i_1>;i_2;Κ_1}:N-S-S");
    const ExprPtr ket_leg = half_transform(L"g{i_2;a_1<i_1>;Κ_1}:N-S-S");
    CAPTURE(toUtf8(serialize(bra_leg)), toUtf8(serialize(ket_leg)));
    CHECK(bra_leg == ket_leg);
  }
  SECTION("slot layout: the PAO basis sorts after every basis of a") {
    const IndexSpace a = ctx.index_basis_registry()->retrieve_space(L"a");
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
