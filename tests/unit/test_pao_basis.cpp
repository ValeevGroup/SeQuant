// PAO as a named basis instance of the particle space is canonically
// equivalent to the μ̃ space.
#include <catch2/catch_test_macros.hpp>
#include "catch2_sequant.hpp"
#include "csv_test_utils.hpp"

#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/domain/mbpt/rules/csv.hpp>
#include <SeQuant/domain/mbpt/rules/df.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace sequant::tests::pao {
using csv::PaoEncoding;
constexpr std::size_t n_pao = 300;

/// the CSV-CCSD registry in the given encoding with the PAO extent set; sized
/// before a Context adopts it, since a registry inside a Context is immutable
std::shared_ptr<IndexSpaceRegistry> sized_registry(PaoEncoding enc) {
  auto isr = csv::csv_cc_registry(enc);
  // either kind of entry: the μ̃ space or the named {a, INT32_MAX}
  isr->approximate_size(L"μ̃", n_pao);
  return isr;
}

/// csv_cc_context over the sized registry, with tmp ordinals allowed in parsed
/// text
Context gate_context(PaoEncoding enc) {
  auto ctx = csv::csv_cc_context(sized_registry(enc));
  ctx.set_first_dummy_index_ordinal(1000000);
  return ctx;
}

/// [E, R1, R2] in mpqc's order: CC{2}.t(2,0) with t granted 0, V1 spintrace,
/// DF, csv_transform(μ̃), flatten
std::vector<ExprPtr> derive(PaoEncoding enc) {
  auto core =
      set_scoped_default_context(csv::csv_cc_context(sized_registry(enc)));
  auto mbpt_ctx = mbpt::set_scoped_default_mbpt_context(
      mbpt::Context({.csv = mbpt::CSV::Yes,
                     .op_registry_ptr = csv::granted_registry({{L"t", 0}})}));
  const auto& isr = get_default_context().index_space_registry();
  Index::reset_tmp_index();
  std::vector<ExprPtr> out;
  for (auto const& eq : mbpt::CC{2}.t(2, 0)) {
    ExprPtr e = mbpt::closed_shell_CC_spintrace(
        eq, {.method = mbpt::BiorthogonalizationMethod::V1});
    e = mbpt::density_fit(e, isr->retrieve(L"Κ"), L"g", L"g");
    e = enc == PaoEncoding::Space
            ? mbpt::csv_transform(e, isr->retrieve(L"μ̃"))
            : mbpt::csv_transform(e, isr->retrieve_basis(L"μ̃"),
                                  /*orthonormal=*/false);
    flatten(e);
    out.push_back(e);
  }
  return out;
}

/// a copy of @p e with every named-PAO index replaced by the μ̃-space index of
/// the same ordinal (copies: no ordinal check), and the number of indices
/// replaced
std::pair<ExprPtr, std::size_t> to_space_encoding(ExprPtr const& e,
                                                  IndexBasis const& pao,
                                                  IndexSpace const& mu_space) {
  container::map<Index, Index> m;
  for (auto const& idx : get_used_indices(e))
    if (idx.basis() == pao) {
      SEQUANT_ASSERT(!idx.has_proto_indices());
      m.emplace(
          idx,
          idx.replace_basis_instance(std::nullopt).replace_space(mu_space));
    }
  return {m.empty() ? e->clone() : transform_expr(e, m), m.size()};
}

/// evidence files next to the ledger when SEQUANT_IBR_WORK is set; UTF-8
/// through a byte stream
void evidence(std::string const& name, std::wstring const& text) {
  if (const char* work = std::getenv("SEQUANT_IBR_WORK"))
    std::ofstream(std::string(work) + "/" + name) << toUtf8(text);
}
}  // namespace sequant::tests::pao

TEST_CASE("pao-named-basis-equivalence", "[mbpt][csv][basis][valgrind_skip]") {
  using namespace sequant;
  using namespace sequant::tests::pao;
  // derived once for all sections; a section clones before it mutates
  static const auto old_eqs = derive(PaoEncoding::Space);
  static const auto new_eqs = derive(PaoEncoding::Basis);
  REQUIRE(old_eqs.size() == new_eqs.size());

  const Context ctx_old = gate_context(PaoEncoding::Space);
  const Context ctx_new = gate_context(PaoEncoding::Basis);
  const IndexBasis pao = ctx_new.index_space_registry()->retrieve_basis(L"μ̃");
  const IndexSpace mu_space = ctx_old.index_space_registry()->retrieve(L"μ̃");

  SECTION("G1 canonical equivalence") {
    auto scoped = set_scoped_default_context(ctx_old);
    for (std::size_t k = 0; k < old_eqs.size(); ++k) {
      INFO("equation " << k);
      auto [mapped, n_mapped] = to_space_encoding(new_eqs[k], pao, mu_space);
      CHECK(n_mapped > 0);
      // (i) before simplify: same tmp numbering, same terms
      CHECK(*old_eqs[k] == *mapped);
      ExprPtr old_c = old_eqs[k]->clone();
      simplify(old_c);
      simplify(mapped);
      // (ii) after simplify
      CHECK(*old_c == *mapped);
      // (iii) the named form canonicalized in its own context, then mapped
      ExprPtr named_c = new_eqs[k]->clone();
      {
        auto s = set_scoped_default_context(ctx_new);
        simplify(named_c);
      }
      auto [mapped_c, n_mapped_c] = to_space_encoding(named_c, pao, mu_space);
      CHECK(n_mapped_c > 0);
      // (iv) evidence, not a gate: the named form's canonical form, mapped,
      // against the μ̃-space one term for term
      {
        const bool same = *old_c == *mapped_c;
        std::wstring text = same ? L"identical\n" : L"differ\n";
        if (!same && old_c->is<Sum>() && mapped_c->is<Sum>()) {
          auto const& x = old_c->as<Sum>().summands();
          auto const& y = mapped_c->as<Sum>().summands();
          for (std::size_t t = 0; t < std::min(x.size(), y.size()); ++t)
            if (*x[t] != *y[t]) {
              text += L"first differing term " + std::to_wstring(t) +
                      L"\nspace: " + serialize(x[t]) + L"\nbasis: " +
                      serialize(y[t]) + L"\n";
              break;
            }
        }
        if (!same)
          text += L"space:\n" + serialize(old_c) + L"\nbasis:\n" +
                  serialize(mapped_c) + L"\n";
        evidence("g1-iv-eq" + std::to_string(k) + ".txt", text);
        INFO("(iv) canonical forms " << (same ? "identical" : "differ"));
        CHECK(true);
      }
      simplify(mapped_c);
      CHECK(*old_c == *mapped_c);
      CHECK(tests::csv::term_count(old_eqs[k]) ==
            tests::csv::term_count(new_eqs[k]));
    }
  }
  SECTION("G2 text (evidence, not a gate)") {
    for (std::size_t k = 0; k < old_eqs.size(); ++k) {
      std::wstring new_text, old_text;
      {
        auto s = set_scoped_default_context(ctx_new);
        new_text = serialize(new_eqs[k]);
      }
      {
        auto s = set_scoped_default_context(ctx_old);
        old_text = serialize(old_eqs[k]);
      }
      evidence("g2-eq" + std::to_string(k) + "-space.txt", old_text);
      evidence("g2-eq" + std::to_string(k) + "-basis.txt", new_text);
      INFO("equation " << k << ": texts "
                       << (old_text == new_text ? "identical"
                                                : "differ (b1-b4)"));
      CHECK(true);
    }
  }
  SECTION("G3 round trip") {
    auto scoped = set_scoped_default_context(ctx_new);
    for (auto const& e : new_eqs)
      CHECK(deserialize<ExprPtr>(serialize(e)) == e);
    CHECK(Index(L"μ̃_3") == Index(pao, 3));
    CHECK(Index(pao, 3).full_label() == L"μ̃_3");
    CHECK(Index(pao, 3).basis_key() == L"μ̃");
    CHECK(Index(pao, 3).space().approximate_size() == n_pao);
    // a temporary minted from the bare (space, instance) pair, as operator
    // grants and the integral projection do, is the basis as registered
    const IndexBasis bare{ctx_new.index_space_registry()->retrieve(L"a"),
                          pao.basis_instance()};
    REQUIRE(bare.space().approximate_size() != n_pao);
    CHECK(Index::make_tmp_index(bare).space().approximate_size() == n_pao);
    CHECK(ctx_new.index_space_registry()->retrieve(L"a").approximate_size() !=
          n_pao);
    CHECK_THROWS_AS(ctx_new.index_space_registry()->retrieve(L"μ̃"),
                    IndexBasisRegistry::not_a_space);
    // the committed fixture: every former μ̃ index is the named basis
    std::ifstream in(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                     "/data/csv_ccsd_doubles_residual_df.txt");
    std::string line;
    REQUIRE(std::getline(in, line));
    const ExprPtr fixture = deserialize<ExprPtr>(toUtf16(line));
    std::size_t n = 0;
    for (auto const& idx : get_used_indices(fixture))
      if (idx.basis_key() == L"μ̃") {
        ++n;
        CHECK(idx.basis() == pao);
        CHECK_FALSE(idx.has_proto_indices());
      }
    CHECK(n > 0);
  }
  SECTION("G5 slot layout") {
    // per canonical slot: (is PAO, has protos)
    using Layout = std::vector<std::pair<bool, bool>>;
    // text -> one layout per encoding
    std::map<std::wstring, std::vector<Layout>> layouts;
    for (auto enc : {PaoEncoding::Space, PaoEncoding::Basis}) {
      auto scoped = set_scoped_default_context(gate_context(enc));
      for (std::wstring text :
           {L"C{μ̃_1;a_1<i_1,i_2;0>}:N-C-S", L"C{μ̃_1;a_1<i_1>}:N-C-S"}) {
        const ExprPtr expr = deserialize<ExprPtr>(text);
        const Tensor t = expr->as<Tensor>();
        const EvalExpr node(t);
        Layout layout;
        for (auto const& ix : node.canon_indices())
          layout.emplace_back(ix.basis_key() == L"μ̃", ix.has_proto_indices());
        layouts[text].push_back(layout);
      }
    }
    for (auto const& [text, per_encoding] : layouts) {
      INFO(toUtf8(text));
      REQUIRE(per_encoding.size() == 2);
      CHECK(per_encoding[0] == per_encoding[1]);
    }
    const IndexSpace a = ctx_new.index_space_registry()->retrieve(L"a");
    CHECK((IndexBasis{a} < pao && IndexBasis{a, 0} < pao &&
           IndexBasis{a, 2} < pao && IndexBasis{a, 10} < pao));
  }
  SECTION("G6 optimizer node counts (evidence, not a gate)") {
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    std::wstring report;
    for (std::size_t k = 0; k < old_eqs.size(); ++k) {
      auto count = [](Context const& ctx, ExprPtr const& e) {
        auto scoped = set_scoped_default_context(ctx);
        const Sum::summands_type terms =
            e->is<Sum>() ? e->as<Sum>().summands() : Sum::summands_type{e};
        // extents: approximate_size; every CSV composite a domain of 12
        const OptimizeOptions opts{
            .inner_pow = [](Index const&, std::size_t) { return 12.0; }};
        std::size_t nodes = 0;
        for (auto const& term : terms)
          binarize(optimize(term, opts)).visit([&nodes](auto const&) {
            ++nodes;
          });
        return nodes;
      };
      const auto n_old = count(ctx_old, old_eqs[k]->clone()),
                 n_new = count(ctx_new, new_eqs[k]->clone());
      report += L"eq " + std::to_wstring(k) + L": space " +
                std::to_wstring(n_old) + L" basis " + std::to_wstring(n_new) +
                L"\n";
      INFO("equation " << k << ": binary nodes space=" << n_old
                       << " basis=" << n_new);
      CHECK(true);
    }
    evidence("g6-node-counts.txt", report);
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  }
}
