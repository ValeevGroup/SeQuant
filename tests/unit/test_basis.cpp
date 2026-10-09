#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/basis.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/options.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/rules/csv.hpp>

#include <atomic>
#include <compare>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

#include <range/v3/view/single.hpp>

using namespace sequant;

namespace {

auto scoped_min_sr_context(CanonicalizeOptions canonicalization_options =
                               CanonicalizeOptions::default_options()) {
  return set_scoped_default_context(
      Context({.index_space_registry_shared_ptr = mbpt::make_min_sr_spaces(),
               .vacuum = Vacuum::SingleProduct,
               .spbasis = SPBasis::Spinor,
               .canonicalization_options = canonicalization_options}));
}

// μ̃ = {a, INT32_MAX} registered by name: the PAO basis of ibr-spec.md §3.3
auto scoped_pao_context() {
  auto isr = mbpt::make_min_sr_spaces();
  const IndexSpace uocc = isr->retrieve(L"a");
  constexpr IndexBasis::instance_type P =
      std::numeric_limits<IndexBasis::instance_type>::max();
  isr->add(L"μ̃", IndexBasis{uocc, P},
           120ul);  // before the Context adopts the registry (#665)
  return set_scoped_default_context(
      Context({.index_space_registry_shared_ptr = std::move(isr),
               .vacuum = Vacuum::SingleProduct}));
}

}  // namespace

// x = a_1<i_1,i_2;1> is the R2 residual leg of CSV-CCSD with t granted 1,
// y = a_7<;1> a CSV::No-granted leg, i_1<;5> an occupied leg under a grant
TEST_CASE("index-basis", "[elements][index][basis]") {
  auto ctx = scoped_min_sr_context();
  const auto& isr = get_default_context().index_space_registry();
  const IndexSpace occ = isr->retrieve(L"i");
  const IndexSpace uocc = isr->retrieve(L"a");
  const Index i1(occ, 1), i2(occ, 2), i3(occ, 3), i4(occ, 4);
  const Index a1(uocc, 1, {i1, i2}), a2(uocc, 2, {i1, i2}),
      a1_34(uocc, 1, {i3, i4});
  const Index x = a1.replace_basis_instance(1);
  const Index y = Index(uocc, 7).replace_basis_instance(1);
  const Index i1g = i1.replace_basis_instance(5),
              i2g = i2.replace_basis_instance(5);
  const Index xg = Index(uocc, 1, {i1g, i2g}).replace_basis_instance(1);
  auto I = [](Index const& idx, IndexBasis::instance_type n) {
    return idx.replace_basis_instance(n);
  };

  // null is the default for every construction from a space or a label
  for (Index const& idx :
       {Index{}, Index(uocc), Index(uocc, 1), a1, Index(L"a_1"),
        Index::make_tmp_index(uocc),
        Index::make_tmp_index(uocc, container::vector<Index>{i1}),
        IndexFactory{}.make(uocc)})
    CHECK_FALSE(idx.basis().has_basis_instance());
  CHECK(hash_value(IndexBasis{uocc}) == hash_value(uocc));

  // identity; the hash differs iff the indices differ
  CHECK(x == a1.replace_basis_instance(1));
  CHECK(x.replace_basis_instance(std::nullopt) == a1);
  for (IndexBasis::optional_instance other :
       std::initializer_list<IndexBasis::optional_instance>{std::nullopt, 0,
                                                            -1}) {
    const Index z = x.replace_basis_instance(other);
    CHECK(z.basis().basis_instance() == other);
    CHECK(z != x);
    CHECK(hash_value(z) != hash_value(x));
  }
  CHECK(hash_value(x) == hash_value(a1.replace_basis_instance(1)));
  CHECK(y != Index(uocc, 7));

  // ordering: the instance is the last key, null first
  using FLC = Index::FullLabelCompare;
  using TC = Index::TypeCompare;
  using TE = Index::TypeEquality;
  using LC = Index::LabelCompare;
  CHECK(I(a1, 9) < a2);
  CHECK(I(a1, 9) < a1_34);
  CHECK((a1 < I(a1, 0) && I(a1, 0) < I(a1, 1)));
  CHECK(I(a1, -1) < I(a1, 0));
  CHECK((FLC{}(I(a1, 9), a2) && !FLC{}(a2, I(a1, 9))));
  CHECK((FLC{}(I(a1, 9), a1_34) && !FLC{}(a1_34, I(a1, 9))));
  CHECK((FLC{}(a1, I(a1, 0)) && FLC{}(I(a1, 0), I(a1, 1)) &&
         !FLC{}(I(a1, 1), I(a1, 0))));
  CHECK((TC{}(I(a1, 9), a1_34) && !TC{}(a1_34, I(a1, 9))));
  CHECK((TC{}(a1, I(a1, 1)) && !TC{}(I(a1, 1), a1)));
  CHECK((!TC{}(I(a2, 1), I(a1, 1)) && !TC{}(I(a1, 1), I(a2, 1))));
  CHECK((a1 < a2 && a1 < a1_34 && FLC{}(a1, a2) && TC{}(a1, a1_34)));
  CHECK((!TE{}(a1, I(a1, 1)) && TE{}(I(a1, 1), I(a2, 1))));
  CHECK((!LC{}(a1, I(a1, 1)) && !LC{}(I(a1, 1), a1)));
  CHECK((IndexBasis{uocc} < IndexBasis{uocc, 0} &&
         IndexBasis{uocc, 0} < IndexBasis{uocc, 1}));

  // colour: sees the instance, of the index and of its proto indices, but not
  // which occupieds the proto indices are
  CHECK(x.color() != a1.color());
  CHECK(x.color() == I(a1_34, 1).color());
  CHECK(xg.color() != x.color());
  CHECK(y.color() != Index(uocc, 7).color());

  // copy paths keep the instance
  {
    Index copied, move_assigned, victim(x), victim2(x), transformed(x);
    copied = x;
    const Index moved(std::move(victim));
    move_assigned = std::move(victim2);
    for (Index const& idx : {Index(x), copied, moved, move_assigned})
      CHECK(idx == x);
    CHECK(Index(x, {i2}) == I(Index(uocc, 1, {i2}), 1));
    CHECK(Index(x, container::vector<Index>{i2}) == I(Index(uocc, 1, {i2}), 1));
    const IndexSpace uocc_alpha = isr->retrieve(L"a↑");
    const Index x_alpha = x.replace_space(uocc_alpha);
    CHECK(x_alpha.space() == uocc_alpha);
    CHECK(x_alpha.basis().basis_instance() == 1);
    CHECK(x_alpha.proto_indices() == x.proto_indices());
    CHECK(x_alpha.symmetric_proto_indices() == x.symmetric_proto_indices());
    CHECK(x.replace_qns(x.space().qns()) == x);

    const Index remade = IndexFactory{}.make(x);
    CHECK(remade.basis() == x.basis());
    CHECK(remade.proto_indices() == x.proto_indices());
    CHECK(remade.symmetric_proto_indices() == x.symmetric_proto_indices());
    CHECK(remade.ordinal() != x.ordinal());

    container::map<Index, Index> map;
    map.emplace(i2, i3);
    REQUIRE(transformed.transform(map));
    CHECK(transformed == I(Index(uocc, 1, {i1, i3}), 1));

    CHECK(Index::make_tmp_index(IndexBasis{uocc, 1}).basis() == x.basis());
    const Index tmp = Index::make_tmp_index(IndexBasis{uocc, 1},
                                            container::vector<Index>{i1, i2});
    CHECK(tmp.basis() == x.basis());
    CHECK(tmp.proto_indices() == x.proto_indices());
  }

  // drop_proto_indices() keeps the instance; the CSV shortcut of a standing
  // CCSD R1 overlap still expands through the null complete set
  CHECK(x.drop_proto_indices() == I(Index(uocc, 1), 1));
  CHECK(x.drop_proto_indices().full_label() == L"a_1<;1>");
  {
    const Index s_bra = I(Index(uocc, 1, {i1}), 1);
    const Index s_ket = I(Index(uocc, 3, {i1, i2}), 1);
    const auto ct = mbpt::csv_transform(make_overlap(s_bra, s_ket), uocc);
    REQUIRE(ct->is<Product>());
    REQUIRE(ct->as<Product>().factors().size() == 2);
    const auto& c1 = ct->as<Product>().factor(0)->as<Tensor>();
    const auto& c2 = ct->as<Product>().factor(1)->as<Tensor>();
    CHECK(c1.label() == L"C");
    CHECK(c2.label() == L"C");
    CHECK(c1.bra().at(0) == s_bra);
    CHECK(c2.ket().at(0) == s_ket);
    const Index dummy = c1.ket().at(0);
    CHECK(dummy == c2.bra().at(0));
    CHECK(dummy.basis() == IndexBasis{uocc});
    CHECK(!dummy.has_proto_indices());
  }
  // ... with a fresh dummy per overlap, also when two overlaps share a leg
  {
    const Index x1(uocc, 1, {i1, i2}), x2(uocc, 2, {i1}), x3(uocc, 3, {i2});
    const auto ct =
        mbpt::csv_transform(make_overlap(x2, x1) * make_overlap(x1, x3), uocc);
    container::map<Index, int> counts;
    ct->visit(
        [&counts](const ExprPtr& e) {
          if (e->is<Tensor>()) {
            for (const auto& idx : e->as<Tensor>().bra()) ++counts[idx];
            for (const auto& idx : e->as<Tensor>().ket()) ++counts[idx];
          }
        },
        /*atoms_only=*/true);
    for (const auto& [idx, n] : counts) {
      CAPTURE(toUtf8(idx.full_label()));
      CHECK(n <= 2);
    }
  }
  // ... and a single C only if the leg without proto indices is in the
  // target basis; else the overlap to that leg stays
  {
    const Index s_bra(uocc, 1, {i1, i2});
    CHECK(mbpt::csv_transform(make_overlap(s_bra, Index(uocc, 2)), uocc)
              ->size() == 1);
    const auto ct =
        mbpt::csv_transform(make_overlap(s_bra, I(Index(uocc, 2), 1)), uocc);
    REQUIRE(ct->is<Product>());
    CHECK(ct->size() == 2);
    CHECK(std::ranges::count_if(ct->as<Product>().factors(), [](auto& f) {
            return f->template as<Tensor>().label() ==
                   reserved::overlap_label();
          }) == 1);
  }
  // ... also between temporary indices, e.g. fresh from Wick
  {
    const Index s_bra = Index::make_tmp_index(IndexBasis{uocc, 1},
                                              container::vector<Index>{i1});
    const Index s_ket = Index::make_tmp_index(IndexBasis{uocc, 1},
                                              container::vector<Index>{i1, i2});
    ExprPtr ct;
    REQUIRE_NOTHROW(ct = mbpt::csv_transform(make_overlap(s_bra, s_ket), uocc));
    REQUIRE(ct->is<Product>());
    REQUIRE(ct->as<Product>().factors().size() == 2);
    const Index dummy = ct->as<Product>().factor(0)->as<Tensor>().ket().at(0);
    CHECK(dummy == ct->as<Product>().factor(1)->as<Tensor>().bra().at(0));
    CHECK(dummy.basis() == IndexBasis{uocc});
    CHECK(!dummy.has_proto_indices());
    CHECK(dummy != s_bra.drop_proto_indices().replace_basis_instance({}));
  }

  // labels: label() never shows the instance
  CHECK(x.label() == L"a_1");
  CHECK(y.label() == L"a_7");
  CHECK(a1.full_label() == L"a_1<i_1, i_2>");
  CHECK(x.full_label() == L"a_1<i_1, i_2;1>");
  CHECK(y.full_label() == L"a_7<;1>");
  CHECK(I(a1, 0).full_label() == L"a_1<i_1, i_2;0>");
  CHECK(xg.full_label() == L"a_1<i_1<;5>, i_2<;5>;1>");
  CHECK(a1.to_latex() == L"{a_1^{{i_1}{i_2}}}");
  CHECK(x.to_latex() == L"{a_1^{{i_1}{i_2};1}}");
  CHECK(I(Index(uocc, 1), 1).to_latex() == L"{a_1^{;1}}");
  {
    Index memoized(a1);
    (void)memoized.full_label();
    CHECK(memoized.replace_basis_instance(1).full_label() == x.full_label());
  }

  // IndexFactory counters are per IndexSpace; the instance is kept
  {
    IndexFactory f;
    CHECK(f.make(uocc).full_label() == L"a_100");
    CHECK(f.make(y).full_label() == L"a_101<;1>");
    CHECK(f.make(uocc).full_label() == L"a_102");
  }
}

TEST_CASE("index-basis-canonicalization", "[algorithms][canonicalize][basis]") {
  SECTION("CCSD R1 term with a standing overlap") {
    auto ctx = scoped_min_sr_context(
        CanonicalizeOptions::default_options().copy_and_set(
            CanonicalizationMethod::Complete));
    const auto& isr = get_default_context().index_space_registry();
    const IndexSpace occ = isr->retrieve(L"i");
    const IndexSpace uocc = isr->retrieve(L"a");
    const Index i1(occ, 1), i2(occ, 2);
    const Index a1 = Index(uocc, 1, {i1}).replace_basis_instance(1);

    // Â{i_1;a_1<i_1;1>} f{i_2;a<i_1,i_2;F>} s{a_1<i_1;1>;b<i_1,i_2;F>}
    // t{a,b;i_1,i_2}: the dummies a, b in amplitude family F
    auto term = [&](std::size_t a, std::size_t b,
                    IndexBasis::instance_type family) {
      const Index ad = Index(uocc, a, {i1, i2}).replace_basis_instance(family);
      const Index bd = Index(uocc, b, {i1, i2}).replace_basis_instance(family);
      return ex<Tensor>(L"Â", bra{i1}, ket{a1}, Symmetry::Antisymm) *
             ex<Tensor>(L"f", bra{i2}, ket{ad}, Symmetry::Nonsymm) *
             make_overlap(a1, bd) *
             ex<Tensor>(L"t", bra{ad, bd}, ket{i1, i2}, Symmetry::Antisymm);
    };

    const auto renumbered = simplify(term(2, 3, 1) - term(7, 8, 1));
    CHECK(renumbered == ex<Constant>(0));
    const auto other_family = simplify(term(2, 3, 1) - term(2, 3, 2));
    REQUIRE(other_family->is<Sum>());
    CHECK(other_family->as<Sum>().size() == 2);
  }

  SECTION("no dummy takes a named index's label in another basis") {
    auto ctx = scoped_min_sr_context();
    const auto& isr = get_default_context().index_space_registry();
    const IndexSpace occ = isr->retrieve(L"i");
    const IndexSpace uocc = isr->retrieve(L"a");
    const Index i1(occ, 1), i2(occ, 2);
    // the tʼ R1 term of CC{2}.tʼ(1,1) (t granted 1, t¹ granted 10) with its
    // symmetrizer stripped: i_1 and a_1<i_1;10> are named
    const Index a1 = Index(uocc, 1, {i1}).replace_basis_instance(10);
    const Index a9 = Index(uocc, 9, {i1}).replace_basis_instance(1);
    const Index a3 = Index(uocc, 3, {i2}).replace_basis_instance(10);
    auto term = ex<Tensor>(L"g", bra{i2, a1}, ket{a9, a3}, Symmetry::Antisymm) *
                ex<Tensor>(L"t", bra{a9}, ket{i1}, Symmetry::Antisymm) *
                ex<Tensor>(L"t¹", bra{a3}, ket{i2}, Symmetry::Antisymm);
    canonicalize(term);

    std::set<Index> indices;
    term->visit(
        [&indices](ExprPtr const& x) {
          if (!x->is<Tensor>()) return;
          for (auto&& idx : x->as<Tensor>().const_braketaux_indices())
            indices.insert(idx);
          if (x->as<Tensor>().label() == L"t")
            CHECK(x->as<Tensor>().bra().at(0).label() != L"a_1");
        },
        /* atoms_only = */ true);
    REQUIRE(indices.contains(a1));
    std::set<std::wstring> labels;
    std::set<std::string> annotations;
    for (Index const& idx : indices) {
      labels.emplace(idx.label());
      annotations.insert(csv_labels(ranges::views::single(idx)));
    }
    CHECK(labels.size() == indices.size());
    CHECK(annotations.size() == indices.size());
  }
}

TEST_CASE("index-basis-annotation-and-hash", "[EvalExpr][basis]") {
  auto ctx = scoped_min_sr_context();
  const auto& isr = get_default_context().index_space_registry();
  const IndexSpace occ = isr->retrieve(L"i");
  const IndexSpace uocc = isr->retrieve(L"a");
  const Index i1(occ, 1), i2(occ, 2);
  const Index a1(uocc, 1, {i1, i2}), a2(uocc, 2, {i1, i2});
  const Index y = Index(uocc, 7).replace_basis_instance(1);
  const Index xg =
      Index(uocc, 1,
            {i1.replace_basis_instance(5), i2.replace_basis_instance(5)})
          .replace_basis_instance(1);

  CHECK(csv_labels(ranges::views::single(a1)) == "a_1i_1i_2");
  CHECK(csv_labels(ranges::views::single(a1.replace_basis_instance(7))) ==
        "a_1i_1i_2#7");
  CHECK(csv_labels(ranges::views::single(
            Index(uocc, 1).replace_basis_instance(3))) == "a_1#3");
  CHECK(csv_labels(ranges::views::single(xg)) == "a_1i_1#5i_2#5#1");
  CHECK(instance_qualified_label(Index(uocc, 1)) == L"a_1");
  CHECK(instance_qualified_label(Index(uocc, 1).replace_basis_instance(3)) ==
        L"a_1#3");
  CHECK(instance_qualified_label(xg) == L"a_1#1");

  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  const auto g0 = ex<Tensor>(L"g", bra{a1}, ket{a2});
  const auto g7 = ex<Tensor>(L"g", bra{a1.replace_basis_instance(7)},
                             ket{a2.replace_basis_instance(7)});
  CHECK(binarize(g0)->hash_value() != binarize(g7)->hash_value());
  CHECK(binarize(g7)->indices_annot().find("#7") != std::string::npos);

  const auto bare0 = ex<Tensor>(L"g", bra{Index(uocc, 7)}, ket{i1});
  const auto bare1 = ex<Tensor>(L"g", bra{y}, ket{i1});
  REQUIRE_NOTHROW(binarize(bare1));
  CHECK(binarize(bare0)->hash_value() != binarize(bare1)->hash_value());
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
}

TEST_CASE("index-basis-serialization", "[serialization][basis]") {
  auto ctx = scoped_min_sr_context();

  SECTION("round trips") {
    const std::vector<std::wstring> expressions = {
        L"t{i_1;a_1<i_1;1>}:N-C-S",
        L"t{i_1;a_1<i_1,i_2;7>}:N-C-S",
        L"t{i_1;a_1<;-2>}:N-C-S",
        L"t{i_1;a_1<i_1<;3>,i_2;1>}:N-C-S",
        L"t{i_1<;3>;a_1<i_1<;3>>}:N-C-S",
        L"t{i_1;a_1<;0>}:N-C-S",
        L"g{i_1,i_2;a_1<;1>,a_2<;2>}:N-C-S",
    };

    for (const std::wstring& current : expressions) {
      ExprPtr expression = deserialize<ExprPtr>(current);
      REQUIRE(serialize(expression, {.annot_symm = true}) == current);
      CHECK(deserialize<ExprPtr>(serialize(expression)) == expression);
    }
  }

  SECTION("';0' is legal and distinct from no instance") {
    auto zero = deserialize<ExprPtr>(L"t{i_1;a_1<;0>}:N-C-S");
    auto null = deserialize<ExprPtr>(L"t{i_1;a_1}:N-C-S");
    const Index& a_zero = zero->as<Tensor>().ket().at(0);
    const Index& a_null = null->as<Tensor>().ket().at(0);
    REQUIRE(a_zero.basis().has_basis_instance());
    CHECK(*a_zero.basis().basis_instance() == 0);
    CHECK(!a_null.basis().has_basis_instance());
    CHECK(a_zero != a_null);
  }

  SECTION("whitespace inside a domain is skipped") {
    auto with_space = deserialize<ExprPtr>(L"t{a_1<i_1 ; 10>;i_1}");
    auto without_space = deserialize<ExprPtr>(L"t{a_1<i_1;10>;i_1}");
    REQUIRE(with_space == without_space);
    CHECK(serialize(with_space) == serialize(without_space));
  }

  SECTION("spin-labelled index with a basis instance") {
    auto spin_ctx = get_default_context();
    spin_ctx.set(mbpt::make_sr_spaces());
    auto resetter = set_scoped_default_context(spin_ctx);

    auto expr = deserialize<ExprPtr>(L"t{a↓_1;i↑_1<;2>}:N-C-S");
    REQUIRE(serialize(expr, {.annot_symm = true}) == L"t{a↓_1;i↑_1<;2>}:N-C-S");
    CHECK(expr->as<Tensor>().ket().at(0).basis().basis_instance() == 2);
  }

  SECTION("malformed or out-of-range domains throw SerializationError") {
    CHECK_THROWS_AS(deserialize<ExprPtr>(L"t{i1;a1<i1;x>}"),
                    io::serialization::SerializationError);
    CHECK_THROWS_AS(deserialize<ExprPtr>(L"t{i1;a1<;2147483648>}"),
                    io::serialization::SerializationError);

    try {
      deserialize<ExprPtr>(L"t{i1;a1<>}");
      FAIL("expected a SerializationError");
    } catch (const io::serialization::SerializationError& e) {
      CHECK(std::string(e.what()).find("empty index domain") !=
            std::string::npos);
    }
  }

  SECTION("a named basis instance round-trips by its label") {
    auto isr = mbpt::make_min_sr_spaces();
    mbpt::add_df_spaces(isr);
    const IndexSpace a = isr->retrieve(L"a"), i = isr->retrieve(L"i");
    constexpr IndexBasis::instance_type P =
        std::numeric_limits<IndexBasis::instance_type>::max();
    isr->add(L"μ̃", IndexBasis{a, P}, 120ul);
    // a named occupied basis (localized occupieds), so a proto can be named
    // too; written with the combining tilde U+0303 the v1 index_name alphabet
    // admits
    isr->add(L"ĩ", IndexBasis{i, 1});
    // a copy shares the registry until set() moves the populated one in
    auto named = get_default_context_snapshot();
    named.set(std::move(isr));
    // fixtures and tmp indices carry ordinals >= 100
    named.set_first_dummy_index_ordinal(1000000);
    auto named_ctx = set_scoped_default_context(std::move(named));
    for (const std::wstring& current :
         {std::wstring(L"C{μ̃_1;a_1<i_1,i_2;0>}:N-C-S"),
          std::wstring(L"g{μ̃_1;i_1;Κ_1}:N-C-S"),
          std::wstring(L"g{μ̃_1;μ̃_2;Κ_1}:N-C-S"),
          std::wstring(L"t{a_1<ĩ_1>;ĩ_1}:N-C-S"),
          std::wstring(L"C{μ̃_1152;a_1<i_1,i_2;0>}:N-C-S")}) {
      ExprPtr e = deserialize<ExprPtr>(current);
      REQUIRE(serialize(e, {.annot_symm = true}) == current);
      CHECK(deserialize<ExprPtr>(serialize(e)) == e);
    }
    const ExprPtr c_named = deserialize<ExprPtr>(L"C{μ̃_1;a_1<i_1,i_2;0>}");
    const Index mu = c_named->as<Tensor>().bra().at(0);
    CHECK(mu.basis() == IndexBasis{a, P});
    CHECK(mu.space().approximate_size() == 120);
    // the generic spelling of the same basis parses to the named entry's
    // metadata and prints by name
    const ExprPtr c_generic =
        deserialize<ExprPtr>(L"C{a_1<;2147483647>;a_2<i_1,i_2;0>}");
    const Index generic = c_generic->as<Tensor>().bra().at(0);
    CHECK(generic == mu);
    CHECK(generic.space().approximate_size() == 120);
    CHECK(generic.full_label() == L"μ̃_1");
    // one spelling per basis: a name with an explicit instance is an error,
    // also as a proto, reported as such and not wrapped as an invalid index
    const auto names_an_instance = [](std::wstring_view input) {
      try {
        deserialize<ExprPtr>(input);
      } catch (const io::serialization::SerializationError& e) {
        const std::string what = e.what();
        return what.find("names a basis instance") != std::string::npos &&
               what.find("Invalid index") == std::string::npos;
      }
      return false;
    };
    CHECK(names_an_instance(L"C{μ̃_1<;3>;a_1<i_1,i_2;0>}"));
    CHECK(names_an_instance(L"C{μ̃_1<;2147483647>;a_1<i_1,i_2;0>}"));
    CHECK(names_an_instance(L"t{a_1<ĩ_1<;3>>;ĩ_1}"));
    // parsed indices carry their names, also under a registry without them
    {
      const ExprPtr named_c = deserialize<ExprPtr>(L"C{μ̃_1;a_1<i_1,i_2;0>}");
      const ExprPtr generic_c =
          deserialize<ExprPtr>(L"C{a_1<;2147483647>;a_2<i_1,i_2;0>}");
      const ExprPtr t = deserialize<ExprPtr>(L"t{a_1<ĩ_1>;ĩ_1}");
      const Index named_mu = named_c->as<Tensor>().bra().at(0),
                  generic_mu = generic_c->as<Tensor>().bra().at(0),
                  with_named_proto = t->as<Tensor>().bra().at(0),
                  named_loc = t->as<Tensor>().ket().at(0);
      auto plain = scoped_min_sr_context();
      CHECK(named_mu.full_label() == L"μ̃_1");
      CHECK(generic_mu.full_label() == L"μ̃_1");
      CHECK(with_named_proto.full_label() == L"a_1<ĩ_1>");
      CHECK(named_loc.full_label() == L"ĩ_1");
    }
    // an unnamed instance still prints as ;N
    CHECK(serialize(deserialize<ExprPtr>(L"t{a_1<;1>;i_1}:N-C-S"),
                    {.annot_symm = true}) == L"t{a_1<;1>;i_1}:N-C-S");
  }
}

// moving an index with a memoized label into a spin space: the label follows
// the new space
TEST_CASE("index-move-into-space-resets-label", "[elements][index]") {
  auto ctx = scoped_min_sr_context();
  const auto& isr = get_default_context().index_space_registry();
  const IndexSpace occ = isr->retrieve(L"i"), uocc = isr->retrieve(L"a"),
                   uocc_alpha = isr->retrieve(L"a↑");
  Index a3(uocc, 3);
  (void)a3.label();
  (void)a3.full_label();
  const Index moved(std::move(a3), uocc_alpha);
  CHECK(moved.space() == uocc_alpha);
  CHECK(moved.label() == L"a↑_3");
  CHECK(moved.full_label() == L"a↑_3");

  // a PNO-style index memoizes full_label too
  Index pno(uocc, 4, {Index(occ, 1)});
  (void)pno.label();
  (void)pno.full_label();
  const Index moved_pno(std::move(pno), uocc_alpha);
  CHECK(moved_pno.label() == L"a↑_4");
  CHECK(moved_pno.full_label() == L"a↑_4<i_1>");
}

TEST_CASE("index-basis-named", "[elements][index][basis]") {
  auto ctx = scoped_pao_context();
  const auto& isr = get_default_context().index_space_registry();
  const IndexSpace uocc = isr->retrieve(L"a");
  const IndexBasis pao = isr->retrieve_basis(L"μ̃");
  const IndexBasis::instance_type P = *pao.basis_instance();

  // the registry's entry carries its name, which is not part of the identity
  CHECK(pao.name() == L"μ̃");
  CHECK_FALSE(IndexBasis(uocc, P).has_name());
  CHECK(IndexBasis(uocc, P) == pao);
  CHECK(hash_value(IndexBasis(uocc, P)) == hash_value(pao));
  CHECK(isr->resolve(IndexBasis(uocc, P)).name() == L"μ̃");

  // labels: the name, no instance suffix anywhere
  const Index m3(pao, 3);
  CHECK(m3.label() == L"μ̃_3");
  CHECK(m3.full_label() == L"μ̃_3");
  CHECK(m3.to_latex() == L"{\\tilde{\\mu}_3}");
  CHECK(m3.to_string() == "μ̃_3");
  CHECK(m3.basis().base_key() == L"μ̃");
  CHECK_FALSE(m3.basis().unnamed_instance());
  CHECK(csv_labels(ranges::views::single(m3)) == "μ̃_3");
  // an unnamed instance (a CSV::No-granted leg) prints as before
  const Index y = Index(uocc, 7).replace_basis_instance(1);
  CHECK(y.label() == L"a_7");
  CHECK(y.full_label() == L"a_7<;1>");
  CHECK(y.to_latex() == L"{a_7^{;1}}");
  CHECK(y.basis().base_key() == L"a");
  CHECK(y.basis().unnamed_instance() == 1);
  CHECK(csv_labels(ranges::views::single(y)) == "a_7#1");
  // a null instance is untouched, and its base key is the space key
  CHECK(Index(uocc, 2).basis().base_key() == L"a");
  CHECK(Index(L"a_2").basis().base_key() == L"a");
  CHECK_FALSE(Index(uocc, 2).basis().unnamed_instance());

  // identity is untouched; the name is not part of it, and a bare instance
  // number carries none
  CHECK(m3 == Index(uocc, 3).replace_basis_instance(P));
  CHECK(Index(uocc, 3).replace_basis_instance(P).full_label() ==
        L"a_3<;" + std::to_wstring(P) + L">");
  CHECK(m3 != Index(uocc, 3));
  CHECK(hash_value(m3) == hash_value(Index(uocc, 3).replace_basis_instance(P)));
  CHECK((Index(uocc, 3) < m3 && Index(uocc, 3).replace_basis_instance(0) < m3));

  // from a label: the named basis with the entry's metadata
  const Index parsed(L"μ̃_3");
  CHECK(parsed == m3);
  CHECK(parsed.basis() == pao);
  CHECK(parsed.space().approximate_size() == 120);
  CHECK(parsed.label() == L"μ̃_3");
  CHECK(m3.space().approximate_size() == 120);
  CHECK(Index(L"a_3").space().approximate_size() == uocc.approximate_size());
  CHECK_THROWS_AS(IndexSpace(L"μ̃"), IndexBasisRegistry::not_a_space);

  // copies of a named index print the name; a basis change drops it
  {
    Index memo(m3);
    (void)memo.full_label();
    Index copied(memo), assigned;
    assigned = memo;
    CHECK(copied.label() == L"μ̃_3");
    CHECK(assigned.full_label() == L"μ̃_3");
    CHECK(memo.replace_basis_instance(std::nullopt).full_label() == L"a_3");
    CHECK(memo.replace_basis_instance(1).full_label() == L"a_3<;1>");
  }

  // the name travels with the basis, whatever registry is current when the
  // index is printed (the nested scope is a thread-local overlay, #655)
  {
    auto plain = set_scoped_default_context(
        Context({.index_space_registry_shared_ptr = mbpt::make_min_sr_spaces(),
                 .vacuum = Vacuum::SingleProduct}));
    CHECK(Index(m3).label() == L"μ̃_3");
    CHECK(Index(pao, 3).full_label() == L"μ̃_3");
    CHECK(Index(pao, 3).basis().base_key() == L"μ̃");
  }
}

// an index minted or renamed in a named basis carries the name: a copy made
// after the context switched to a registry without the name still prints it
TEST_CASE("index-basis-named-at-minting", "[elements][index][basis]") {
  auto ctx = scoped_pao_context();
  const auto& isr = get_default_context().index_space_registry();
  const IndexSpace occ = isr->retrieve(L"i"), uocc = isr->retrieve(L"a");
  const IndexBasis pao = isr->retrieve_basis(L"μ̃");
  const Index i1(occ, 1), i2(occ, 2);

  // minted, never labelled here
  IndexFactory factory;
  const Index from_basis = factory.make(pao);
  const Index from_index = factory.make(Index(pao, 3));
  const Index tmp = Index::make_tmp_index(pao);
  // from the bare instance number: the name is the registry's at minting
  const Index from_number =
      Index::make_tmp_index(IndexBasis{uocc, *pao.basis_instance()});

  // renamed by the canonicalizer: the PAO Fock coupling of a CSV R2 term,
  // C{a<i_1,i_2;0>;μ̃} f{μ̃;μ̃} C{μ̃;c<i_1,i_2;0>} t{c,b;i_1,i_2}, with the μ̃
  // dummies numbered off the canonical order
  const Index a1 = Index(uocc, 1, {i1, i2}).replace_basis_instance(0);
  const Index c2 = Index(uocc, 2, {i1, i2}).replace_basis_instance(0);
  const Index b3 = Index(uocc, 3, {i1, i2}).replace_basis_instance(0);
  const Index m7(pao, 7), m8(pao, 8);
  ExprPtr term = ex<Tensor>(L"C", bra{a1}, ket{m7}, Symmetry::Nonsymm) *
                 ex<Tensor>(L"f", bra{m7}, ket{m8}, Symmetry::Nonsymm) *
                 ex<Tensor>(L"C", bra{m8}, ket{c2}, Symmetry::Nonsymm) *
                 ex<Tensor>(L"t", bra{c2, b3}, ket{i1, i2}, Symmetry::Antisymm);
  simplify(term);
  std::vector<Index> renamed;
  for (const Index& idx : get_used_indices(term))
    if (idx.basis() == pao) renamed.push_back(idx);
  REQUIRE(renamed.size() == 2);
  REQUIRE((renamed[0].ordinal() != 7 || renamed[1].ordinal() != 8));

  auto plain = scoped_min_sr_context();
  auto named = [](Index copy) { return copy.label().starts_with(L"μ̃_"); };
  CHECK(named(from_basis));
  CHECK(named(from_index));
  CHECK(named(tmp));
  for (const Index& idx : renamed) CHECK(named(idx));
  CHECK(named(from_number));
}

// several threads copy one index and read its label at once, each labelling its
// own copy; factory-minted instance-bearing indices carry their label memo, so
// this only reads a settled value (ibr-spec.md §6.1); the CSV composite
// a<i_1,i_2;0> is what OpMaker mints, μ̃ the PAO basis. The std::threads below
// see the process-wide context, not this test's scoped overlay (#655), whose
// registry lacks μ̃: a probe whose memo was not settled at minting prints a_N
// there
TEST_CASE("index-basis-copy-while-labelling", "[elements][index][basis]") {
  auto ctx = scoped_pao_context();
  const auto& isr = get_default_context().index_space_registry();
  const IndexSpace occ = isr->retrieve(L"i"), uocc = isr->retrieve(L"a");
  const IndexBasis pao = isr->retrieve_basis(L"μ̃");
  const Index i1(occ, 1), i2(occ, 2);
  IndexFactory factory;
  const std::vector<Index> minted = {
      Index::make_tmp_index(IndexBasis{uocc, 0},
                            container::vector<Index>{i1, i2}),
      Index::make_tmp_index(IndexBasis{uocc, 0}),
      factory.make(IndexBasis{uocc, 0}),
      factory.make(Index(uocc, 1, {i1, i2}).replace_basis_instance(0)),
      Index::make_tmp_index(pao),
      factory.make(pao),
      factory.make(Index(pao, 3))};
  // an index renamed by the canonicalizer (Index::transform re-populates the
  // memo it resets): the CSV R1 term of the "index-basis-canonicalization"
  // case, simplified
  const Index a1 = Index(uocc, 1, {i1}).replace_basis_instance(1);
  const Index ad = Index(uocc, 2, {i1, i2}).replace_basis_instance(1);
  const Index bd = Index(uocc, 3, {i1, i2}).replace_basis_instance(1);
  ExprPtr term = ex<Tensor>(L"Â", bra{i1}, ket{a1}, Symmetry::Antisymm) *
                 ex<Tensor>(L"f", bra{i2}, ket{ad}, Symmetry::Nonsymm) *
                 make_overlap(a1, bd) *
                 ex<Tensor>(L"t", bra{ad, bd}, ket{i1, i2}, Symmetry::Antisymm);
  simplify(term);
  std::vector<Index> canonicalized;
  for (const Index& idx : get_used_indices(term))
    if (idx.basis().has_basis_instance()) canonicalized.push_back(idx);
  REQUIRE_FALSE(canonicalized.empty());

  std::vector<Index> probes = minted;
  probes.insert(probes.end(), canonicalized.begin(), canonicalized.end());
  for (const Index& m : probes) {
    // expected labels from an index rebuilt without m's memo, so that m is
    // never labelled on this thread
    const Index rebuilt(m.drop_proto_indices(), m.proto_indices(),
                        m.symmetric_proto_indices());
    REQUIRE(rebuilt == m);
    const std::wstring expected_full(rebuilt.full_label()),
        expected(rebuilt.label());
    std::vector<std::thread> threads;
    std::atomic<int> mismatches{0};
    for (int t = 0; t < 8; ++t)
      threads.emplace_back([&m, &expected_full, &expected, &mismatches] {
        for (int k = 0; k < 1000; ++k) {
          Index copy(m);
          if (copy.full_label() != expected_full || copy.label() != expected ||
              m.label() != expected)
            ++mismatches;
        }
      });
    for (auto& th : threads) th.join();
    CHECK(mismatches == 0);
  }
}

namespace {

template <typename Basis, typename... Args>
concept csv_transform_callable = requires(ExprPtr e, Basis b, Args... args) {
  mbpt::csv_transform(e, b, args...);
};

template <typename Basis>
concept csv_transform_takes_label_literal =
    requires(ExprPtr e, Basis b) { mbpt::csv_transform(e, b, L"C"); };

// a label in the orthonormal slot would convert to true
static_assert(!csv_transform_takes_label_literal<IndexBasis>);
static_assert(csv_transform_takes_label_literal<IndexSpace>);
static_assert(!csv_transform_callable<IndexBasis, const wchar_t*>);
static_assert(!csv_transform_callable<IndexBasis, wchar_t*>);
static_assert(!csv_transform_callable<IndexBasis, const char*>);
static_assert(csv_transform_callable<IndexBasis, bool>);
static_assert(csv_transform_callable<IndexBasis, bool, const wchar_t*>);
static_assert(csv_transform_callable<IndexBasis, bool, wchar_t*>);
static_assert(csv_transform_callable<IndexBasis, bool, std::wstring,
                                     container::svector<std::wstring>>);
static_assert(csv_transform_callable<IndexSpace>);
static_assert(csv_transform_callable<IndexSpace, const wchar_t*>);
static_assert(csv_transform_callable<IndexSpace, wchar_t*>);
static_assert(csv_transform_callable<IndexSpace, const wchar_t*,
                                     container::svector<std::wstring>>);

}  // namespace

TEST_CASE("csv-transform-named-basis", "[mbpt][csv][basis]") {
  auto isr = mbpt::make_min_sr_spaces();
  const IndexSpace occ = isr->retrieve(L"i"), uocc = isr->retrieve(L"a");
  constexpr IndexBasis::instance_type P =
      std::numeric_limits<IndexBasis::instance_type>::max();
  isr->add(L"μ̃", IndexBasis{uocc, P},
           120ul);  // before the Context adopts the registry (#665)
  // an orthonormal unoccupied basis other than the canonical one, e.g.
  // localized virtuals
  isr->add(L"ã", IndexBasis{uocc, 2});
  auto ctx = set_scoped_default_context(
      Context({.index_space_registry_shared_ptr = std::move(isr),
               .vacuum = Vacuum::SingleProduct}));
  const auto& registry = *get_default_context().index_space_registry();
  const Index i1(occ, 1), i2(occ, 2);
  const Index x =
      Index(uocc, 1, {i1, i2})
          .replace_basis_instance(0);  // the R2 leg of CSV-CCSD, t granted 0
  const auto f = ex<Tensor>(
      L"f", bra{x}, ket{Index(uocc, 2, {i1, i2}).replace_basis_instance(0)});

  SECTION("PAO target: non-orthonormal, minted from the registry entry") {
    Index::reset_tmp_index();
    const ExprPtr out = mbpt::csv_transform(f, registry.retrieve_basis(L"μ̃"),
                                            /*orthonormal=*/false);
    REQUIRE(out->is<Product>());
    const auto& prod = out->as<Product>();
    REQUIRE(prod.factors().size() == 3);  // f{μ̃;μ̃} C C
    const Tensor ft = prod.factor(0)->as<Tensor>();
    for (const Index& idx : ft.const_braket_indices()) {
      CHECK(idx.basis() == IndexBasis{uocc, P});
      CHECK(idx.basis().base_key() == L"μ̃");
      CHECK(idx.space().approximate_size() == 120);
      CHECK(idx.full_label().find(L'<') == std::wstring::npos);
    }
    // a hand-built basis equal to the entry is resolved to the entry too
    Index::reset_tmp_index();
    const ExprPtr out2 = mbpt::csv_transform(f, IndexBasis{uocc, P}, false);
    CHECK(out2->as<Product>()
              .factor(0)
              ->as<Tensor>()
              .bra()
              .at(0)
              .space()
              .approximate_size() == 120);
    // the overlap stays (no orthonormal shortcut)
    const ExprPtr s = mbpt::csv_transform(
        make_overlap(x, Index(uocc, 3, {i1}).replace_basis_instance(0)),
        registry.retrieve_basis(L"μ̃"), false);
    REQUIRE(s->is<Product>());
    CHECK(s->as<Product>().factors().size() == 3);
  }
  SECTION("named orthonormal target: the overlap's dummy is in that basis") {
    const ExprPtr s = mbpt::csv_transform(
        make_overlap(x, Index(uocc, 3, {i1}).replace_basis_instance(0)),
        registry.retrieve_basis(L"ã"), /*orthonormal=*/true);
    REQUIRE(s->is<Product>());
    REQUIRE(s->as<Product>().factors().size() == 2);  // C C
    const Index dummy = s->as<Product>().factor(0)->as<Tensor>().ket().at(0);
    CHECK(dummy == s->as<Product>().factor(1)->as<Tensor>().bra().at(0));
    CHECK(dummy.basis() == IndexBasis{uocc, 2});
    CHECK(dummy.basis().base_key() == L"ã");
    CHECK(!dummy.has_proto_indices());
  }
  SECTION("spin-resolved legs, as an open-shell spintrace leaves them") {
    const IndexSpace occ_a = registry.retrieve(L"i↑"),
                     uocc_a = registry.retrieve(L"a↑");
    const Index ia1(occ_a, 1), ia2(occ_a, 2);
    const ExprPtr s_a =
        make_overlap(Index(uocc_a, 1, {ia1, ia2}).replace_basis_instance(0),
                     Index(uocc_a, 3, {ia1}).replace_basis_instance(0));
    // a space target: the dummy keeps the legs' spin
    const ExprPtr out = mbpt::csv_transform(s_a, uocc);
    REQUIRE(out->is<Product>());
    REQUIRE(out->as<Product>().factors().size() == 2);  // C C
    const Index dummy = out->as<Product>().factor(0)->as<Tensor>().ket().at(0);
    CHECK(dummy == out->as<Product>().factor(1)->as<Tensor>().bra().at(0));
    CHECK(dummy.space() == uocc_a);
    CHECK(dummy.basis() == IndexBasis{uocc_a});
    // a named target in the spin-free space is not supported
    CHECK_THROWS_AS(
        mbpt::csv_transform(s_a, registry.retrieve_basis(L"ã"), true),
        Exception);
    // ... on the general path, too
    CHECK_THROWS_AS(
        mbpt::csv_transform(s_a, registry.retrieve_basis(L"μ̃"), false),
        Exception);
    CHECK_THROWS_AS(
        mbpt::csv_transform(
            ex<Tensor>(L"f", bra{Index(uocc_a, 1, {ia1, ia2})}, ket{ia1}),
            registry.retrieve_basis(L"μ̃"), false),
        Exception);
  }
  SECTION("an unnamed instance basis as target throws") {
    CHECK_THROWS_AS(mbpt::csv_transform(f, IndexBasis{uocc, 5}, false),
                    Exception);
  }
  SECTION(
      "the IndexSpace overload forwards with orthonormality read from the qns "
      "bits") {
    auto isr2 =
        mbpt::make_min_sr_spaces();  // a fresh registry: μ̃ is a space here
    mbpt::add_pao_spaces(isr2, IndexSpace::QuantumNumbers{mbpt::Spin::any});
    const IndexSpace occ2 = isr2->retrieve(L"i"), uocc2 = isr2->retrieve(L"a"),
                     mu2 = isr2->retrieve(L"μ̃");
    auto ctx2 = set_scoped_default_context(
        Context({.index_space_registry_shared_ptr = std::move(isr2),
                 .vacuum = Vacuum::SingleProduct}));
    const Index p(occ2, 1), q(occ2, 2);
    const Index xb = Index(uocc2, 1, {p, q}).replace_basis_instance(0);
    const Index xk = Index(uocc2, 3, {p}).replace_basis_instance(0);
    // PAO space: not orthonormal, the overlap stays (C s C); unoccupied MOs:
    // the shortcut (C C)
    const ExprPtr pao_out = mbpt::csv_transform(make_overlap(xb, xk), mu2);
    const ExprPtr mo_out = mbpt::csv_transform(make_overlap(xb, xk), uocc2);
    CHECK(pao_out->as<Product>().factors().size() == 3);
    CHECK(mo_out->as<Product>().factors().size() == 2);
  }
}
