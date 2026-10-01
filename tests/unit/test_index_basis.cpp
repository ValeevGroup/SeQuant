#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/index_basis.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/options.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/rules/csv.hpp>

#include <compare>
#include <cstddef>
#include <initializer_list>
#include <optional>
#include <set>
#include <string>
#include <utility>

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
    CHECK(dummy.full_label() == L"a_1");
    CHECK(dummy.basis() == IndexBasis{uocc});
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
