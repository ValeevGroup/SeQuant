//
// Created by Eduard Valeyev on 10/12/22.
//

#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <iostream>
#include <memory>
#include <string>

TEST_CASE("context", "[runtime]") {
  using namespace sequant;

  SECTION("constructors") { REQUIRE_NOTHROW(Context{}); }

  SECTION("default context") {
    CHECK_NOTHROW(get_default_context());
    auto initial_ctx = get_default_context();

    // basic set_default_context test
    CHECK_NOTHROW(set_default_context(
        {.index_space_registry_shared_ptr = mbpt::make_sr_spaces(),
         .vacuum = Vacuum::SingleProduct,
         .metric = IndexSpaceMetric::Unit,
         .spbasis = SPBasis::Spinfree}));
    CHECK(get_default_context().vacuum() == Vacuum::SingleProduct);
    CHECK(get_default_context().metric() == IndexSpaceMetric::Unit);
    CHECK(get_default_context().spbasis() == SPBasis::Spinfree);

    // set distinct contexts for fermi and bose statistics
    auto [fermi_isr, bose_isr] = mbpt::make_fermi_and_bose_spaces();
    CHECK(fermi_isr->spaces() ==
          bose_isr->spaces());  // fermi_isr and bose_isr share the space set
    CHECK_NOTHROW(set_default_context(
        {{Statistics::FermiDirac,
          Context({.index_space_registry_shared_ptr = fermi_isr,
                   .vacuum = Vacuum::SingleProduct})},
         {Statistics::BoseEinstein,
          Context({.index_space_registry_shared_ptr = bose_isr,
                   .vacuum = Vacuum::Physical})}}));
    CHECK(get_default_context(Statistics::Arbitrary).vacuum() ==
          Vacuum::SingleProduct);
    CHECK(get_default_context(Statistics::FermiDirac).vacuum() ==
          Vacuum::SingleProduct);
    CHECK(get_default_context(Statistics::BoseEinstein).vacuum() ==
          Vacuum::Physical);

    // reset back to default
    CHECK_NOTHROW(reset_default_context());
    CHECK(get_default_context().vacuum() == Vacuum::Physical);
    CHECK(get_default_context().metric() == IndexSpaceMetric::Unit);
    CHECK(get_default_context().spbasis() == SPBasis::Spinor);

    // reset back to initial context
    CHECK_NOTHROW(set_default_context(initial_ctx));
    CHECK(get_default_context() == initial_ctx);
    CHECK(!(get_default_context() != initial_ctx));

    // scoped changes to default context
    {
      // if we do not save the resetter context is reset back immediately
      CHECK_NOTHROW(set_scoped_default_context(
          {.index_space_registry_shared_ptr = mbpt::make_sr_spaces(),
           .vacuum = Vacuum::SingleProduct,
           .metric = IndexSpaceMetric::Unit,
           .spbasis = SPBasis::Spinfree}));
      CHECK(get_default_context() == initial_ctx);

      auto ctx = get_default_context();
      ctx.set(mbpt::make_sr_spaces());
      ctx.set(SPBasis::Spinfree);
      const auto ctx_copy = ctx;
      auto resetter = set_scoped_default_context(ctx);
      CHECK(get_default_context() == ctx_copy);
    }
    // leaving scope resets the context back
    CHECK(get_default_context() == initial_ctx);
  }

  SECTION("tensor canonicalizers") {
    const Context ctx(
        {.index_space_registry_shared_ptr = mbpt::make_sr_spaces()});

    // defaults
    REQUIRE(ctx.tensor_canonicalizer_ptr(L""));
    CHECK(std::dynamic_pointer_cast<DefaultTensorCanonicalizer>(
        ctx.tensor_canonicalizer_ptr(L"")));
    CHECK(ctx.tensor_canonicalizer_ptr(L"anything") ==
          ctx.tensor_canonicalizer_ptr(L""));
    CHECK(&ctx.tensor_canonicalizer(L"anything") ==
          ctx.tensor_canonicalizer_ptr(L"").get());
    CHECK(ctx.nondefault_tensor_canonicalizer_ptr(L"anything") == nullptr);
    CHECK(ctx.cardinal_tensor_labels() ==
          container::vector<std::wstring>{reserved::antisymm_label(),
                                          reserved::symm_label(),
                                          reserved::transposition_label()});
    {
      const Index i1(L"i_1"), i2(L"i_2"), a1(L"a_1");
      const auto& cmp = ctx.index_comparer();
      REQUIRE(cmp);
      CHECK(cmp(i1, i2));
      CHECK(!cmp(i2, i1));
      CHECK(!cmp(i1, i1));
      CHECK(cmp(i1, a1) == (i1.space() < a1.space()));
      const auto& paircmp = ctx.index_pair_comparer();
      REQUIRE(paircmp);
      CHECK(paircmp({i1, a1}, {i2, a1}));
      CHECK(!paircmp({i2, a1}, {i1, a1}));
    }

    // two default-constructed states compare equal
    const Context ctx2({.index_space_registry_shared_ptr =
                            ctx.mutable_index_space_registry()});
    CHECK(ctx == ctx2);

    // setters act on a copy only
    auto ctx_q = ctx;
    CHECK(ctx_q == ctx);
    const auto null_canon = std::make_shared<NullTensorCanonicalizer>();
    ctx_q.set_tensor_canonicalizer(L"Q", null_canon);
    CHECK(ctx_q.nondefault_tensor_canonicalizer_ptr(L"Q") == null_canon);
    CHECK(ctx_q.tensor_canonicalizer_ptr(L"Q") == null_canon);
    CHECK(ctx_q.tensor_canonicalizer_ptr(L"R") ==
          ctx.tensor_canonicalizer_ptr(L""));
    CHECK(ctx.nondefault_tensor_canonicalizer_ptr(L"Q") == nullptr);
    CHECK(ctx_q != ctx);
    ctx_q.unset_tensor_canonicalizer(L"Q");
    CHECK(ctx_q.nondefault_tensor_canonicalizer_ptr(L"Q") == nullptr);
    CHECK(ctx_q == ctx);

    // no canonicalizer at all
    auto ctx_none = ctx;
    ctx_none.unset_tensor_canonicalizer(L"");
    CHECK(ctx_none.tensor_canonicalizer_ptr(L"anything") == nullptr);
    CHECK_THROWS_AS(ctx_none.tensor_canonicalizer(L"anything"), Exception);
    CHECK(ctx_none != ctx);

    // comparers compare by identity: a behaviourally identical replacement is
    // unequal, while a copy is equal
    auto ctx_cmp = ctx;
    ctx_cmp.set_index_comparer(TensorCanonicalizer::index_comparer_t(
        [](const Index& a, const Index& b) { return a < b; }));
    CHECK(ctx_cmp != ctx);
    const auto ctx_cmp_copy = ctx_cmp;
    CHECK(ctx_cmp_copy == ctx_cmp);
    CHECK(ctx_cmp.index_comparer()(Index(L"i_1"), Index(L"i_2")));
    auto ctx_paircmp = ctx;
    ctx_paircmp.set_index_pair_comparer(
        TensorCanonicalizer::index_pair_comparer_t(
            [](const TensorCanonicalizer::index_pair_t& a,
               const TensorCanonicalizer::index_pair_t b) {
              return a.first < b.first;
            }));
    CHECK(ctx_paircmp != ctx);

    // cardinal labels
    auto ctx_card = ctx;
    ctx_card.set_cardinal_tensor_labels({L"X", L"Y"});
    CHECK(ctx_card.cardinal_tensor_labels() ==
          container::vector<std::wstring>{L"X", L"Y"});
    CHECK(ctx_card != ctx);
    ctx_card.set_cardinal_tensor_labels(ctx.cardinal_tensor_labels());
    CHECK(ctx_card == ctx);

    // chaining
    auto ctx_chain = ctx;
    CHECK(&ctx_chain.set_tensor_canonicalizer(L"Q", null_canon)
               .unset_tensor_canonicalizer(L"Q")
               .set_cardinal_tensor_labels(ctx.cardinal_tensor_labels())
               .set_index_comparer(ctx.index_comparer())
               .set_index_pair_comparer(ctx.index_pair_comparer()) ==
          &ctx_chain);

    // named-parameter construction
    const Context ctx_opts(
        {.index_space_registry_shared_ptr = ctx.mutable_index_space_registry(),
         .tensor_canonicalizers =
             container::map<std::wstring, std::shared_ptr<TensorCanonicalizer>>{
                 {L"Q", null_canon}},
         .index_comparer = ctx_cmp.index_comparer(),
         .cardinal_tensor_labels = container::vector<std::wstring>{L"Z"}});
    CHECK(ctx_opts.tensor_canonicalizer_ptr(L"Q") == null_canon);
    CHECK(ctx_opts.tensor_canonicalizer_ptr(L"R") == nullptr);
    CHECK(ctx_opts.cardinal_tensor_labels() ==
          container::vector<std::wstring>{L"Z"});
    CHECK(ctx_opts.index_comparer());
    CHECK(ctx_opts.index_pair_comparer());
    CHECK(ctx_opts != ctx);
  }
}
