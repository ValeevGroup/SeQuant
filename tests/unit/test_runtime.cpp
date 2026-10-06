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
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/scope.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <atomic>
#include <cstdint>
#include <exception>
#include <functional>
#include <iostream>
#include <latch>
#include <memory>
#include <mutex>
#include <numeric>
#include <set>
#include <string>
#include <thread>
#include <vector>

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

    // the version of the installed default context is the current version
    {
      Context modified({.vacuum = Vacuum::SingleProduct});
      modified.set(SPBasis::Spinfree);
      set_default_context(modified);
      CHECK(current_context_version() == modified.version());
    }

    // set distinct contexts for fermi and bose statistics
    auto [fermi_isr, bose_isr] = mbpt::make_fermi_and_bose_spaces();
    CHECK(fermi_isr->spaces() ==
          bose_isr->spaces());  // fermi_isr and bose_isr have the same spaces
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
    CHECK(current_context_version() == initial_ctx.version());

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
    const Context ctx2({.index_space_registry = *ctx.index_space_registry()});
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

    // Context is copy-only: moving from it leaves it intact
    {
      auto source = ctx;
      const auto moved_to = std::move(source);
      CHECK(source == moved_to);
      CHECK(source.index_comparer());
      CHECK(source.tensor_canonicalizer_ptr(L""));
    }

    // a comparer re-installed through its shared pointer keeps equality
    auto ctx_same_cmp = ctx;
    ctx_same_cmp.set_index_comparer(ctx.index_comparer_ptr())
        .set_index_pair_comparer(ctx.index_pair_comparer_ptr());
    CHECK(ctx_same_cmp == ctx);
    CHECK(ctx_same_cmp.version() == ctx.version());
    CHECK(Context(ctx).set_index_comparer(ctx.index_comparer()) != ctx);

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
        {.index_space_registry = *ctx.index_space_registry(),
         .tensor_canonicalizers =
             container::map<std::wstring, std::shared_ptr<TensorCanonicalizer>>{
                 {L"Q", null_canon}},
         .index_comparer = ctx_cmp.index_comparer(),
         .cardinal_tensor_labels = container::vector<std::wstring>{L"Z"}});
    CHECK(ctx_opts.tensor_canonicalizer_ptr(L"Q") == null_canon);
    // the empty label maps to the default canonicalizer unless given
    CHECK(ctx_opts.tensor_canonicalizer_ptr(L"R") ==
          Context{}.tensor_canonicalizer_ptr(L""));
    CHECK(Context({.tensor_canonicalizers =
                       container::map<std::wstring,
                                      std::shared_ptr<TensorCanonicalizer>>{
                           {L"", null_canon}}})
              .tensor_canonicalizer_ptr(L"R") == null_canon);
    CHECK(ctx_opts.cardinal_tensor_labels() ==
          container::vector<std::wstring>{L"Z"});
    CHECK(ctx_opts.index_comparer());
    CHECK(ctx_opts.index_pair_comparer());
    CHECK(ctx_opts != ctx);

    // null canonicalizers are rejected
    CHECK_THROWS_AS(
        Context({.tensor_canonicalizers =
                     container::map<std::wstring,
                                    std::shared_ptr<TensorCanonicalizer>>{
                         {L"Q", nullptr}}}),
        Exception);
    CHECK_THROWS_AS(Context{}.set_tensor_canonicalizer(L"Q", nullptr),
                    Exception);
  }

  SECTION("index space registry is owned") {
    // a shared_ptr whose object has other owners is copied, so modifying the
    // object does not affect the context
    auto shared = mbpt::make_sr_spaces();
    Context ctx({.index_space_registry_shared_ptr = shared});
    CHECK(ctx.index_space_registry().get() != shared.get());
    CHECK(*ctx.index_space_registry() == *shared);
    shared->add(L"q", 0b10000);
    CHECK(!ctx.index_space_registry()->contains(L"q"));
    ctx.set(shared);
    CHECK(ctx.index_space_registry().get() != shared.get());
    CHECK(ctx.index_space_registry()->contains(L"q"));

    // the only owner of its object is adopted, without a copy
    {
      auto unique = mbpt::make_sr_spaces();
      const auto* object = unique.get();
      const Context adopted(
          {.index_space_registry_shared_ptr = std::move(unique)});
      CHECK(adopted.index_space_registry().get() == object);
      auto set_unique = mbpt::make_sr_spaces();
      const auto* set_object = set_unique.get();
      ctx.set(std::move(set_unique));
      CHECK(ctx.index_space_registry().get() == set_object);
    }

    // a registry given by value is moved in, keeping its storage
    {
      IndexSpaceRegistry by_value = *mbpt::make_sr_spaces();
      const auto* storage = &*by_value.begin();
      const Context from_value({.index_space_registry = std::move(by_value)});
      CHECK(&*from_value.index_space_registry()->begin() == storage);
      IndexSpaceRegistry set_value = *mbpt::make_sr_spaces();
      const auto* set_storage = &*set_value.begin();
      ctx.set(std::move(set_value));
      CHECK(&*ctx.index_space_registry()->begin() == set_storage);
    }

    // the Options given to set_default_context() and
    // set_scoped_default_context() hand over their registry the same way
    {
      const Context initial_ctx = get_default_context_snapshot();
      auto unique = mbpt::make_sr_spaces();
      const auto* object = unique.get();
      set_default_context(
          {.index_space_registry_shared_ptr = std::move(unique)});
      CHECK(get_default_context().index_space_registry().get() == object);
      IndexSpaceRegistry by_value = *mbpt::make_sr_spaces();
      const auto* storage = &*by_value.begin();
      set_default_context({.index_space_registry = std::move(by_value)});
      CHECK(&*get_default_context().index_space_registry()->begin() == storage);
      set_default_context(initial_ctx);

      auto scoped_unique = mbpt::make_sr_spaces();
      const auto* scoped_object = scoped_unique.get();
      const auto resetter = set_scoped_default_context(
          {.index_space_registry_shared_ptr = std::move(scoped_unique)});
      CHECK(get_default_context().index_space_registry().get() ==
            scoped_object);
    }

    // copies of a context share its registry
    const Context copy(ctx);
    CHECK(copy.index_space_registry() == ctx.index_space_registry());

    // to modify a context's registry, modify a copy and set it
    IndexSpaceRegistry modified = *ctx.index_space_registry();
    modified.add(L"q", 0b10000);
    ctx.set(std::move(modified));
    CHECK(ctx.index_space_registry()->contains(L"q"));
  }

  SECTION("version") {
    Context ctx;
    const auto v0 = ctx.version();
    CHECK(v0 != 0);
    // the version identifies the canonicalization configuration
    CHECK(Context{}.version() == v0);
    CHECK(Context({.vacuum = Vacuum::SingleProduct}).version() == v0);

    // copies keep the version
    const Context copy(ctx);
    CHECK(copy.version() == v0);
    Context assigned;
    assigned = ctx;
    CHECK(assigned.version() == v0);
    const Context with_registry({.index_space_registry = IndexSpaceRegistry{}});
    CHECK(with_registry.version() != v0);
    // contexts that share a registry, such as a copy with another vacuum, share
    // the version
    Context same_registry(with_registry);
    same_registry.set(Vacuum::SingleProduct);
    CHECK(same_registry.version() == with_registry.version());
    // registries are compared by value, the version keys on the registry
    // object; a context need not have one
    CHECK(Context({.index_space_registry = IndexSpaceRegistry{}}) ==
          with_registry);
    CHECK(Context({.index_space_registry = IndexSpaceRegistry{}}).version() !=
          with_registry.version());
    {
      IndexSpaceRegistry other;
      other.add(L"q", 0b01);
      CHECK(Context({.index_space_registry = other}) != with_registry);
      // including the approximate sizes of their spaces
      IndexSpaceRegistry resized = other;
      resized.retrieve_ptr(L"q")->approximate_size(
          other.retrieve(L"q").approximate_size() + 1);
      CHECK(Context({.index_space_registry = std::move(resized)}) !=
            Context({.index_space_registry = std::move(other)}));
    }
    CHECK(Context{} == Context{});
    CHECK(Context{} != with_registry);

    // the version of a configuration whose objects are gone is not reused,
    // even if a new object takes the address of a dead one
    {
      std::uint64_t dead_version = 0;
      {
        const Context dead({.index_space_registry = IndexSpaceRegistry{}});
        dead_version = dead.version();
      }
      for (int i = 0; i != 8; ++i)
        CHECK(
            Context({.index_space_registry = IndexSpaceRegistry{}}).version() !=
            dead_version);
    }

    // the canonicalization settings change the version, the others do not
    auto changes = [&ctx](auto&& set) {
      const auto before = ctx.version();
      set(ctx);
      return ctx.version() != before;
    };
    CHECK(!changes([](Context& c) { c.set(Vacuum::SingleProduct); }));
    CHECK(!changes([](Context& c) { c.set(IndexSpaceMetric::General); }));
    CHECK(!changes([](Context& c) { c.set(AssertStrictBraKetSymmetry::No); }));
    CHECK(!changes([](Context& c) { c.set_first_dummy_index_ordinal(200); }));
    CHECK(!changes([](Context& c) { c.set(BraKetTypesetting::KetSub); }));
    CHECK(!changes([](Context& c) { c.set(BraKetSlotTypesetting::Naive); }));
    CHECK(!changes([](Context& c) { c.set(Symmetry::Symm); }));
    CHECK(!changes([](Context& c) { c.set(Hermiticity::Hermitian); }));
    CHECK(!changes([](Context& c) { c.set(ColumnSymmetry::Symm); }));
    CHECK(changes([](Context& c) { c.set(SPBasis::Spinfree); }));
    CHECK(changes([](Context& c) { c.set(IndexSpaceRegistry{}); }));
    CHECK(changes(
        [](Context& c) { c.set(std::make_shared<IndexSpaceRegistry>()); }));
    CHECK(changes([](Context& c) { c.set(CanonicalizeOptions{}); }));
    CHECK(changes([](Context& c) {
      c.set(CanonicalizeOptions::default_options().copy_and_set(
          CanonicalizeOptions::IgnoreNamedIndexLabel::No));
    }));
    {
      const auto before = ctx.version();
      CHECK(changes([](Context& c) {
        c.set_tensor_canonicalizer(L"Q",
                                   std::make_shared<NullTensorCanonicalizer>());
      }));
      CHECK(changes([](Context& c) { c.unset_tensor_canonicalizer(L"Q"); }));
      // undoing a change restores the version
      CHECK(ctx.version() == before);
    }
    CHECK(changes([](Context& c) {
      c.set_index_comparer(TensorCanonicalizer::default_index_comparer());
    }));
    CHECK(changes([](Context& c) {
      c.set_index_pair_comparer(
          TensorCanonicalizer::default_index_pair_comparer());
    }));
    CHECK(!changes(
        [](Context& c) { c.set_index_comparer(c.index_comparer_ptr()); }));
    CHECK(changes([](Context& c) {
      c.set_cardinal_tensor_labels(container::vector<std::wstring>{L"Z"});
    }));

    // the copy was not affected
    CHECK(copy.version() == v0);
  }

  SECTION("snapshot") {
    const auto initial_ctx = get_default_context();
    auto restore = sequant::detail::make_scope_exit(
        [&initial_ctx] { set_default_context(initial_ctx); });

    // the process-wide context is the sole owner of its labels and comparer
    set_default_context(
        Context(initial_ctx)
            .set_cardinal_tensor_labels(container::vector<std::wstring>{L"X"})
            .set_index_comparer(TensorCanonicalizer::index_comparer_t(
                [](const Index& i1, const Index& i2) {
                  return i2.label() < i1.label();
                })));
    const auto snapshot = get_default_context_snapshot();
    CHECK(snapshot == get_default_context());
    CHECK(snapshot.version() == current_context_version());
    const auto& labels = snapshot.cardinal_tensor_labels();
    const auto& comparer = snapshot.index_comparer();

    // replacing the process-wide context leaves the snapshot intact
    set_default_context(initial_ctx);
    CHECK(labels == container::vector<std::wstring>{L"X"});
    CHECK(comparer(Index(L"i_2"), Index(L"i_1")));

    // the next snapshot sees every change of the process-wide context,
    // including one made on another thread
    std::thread([&initial_ctx] {
      set_default_context(Context(initial_ctx)
                              .set_cardinal_tensor_labels(
                                  container::vector<std::wstring>{L"Y"}));
    }).join();
    CHECK(get_default_context_snapshot().cardinal_tensor_labels() ==
          container::vector<std::wstring>{L"Y"});
    std::thread([] { reset_default_context(); }).join();
    CHECK(get_default_context_snapshot().cardinal_tensor_labels() ==
          Context{}.cardinal_tensor_labels());
    set_default_context(initial_ctx);

    // on a thread with a scoped context, a copy of that context
    auto scoped = set_scoped_default_context(snapshot);
    CHECK(get_default_context_snapshot() == snapshot);
  }
}

TEST_CASE("scoped contexts", "[runtime]") {
  using namespace sequant;

  auto q_canonicalizer = [] {
    return get_default_context().nondefault_tensor_canonicalizer_ptr(L"Q");
  };
  auto with_q_canonicalizer =
      [](std::shared_ptr<TensorCanonicalizer> canonicalizer) {
        return Context(get_default_context())
            .set_tensor_canonicalizer(L"Q", std::move(canonicalizer));
      };
  REQUIRE(!q_canonicalizer());

  SECTION("each thread sees only its own") {
    constexpr int nthreads = 6;
    std::latch all_scoped(nthreads);
    std::atomic<int> mismatches = 0;
    std::vector<std::thread> threads;
    for (int t = 0; t != nthreads; ++t) {
      threads.emplace_back([&] {
        const auto mine = std::make_shared<NullTensorCanonicalizer>();
        {
          auto scoped = set_scoped_default_context(with_q_canonicalizer(mine));
          all_scoped.arrive_and_wait();
          for (int i = 0; i != 1000; ++i)
            if (q_canonicalizer() != mine) ++mismatches;
        }
        if (q_canonicalizer()) ++mismatches;
      });
    }
    for (auto& thread : threads) thread.join();
    CHECK(mismatches == 0);
    CHECK(!q_canonicalizer());
  }

  SECTION("parallel workers see the scoped context of the caller") {
    const auto mine = std::make_shared<NullTensorCanonicalizer>();
    std::vector<int> items(64);
    std::atomic<int> mismatches = 0;
    {
      auto scoped = set_scoped_default_context(with_q_canonicalizer(mine));
      auto scoped_mbpt = mbpt::set_scoped_default_mbpt_context(
          mbpt::Context({.csv = mbpt::CSV::Yes}));
      sequant::for_each(items, [&](int&) {
        if (q_canonicalizer() != mine) ++mismatches;
        if (mbpt::get_default_mbpt_context().csv() != mbpt::CSV::Yes)
          ++mismatches;
        // nested scopes compose
        const auto inner = std::make_shared<NullTensorCanonicalizer>();
        {
          auto scoped = set_scoped_default_context(with_q_canonicalizer(inner));
          if (q_canonicalizer() != inner) ++mismatches;
        }
        if (q_canonicalizer() != mine) ++mismatches;
      });
      CHECK(sequant::transform_reduce(items, 0, std::plus<int>{}, [&](int) {
              return q_canonicalizer() == mine ? 0 : 1;
            }) == 0);
      CHECK(q_canonicalizer() == mine);
    }
    CHECK(mismatches == 0);

    // after the scope workers see the process-wide context again
    sequant::for_each(items, [&](int&) {
      if (q_canonicalizer()) ++mismatches;
    });
    CHECK(mismatches == 0);
  }

  // unlike for_each, whose execution-policy backend may run serially,
  // parallel_do always runs on num_threads() distinct threads
  SECTION("parallel_do threads see the scoped context of the caller") {
    const auto mine = std::make_shared<NullTensorCanonicalizer>();
    std::atomic<int> invocations = 0;
    std::atomic<int> mismatches = 0;
    std::mutex thread_ids_mtx;
    std::set<std::thread::id> thread_ids;
    {
      auto scoped = set_scoped_default_context(with_q_canonicalizer(mine));
      sequant::parallel_do([&](int) {
        ++invocations;
        {
          std::scoped_lock lock(thread_ids_mtx);
          thread_ids.insert(std::this_thread::get_id());
        }
        if (q_canonicalizer() != mine) ++mismatches;
      });
      CHECK(q_canonicalizer() == mine);
    }
    CHECK(invocations == num_threads());
    CHECK(mismatches == 0);
    if (num_threads() > 1) CHECK(thread_ids.size() > 1);

    sequant::parallel_do([&](int) {
      if (q_canonicalizer()) ++mismatches;
    });
    CHECK(mismatches == 0);
  }

  SECTION("a modification applies to the context of every statistics") {
    const auto labels = [](Statistics s) {
      return get_default_context(s).cardinal_tensor_labels();
    };
    const container::vector<std::wstring> arbitrary_labels{L"A"};
    const container::vector<std::wstring> fermi_labels{L"F"};
    auto scoped =
        set_scoped_default_context(container::map<Statistics, Context>{
            {Statistics::Arbitrary,
             Context(get_default_context())
                 .set_cardinal_tensor_labels(arbitrary_labels)},
            {Statistics::FermiDirac,
             Context(get_default_context())
                 .set_cardinal_tensor_labels(fermi_labels)}});
    const auto mine = std::make_shared<NullTensorCanonicalizer>();
    auto q_canonicalizer_of = [](Statistics s) {
      return get_default_context(s).nondefault_tensor_canonicalizer_ptr(L"Q");
    };
    {
      auto modified = set_scoped_modified_default_context(
          [&mine](Context& ctx) { ctx.set_tensor_canonicalizer(L"Q", mine); });
      CHECK(labels(Statistics::Arbitrary) == arbitrary_labels);
      CHECK(labels(Statistics::FermiDirac) == fermi_labels);
      CHECK(labels(Statistics::BoseEinstein) == arbitrary_labels);
      CHECK(q_canonicalizer_of(Statistics::Arbitrary) == mine);
      CHECK(q_canonicalizer_of(Statistics::FermiDirac) == mine);
    }
    CHECK(labels(Statistics::Arbitrary) == arbitrary_labels);
    CHECK(labels(Statistics::FermiDirac) == fermi_labels);
    CHECK(!q_canonicalizer_of(Statistics::Arbitrary));
    CHECK(!q_canonicalizer_of(Statistics::FermiDirac));
  }

  SECTION("mbpt::load refuses to run under a scoped context") {
    const auto process_wide_version = current_context_version();
    {
      auto scoped = set_scoped_default_context(
          with_q_canonicalizer(std::make_shared<NullTensorCanonicalizer>()));
      CHECK_THROWS_AS(mbpt::load(), Exception);
    }
    CHECK(current_context_version() == process_wide_version);
    CHECK(!q_canonicalizer());
  }

  // e.g. the scopes of successive top-level WickTheorems, which map the
  // normal operator labels to one shared NullTensorCanonicalizer
  SECTION("equal modifications of a context scope equal versions") {
    const auto& shared = NullTensorCanonicalizer::instance();
    REQUIRE(shared == NullTensorCanonicalizer::instance());
    auto scoped_version = [](std::shared_ptr<TensorCanonicalizer> c) {
      auto modified = set_scoped_modified_default_context(
          [&c](Context& ctx) { ctx.set_tensor_canonicalizer(L"Q", c); });
      return current_context_version();
    };
    const auto v1 = scoped_version(shared);
    CHECK(v1 != current_context_version());
    CHECK(scoped_version(shared) == v1);
    CHECK(scoped_version(std::make_shared<NullTensorCanonicalizer>()) != v1);
  }

  SECTION("the current version follows the effective context") {
    const auto v0 = current_context_version();
    CHECK(v0 == get_default_context().version());

    const Context other;
    {
      auto scoped = set_scoped_default_context(other);
      CHECK(current_context_version() == other.version());
      CHECK(current_context_version() != v0);
      {
        auto modified = set_scoped_modified_default_context([](Context& ctx) {
          ctx.set_cardinal_tensor_labels(container::vector<std::wstring>{L"Z"});
        });
        CHECK(current_context_version() != other.version());
        CHECK(current_context_version() != v0);
        CHECK(current_context_version(Statistics::FermiDirac) ==
              current_context_version());
      }
      CHECK(current_context_version() == other.version());

      // a thread without the scope sees the process-wide default
      std::uint64_t seen = 0;
      std::thread([&seen] { seen = current_context_version(); }).join();
      CHECK(seen == v0);
    }
    CHECK(current_context_version() == v0);
  }
}

TEST_CASE("parallel exceptions", "[runtime]") {
  using namespace sequant;

  struct ItemError : Exception {
    explicit ItemError(int item)
        : Exception("item " + std::to_string(item)), item(item) {}
    int item;
  };
  auto rethrown_item = [](const std::exception_ptr& e) {
    try {
      std::rethrow_exception(e);
    } catch (const ItemError& error) {
      return error.item;
    } catch (...) {
    }
    return -1;
  };

  const auto nthreads = num_threads();
  set_num_threads(4);
  auto restore_nthreads = sequant::detail::make_scope_exit(
      [nthreads] { set_num_threads(nthreads); });

  std::vector<int> items(16);
  std::iota(items.begin(), items.end(), 0);

  SECTION("for_each rethrows a single exception as is") {
    std::atomic<int> ran = 0;
    CHECK_THROWS_AS(sequant::for_each(items,
                                      [&ran](int& i) {
                                        ++ran;
                                        if (i == 5) throw ItemError(i);
                                      }),
                    ItemError);
    // the other items were not abandoned
    CHECK(ran == 16);
  }

  SECTION("for_each collects several exceptions, by item") {
    try {
      sequant::for_each(items, [](int& i) {
        if (i % 5 == 1) throw ItemError(i);
      });
      FAIL("no exception");
    } catch (const ParallelExceptions& e) {
      REQUIRE(e.exceptions().size() == 3);
      for (std::size_t k = 0; k != 3; ++k) {
        CHECK(e.exceptions()[k].first == 5 * k + 1);
        CHECK(rethrown_item(e.exceptions()[k].second) ==
              static_cast<int>(5 * k + 1));
      }
      CHECK(std::string(e.what()).find("item 1") != std::string::npos);
    }
  }

  SECTION("nested exceptions are attributed to the outer item") {
    std::vector<int> outer{0, 1, 2};
    try {
      sequant::for_each(outer, [&items](int& o) {
        if (o == 0) return;
        sequant::for_each(items, [o](int& i) {
          if (i < o) throw ItemError(i);
        });
      });
      FAIL("no exception");
    } catch (const ParallelExceptions& e) {
      // item 1 throws once (rethrown as is), item 2 twice (ParallelExceptions)
      REQUIRE(e.exceptions().size() == 3);
      CHECK(e.exceptions()[0].first == 1);
      CHECK(e.exceptions()[1].first == 2);
      CHECK(e.exceptions()[2].first == 2);
    }
  }

  SECTION("transform_reduce") {
    CHECK(sequant::transform_reduce(items, 0, std::plus<int>{},
                                    [](int i) { return i; }) == 120);
    CHECK_THROWS_AS(sequant::transform_reduce(items, 0, std::plus<int>{},
                                              [](int i) {
                                                if (i == 7) throw ItemError(i);
                                                return i;
                                              }),
                    ItemError);
    CHECK_THROWS_AS(sequant::transform_reduce(items, 0, std::plus<int>{},
                                              [](int i) {
                                                if (i > 13) throw ItemError(i);
                                                return i;
                                              }),
                    ParallelExceptions);
  }

  SECTION("parallel_do") {
    CHECK_THROWS_AS(sequant::parallel_do([](int thread_id) {
                      if (thread_id == 0) throw ItemError(thread_id);
                    }),
                    ItemError);
    try {
      sequant::parallel_do([](int thread_id) { throw ItemError(thread_id); });
      FAIL("no exception");
    } catch (const ParallelExceptions& e) {
      REQUIRE(e.exceptions().size() == 4);
      for (std::size_t k = 0; k != 4; ++k)
        CHECK(rethrown_item(e.exceptions()[k].second) == static_cast<int>(k));
    }
  }
}
