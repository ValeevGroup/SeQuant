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
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/scope.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <iostream>
#include <latch>
#include <memory>
#include <mutex>
#include <new>
#include <numeric>
#include <set>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>

TEST_CASE("context", "[runtime]") {
  using namespace sequant;

  SECTION("constructors") { REQUIRE_NOTHROW(Context{}); }

  SECTION("default context") {
    CHECK_NOTHROW(get_default_context());
    auto initial_ctx = get_default_context();

    // basic set_default_context test
    CHECK_NOTHROW(set_default_context(
        {.index_basis_registry_shared_ptr = mbpt::make_sr_spaces(),
         .vacuum = Vacuum::SingleProduct,
         .spbasis = SPBasis::Spinfree}));
    CHECK(get_default_context().vacuum() == Vacuum::SingleProduct);
    CHECK(get_default_context().spbasis() == SPBasis::Spinfree);

    // the version of the installed default context is the current version
    {
      Context modified({.vacuum = Vacuum::SingleProduct});
      modified.set(SPBasis::Spinfree);
      set_default_context(modified);
      CHECK(current_context_version() == modified.version());
    }

    // set distinct contexts for fermi and bose statistics
    auto [fermi_ibr, bose_ibr] = mbpt::make_fermi_and_bose_spaces();
    CHECK(fermi_ibr->bases() ==
          bose_ibr->bases());  // fermi_ibr and bose_ibr have the same spaces
    CHECK_NOTHROW(set_default_context(
        {{Statistics::FermiDirac,
          Context({.index_basis_registry_shared_ptr = fermi_ibr,
                   .vacuum = Vacuum::SingleProduct})},
         {Statistics::BoseEinstein,
          Context({.index_basis_registry_shared_ptr = bose_ibr,
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
          {.index_basis_registry_shared_ptr = mbpt::make_sr_spaces(),
           .vacuum = Vacuum::SingleProduct,
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

  SECTION("deprecated index space registry spellings forward") {
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    static_assert(std::is_same_v<IndexSpaceRegistry, IndexBasisRegistry>);
    auto shared = mbpt::make_sr_spaces();
    const Context from_shared({.index_space_registry_shared_ptr = shared,
                               .vacuum = Vacuum::Physical});
    REQUIRE(from_shared.index_space_registry());
    CHECK(from_shared.index_space_registry() ==
          from_shared.index_basis_registry());
    CHECK(*from_shared.index_basis_registry() == *shared);
    CHECK(from_shared.vacuum() == Vacuum::Physical);
    const Context from_object(
        {.index_space_registry = IndexBasisRegistry(*shared)});
    REQUIRE(from_object.index_basis_registry());
    CHECK(*from_object.index_basis_registry() == *shared);
    auto scoped = set_scoped_default_context(from_shared);
    CHECK(get_default_index_space_registry() ==
          get_default_index_basis_registry());
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  }

  SECTION("index basis registry is owned") {
    // a shared_ptr whose object has other owners is copied, so modifying the
    // object does not affect the context
    auto shared = mbpt::make_sr_spaces();
    Context ctx({.index_basis_registry_shared_ptr = shared});
    CHECK(ctx.index_basis_registry().get() != shared.get());
    CHECK(*ctx.index_basis_registry() == *shared);
    shared->add(L"q", 0b10000);
    CHECK(!ctx.index_basis_registry()->contains(L"q"));
    ctx.set(shared);
    CHECK(ctx.index_basis_registry().get() != shared.get());
    CHECK(ctx.index_basis_registry()->contains(L"q"));

    // the only owner of its object is adopted, without a copy
    {
      auto unique = mbpt::make_sr_spaces();
      const auto* object = unique.get();
      const Context adopted(
          {.index_basis_registry_shared_ptr = std::move(unique)});
      CHECK(adopted.index_basis_registry().get() == object);
      auto set_unique = mbpt::make_sr_spaces();
      const auto* set_object = set_unique.get();
      ctx.set(std::move(set_unique));
      CHECK(ctx.index_basis_registry().get() == set_object);
    }

    // a registry given by value is moved in, keeping its storage
    {
      IndexBasisRegistry by_value = *mbpt::make_sr_spaces();
      const auto* storage = &*by_value.begin();
      const Context from_value({.index_basis_registry = std::move(by_value)});
      CHECK(&*from_value.index_basis_registry()->begin() == storage);
      IndexBasisRegistry set_value = *mbpt::make_sr_spaces();
      const auto* set_storage = &*set_value.begin();
      ctx.set(std::move(set_value));
      CHECK(&*ctx.index_basis_registry()->begin() == set_storage);
    }

    // the Options given to set_default_context() and
    // set_scoped_default_context() hand over their registry the same way
    {
      const Context initial_ctx = get_default_context_snapshot();
      auto unique = mbpt::make_sr_spaces();
      const auto* object = unique.get();
      set_default_context(
          {.index_basis_registry_shared_ptr = std::move(unique)});
      CHECK(get_default_context().index_basis_registry().get() == object);
      IndexBasisRegistry by_value = *mbpt::make_sr_spaces();
      const auto* storage = &*by_value.begin();
      set_default_context({.index_basis_registry = std::move(by_value)});
      CHECK(&*get_default_context().index_basis_registry()->begin() == storage);
      set_default_context(initial_ctx);

      auto scoped_unique = mbpt::make_sr_spaces();
      const auto* scoped_object = scoped_unique.get();
      const auto resetter = set_scoped_default_context(
          {.index_basis_registry_shared_ptr = std::move(scoped_unique)});
      CHECK(get_default_context().index_basis_registry().get() ==
            scoped_object);
    }

    // copies of a context share its registry
    const Context copy(ctx);
    CHECK(copy.index_basis_registry() == ctx.index_basis_registry());

    // to modify a context's registry, modify a copy and set it
    IndexBasisRegistry modified = *ctx.index_basis_registry();
    modified.add(L"q", 0b10000);
    ctx.set(std::move(modified));
    CHECK(ctx.index_basis_registry()->contains(L"q"));
  }

  SECTION("version") {
    Context ctx;
    const auto v0 = ctx.version();
    CHECK(v0 != 0);
    // the version identifies the configuration canonicalization depends on:
    // the index basis registry and the single-particle basis
    CHECK(Context{}.version() == v0);
    CHECK(Context({.vacuum = Vacuum::SingleProduct}).version() == v0);

    // copies keep the version
    const Context copy(ctx);
    CHECK(copy.version() == v0);
    Context assigned;
    assigned = ctx;
    CHECK(assigned.version() == v0);
    const Context with_registry({.index_basis_registry = IndexBasisRegistry{}});
    CHECK(with_registry.version() != v0);
    // contexts that share a registry, such as a copy with another vacuum, share
    // the version
    Context same_registry(with_registry);
    same_registry.set(Vacuum::SingleProduct);
    CHECK(same_registry.version() == with_registry.version());
    // registries are compared by value, the version keys on the registry
    // object; a context need not have one
    CHECK(Context({.index_basis_registry = IndexBasisRegistry{}}) ==
          with_registry);
    CHECK(Context({.index_basis_registry = IndexBasisRegistry{}}).version() !=
          with_registry.version());
    {
      IndexBasisRegistry other;
      other.add(L"q", 0b01);
      CHECK(Context({.index_basis_registry = other}) != with_registry);
      // including the approximate sizes of their spaces
      IndexBasisRegistry resized = other;
      resized.extent(L"q", other.retrieve_space(L"q").dimension() + 1);
      CHECK(Context({.index_basis_registry = std::move(resized)}) !=
            Context({.index_basis_registry = std::move(other)}));
    }
    CHECK(Context{} == Context{});
    CHECK(Context{} != with_registry);

    // the version of a configuration whose objects are gone is not reused,
    // even if a new object takes the address of a dead one: the registries
    // are constructed in one buffer, so the second has the address of the
    // first; a sole-owner shared_ptr is adopted rather than copied
    {
      alignas(IndexBasisRegistry) std::byte buffer[sizeof(IndexBasisRegistry)];
      auto registry_in_buffer = [&buffer] {
        return std::shared_ptr<const IndexBasisRegistry>(
            new (buffer) IndexBasisRegistry{},
            [](const IndexBasisRegistry* r) { r->~IndexBasisRegistry(); });
      };
      std::uint64_t dead_version = 0;
      {
        const Context dead(
            {.index_basis_registry_shared_ptr = registry_in_buffer()});
        REQUIRE(static_cast<const void*>(dead.index_basis_registry().get()) ==
                static_cast<const void*>(buffer));
        dead_version = dead.version();
      }
      const Context reborn(
          {.index_basis_registry_shared_ptr = registry_in_buffer()});
      REQUIRE(static_cast<const void*>(reborn.index_basis_registry().get()) ==
              static_cast<const void*>(buffer));
      CHECK(reborn.version() != dead_version);
    }

    // the registry and the single-particle basis change the version, the
    // other settings do not
    auto changes = [&ctx](auto&& set) {
      const auto before = ctx.version();
      set(ctx);
      return ctx.version() != before;
    };
    CHECK(!changes([](Context& c) { c.set(Vacuum::SingleProduct); }));
    CHECK(!changes([](Context& c) { c.set(AssertStrictBraKetSymmetry::No); }));
    CHECK(!changes([](Context& c) { c.set_first_dummy_index_ordinal(200); }));
    CHECK(!changes([](Context& c) { c.set(BraKetTypesetting::KetSub); }));
    CHECK(!changes([](Context& c) { c.set(BraKetSlotTypesetting::Naive); }));
    CHECK(!changes([](Context& c) {
      c.set(PrintOptions{.ordering = LexicographicCompareOptions{}});
    }));
    CHECK(!changes([](Context& c) { c.set(Symmetry::Symm); }));
    CHECK(!changes([](Context& c) { c.set(Hermiticity::Hermitian); }));
    CHECK(!changes([](Context& c) { c.set(ColumnSymmetry::Symm); }));
    CHECK(changes([](Context& c) { c.set(SPBasis::Spinfree); }));
    CHECK(changes([](Context& c) { c.set(IndexBasisRegistry{}); }));
    CHECK(changes(
        [](Context& c) { c.set(std::make_shared<IndexBasisRegistry>()); }));
    {
      const auto before = ctx.version();
      CHECK(changes([](Context& c) { c.set(SPBasis::Spinor); }));
      CHECK(changes([](Context& c) { c.set(SPBasis::Spinfree); }));
      // undoing a change restores the version
      CHECK(ctx.version() == before);
    }
    // the version is assigned when read, so the intermediate configurations
    // of a chain of setters get none: versions are assigned in increasing
    // order, and the intermediate one, read after the final one, is newer
    {
      const auto registry = std::make_shared<IndexBasisRegistry>();
      Context chained(ctx);
      chained.set(registry).set(SPBasis::Spinor);
      const auto final_version = chained.version();
      Context intermediate(ctx);
      intermediate.set(registry);
      CHECK(intermediate.version() > final_version);
    }

    // the copy was not affected
    CHECK(copy.version() == v0);
  }

  SECTION("snapshot") {
    const auto initial_ctx = get_default_context();
    auto restore = sequant::detail::make_scope_exit(
        [&initial_ctx] { set_default_context(initial_ctx); });

    // the process-wide context is the sole owner of its registry
    {
      IndexBasisRegistry registry;
      registry.add(L"q", 0b01);
      set_default_context(Context(initial_ctx)
                              .set(std::move(registry))
                              .set_first_dummy_index_ordinal(200));
    }
    const auto snapshot = get_default_context_snapshot();
    CHECK(snapshot == get_default_context());
    CHECK(snapshot.version() == current_context_version());
    const auto& registry = *snapshot.index_basis_registry();

    // replacing the process-wide context leaves the snapshot intact
    set_default_context(initial_ctx);
    CHECK(registry.contains(L"q"));
    CHECK(snapshot.first_dummy_index_ordinal() == 200);

    // the next snapshot sees every change of the process-wide context,
    // including one made on another thread
    std::thread([&initial_ctx] {
      set_default_context(
          Context(initial_ctx).set_first_dummy_index_ordinal(300));
    }).join();
    CHECK(get_default_context_snapshot().first_dummy_index_ordinal() == 300);
    std::thread([] { reset_default_context(); }).join();
    CHECK(get_default_context_snapshot().first_dummy_index_ordinal() ==
          Context{}.first_dummy_index_ordinal());
    set_default_context(initial_ctx);

    // on a thread with a scoped context, a copy of that context
    auto scoped = set_scoped_default_context(snapshot);
    CHECK(get_default_context_snapshot() == snapshot);
  }
}

TEST_CASE("scoped contexts", "[runtime]") {
  using namespace sequant;

  // the first dummy index ordinal serves as a marker of the context in effect
  auto marker = [] {
    return get_default_context().first_dummy_index_ordinal();
  };
  auto with_marker = [](std::size_t m) {
    return Context(get_default_context()).set_first_dummy_index_ordinal(m);
  };
  const auto initial_marker = marker();

  SECTION("each thread sees only its own") {
    constexpr int nthreads = 6;
    std::latch all_scoped(nthreads);
    std::atomic<int> mismatches = 0;
    std::vector<std::thread> threads;
    for (int t = 0; t != nthreads; ++t) {
      threads.emplace_back([&, t] {
        const std::size_t mine = initial_marker + 1 + t;
        {
          auto scoped = set_scoped_default_context(with_marker(mine));
          all_scoped.arrive_and_wait();
          for (int i = 0; i != 1000; ++i)
            if (marker() != mine) ++mismatches;
        }
        if (marker() != initial_marker) ++mismatches;
      });
    }
    for (auto& thread : threads) thread.join();
    CHECK(mismatches == 0);
    CHECK(marker() == initial_marker);
  }

  SECTION("parallel workers see the scoped context of the caller") {
    const std::size_t mine = initial_marker + 1;
    std::vector<int> items(64);
    std::atomic<int> mismatches = 0;
    {
      auto scoped = set_scoped_default_context(with_marker(mine));
      auto scoped_mbpt = mbpt::set_scoped_default_mbpt_context(
          mbpt::Context({.csv = mbpt::CSV::Yes}));
      sequant::for_each(items, [&](int&) {
        if (marker() != mine) ++mismatches;
        if (mbpt::get_default_mbpt_context().csv() != mbpt::CSV::Yes)
          ++mismatches;
        // nested scopes compose
        const std::size_t inner = initial_marker + 2;
        {
          auto scoped = set_scoped_default_context(with_marker(inner));
          if (marker() != inner) ++mismatches;
        }
        if (marker() != mine) ++mismatches;
      });
      CHECK(sequant::transform_reduce(items, 0, std::plus<int>{}, [&](int) {
              return marker() == mine ? 0 : 1;
            }) == 0);
      CHECK(marker() == mine);
    }
    CHECK(mismatches == 0);

    // after the scope workers see the process-wide context again
    sequant::for_each(items, [&](int&) {
      if (marker() != initial_marker) ++mismatches;
    });
    CHECK(mismatches == 0);
  }

  // unlike for_each, whose execution-policy backend may run serially,
  // parallel_do always runs on num_threads() distinct threads
  SECTION("parallel_do threads see the scoped context of the caller") {
    const std::size_t mine = initial_marker + 1;
    std::atomic<int> invocations = 0;
    std::atomic<int> mismatches = 0;
    std::mutex thread_ids_mtx;
    std::set<std::thread::id> thread_ids;
    {
      auto scoped = set_scoped_default_context(with_marker(mine));
      sequant::parallel_do([&](int) {
        ++invocations;
        {
          std::scoped_lock lock(thread_ids_mtx);
          thread_ids.insert(std::this_thread::get_id());
        }
        if (marker() != mine) ++mismatches;
      });
      CHECK(marker() == mine);
    }
    CHECK(invocations == num_threads());
    CHECK(mismatches == 0);
    if (num_threads() > 1) CHECK(thread_ids.size() > 1);

    sequant::parallel_do([&](int) {
      if (marker() != initial_marker) ++mismatches;
    });
    CHECK(mismatches == 0);
  }

  SECTION("a modification applies to the context of every statistics") {
    const auto vacuum = [](Statistics s) {
      return get_default_context(s).vacuum();
    };
    auto scoped =
        set_scoped_default_context(container::map<Statistics, Context>{
            {Statistics::Arbitrary,
             Context(get_default_context()).set(Vacuum::SingleProduct)},
            {Statistics::FermiDirac,
             Context(get_default_context()).set(Vacuum::Physical)}});
    const std::size_t mine = initial_marker + 1;
    auto marker_of = [](Statistics s) {
      return get_default_context(s).first_dummy_index_ordinal();
    };
    {
      auto modified = set_scoped_modified_default_context(
          [mine](Context& ctx) { ctx.set_first_dummy_index_ordinal(mine); });
      CHECK(vacuum(Statistics::Arbitrary) == Vacuum::SingleProduct);
      CHECK(vacuum(Statistics::FermiDirac) == Vacuum::Physical);
      CHECK(vacuum(Statistics::BoseEinstein) == Vacuum::SingleProduct);
      CHECK(marker_of(Statistics::Arbitrary) == mine);
      CHECK(marker_of(Statistics::FermiDirac) == mine);
    }
    CHECK(vacuum(Statistics::Arbitrary) == Vacuum::SingleProduct);
    CHECK(vacuum(Statistics::FermiDirac) == Vacuum::Physical);
    CHECK(marker_of(Statistics::Arbitrary) == initial_marker);
    CHECK(marker_of(Statistics::FermiDirac) == initial_marker);
  }

  SECTION("mbpt::load refuses to run under a scoped context") {
    const auto process_wide_version = current_context_version();
    {
      auto scoped = set_scoped_default_context(with_marker(initial_marker + 1));
      CHECK_THROWS_AS(mbpt::load(), Exception);
    }
    CHECK(current_context_version() == process_wide_version);
    CHECK(marker() == initial_marker);
  }

  SECTION("equal modifications of a context scope equal versions") {
    const auto other_spbasis =
        get_default_context().spbasis() == SPBasis::Spinor ? SPBasis::Spinfree
                                                           : SPBasis::Spinor;
    auto scoped_version = [](auto&& modify) {
      auto modified = set_scoped_modified_default_context(modify);
      return current_context_version();
    };
    auto set_spbasis = [other_spbasis](Context& ctx) {
      ctx.set(other_spbasis);
    };
    const auto v1 = scoped_version(set_spbasis);
    CHECK(v1 != current_context_version());
    CHECK(scoped_version(set_spbasis) == v1);
    // a (copied) registry of its own is a configuration of its own
    CHECK(scoped_version([other_spbasis](Context& ctx) {
            ctx.set(other_spbasis)
                .set(IndexBasisRegistry(*ctx.index_basis_registry()));
          }) != v1);
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
        auto modified = set_scoped_modified_default_context(
            [](Context& ctx) { ctx.set(SPBasis::Spinfree); });
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

// a pin scopes a copy of the contexts in effect, with their versions, that a
// change of the process-wide contexts made meanwhile does not reach; under a
// scoped context it installs nothing
TEST_CASE("pin_default_contexts", "[runtime]") {
  using namespace sequant;
  using Contexts = container::map<Statistics, Context>;
  REQUIRE(detail::implicit_context_overlay<Contexts>() == nullptr);
  const auto version = current_contexts_version();
  const auto vacuum = get_default_context().vacuum();
  const auto other_vacuum =
      vacuum == Vacuum::Physical ? Vacuum::SingleProduct : Vacuum::Physical;
  {
    auto pin = pin_default_contexts();
    const auto* overlay = detail::implicit_context_overlay<Contexts>();
    REQUIRE(overlay != nullptr);
    CHECK(current_contexts_version() == version);
    auto changed = get_default_context();
    changed.set(other_vacuum);
    set_default_context(changed);
    CHECK(get_default_context().vacuum() == vacuum);
    auto inner = pin_default_contexts();
    CHECK(detail::implicit_context_overlay<Contexts>() == overlay);
  }
  CHECK(get_default_context().vacuum() == other_vacuum);
  auto restored = get_default_context();
  restored.set(vacuum);
  set_default_context(restored);
  CHECK(current_contexts_version() == version);
}
