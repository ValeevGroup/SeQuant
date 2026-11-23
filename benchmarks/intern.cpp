#include <benchmark/benchmark.h>

#include <SeQuant/core/basis.hpp>
#include <SeQuant/core/space.hpp>

#include <cstddef>
#include <string>
#include <vector>

using namespace sequant;

namespace {

// @p n distinct bases of distinct spaces, as values that are not themselves
// interned (so intern() cannot short-circuit on their address)
std::vector<IndexBasis> make_bases(std::size_t n) {
  std::vector<IndexBasis> bases;
  bases.reserve(n);
  for (std::size_t k = 0; k != n; ++k)
    bases.emplace_back(
        IndexSpace(L"s" + std::to_wstring(k), TypeAttr(1u << (k % 30))));
  return bases;
}

}  // namespace

// Every thread cycles through the same working set of state.range(0) bases;
// the time per iteration is that of one intern() call
static void intern_basis(benchmark::State &state) {
  const auto bases = make_bases(state.range(0));
  for (const auto &b : bases) intern(b);
  std::size_t i = state.thread_index() % bases.size();
  for (auto _ : state) {
    benchmark::DoNotOptimize(&intern(bases[i]));
    if (++i == bases.size()) i = 0;
  }
  state.SetItemsProcessed(state.iterations());
}
BENCHMARK(intern_basis)
    ->ArgName("bases")
    ->Arg(1)
    ->Arg(4)
    ->Arg(8)
    ->Arg(16)
    ->ThreadRange(1, 16)
    ->UseRealTime();

// As intern_basis, but with the interned bases themselves as arguments, as when
// an Index is built from another's basis()
static void intern_interned_basis(benchmark::State &state) {
  std::vector<const IndexBasis *> bases;
  for (const auto &b : make_bases(state.range(0))) bases.push_back(&intern(b));
  std::size_t i = state.thread_index() % bases.size();
  for (auto _ : state) {
    benchmark::DoNotOptimize(&intern(*bases[i]));
    if (++i == bases.size()) i = 0;
  }
  state.SetItemsProcessed(state.iterations());
}
BENCHMARK(intern_interned_basis)
    ->ArgName("bases")
    ->Arg(1)
    ->Arg(4)
    ->Arg(8)
    ->ThreadRange(1, 16)
    ->UseRealTime();
