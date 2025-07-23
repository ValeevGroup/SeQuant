#include <benchmark/benchmark.h>

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network_v2.hpp>

#include <range/v3/all.hpp>

#include <cassert>
#include <format>
#include <random>

using namespace sequant;

ProductPtr create_random_network(const std::size_t testcase,
                                 const std::size_t num_indices) {
  int n = 1;
  switch (testcase) {
    case 0:
      n = 1;
      break;
    case 1:
      n = num_indices / 2;
      break;
    case 2:
      n = num_indices;
      break;
    case 3:
      n = 1;
      break;
    default:
      abort();
  }

  if (n == 0 || n > num_indices) {
    return nullptr;
  }

  container::vector<Index> covariant_indices;
  for (auto i = 0; i != num_indices; ++i) {
    covariant_indices.emplace_back(std::format(L"i_{}", i));
  }

  auto contravariant_indices = covariant_indices;

  std::random_device rd;

  // Randomize connectivity by randomizing index sets
  std::shuffle(covariant_indices.begin(), covariant_indices.end(),
               std::mt19937{rd()});
  std::shuffle(contravariant_indices.begin(), contravariant_indices.end(),
               std::mt19937{rd()});

  auto utensors =
      covariant_indices | ranges::views::chunk(num_indices / n) |
      ranges::views::transform([&](const auto& idxs) {
        return ex<Tensor>(
            L"u", bra(idxs), ket{},
            (testcase == 3 ? Symmetry::nonsymm
                           : ((n == 1) ? Symmetry::symm : Symmetry::nonsymm)));
      }) |
      ranges::to_vector;

  assert(utensors.size() == static_cast<std::size_t>(n));

  auto dtensors =
      contravariant_indices | ranges::views::chunk(num_indices) |
      ranges::views::transform([&](const auto& idxs) {
        return ex<Tensor>(L"d", bra{}, ket(idxs), Symmetry::nonsymm);
      }) |
      ranges::to_vector;
  assert(dtensors.size() == 1);

  ExprPtr expr;
  for (int g = 0; g < n; ++g) {
    if (g == 0)
      expr = utensors[0] * dtensors[0];
    else
      expr = expr * utensors[g];
  }

  return std::dynamic_pointer_cast<Product>(expr);
}

static void random_tensor_network(benchmark::State& state) {
  const std::size_t testcase = state.range(0);
  const std::size_t num_indices = state.range(1);

  auto ctx_resetter = set_scoped_default_context(
      Context(get_default_context())
          .set_first_dummy_index_ordinal(num_indices + 1));

  const ProductPtr& prod = create_random_network(testcase, num_indices);

  if (!prod) {
    state.SkipWithMessage("Invalid");
  }

  for (auto _ : state) {
    // Need to clone in order to avoid mutating original expression
    TensorNetworkV2 tn(prod->clone()->as<Product>().factors());

    const bool fast = false;
    ExprPtr expr =
        tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), fast);

    // Prevent the compiler from optimizing the canonicalization away
    benchmark::DoNotOptimize(expr);
  }
}

// 10.1016/j.cpc.2018.02.014
// - testcase=0,2 are "equivalent" and correspond to the "frustrated"
//   case in Section 5.3 of the reference
// - testcase=1 corresponds to the "frustrated" case in Section 5.4 of
//   the reference
// - testcase=3 corresponds to the "No symmetry dummy"
//   case in Section 5.1 of the reference
BENCHMARK(random_tensor_network)
    ->ArgsProduct({benchmark::CreateDenseRange(0, 3, 1),
                   benchmark::CreateRange(1, 256, 2)});
