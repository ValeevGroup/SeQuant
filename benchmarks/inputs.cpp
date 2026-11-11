#include "inputs.hpp"

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/abstract_tensor.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/index_basis_registry.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <array>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <utility>
#include <vector>

namespace sequant::bench {

const ExprPtr& cc_residual(std::size_t cc_rank, std::size_t projection) {
  static std::map<std::size_t, std::vector<ExprPtr>> residuals;

  auto it = residuals.find(cc_rank);
  if (it == residuals.end()) {
    it = residuals.emplace(cc_rank, mbpt::CC(cc_rank).t()).first;
  }

  return it->second.at(projection);
}

const ExprPtr& unsimplified_cc_residual(std::size_t cc_rank,
                                        std::size_t projection) {
  static std::map<std::pair<std::size_t, std::size_t>, ExprPtr> cache;

  const auto key = std::make_pair(cc_rank, projection);
  auto it = cache.find(key);
  if (it == cache.end()) {
    // Mirrors the steps closed_shell_spintrace takes before simplifying
    ExprPtr expr = cc_residual(cc_rank, projection)->clone();
    expand(expr);
    if (has_tensor(expr, reserved::antisymm_label())) {
      expr = mbpt::symmetrize_expr(expr);
    }
    expr = mbpt::expand_antisymm(expr);

    it = cache.emplace(key, std::move(expr)).first;
  }

  return it->second;
}

const ExprPtr& spinfree_cc_residual(std::size_t cc_rank,
                                    std::size_t projection) {
  static std::map<std::pair<std::size_t, std::size_t>, ExprPtr> cache;

  const auto key = std::make_pair(cc_rank, projection);
  auto it = cache.find(key);
  if (it == cache.end()) {
    const ExprPtr& residual = cc_residual(cc_rank, projection);
    it = cache
             .emplace(key,
                      mbpt::closed_shell_spintrace(
                          residual->clone(), external_indices(residual->at(0))))
             .first;
  }

  return it->second;
}

namespace {

constexpr std::array<const char*, 12> nevpt2_files = {
    "nevpt2_en0.inp",
    "nevpt2_en.inp",
    "nevpt2_res1_i1.inp",
    "nevpt2_res1_s0.inp",
    "nevpt2_res1_s1.inp",
    "nevpt2_res2_i2.inp",
    "nevpt2_res2_p0.inp",
    "nevpt2_res2_p1.inp",
    "nevpt2_res2_p2.inp",
    "nevpt2_res2_s1.inp",
    "nevpt2_res2_s1_singles.inp",
    "nevpt2_res2_s2.inp"};

}  // namespace

struct Nevpt2::Impl {
  using Resetter = decltype(set_scoped_default_context(Context{}));

  Resetter resetter;
  std::vector<ExprPtr> residuals;

  static Context make_context() {
    // Same set of spaces as the external interface uses for these inputs
    IndexBasisRegistry registry;
    const std::array<std::pair<const wchar_t*, int>, 3> spaces = {
        {{L"a", 1000}, {L"u", 5}, {L"i", 80}}};
    for (std::size_t i = 0; i < spaces.size(); ++i) {
      registry.add(spaces[i].first, IndexSpace::Type(1 << i),
                   {.dimension = spaces[i].second,
                    .qns = mbpt::Spin::any,
                    .field = Field::Real});
    }
    mbpt::add_fermi_spin(registry);

    Context ctx({.index_basis_registry = std::move(registry),
                 .vacuum = Vacuum::SingleProduct});
    ctx.set(CanonicalizeOptions{.method = CanonicalizationMethod::Complete});
    ctx.set_cardinal_tensor_labels(mbpt::cardinal_tensor_labels());
    return ctx;
  }

  Impl() : resetter(set_scoped_default_context(make_context())) {
    const std::filesystem::path dir(SEQUANT_BENCH_NEVPT2_DIR);
    const io::serialization::DeserializationOptions options{
        .def_perm_symm = Symmetry::Antisymm,
        .def_braket_symm = BraKetSymmetry::Nonsymm,
        .def_col_symm = ColumnSymmetry::Symm};

    for (const char* file : nevpt2_files) {
      std::ifstream in(dir / file);
      if (!in) {
        throw Exception(std::string("Unable to read ") + (dir / file).string());
      }
      const std::string contents(std::istreambuf_iterator<char>(in), {});
      residuals.push_back(
          io::serialization::from_string<ResultExpr>(contents, options)
              .expression());
    }
  }
};

Nevpt2::Nevpt2() : impl_(std::make_unique<Impl>()) {}

Nevpt2::~Nevpt2() = default;

std::size_t Nevpt2::size() { return nevpt2_files.size(); }

const std::string& Nevpt2::name(std::size_t i) {
  static const std::vector<std::string> names(nevpt2_files.begin(),
                                              nevpt2_files.end());
  return names.at(i);
}

const ExprPtr& Nevpt2::residual(std::size_t i) const {
  return impl_->residuals.at(i);
}

void cc_residual_args(benchmark::Benchmark* benchmark) {
  for (std::size_t n = 1; n <= max_cc_rank; ++n) {
    for (std::size_t p = 1; p <= n; ++p) {
      benchmark->Args({static_cast<int64_t>(n), static_cast<int64_t>(p)});
    }
  }
}

}  // namespace sequant::bench
