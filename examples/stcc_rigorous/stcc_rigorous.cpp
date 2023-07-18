#include <SeQuant/core/math.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/timer.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <clocale>

#include <Eigen/Eigenvalues>

using namespace sequant;

ExprPtr symmetrize_expr(
    ExprPtr& expr,
    const container::vector<container::vector<Index>>& ext_index_groups = {{}});

#define runtime_assert(tf)                                         \
  if (!(tf)) {                                                     \
    std::ostringstream oss;                                        \
    oss << "failed assert at line " << __LINE__                    \
        << " in closed-shell spin-traced coupled cluster example"; \
    throw std::runtime_error(oss.str().c_str());                   \
  }

int main(int argc, char* argv[]) {
  std::setlocale(LC_ALL, "en_US.UTF-8");
  std::wcout.precision(std::numeric_limits<double>::max_digits10);
  std::wcerr.precision(std::numeric_limits<double>::max_digits10);
  std::wcout.sync_with_stdio(false);
  std::wcerr.sync_with_stdio(false);
  std::wcout.imbue(std::locale("en_US.UTF-8"));
  std::wcerr.imbue(std::locale("en_US.UTF-8"));
  std::wcout.sync_with_stdio(true);
  std::wcerr.sync_with_stdio(true);

  sequant::set_default_context(
      Context(Vacuum::SingleProduct, IndexSpaceMetric::Unit,
              BraKetSymmetry::conjugate, SPBasis::spinorbital));
  mbpt::set_default_convention();
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

#ifndef NDEBUG
  const size_t DEFAULT_NMAX = 2;
#else
  const size_t DEFAULT_NMAX = 3;
#endif
  const size_t NMAX = argc > 1 ? std::atoi(argv[1]) : DEFAULT_NMAX;

  /// Make external index
  auto ext_idx_list = [](const int i_max) {
    container::vector<container::vector<Index>> ext_idx_list;

    for (size_t i = 1; i <= i_max; ++i) {
      auto label = std::to_wstring(i);
      auto occ_i = Index::make_label_index(
          IndexSpace::instance(IndexSpace::active_occupied), label);
      auto virt_i = Index::make_label_index(
          IndexSpace::instance(IndexSpace::active_unoccupied), label);
      container::vector<Index> pair = {occ_i, virt_i};
      ext_idx_list.push_back(pair);
    }
    return ext_idx_list;
  };

  // Spin-orbital coupled cluster
  auto cc_r = sequant::mbpt::sr::cceqs{NMAX}.t();
  for (auto i = 1; i < cc_r.size(); ++i) {
    std::cout << "Spin-orbital CC R" << i << " size: " << cc_r[i]->size()
              << "\n";
  }

  //
  // Closed-shell spintrace (fast)
  //
  std::cout << "\nClosed-shell coupled cluster:\n";
  std::vector<ExprPtr> cc_st_r(cc_r.size());
  for (auto i = 1; i < cc_r.size(); ++i) {
    const auto tstart = std::chrono::high_resolution_clock::now();
    auto ext_idx = ext_idx_list(i);
    cc_st_r[i] = sequant::spintrace(cc_r[i], ext_idx_list(i));
    canonicalize(cc_st_r[i]);

    // Remove S operator
    for (auto& term : *cc_st_r[i]) {
      if (term->is<Product>()) term = remove_tensor(term->as<Product>(), L"S");
    }

    // Biorthogonal transformation
    cc_st_r[i] = biorthogonal_transform(cc_st_r[i], i, ext_idx);

    // The symmetrizer operator is required for canonicalizer to give the
    // correct result
    if (i != 1) cc_st_r[i] = symmetrize_expr(cc_st_r[i], ext_idx);
    simplify(cc_st_r[i]);

    // Remove S operator
    for (auto& term : *cc_st_r[i]) {
      if (term->is<Product>()) term = remove_tensor(term->as<Product>(), L"S");
    }

    auto tstop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_elapsed = tstop - tstart;
    printf("CC R%d size: %lu time: %5.3f sec.\n", i, cc_st_r[i]->size(),
           time_elapsed.count());
  }

  if (NMAX == 2) {
    runtime_assert(cc_st_r.size() == 3)
        runtime_assert(cc_st_r.at(1)->size() == 26)  // T1
        runtime_assert(cc_st_r.at(2)->size() == 55)  // T2
  } else if (NMAX == 3) {
    runtime_assert(cc_st_r.size() == 4)
        runtime_assert(cc_st_r.at(1)->size() == 30)   // T1
        runtime_assert(cc_st_r.at(2)->size() == 73)   // T2
        runtime_assert(cc_st_r.at(3)->size() == 490)  // T3
  }
}

// Generate S operator from external index list
ExprPtr symmetrize_expr(
    ExprPtr& expr,
    const container::vector<container::vector<Index>>& ext_index_groups) {
  container::vector<Index> bra_list, ket_list;
  for (auto&& idx_group : ext_index_groups) {
    bra_list.push_back(*idx_group.begin());
    ket_list.push_back(*(idx_group.begin() + 1));
  }

  assert(bra_list.size() == ket_list.size());
  auto S = Tensor(L"S", bra_list, ket_list, Symmetry::nonsymm);
  return ex<Tensor>(S) * expr;
}
