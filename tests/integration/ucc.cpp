#include <SeQuant/version.hpp>

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/conversion.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/core/utility/timer.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>

#include <cstddef>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Unitary CC (UCC) equation derivation: the srcc.cpp analogue for the unitary
// ansatz, covering both H̄ expansions: the standard BCH commutator series and
// the Bernoulli expansion of 10.1063/1.5030344.
//
// CC::t() yields the whole equation set in one derivation: element 0 is the
// energy <0|H̄|0>, element R>0 the residual <R|H̄|0>. Term counts are pinned
// below.
//
// Usage: ucc [N] [bch|bernoulli] [RANK] [print]
//        N     cluster/excitation rank of T                       (default 2)
//        RANK  commutator truncation rank of H̄                    (default 2)

using namespace sequant;
using namespace sequant::mbpt;

namespace {

#define runtime_assert(tf)                                         \
  if (!(tf)) {                                                     \
    std::ostringstream oss;                                        \
    oss << "failed assert at line " << __LINE__ << " in function " \
        << __func__;                                               \
    throw std::runtime_error(oss.str().c_str());                   \
  }

TimerPool<32> tpool;

using Hbar = CC::HbarExpansion;

const std::map<std::string, Hbar> str2expansion = {
    {"bch", Hbar::BCH}, {"bernoulli", Hbar::Bernoulli}};

/// pinned term count of one equation
struct TermCounts {
  Hbar expansion;
  std::size_t n;       ///< cluster rank
  std::size_t rank;    ///< H̄ commutator truncation rank
  std::size_t r;       ///< projection manifold rank; 0 = energy
  std::size_t nterms;  ///< expected number of terms
};

// Regression pins, not independent references
const std::vector<TermCounts> pins = {
    // clang-format off
    // expansion,        N, rank, R, terms
    {Hbar::BCH,          2, 2,    0,   20}, {Hbar::BCH,       2, 2, 1,   44}, {Hbar::BCH,       2, 2, 2,   42},
    {Hbar::BCH,          2, 3,    0,   74}, {Hbar::BCH,       2, 3, 1,  219}, {Hbar::BCH,       2, 3, 2,  267},
    {Hbar::BCH,          2, 4,    0,  307}, {Hbar::BCH,       2, 4, 1, 1100}, {Hbar::BCH,       2, 4, 2, 1433},
    {Hbar::Bernoulli,    2, 2,    0,    6}, {Hbar::Bernoulli, 2, 2, 1,   32}, {Hbar::Bernoulli, 2, 2, 2,   38},
    {Hbar::Bernoulli,    2, 3,    0,   46}, {Hbar::Bernoulli, 2, 3, 1,  141}, {Hbar::Bernoulli, 2, 3, 2,  191},
    {Hbar::Bernoulli,    2, 4,    0,  203}, {Hbar::Bernoulli, 2, 4, 1,  722}, {Hbar::Bernoulli, 2, 4, 2, 1044},
    // clang-format on
};

void check(Hbar expansion, std::size_t n, std::size_t rank, std::size_t r,
           std::size_t nterms) {
  for (const auto& p : pins)
    if (expansion == p.expansion && n == p.n && rank == p.rank && r == p.r) {
      if (nterms != p.nterms)
        std::wcout << "MISMATCH: expected " << p.nterms << " terms, got "
                   << nterms << std::endl;
      runtime_assert(nterms == p.nterms);
      return;
    }
}

}  // namespace

int main(int argc, char* argv[]) {
  std::wcout.precision(std::numeric_limits<double>::max_digits10);
  sequant::set_locale();

  const std::size_t N = argc > 1 ? string_to<std::size_t>(argv[1]) : 2;
  const std::string expansion_str = argc > 2 ? argv[2] : "bch";
  const auto expansion = str2expansion.at(expansion_str);
  const std::size_t RANK = argc > 3 ? string_to<std::size_t>(argv[3]) : 2;
  const bool print = argc > 4 && std::string(argv[4]) == "print";

  sequant::detail::OpIdRegistrar op_id_registrar;
  set_default_context({.index_space_registry_shared_ptr = make_sr_spaces(),
                       .vacuum = Vacuum::SingleProduct,
                       .metric = IndexSpaceMetric::Unit,
                       .spbasis = SPBasis::Spinor,
                       .first_dummy_index_ordinal = 100});
  TensorCanonicalizer::set_cardinal_tensor_labels(cardinal_tensor_labels());
  set_default_mbpt_context(
      {.csv = mbpt::CSV::No, .op_registry_ptr = make_legacy_registry()});

  std::cout << "SeQuant revision: " << sequant::git_revision() << "\n";
  std::cout << "Number of threads: " << sequant::num_threads() << "\n";

  const CC cc(N, {.ansatz = CC::Ansatz::U,
                  .hbar_comm_rank = RANK,
                  .hbar_expansion = expansion});

  tpool.clear();
  tpool.start(0);
  const auto eqvec = cc.t();
  tpool.stop(0);

  std::wcout << "UCC equations [rank=" << N
             << ",expansion=" << sequant::toUtf16(expansion_str)
             << ",hbar_comm_rank=" << RANK << "] computed in " << tpool.read(0)
             << " seconds" << std::endl;

  for (std::size_t R = 0; R < eqvec.size(); ++R) {
    std::wcout << (R == 0 ? "E" : "R") << (R == 0 ? L"" : std::to_wstring(R))
               << "(expU" << N << ") has " << eqvec[R]->size()
               << " terms:" << std::endl;
    if (print) std::wcout << to_latex_align(eqvec[R], 20, 1) << std::endl;
    check(expansion, N, RANK, R, eqvec[R]->size());
  }

  return 0;
}
