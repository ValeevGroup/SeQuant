//
// Created by Ajay Melekamburath on 2/3/25.
//
#include <SeQuant/version.hpp>

#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/timer.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>

using namespace sequant;
using namespace sequant::mbpt;

namespace {
#define runtime_assert(tf)                                   \
  if (!(tf)) {                                               \
    std::ostringstream oss;                                  \
    oss << "failed assert at line " << __LINE__              \
        << " in equation-of-motion coupled cluster example"; \
    throw std::runtime_error(oss.str().c_str());             \
  }

TimerPool<32> timer_pool;

std::pair<size_t, size_t> parse_excitation_manifold(std::string& str) {
  std::pair<size_t, size_t> result;

  ranges::transform(str, str.begin(), ::tolower);
  const auto h_pos = str.find('h');
  const auto p_pos = str.find('p');

  if (h_pos == std::string::npos && p_pos == std::string::npos) {
    throw std::runtime_error(
        "Invalid excitation manifold string: must contain 'h' or 'p'");
  }

  if (h_pos != std::string::npos && p_pos == std::string::npos) {
    // holes only
    result.first = std::stoi(str.substr(0, h_pos));
    result.second = 0;
  } else if (p_pos != std::string::npos && h_pos == std::string::npos) {
    // particles only
    result.first = 0;
    result.second = std::stoi(str.substr(0, p_pos));
  } else {
    // hp/ph cases
    if (h_pos < p_pos) {
      result.first = std::stoi(str.substr(0, h_pos));
      result.second = std::stoi(str.substr(h_pos + 1, p_pos - h_pos - 1));
    } else {
      result.first = std::stoi(str.substr(p_pos + 1, h_pos - p_pos - 1));
      result.second = std::stoi(str.substr(0, p_pos));
    }
  }

  if (result.first == 0 && result.second == 0)
    throw std::runtime_error(
        "Invalid excitation manifold: both particle and hole ranks cannot be "
        "zero");

  return result;
}

enum class EqnType { left, right };

inline const container::map<std::string, EqnType> str2type = {
    {"L", EqnType::left}, {"R", EqnType::right}};

inline const container::map<EqnType, std::wstring> type2wstr = {
    {EqnType::left, L"L"}, {EqnType::right, L"R"}};

class compute_eomcc {
  size_t N, np, nh;
  std::string manifold;
  EqnType type;

 public:
  compute_eomcc(size_t n, const std::string& exc_manifold,
                EqnType t = EqnType::right)
      : N(n), manifold(exc_manifold), type(t) {
    std::tie(nh, np) = parse_excitation_manifold(manifold);
  }

  void operator()(bool print) {
    SEQUANT_ASSERT(get_default_context().spbasis() == SPBasis::Spinor);
    timer_pool.start(N);
    std::vector<ExprPtr> eqvec;
    switch (type) {
      case EqnType::right:
        eqvec = CC{N}.eom_r(nₚ(np), nₕ(nh));
        break;
      case EqnType::left:
        eqvec = CC{N}.eom_l(nₚ(np), nₕ(nh));
        break;
    }
    timer_pool.stop(N);

    std::wcout << std::boolalpha
               << "EOM-CC Equations [type=" << type2wstr.at(type)
               << ", CC rank=" << N
               << ", manifold=" << sequant::toUtf16(manifold) << "]"
               << " computed in " << timer_pool.read(N) << " s\n";
    for (auto i = 0; i < eqvec.size(); i++) {
      if (eqvec[i] == nullptr) continue;
      std::wcout << "R[" << i << "] has " << eqvec[i].size() << " terms\n";
      if (print) std::wcout << to_latex_align(eqvec[i], 20, 1) << "\n";
      if (N == 2 && type == EqnType::right) {
        if (np == 2 && nh == 2) {  // EOM-CCSD(2h2p)
          if (i == 1) runtime_assert(eqvec[i].size() == 21);
          if (i == 2) runtime_assert(eqvec[i].size() == 53);
        }
        if (np == 2 && nh == 1) {  // EA-EOM-CCSD(1h2p)
          if (i == 0) runtime_assert(eqvec[i].size() == 9);
          if (i == 1) runtime_assert(eqvec[i].size() == 32);
        }
        if (np == 1 && nh == 2) {  // IP-EOM-CCSD(2h1p)
          if (i == 0) runtime_assert(eqvec[i].size() == 9);
          if (i == 1) runtime_assert(eqvec[i].size() == 32);
        }
        if (np == 1 && nh == 3) {  // DIP-EOM-CCSD(3h1p)
          if (i == 0) runtime_assert(eqvec[i].size() == 13);
          if (i == 1) runtime_assert(eqvec[i].size() == 34);
        }
        if (np == 3 && nh == 1) {  // DEA-EOM-CCSD(1h3p)
          if (i == 0) runtime_assert(eqvec[i].size() == 13);
          if (i == 1) runtime_assert(eqvec[i].size() == 34);
        }
        if (np == 2 && nh == 4) {  // DIP-EOM-CCSD(4h2p)
          if (i == 0) runtime_assert(eqvec[i].size() == 14);
          if (i == 1) runtime_assert(eqvec[i].size() == 40);
          if (i == 2) runtime_assert(eqvec[i].size() == 65);
        }
      }
      if (N == 3 && type == EqnType::right) {
        if (np == 3 && nh == 3) {  // EOM-CCSDT(3h3p)
          if (i == 1) runtime_assert(eqvec[i].size() == 22);
          if (i == 2) runtime_assert(eqvec[i].size() == 62);
          if (i == 3) runtime_assert(eqvec[i].size() == 99);
        }
      }
    }
  }
};  // class compute_eomcc

}  // namespace

int main(int argc, char* argv[]) {
  std::wcout.precision(std::numeric_limits<double>::max_digits10);
  std::wcerr.precision(std::numeric_limits<double>::max_digits10);
  sequant::set_locale();

  std::cout << "SeQuant revision: " << sequant::git_revision() << "\n";
  std::cout << "Number of threads: " << sequant::num_threads() << "\n\n";

  sequant::detail::OpIdRegistrar op_id_registrar;
  sequant::set_default_context(
      sequant::Context({.index_space_registry_shared_ptr = make_min_sr_spaces(),
                        .vacuum = Vacuum::SingleProduct}));
  mbpt::set_default_mbpt_context(
      {.op_registry_ptr = mbpt::make_minimal_registry()});
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  // change to true to print stats
  Logger::instance().wick_stats = false;

  if (argc > 1) {
    // Usage: eomcc <N> <manifold> [R|L] [print]
    const size_t N = std::stoi(argv[1]);
    SEQUANT_ASSERT(N > 0 && "Invalid N");
    std::string exc_manifold =
        argc > 2 ? argv[2]
                 : (std::to_string(N) + "h" + std::to_string(N) + "p");
    SEQUANT_ASSERT(!exc_manifold.empty() && "Invalid excitation manifold");
    const std::string eqn_type = argc > 3 ? argv[3] : "R";
    const std::string print_str = argc > 4 ? argv[4] : "noprint";
    const bool print = print_str == "print";

    compute_eomcc{N, exc_manifold, str2type.at(eqn_type)}(print);
  } else {
    // default: EE-EOM-CCSD right
    compute_eomcc{2, "2h2p", EqnType::right}(false);
  }
}
