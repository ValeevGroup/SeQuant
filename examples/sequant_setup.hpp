#include <clocale>
#include <iostream>

#include <boost/numeric/interval.hpp>
#include <boost/math/special_functions/factorials.hpp>

#include <SeQuant/core/timer.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>
#include <SeQuant/domain/mbpt/sr/sr.hpp>
//#include <SeQuant/domain/evaluate/factorizer.hpp>

// TODO
  using namespace sequant;
  using namespace sequant::mbpt::sr::so;

    using ispace_pair = std::pair<sequant::IndexSpace::Type, size_t>;
    using ispace_map = sequant::container::map<ispace_pair::first_type,
    ispace_pair::second_type>;
// 
//using namespace sequant::mbpt::sr::so::csv;

// print_expr adds extra parentheses
// to make product of products more visible
// however, the scalars of a product are dropped.
void print_expr(const ExprPtr& expr_ptr){
  if (expr_ptr->is<sequant::Product>()){
    // Note: the scalar factor is not printed.
    for (auto &&fac: expr_ptr->as<sequant::Product>()){
      if (fac->is<sequant::Product>())
        std::wcout << LR"(\left()" << fac->to_latex() << LR"(\right))";
      else std::wcout << fac->to_latex();
    }
  } else if (expr_ptr->is<sequant::Sum>()) {
    for (auto &&sumand: expr_ptr->as<sequant::Sum>()) {
      print_expr(sumand);
      std::wcout << " + ";
    }
  } else {
    std::wcout << expr_ptr->to_latex(); }
}

namespace {

/// computes VEV for A(P)*H*T(N)^K using excitation level screening (unless @c
/// screen is set) + computes only canonical (with T ranks increasing) terms

class screened_vac_av {
 private:
  size_t K;

 public:

  screened_vac_av(size_t k) : K(k) {}

  ExprPtr operator()(const ExprPtr& expr,
                     std::initializer_list<std::pair<int, int>> op_connections,
                     bool screen = true, bool use_topology = true,
                     bool canonical_only = true) {
    if (!screen) return sequant::mbpt::sr::so::vac_av(expr, op_connections, use_topology);

    ExprPtr input = expr;
    // expand, if possible
    if (input->is<Product>()) {
      expand(input);
      if (input->is<Product>()) input = ex<Sum>(ExprPtrList{input});
    }
    assert(input->is<Sum>());
    auto input_sum = input->as<Sum>();

    // this will collect all canonical nonzero terms
    SumPtr screened_input = std::make_shared<Sum>();
    for (auto&& term : input_sum.summands()) {
      assert(term->is<Product>());
      auto& term_prod = term->as<Product>();
      assert(term_prod.factors().size() == 4 + 2 * K);

      // locate projector
      assert(term_prod.factor(0)->is<Tensor>());
      assert(term_prod.factor(0)->as<Tensor>().label() == L"A");
      const int P = term_prod.factor(0)->as<Tensor>().rank();

      // locate Hamiltonian
      assert(term_prod.factor(2)->is<Tensor>());
      auto hlabel = term_prod.factor(2)->as<Tensor>().label();
      assert(hlabel == L"f" || hlabel == L"g");
      const int R = term_prod.factor(2)->as<Tensor>().rank();
      const int max_exlev_R = R - K;  // at least K lines must point down

      auto exlev = -P;

      bool canonical = true;
      // number of possible permutations within same-rank partitions of T
      // number of possible permutations other than these
      // degeneracy = M1! M2! .. where M1, M2 ... are sizes of each partition
      double degeneracy =
          canonical_only ? boost::math::factorial<double>(K) : 1;
      int total_T_rank = 0;
      int prev_rank = 0;
      int current_partition_size =
          1;  // size of current same-rank partition of T
      for (size_t k = 0; k != K && canonical; ++k) {
        auto p = 4 + k * 2;
        assert(term_prod.factor(p)->is<Tensor>());
        assert(term_prod.factor(p)->as<Tensor>().label() == L"t");
        const auto current_rank = term_prod.factor(p)->as<Tensor>().rank();
        exlev += current_rank;
        total_T_rank += current_rank;
        // screen out the noncanonical terms, if needed
        if (canonical_only) {
          if (current_rank < prev_rank)  // if T ranks are not increasing, omit
            canonical = false;
          else {  // else keep track of degeneracy
            assert(current_rank != 0);
            if (current_rank == prev_rank) {
              ++current_partition_size;
            } else {
              if (current_partition_size > 1)
                degeneracy /=
                    boost::math::factorial<double>(current_partition_size);
              current_partition_size = 1;
              prev_rank = current_rank;
            }
          }
        }
      }
      if (canonical_only)
        degeneracy /= boost::math::factorial<double>(
            current_partition_size);  // account for the last partition
      const int min_exlev_R = std::max(
          -R,
          R - 2 * total_T_rank);  // at most 2*total_T_rank lines can point down

      if (canonical || !canonical_only) {
        using interval = boost::numeric::interval<int>;
        if (exlev + min_exlev_R <= 0 && 0 <= exlev + max_exlev_R) {  // VEV != 0
          assert(min_exlev_R <= max_exlev_R);
          screened_input->append(
              degeneracy == 1 ? term : ex<Constant>(degeneracy) * term);
        }
      }
    }  // term loop

    if (screened_input->size() == 0)
      return ex<Constant>(0);
    else {
      return sequant::mbpt::sr::so::vac_av(screened_input, op_connections, use_topology);
    }
  }
};  // screened_vac_av


class ccresidual {
  size_t P, N;
 public:
  ccresidual(size_t p, size_t n) : P(p), N(n) {}

  ExprPtr operator()(bool screen, bool use_topology, bool use_connectivity, bool canonical_only) {
    auto ahbar = [=](const bool screen) {
      auto connect = [=](std::initializer_list<std::pair<int, int>> connlist) {
        if (use_connectivity)
          return connlist;
        else
          return std::initializer_list<std::pair<int, int>>{};
      };
      auto result =
          screened_vac_av{0}(A(P) * H(), connect({}), screen, use_topology,
                             canonical_only) +
          screened_vac_av{1}(A(P) * H() * T(N), connect({{1, 2}}), screen,
                             use_topology, canonical_only) +
          ex<Constant>(1. / 2) * screened_vac_av{2}(A(P) * H() * T(N) * T(N),
                                                    connect({{1, 2}, {1, 3}}),
                                                    screen, use_topology,
                                                    canonical_only) +
          ex<Constant>(1. / 6) *
              screened_vac_av{3}(A(P) * H() * T(N) * T(N) * T(N),
                                 connect({{1, 2}, {1, 3}, {1, 4}}), screen,
                                 use_topology, canonical_only) +
          ex<Constant>(1. / 24) *
              screened_vac_av{4}(A(P) * H() * T(N) * T(N) * T(N) * T(N),
                                 connect({{1, 2}, {1, 3}, {1, 4}, {1, 5}}),
                                 screen, use_topology, canonical_only);
      simplify(result);

      return result;
    };

    return ahbar(screen);
  }
};  // class ccresidual

class ccresidual_vec {
  size_t P, PMIN, N;

 public:
  ccresidual_vec(size_t p, size_t pmin, size_t n) : P(p), PMIN(pmin), N(n) {}

  void operator()(std::vector<ExprPtr>& result, bool screen,
                      bool use_topology, bool use_connectivity,
                      bool canonical_only) {
    result[P] = ccresidual{P, N}(screen, use_topology, use_connectivity,
                                 canonical_only);
    rapid_simplify(result[P]);
    if (P > PMIN)
      ccresidual_vec{P - 1, PMIN, N}(result, screen, use_topology,
                                     use_connectivity, canonical_only);
  }
};  // class ccresidual_vec

}  // namespace

class cceqvec {
  size_t N, P, PMIN;

 public:
  cceqvec(size_t n, size_t p = std::numeric_limits<size_t>::max(),
          size_t pmin = 1)
      : N(n), P(p == std::numeric_limits<size_t>::max() ? n : p), PMIN(pmin) {}
  std::vector<ExprPtr> operator()(bool screen, bool use_topology,
                                  bool use_connectivity, bool canonical_only) {
    std::vector<ExprPtr> result(P + 1);
    ccresidual_vec{P, PMIN, N}(result, screen, use_topology, use_connectivity,
                               canonical_only);
    return result;
  }
};  // class cceqvec

#define runtime_assert(tf) if (!(tf)) { \
                             std::ostringstream oss; \
                             oss << "failed assert at line " << __LINE__ << " in function " << __func__; \
                             throw std::runtime_error(oss.str().c_str()); \
                           }

TimerPool<32> tpool;

class compute_cceqvec {
  size_t P, PMIN, N;

 public:
  compute_cceqvec(size_t p, size_t pmin, size_t n) : P(p), PMIN(pmin), N(n) {}

  void operator()(bool print, bool screen, bool use_topology,
                       bool use_connectivity, bool canonical_only) {
    tpool.start(N);
    auto eqvec =
        cceqvec{N, P}(screen, use_topology, use_connectivity, canonical_only);
    tpool.stop(N);
    std::wcout << std::boolalpha << "expS" << N << "[screen=" << screen
               << ",use_topology=" << use_topology
               << ",use_connectivity=" << use_connectivity
               << ",canonical_only=" << canonical_only << "] computed in "
               << tpool.read(N) << " seconds" << std::endl;
    for (size_t R = PMIN; R <= P; ++R) {
      std::wcout << "R" << R << "(expS" << N << ") has " << eqvec[R]->size()
                 << " terms:" << std::endl;
      if (print) std::wcout << to_latex_align(eqvec[R], 20, 3) << std::endl;

      // validate known sizes of some CC residuals
      if (R == 1 && N == 2) runtime_assert(eqvec[R]->size() == 14);
      if (R == 2 && N == 2) runtime_assert(eqvec[R]->size() == 31);
      if (R == 3 && N == 3) runtime_assert(eqvec[R]->size() == 47);
      if (R == 4 && N == 4) runtime_assert(eqvec[R]->size() == 74);
      if (R == 5 && N == 5) runtime_assert(eqvec[R]->size() == 99);
    }
  }
};  // class compute_cceqvec

class compute_all {
  size_t NMAX;

 public:
  compute_all(size_t nmax) : NMAX(nmax) {}

  void operator()(bool print = true, bool screen = true,
                  bool use_topology = true, bool use_connectivity = true,
                  bool canonical_only = true) {
    for (size_t N = 2; N <= NMAX; ++N)
      compute_cceqvec{N, 1, N}(print, screen, use_topology,
                               use_connectivity, canonical_only);
  }
};  // class compute_all

