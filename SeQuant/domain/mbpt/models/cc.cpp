#include <SeQuant/domain/mbpt/formalism.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>

#include <clocale>
#include <iostream>

#include <boost/math/special_functions/factorials.hpp>
// boost/numeric/interval does not know about arm rounding .. on arm64/macos use
// c99 rounding
#if defined(__arm64__) && defined(__APPLE__) && !defined(__USE_ISOC99)
#define __USE_ISOC99 1
#include <boost/numeric/interval.hpp>
#undef __USE_ISOC99
#else
#include <boost/numeric/interval.hpp>
#endif

#include <SeQuant/core/op.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <SeQuant/domain/mbpt/sr/sr.hpp>

namespace sequant::mbpt::sr::so {

namespace {

/// computes VEV using excitation level screening (unless @c
/// screen is set) + computes only canonical (with T ranks increasing) terms
class screened_vac_av {
 private:
  size_t K;

 public:
  /// \param k power of T
  screened_vac_av(size_t k) : K(k) {}

  //  ExprPtr operator()(const ExprPtr& expr,
  //                     std::initializer_list<std::pair<int, int>>
  //                     op_connections, bool screen = true, bool use_topology =
  //                     true, bool canonical_only = true) {}

  /// computes VEV for A(P)*H*T(N)^K using excitation-level screening,
  /// including only canonical terms
  /// and
  /// \param expr input expression, must contain `A`, `H` (`f` or `g`),
  ///        and `K` `T`'s
  /// \param op_connections specifies the connectivity
  /// \param screen if false, will use brute-force evaluation
  /// \param use_topology if true, forces topological optimization
  /// \param canonical_only if true AND \p screen is true then optimize
  ///        evaluation by combining equivalent terms such that VEV evaluation
  ///        only involves canonical products (e.g., evaluate only H*T1*T2 and
  ///        not H*T2*T1)
  /// \return the resulting VEV
  ExprPtr t(const ExprPtr& expr,
            std::initializer_list<std::pair<int, int>> op_connections,
            bool screen = true, bool use_topology = true,
            bool canonical_only = true) {
    // TODO: Implement antisymm here
    if (!screen)
      return sequant::mbpt::sr::so::vac_av(expr, op_connections, use_topology);

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
        auto p = 4 + 2 * k;
        assert(term_prod.factor(p)->is<Tensor>());
        assert(term_prod.factor(p)->as<Tensor>().label() == L"t");
        const auto current_rank = term_prod.factor(p)->as<Tensor>().rank();
        exlev += current_rank;
        total_T_rank += current_rank;
        // screen out the noncanonical terms, if needed
        if (canonical_only) {
          if (current_rank < prev_rank)  // if T ranks are not increasing, omit
            canonical = false;
          else {                         // else keep track of degeneracy
            assert(current_rank != 0);
            if (current_rank == prev_rank)
              ++current_partition_size;
            else {
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
      return sequant::mbpt::sr::so::vac_av(screened_input, op_connections,
                                           use_topology);
    }
  }  // screened_vac_av_t

  /// computes VEV for (1+Λ)*H*T(N)^K*adj(A(P)) using excitation-level
  /// screening, including only canonical terms
  /// and
  /// \param expr input expression, must contain `A`, `H` (`f` or `g`), `Λ`,
  ///        and `K` `T`'s
  /// \param op_connections specifies the connectivity
  /// \param screen if false, will use brute-force evaluation
  /// \param use_topology if true, forces topological optimization
  /// \param canonical_only if true AND \p screen is true then optimize
  ///        evaluation by combining equivalent terms such that VEV evaluation
  ///        only involves canonical products (e.g., evaluate only H*T1*T2 and
  ///        not H*T2*T1)
  /// \return the resulting VEV
  ExprPtr lambda(const ExprPtr& expr,
                 std::initializer_list<std::pair<int, int>> op_connections,
                 bool screen = true, bool use_topology = true,
                 bool canonical_only = true) {
    if (!screen)
      return sequant::mbpt::sr::so::vac_av(expr, op_connections, use_topology);

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

      // locate projector - will be the right hand end factor in this case
      assert(term_prod.factor(term_prod.factors().size() - 1)->is<Tensor>());
      assert(term_prod.factor(term_prod.factors().size() - 1)
                 ->as<Tensor>()
                 .label() == L"A");
      const int P =
          term_prod.factor(term_prod.factors().size() - 1)->as<Tensor>().rank();

      // locate hamiltonian
      // can be the first or the third factor because of (1 + Λ)
      assert(term_prod.factor(0)->is<Tensor>());
      auto h_term = term_prod.factor(0)->as<Tensor>();
      if (h_term.label() != L"f" &&
          h_term.label() != L"g") {  // if term has λ amplitudes
        h_term = term_prod.factor(2)->as<Tensor>();
        assert(term_prod.factors().size() == 6 + 2 * K);
      } else {  // if term has no λ amplitudes
        assert(term_prod.factors().size() == 4 + 2 * K);
      }

      assert(h_term.label() == L"f" || h_term.label() == L"g");
      const int R = h_term.rank();
      const int max_exlev_R = R - K;

      auto exlev = -P;

      bool canonical = true;
      double degeneracy =
          canonical_only ? boost::math::factorial<double>(K) : 1;
      int total_T_rank = 0;
      int prev_rank = 0;
      int current_partition_size = 1;
      for (size_t k = 0; k < K && canonical; ++k) {
        // k goes from 0 to K-1, since the last two terms are from A
        // if λ amplitudes are present, p = 4 + 2k; else, p = 2 + 2k
        auto p = term_prod.factor(0)->as<Tensor>().label() == L"λ"
                     ? (4 + 2 * k)
                     : (2 + 2 * k);

        assert(term_prod.factor(p)->is<Tensor>());
        assert(term_prod.factor(p)->as<Tensor>().label() == L"t");
        const auto current_rank = term_prod.factor(p)->as<Tensor>().rank();
        exlev += current_rank;
        total_T_rank += current_rank;
        // screen out noncanonical terms, if needed
        if (canonical_only) {
          if (current_rank < prev_rank)  // if T ranks are not increasing, omit
            canonical = false;
          else {                         // else keep track of degeneracy
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
      const int min_exlev_R = std::max(-R, R - 2 * total_T_rank);
      if (canonical || !canonical_only) {
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
      return sequant::mbpt::sr::so::vac_av(screened_input, op_connections,
                                           use_topology);
    }
  }  // screened_vac_av_lambda

};   // screened_vac_av

/// Evaluates coupled-cluster amplitude equation, `<P|(H exp(T(N))_c|0>`,
/// for particular `P` and `N`
class cceqs_t {
  size_t P, N;

 public:
  cceqs_t(size_t p, size_t n) : P(p), N(n) {}

  /// Evalaute the coupled-cluster amplitude equations, <P|(H exp(T(N))_c|0>
  /// \param screen if true, will use screening (see `screened_vac_av`)
  /// \param use_topology if true, forces topological optimization
  /// \param use_connectivity if true, will tell Wick engine to prune search
  ///                         tree using known connectivity information
  /// \param canonical_only if true AND \p screen is true then combine
  ///                       equivalent terms (see `screened_vac_av`)
  /// \return the result
  ExprPtr operator()(bool screen, bool use_topology, bool use_connectivity,
                     bool canonical_only) {
    // currently topological equivalence of indices within a normal operator is
    // not detected, assumed based on use_topology ... so turn off use of
    // topology if antisymm=false
    if (get_default_formalism().two_body_interaction() !=
        TwoBodyInteraction::Antisymm)
      use_topology = false;

    auto ahbar = [=](const bool screen) {
      auto connect = [=](std::initializer_list<std::pair<int, int>> connlist) {
        if (use_connectivity)
          return connlist;
        else
          return std::initializer_list<std::pair<int, int>>{};
      };
      auto result =
          screened_vac_av{0}.t(A(P) * H(), connect({}), screen, use_topology,
                               canonical_only) +
          screened_vac_av{1}.t(A(P) * H() * T(N, N), connect({{1, 2}}), screen,
                               use_topology, canonical_only) +
          ex<Constant>(1. / 2) *
              screened_vac_av{2}.t(A(P) * H() * T(N, N) * T(N, N),
                                   connect({{1, 2}, {1, 3}}), screen,
                                   use_topology, canonical_only) +
          ex<Constant>(1. / 6) *
              screened_vac_av{3}.t(A(P) * H() * T(N, N) * T(N, N) * T(N, N),
                                   connect({{1, 2}, {1, 3}, {1, 4}}), screen,
                                   use_topology, canonical_only) +
          ex<Constant>(1. / 24) *
              screened_vac_av{4}.t(
                  A(P) * H() * T(N, N) * T(N, N) * T(N, N) * T(N, N),
                  connect({{1, 2}, {1, 3}, {1, 4}, {1, 5}}), screen,
                  use_topology, canonical_only);
      simplify(result);

      return result;
    };

    return ahbar(screen);
  }
};  // class cceqs_t

/// Evaluates coupled-cluster λ amplitude equation, `<0|(1+Λ)(H exp(T(N))_c|P>`,
/// for particular `P` and `N`
class cceqs_lambda {
  size_t P, N;

 public:
  cceqs_lambda(size_t p, size_t n) : P(p), N(n) {}

  /// Evaluate the coupled-cluster λ amplitude equations,
  /// `<0|(1+Λ)(H exp(T(N))_c|P>` for particular `P` and `N`
  /// \param screen if true, will use screening (see `screened_vac_av`)
  /// \param use_topology if true, will tell Wick engine to prune search
  ///                         tree using known connectivity information
  /// \param canonical_only if true AND \p screen is true then combine
  ///                       equivalent terms (see `screened_vac_av`)
  /// \return the result

  ExprPtr operator()(bool screen, bool use_topology, bool use_connectivity,
                     bool canonical_only) {
    if (get_default_formalism().two_body_interaction() !=
        TwoBodyInteraction::Antisymm)
      use_topology = false;

    auto ahbar = [=](const bool screen) {
      auto connect = [=](std::initializer_list<std::pair<int, int>> connlist) {
        if (use_connectivity)
          return connlist;
        else
          return std::initializer_list<std::pair<int, int>>{};
      };
      const auto One = ex<Constant>(1);
      auto result =
          screened_vac_av{0}.lambda((One + Lambda(N, N)) * H() * adjoint(A(P)),
                                    connect({}), screen, use_topology,
                                    canonical_only) +
          screened_vac_av{1}.lambda(
              (One + Lambda(N, N)) * H() * T(N, N) * adjoint(A(P)),
              connect({{1, 2}}), screen, use_topology, canonical_only) +
          ex<Constant>(1. / 2) *
              screened_vac_av{2}.lambda((One + Lambda(N, N)) * H() * T(N, N) *
                                            T(N, N) * adjoint(A(P)),
                                        connect({{1, 2}, {1, 3}}), screen,
                                        use_topology, canonical_only) +
          ex<Constant>(1. / 6) *
              screened_vac_av{3}.lambda((One + Lambda(N, N)) * H() * T(N, N) *
                                            T(N, N) * T(N, N) * adjoint(A(P)),
                                        connect({{1, 2}, {1, 3}, {1, 4}}),
                                        screen, use_topology, canonical_only);
      simplify(result);
      return result;
    };
    return ahbar(screen);
  }
};  // class cceqs_λ

}  // namespace

cceqs::cceqs(size_t n, size_t p, size_t pmin)
    : N(n), P(p == std::numeric_limits<size_t>::max() ? n : p), PMIN(pmin) {}

std::vector<ExprPtr> cceqs::t(bool screen, bool use_topology,
                              bool use_connectivity, bool canonical_only) {
  std::vector<ExprPtr> result(P + 1);
  for (auto p = P; p >= PMIN; --p) {
    result.at(p) =
        cceqs_t{p, N}(screen, use_topology, use_connectivity, canonical_only);
    rapid_simplify(result[p]);
  }
  return result;
}

std::vector<ExprPtr> cceqs::lambda(bool screen, bool use_topology,
                                   bool use_connectivity, bool canonical_only) {
  std::vector<ExprPtr> result(P + 1);
  for (auto p = P; p >= PMIN; --p) {
    result.at(p) = cceqs_lambda{p, N}(screen, use_topology, use_connectivity,
                                      canonical_only);
    rapid_simplify(result[p]);
  }
  return result;
}

}  // namespace sequant::mbpt::sr::so
