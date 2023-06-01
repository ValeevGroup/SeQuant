#include <SeQuant/domain/mbpt/formalism.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>

#include <clocale>
#include <iostream>

#include <SeQuant/core/math.hpp>

#include <SeQuant/core/op.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <SeQuant/domain/mbpt/sr/sr.hpp>

namespace sequant::mbpt::sr {

namespace {

// TODO: generalize screening logic

/// computes VEV using excitation level screening (unless @c
/// screen is set) + computes only canonical (with T ranks increasing) terms
class screened_vac_av {
 private:
  size_t K;

 public:
  /// \param k power of T
  screened_vac_av(size_t k) : K(k) {}

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
      return sequant::mbpt::sr::vac_av(expr, op_connections, use_topology);

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
      rational degeneracy = canonical_only ? sequant::factorial(K) : 1;
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
                degeneracy /= sequant::factorial(current_partition_size);
              current_partition_size = 1;
              prev_rank = current_rank;
            }
          }
        }
      }
      if (canonical_only)
        degeneracy /= sequant::factorial(
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
      return sequant::mbpt::sr::vac_av(screened_input, op_connections,
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
    // screening for lambda is not available now
    assert(!screen &&
           "screening for λ residual equations is not available now");
    if (!screen)
      return sequant::mbpt::sr::vac_av(expr, op_connections, use_topology);

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
      // obtaining rank from the projector
      const int P =
          term_prod.factor(term_prod.factors().size() - 1)->as<Tensor>().rank();

      // initialize h_factor and l_factor
      sequant::Tensor h_factor;  // just Tensor should be enough
      sequant::Tensor l_factor;

      // locate hamiltonian and λ factors - the hamiltonian factor could be
      // first or third depending on if λ amplitudes are present
      assert(term_prod.factor(0)->is<Tensor>());
      if (term_prod.factor(0)->as<Tensor>().label() == L"λ") {
        // if first factor is λ, third factor should be f or g
        assert(term_prod.factor(2)->as<Tensor>().label() == L"f" ||
               term_prod.factor(2)->as<Tensor>().label() == L"g");
        // assert term size
        assert(term_prod.factors().size() == 6 + 2 * K);
        h_factor = term_prod.factor(2)->as<Tensor>();
        l_factor = term_prod.factor(0)->as<Tensor>();
      } else {
        // else, the first factor should be from hamiltonian
        assert(term_prod.factor(0)->as<Tensor>().label() == L"f" ||
               term_prod.factor(0)->as<Tensor>().label() == L"g");
        // assert term size
        assert(term_prod.factors().size() == 4 + 2 * K);
        h_factor = term_prod.factor(0)->as<Tensor>();
        // no λ amplitudes will be present
      }

      const int R = h_factor.rank();
      const int L = l_factor.rank();

      const int max_exlev_R = R - K;

      auto exlev = -P;

      bool canonical = true;
      rational degeneracy = canonical_only ? sequant::factorial(K) : 1;
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
                degeneracy /= sequant::factorial(current_partition_size);
              current_partition_size = 1;
              prev_rank = current_rank;
            }
          }
        }
      }
      if (canonical_only)
        degeneracy /= sequant::factorial(
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
      return sequant::mbpt::sr::vac_av(screened_input, op_connections,
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
          ex<Constant>(rational{1, 2}) *
              screened_vac_av{2}.t(A(P) * H() * T(N, N) * T(N, N),
                                   connect({{1, 2}, {1, 3}}), screen,
                                   use_topology, canonical_only) +
          ex<Constant>(rational{1, 6}) *
              screened_vac_av{3}.t(A(P) * H() * T(N, N) * T(N, N) * T(N, N),
                                   connect({{1, 2}, {1, 3}, {1, 4}}), screen,
                                   use_topology, canonical_only) +
          ex<Constant>(rational{1, 24}) *
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
      auto result =
          screened_vac_av{0}.lambda(H() * adjoint(A(P)), connect({{0, 1}}),
                                    screen, use_topology, canonical_only) +
          screened_vac_av{0}.lambda(Lambda(N, N) * H() * adjoint(A(P)),
                                    connect({{1, 2}}), screen, use_topology,
                                    canonical_only) +
          screened_vac_av{1}.lambda(H() * T(N, N) * adjoint(A(P)),
                                    connect({{0, 1}, {0, 2}}), screen,
                                    use_topology, canonical_only) +
          screened_vac_av{1}.lambda(
              Lambda(N, N) * H() * T(N, N) * adjoint(A(P)),
              connect({{1, 2}, {1, 3}}), screen, use_topology, canonical_only) +
          ex<Constant>(rational{1, 2}) *
              screened_vac_av{2}.lambda(H() * T(N, N) * T(N, N) * adjoint(A(P)),
                                        connect({{0, 1}, {0, 2}, {0, 3}}),
                                        screen, use_topology, canonical_only) +
          ex<Constant>(rational{1, 2}) *
              screened_vac_av{2}.lambda(
                  Lambda(N, N) * H() * T(N, N) * T(N, N) * adjoint(A(P)),
                  connect({{1, 2}, {1, 3}, {1, 4}}), screen, use_topology,
                  canonical_only) +
          ex<Constant>(rational{1, 6}) *
              screened_vac_av{3}.lambda(
                  H() * T(N, N) * T(N, N) * T(N, N) * adjoint(A(P)),
                  connect({{0, 1}, {0, 2}, {0, 3}, {0, 4}}), screen,
                  use_topology, canonical_only) +
          ex<Constant>(rational{1, 6}) *
              screened_vac_av{3}.lambda(
                  Lambda(N, N) * H() * T(N, N) * T(N, N) * T(N, N) *
                      adjoint(A(P)),
                  connect({{1, 2}, {1, 3}, {1, 4}, {1, 5}}), screen,
                  use_topology, canonical_only);

      simplify(result);
      return result;
    };
    return ahbar(screen);
  }
};  // class cceqs_lambda

}  // namespace

cceqs::cceqs(size_t n, size_t p, size_t pmin)
    : N(n), P(p == std::numeric_limits<size_t>::max() ? n : p), PMIN(pmin) {}

std::vector<ExprPtr> cceqs::t(bool screen, bool use_topology,
                              bool use_connectivity, bool canonical_only) {
  constexpr bool use_ops = true;
  if (use_ops) {
    // 1. construct hbar(op) in canonical form
    auto hbar = op::H();
    auto H_Tk = hbar;
    for (int64_t k = 1; k <= 4; ++k) {
      H_Tk = simplify(ex<Constant>(rational{1, k}) * H_Tk * op::T(N));
      hbar += H_Tk;
    }

    // 2. project onto each manifold, screen, lower to tensor form and wick it
    std::vector<ExprPtr> result(P + 1);
    for (auto p = P; p >= PMIN; --p) {
      // 2.a. screen out terms that cannot give nonzero after projection onto
      // <P|
      std::shared_ptr<Sum>
          hbar_p;     // products that can produce excitations of rank p
      std::shared_ptr<Sum>
          hbar_le_p;  // keeps products that can produce excitations rank <=p
      for (auto& term : *hbar) {
        assert(term->is<Product>() || term->is<op_t>());

        if (op::contains_up_to_rank(term, p)) {
          if (!hbar_le_p)
            hbar_le_p = std::make_shared<Sum>(ExprPtrList{term});
          else
            hbar_le_p->append(term);
          if (op::contains_rank(term, p)) {
            if (!hbar_p)
              hbar_p = std::make_shared<Sum>(ExprPtrList{term});
            else
              hbar_p->append(term);
          }
        }
      }
      hbar = hbar_le_p;

      // 2.b multiply by A(P)
      auto A_hbar = simplify(op::A(p) * hbar_p);

      // 2.c compute vacuum average
      result.at(p) = op::vac_av(A_hbar);
      simplify(result.at(p));
    }

    return result;
  } else {
    std::vector<ExprPtr> result(P + 1);
    for (auto p = P; p >= PMIN; --p) {
      result.at(p) =
          cceqs_t{p, N}(screen, use_topology, use_connectivity, canonical_only);
    }
    return result;
  }
}

std::vector<ExprPtr> cceqs::lambda(bool screen, bool use_topology,
                                   bool use_connectivity, bool canonical_only) {
  constexpr bool use_ops = true;
  if (use_ops) { // in development
    // construct hbar
    auto hbar = op::H();
    auto H_Tk = hbar;
    for (int64_t k = 1; k <= 3; ++k) {
      H_Tk = simplify(ex<Constant>(rational{1, k}) * H_Tk * op::T(N));
      hbar += H_Tk;
    }

//    std::wcout << "hbar: \n" << to_latex_align(hbar, 0, 4) << std::endl;
    // multiply with (1 + Λ)
    const auto One = ex<Constant>(1);
    auto lhbar = simplify((One + op::Lambda(N)) * hbar);

//    std::wcout << "lhbar: \n" << to_latex_align(lhbar, 0, 4) << std::endl;

    // 2. project onto each manifold, screen, lower to tensor form and wick it
    std::vector<ExprPtr> result(P + 1);
    for (auto p = P; p >= PMIN; --p) {
      // 2.a. screen out terms that cannot give nonzero after projection onto
      // <P|
      std::shared_ptr<Sum>
          hbar_p;     // products that can produce excitations of rank p
      std::shared_ptr<Sum>
          hbar_le_p;  // keeps products that can produce excitations rank <=p
      for (auto& term : *lhbar) {  // pick terms from lhbar
        assert(term->is<Product>() || term->is<op_t>());

        if (op::contains_up_to_rank(term, p)) {
          if (!hbar_le_p)
            hbar_le_p = std::make_shared<Sum>(ExprPtrList{term});
          else
            hbar_le_p->append(term);
          if (op::contains_rank(term, p)) {
            if (!hbar_p)
              hbar_p = std::make_shared<Sum>(ExprPtrList{term});
            else
              hbar_p->append(term);
          }
        }
      }
      lhbar = hbar_le_p; // not needed

//      std::wcout << "p = " << p << std::endl;
//      std::wcout << "hbar-le-p: \n"
//                 << to_latex_align(hbar_le_p, 0, 4) << std::endl;
//      std::wcout << "hbar-p: \n" << to_latex_align(hbar_p, 0, 4) << std::endl;

      // 2.b multiply by adjoint of A(P) on the right side

      auto A_hbar = simplify(hbar_p * adjoint(op::A(p)));
//      std::wcout << "Ahbar: \n" << to_latex_align(A_hbar, 0, 4) << std::endl;

      // temp
      std::vector<std::pair<std::wstring, std::wstring>> new_op_connect = {
          {L"h", L"t"}, {L"f", L"t"}, {L"g", L"t"},
          {L"h", L"A"}, {L"f", L"A"}, {L"g", L"A"}};

      // 2.c compute vacuum average
      result.at(p) = op::vac_av(A_hbar, new_op_connect);
      simplify(result.at(p));
//      std::wcout << "result.at(p): \n"
//                 << to_latex_align(result.at(p), 0, 4) << std::endl;
    }
    return result;
  } else {
    std::vector<ExprPtr> result(P + 1);
    for (auto p = P; p >= PMIN; --p) {
      result.at(p) = cceqs_lambda{p, N}(screen, use_topology, use_connectivity,
                                        canonical_only);
    }
    return result;
  }
}

}  // namespace sequant::mbpt::sr
