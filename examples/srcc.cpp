#include <clocale>
#include <iostream>

#include <boost/numeric/interval.hpp>
#include <boost/math/special_functions/factorials.hpp>

#include "../src/domain/mbpt/spin.hpp"
#include "../src/domain/mbpt/sr/sr.hpp"

using namespace sequant;

using namespace sequant::mbpt::sr::so;
// using namespace sequant::mbpt::sr::so::pno;

namespace {

/// computes VEE for A(P)*H*T(N)^K using excitation level screening (unless @c
/// screen is set) + computes only canonical (with T ranks increasing) terms
template <size_t K>
ExprPtr screened_vac_av(
    const ExprPtr& expr,
    std::initializer_list<std::pair<int, int>> op_connections,
    bool screen = true) {
  if (!screen) return vac_av(expr, op_connections);

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
    // number of possible permutations excluding permutations within same-rank partitions of T
    // degeneracy = K! / M1! M2! .. where M1, M2 ... are sizes of each partition
    double degeneracy = boost::math::factorial<double>(K);
    int total_T_rank = 0;
    int prev_rank = 0;
    int current_partition_size = 1;  // size of current same-rank partition oF T
    for (size_t k = 0; k != K && canonical; ++k) {
      auto p = 4 + k * 2;
      assert(term_prod.factor(p)->is<Tensor>());
      assert(term_prod.factor(p)->as<Tensor>().label() == L"t");
      const auto current_rank = term_prod.factor(p)->as<Tensor>().rank();
      if (current_rank < prev_rank)  // if T ranks are not increasing, omit
        canonical = false;
      else {  // else keep track of degeneracy
        assert(current_rank != 0);
        if (current_rank == prev_rank) {
          ++current_partition_size;
        }
        else {
          if (current_partition_size > 1)
            degeneracy /= boost::math::factorial<double>(current_partition_size);
          current_partition_size = 1;
          prev_rank = current_rank;
        }
        exlev += current_rank;
        total_T_rank += current_rank;
      }
    }
    const int min_exlev_R = std::max(-R, R-2*total_T_rank);  // at most 2*total_T_rank lines can point down

    if (canonical) {
      using interval = boost::numeric::interval<int>;
      assert(min_exlev_R <= max_exlev_R);
      if (exlev + min_exlev_R <= 0 && 0 <= exlev + max_exlev_R) {  // VEE != 0
        screened_input->append(
            degeneracy == 1 ? term : ex<Constant>(degeneracy) * term);
      }
    }
  }  // term loop

  if (screened_input->size() == 0)
    return ex<Constant>(0);
  else {
    return vac_av(screened_input, op_connections);
  }
}

template <size_t P, size_t N>
auto ccresidual() {
  auto ahbar = [](const bool screen) {
    auto result =
        screened_vac_av<0>(A<P>() * H(), {}, screen) +
        screened_vac_av<1>(A<P>() * H() * T<N>(), {{1, 2}}, screen) +
        ex<Constant>(1. / 2) *
            screened_vac_av<2>(A<P>() * H() * T<N>() * T<N>(), {{1, 2}, {1, 3}},
                               screen) +
        ex<Constant>(1. / 6) *
            screened_vac_av<3>(A<P>() * H() * T<N>() * T<N>() * T<N>(),
                               {{1, 2}, {1, 3}, {1, 4}}, screen) +
        ex<Constant>(1. / 24) *
            screened_vac_av<4>(A<P>() * H() * T<N>() * T<N>() * T<N>() * T<N>(),
                               {{1, 2}, {1, 3}, {1, 4}, {1, 5}}, screen);
    simplify(result);

    return result;
  };

  return ahbar(true);
}

template <size_t P, size_t PMIN, size_t N>
void ccresidual_rec(std::vector<ExprPtr>& result) {
  result[P] = ccresidual<P, N>();
  rapid_simplify(result[P]);
  if constexpr (P > PMIN) ccresidual_rec<P - 1, PMIN, N>(result);
}

}  // namespace

template <size_t N, size_t P = N, size_t PMIN = 1>
std::vector<ExprPtr> cceqvec() {
  std::vector<ExprPtr> result(P + 1);
  ccresidual_rec<P, PMIN, N>(result);
  return result;
}

int main(int argc, char* argv[]) {
  std::setlocale(LC_ALL, "en_US.UTF-8");
  std::cout.precision(std::numeric_limits<double>::max_digits10);
  sequant::IndexSpace::register_standard_instances();
  sequant::detail::OpIdRegistrar op_id_registrar;
  TensorCanonicalizer::set_cardinal_tensor_labels({L"A", L"f", L"g", L"t"});
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  {  // CC amplitude eqs
    constexpr size_t N = 6;
    constexpr size_t P = 6;
    constexpr size_t PMIN = 1;
    auto eqvec = cceqvec<N, P>();
    for (size_t R = PMIN; R <= P; ++R) {
      std::wcout << "R" << R << "(expS" << N << ") has " << eqvec[R]->size()
                 << " terms:\n"
                 << to_latex_align(eqvec[R], 20, 5) << std::endl;

      // validate known sizes of some CC residuals
      if (R == 1 && N == 2) assert(eqvec[R]->size() == 14);
      if (R == 2 && N == 2) assert(eqvec[R]->size() == 31);
      if (R == 3 && N == 3) assert(eqvec[R]->size() == 47);
      if (R == 4 && N == 4) assert(eqvec[R]->size() == 74);
      if (R == 5 && N == 5) assert(eqvec[R]->size() == 99);
    }
  }

  return 0;
}
