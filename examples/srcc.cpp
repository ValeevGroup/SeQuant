#include <clocale>
#include <iostream>
#include "../src/SeQuant2/mbpt/spin.hpp"
#include "../src/SeQuant2/mbpt/sr/sr.hpp"

using namespace sequant2;

namespace {

template <size_t P, size_t N>
auto ccresidual() {
  using namespace sequant2::mbpt::sr::so;
  // using namespace sequant2::mbpt::sr::so::pno;

  // this is clearly very wasteful since we are including many terms that are 0
  // TODO screen out zeroes
  auto result =
      vac_av(A<P>() * H()) + vac_av(A<P>() * H() * T<N>(), {{1, 2}}) +
      ex<Constant>(1. / 2) *
          vac_av(A<P>() * H() * T<N>() * T<N>(), {{1, 2}, {1, 3}}) +
      ex<Constant>(1. / 6) * vac_av(A<P>() * H() * T<N>() * T<N>() * T<N>(),
                                    {{1, 2}, {1, 3}, {1, 4}}) +
      ex<Constant>(1. / 24) *
          vac_av(A<P>() * H() * T<N>() * T<N>() * T<N>() * T<N>(),
                 {{1, 2}, {1, 3}, {1, 4}, {1, 5}});
  simplify(result);
  expand(result);
  canonicalize(result);
  return result;
}

template <size_t P, size_t PMIN, size_t N>
void ccresidual_rec(std::vector<ExprPtr>& result) {
  result[P] = ccresidual<P, N>();
  simplify(result[P]);
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
  sequant2::IndexSpace::register_standard_instances();
  sequant2::detail::OpIdRegistrar op_id_registrar;
  TensorCanonicalizer::set_cardinal_tensor_labels({L"A", L"f", L"g", L"t"});
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  {  // CC amplitude eqs
    constexpr size_t N = 1;
    constexpr size_t P = 1;
    constexpr size_t PMIN = 1;
    auto eqvec = cceqvec<N, P>();
    for (size_t R = PMIN; R <= P; ++R) {
      std::wcout << "R" << R << "(expS" << N << ") has " << eqvec[R]->size()
                 << " terms:\n"
                 << to_latex_align(eqvec[R], 20, 5) << std::endl;
      if (R == 2 && N == 2) assert(eqvec[R]->size() == 31);
    }
  }

  return 0;
}
