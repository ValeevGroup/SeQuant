#include <SeQuant/core/expr.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/op.hpp>

// SEQUANT_DOC_EXAMPLE: DO NOT RUN

int main() {
  using namespace sequant;
  using namespace sequant::mbpt;
  set_default_context({.index_space_registry_shared_ptr = make_min_sr_spaces(),
                       .vacuum = Vacuum::SingleProduct});

  // start-snippet-1
  auto t2 = ex<Constant>(rational(1, 4)) *
            ex<Tensor>(L"t", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"},
                       Symmetry::antisymm) *
            ex<FNOperator>(cre{L"a_1", L"a_2"}, ann{L"i_1", L"i_2"});
  // end-snippet-1

  // start-snippet-2
  auto H = ex<Tensor>(L"f", bra{L"p_1"}, ket{L"p_2"}, Symmetry::nonsymm) *
               ex<FNOperator>(cre{L"p_1"}, ann{L"p_2"}) +
           ex<Constant>(rational(1, 4)) *
               ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"},
                          Symmetry::antisymm) *
               ex<FNOperator>(cre{L"p_1", L"p_2"}, ann{L"p_3", L"p_4"});
  // end-snippet-2

  // start-snippet-3
  auto commutator = [](auto op1, auto op2) {
    return simplify(op1 * op2 - op2 * op1);
  };

  auto c_ht = commutator(H, t2);
  // end-snippet-3

  // start-snippet-4
  auto c_htt_v1 =
      ex<Constant>(rational(1, 2)) * commutator(commutator(H, t2), t2);
  // end-snippet-4

  // start-snippet-5
  auto t2_0 = ex<Constant>(rational(1, 4)) *
              ex<Tensor>(L"t", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"},
                         Symmetry::antisymm) *
              ex<FNOperator>(cre{L"a_1", L"a_2"}, ann{L"i_1", L"i_2"});

  auto t2_1 = ex<Constant>(rational(1, 4)) *
              ex<Tensor>(L"t", bra{L"a_3", L"a_4"}, ket{L"i_3", L"i_4"},
                         Symmetry::antisymm) *
              ex<FNOperator>(cre{L"a_3", L"a_4"}, ann{L"i_3", L"i_4"});

  auto c_htt =
      ex<Constant>(rational(1, 4)) * commutator(commutator(H, t2_0), t2_1);
  // end-snippet-5

  return 0;
}
