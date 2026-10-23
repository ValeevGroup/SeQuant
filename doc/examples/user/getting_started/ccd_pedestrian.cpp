#include <SeQuant/core/expr.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/op.hpp>

// SEQUANT_DOC_EXAMPLE: DO NOT RUN

int main() {
  using namespace sequant;
  using namespace sequant::mbpt;
  set_default_context({.index_basis_registry_shared_ptr = make_min_sr_spaces(),
                       .vacuum = Vacuum::SingleProduct});
  set_default_mbpt_context({.op_registry_ptr = make_legacy_registry()});

  // start-snippet-1
  auto t2 = ex<Constant>(rational(1, 4)) *
            ex<Tensor>("t", bra{"a_1", "a_2"}, ket{"i_1", "i_2"},
                       Symmetry::Antisymm) *
            ex<FNOperator>(cre{"a_1", "a_2"}, ann{"i_1", "i_2"});
  // end-snippet-1

  // start-snippet-2
  auto H = ex<Tensor>("f", bra{"p_1"}, ket{"p_2"}, Symmetry::Nonsymm) *
               ex<FNOperator>(cre{"p_1"}, ann{"p_2"}) +
           ex<Constant>(rational(1, 4)) *
               ex<Tensor>("g", bra{"p_1", "p_2"}, ket{"p_3", "p_4"},
                          Symmetry::Antisymm) *
               ex<FNOperator>(cre{"p_1", "p_2"}, ann{"p_3", "p_4"});
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
              ex<Tensor>("t", bra{"a_1", "a_2"}, ket{"i_1", "i_2"},
                         Symmetry::Antisymm) *
              ex<FNOperator>(cre{"a_1", "a_2"}, ann{"i_1", "i_2"});

  auto t2_1 = ex<Constant>(rational(1, 4)) *
              ex<Tensor>("t", bra{"a_3", "a_4"}, ket{"i_3", "i_4"},
                         Symmetry::Antisymm) *
              ex<FNOperator>(cre{"a_3", "a_4"}, ann{"i_3", "i_4"});

  auto c_htt =
      ex<Constant>(rational(1, 4)) * commutator(commutator(H, t2_0), t2_1);
  // end-snippet-5

  return 0;
}
