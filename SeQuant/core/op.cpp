//
// Created by Eduard Valeyev on 2019-02-18.
//

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/op.hpp>

namespace sequant {
template <Statistics S>
const container::svector<std::wstring>& NormalOperator<S>::labels() {
  using namespace std::literals;
  static_assert(S == Statistics::FermiDirac || S == Statistics::BoseEinstein);
  static const container::svector<std::wstring> labels_ =
      S == Statistics::FermiDirac
          ? container::svector<std::wstring>{L"a"s, L"ã"s}
          : container::svector<std::wstring>{L"b"s, L"b̃"s};
  return labels_;
}

template class SEQUANT_SYMB_TEMPLATE_DEF_API Operator<Statistics::FermiDirac>;
template class SEQUANT_SYMB_TEMPLATE_DEF_API Operator<Statistics::BoseEinstein>;
template class SEQUANT_SYMB_TEMPLATE_DEF_API
    NormalOperator<Statistics::FermiDirac>;
template class SEQUANT_SYMB_TEMPLATE_DEF_API
    NormalOperator<Statistics::BoseEinstein>;
template class SEQUANT_SYMB_TEMPLATE_DEF_API
    NormalOperatorSequence<Statistics::FermiDirac>;
template class SEQUANT_SYMB_TEMPLATE_DEF_API
    NormalOperatorSequence<Statistics::BoseEinstein>;

}  // namespace sequant
