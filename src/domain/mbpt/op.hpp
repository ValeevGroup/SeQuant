//
// Created by Eduard Valeyev on 2019-03-26.
//

#ifndef SEQUANT_OP_HPP
#define SEQUANT_OP_HPP

namespace sequant {
namespace mbpt {

enum class OpType { f, g, t, l, A, L, R };

const static std::vector<std::wstring> cardinal_tensor_labels = {
    L"A", L"L", L"λ", L"f", L"g", L"t", L"R", L"S", L"a", L"ã", L"b", L"ᵬ"};

}  // namespace mbpt
}  // namespace sequant

#endif //SEQUANT_OP_HPP
