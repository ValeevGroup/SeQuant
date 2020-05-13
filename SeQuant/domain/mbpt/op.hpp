//
// Created by Eduard Valeyev on 2019-03-26.
//

#ifndef SEQUANT_DOMAIN_MBPT_OP_HPP
#define SEQUANT_DOMAIN_MBPT_OP_HPP

#include <string>
#include <vector>

namespace sequant {
namespace mbpt {

enum class OpType { h, f, g, t, l, A, L, R, R12 };
enum class OpClass { ex, deex, gen };

/// @return the tensor labels in the cardinal order
std::vector<std::wstring> cardinal_tensor_labels();

std::wstring to_wstring(OpType op);
OpClass to_class(OpType op);

}  // namespace mbpt
}  // namespace sequant

#endif //SEQUANT_DOMAIN_MBPT_OP_HPP
