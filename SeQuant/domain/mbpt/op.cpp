#include "SeQuant/domain/mbpt/op.hpp"

#include "SeQuant/core/tensor.hpp"

#include <stdexcept>

namespace sequant {
namespace mbpt {

std::vector<std::wstring>
    cardinal_tensor_labels() {
  return {L"\\lambda",L"\\gamma",L"\\Gamma", L"A", L"S", L"P", L"L", L"λ", L"h", L"f", L"g",
          L"t", L"R", L"F", overlap_label(), L"a", L"ã", L"b", L"ᵬ", L"E"};
}

std::wstring to_wstring(OpType op) {
  switch (op) {
    case OpType::h:
      return L"h";
    case OpType::f:
      return L"f";
    case OpType::g:
      return L"g";
    case OpType::t:
      return L"t";
    case OpType::l:
      return L"λ";
    case OpType::A:
      return L"A";
    case OpType::L:
      return L"L";
    case OpType::R:
      return L"R";
    case OpType::R12:
      return L"F";
    default:
      throw std::invalid_argument("to_wstring(OpType op): invalid op");
  }
}

OpClass to_class(OpType op) {
  switch (op) {
    case OpType::h:
    case OpType::f:
    case OpType::g:
      return OpClass::gen;
    case OpType::t:
    case OpType::R:
    case OpType::R12:
      return OpClass::ex;
    case OpType::l:
    case OpType::A:
    case OpType::L: return OpClass::deex;
    default:
      throw std::invalid_argument("to_class(OpType op): invalid op");
  }
}

}  // namespace mbpt
}  // namespace sequant
