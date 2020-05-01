#include "SeQuant/domain/mbpt/op.hpp"
#include <stdexcept>

namespace sequant {
namespace mbpt {

std::wstring to_wstring(OpType op) {
  switch (op) {
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
