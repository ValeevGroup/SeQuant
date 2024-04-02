//
// Created by Eduard Valeyev on 3/27/18.
//

#include <SeQuant/core/space.hpp>
#include <SeQuant/core/container.hpp>

#include <bitset>

int32_t sequant::TypeAttr::used_bits = 0xffff;

sequant::container::map<sequant::IndexSpace::Attr, std::wstring>
    sequant::IndexSpace::attr2basekey_{};
sequant::container::map<std::wstring, sequant::IndexSpace::Attr,
                        sequant::IndexSpace::KeyCompare>
    sequant::IndexSpace::key2attr_{};
sequant::container::map<sequant::IndexSpace::Attr, sequant::IndexSpace>
    sequant::IndexSpace::instances_{};
sequant::IndexSpace sequant::IndexSpace::null_instance_{
    sequant::IndexSpace::Attr::null()};

namespace sequant {

std::wstring to_wolfram(const IndexSpace& space) {
  std::wstring result = L"particleSpace[";

  // this is a hack due to partial representation of spaces in SeQuant
  if (space.type() == IndexSpace::active_occupied)
    result += L"occupied";
  else if (space.type() == IndexSpace::active_unoccupied)
    result += L"virtual";
  else if (space.type() == IndexSpace::all)
    result += L"occupied,virtual";
  else if (space.type() == IndexSpace::other_unoccupied)
    result += L"othervirtual";
  else if (space.type() == IndexSpace::complete_unoccupied)
    result += L"virtual,othervirtual";
  else if (space.type() == IndexSpace::complete)
    result += L"occupied,virtual,othervirtual";
  else
    throw std::invalid_argument(
        "to_wolfram(IndexSpace) received a nonstandard space");

  result += L"]";

  return result;
}

std::wstring to_wstring(TypeAttr type) {
  using TypeBitset =
      std::bitset<sizeof(decltype(IndexSpace::TypeAttr::bitset))>;
  if (type == IndexSpace::frozen_occupied) {
    return L"frozen_occupied";
  } else if (type == IndexSpace::inactive_occupied) {
    return L"inactive_occupied";
  } else if (type == IndexSpace::active_occupied) {
    return L"active_occupied";
  } else if (type == IndexSpace::occupied) {
    return L"occupied";
  } else if (type == IndexSpace::active) {
    return L"active";
  } else if (type == IndexSpace::maybe_occupied) {
    return L"maybe_occupied";
  } else if (type == IndexSpace::active_maybe_occupied) {
    return L"active_maybe_occupied";
  } else if (type == IndexSpace::active_unoccupied) {
    return L"active_unoccupied";
  } else if (type == IndexSpace::inactive_unoccupied) {
    return L"inactive_unoccupied";
  } else if (type == IndexSpace::unoccupied) {
    return L"unoccupied";
  } else if (type == IndexSpace::maybe_unoccupied) {
    return L"maybe_unoccupied";
  } else if (type == IndexSpace::active_maybe_unoccupied) {
    return L"active_maybe_unoccupied";
  } else if (type == IndexSpace::all_active) {
    return L"all_active";
  } else if (type == IndexSpace::all) {
    return L"all";
  } else if (type == IndexSpace::other_unoccupied) {
    return L"other_unoccupied";
  } else if (type == IndexSpace::complete_unoccupied) {
    return L"complete_unoccupied";
  } else if (type == IndexSpace::complete_maybe_unoccupied) {
    return L"complete_maybe_unoccupied";
  } else if (type == IndexSpace::complete) {
    return L"complete";
  } else {
    std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;
    return std::wstring(L"Custom (") +
           converter.from_bytes(TypeBitset(type.bitset).to_string()) + L")";
  }
}

std::wstring to_wstring(QuantumNumbersAttr qns) {
  using QNBitset =
      std::bitset<sizeof(decltype(IndexSpace::QuantumNumbersAttr::bitset))>;

  if (qns == IndexSpace::alpha) {
    return L"↑";
  } else if (qns == IndexSpace::beta) {
    return L"↓";
  } else if (qns == IndexSpace::nullqns) {
    return L"";
  } else {
    std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;
    return L"Custom (" +
           converter.from_bytes(QNBitset(qns.bitset).to_string()) + L")";
  }
}

}  // namespace sequant
